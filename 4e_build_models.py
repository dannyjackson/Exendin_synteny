#!/usr/bin/env python3
"""
Protein and CDS for every curated locus (4d curated_loci.tsv), for the gene tree (4g/4h).

Default method 'hsp': the model is built from the forward tblastn alignments that
defined the locus (4c keeps them in loci/scan/<code>/<contig>.blast.tsv). A model is
built from EACH query of the locus's own name group that hit it (any query, if none of
its group did), and the best is kept: hormone-domain motif first, then coverage >=
--min-cover, then 60-150 aa, then collinear, then coverage. For each query:
  - its HSPs inside the locus are taken best-first, skipping any that would add nothing;
  - they are put in query order, and where two HSPs overlap in query coordinates at an
    exon junction, the later one is trimmed so each query position is used once;
  - the protein is the concatenated aligned SUBJECT residues (gaps removed); the CDS is
    the codon under each of those residues, cut from the genome with samtools. Each HSP
    is in frame on its own, so the CDS translates to the protein exactly (checked).
This does not depend on finding splice sites or an ATG: tblastn HSPs end a residue or two
off the real exon boundaries, and splicing them and searching for the longest ORF (the
old 2e approach, which is built for cDNA queries) shifts the frame and gives fragments.
Residues the query does not align to (diverged ends) are not in the model; coverage of
the query is reported.

Method 'rescue' runs lib_rescue_locus.py (the old repo's 2e) in protein mode, as before.
Use it only to compare.

Status per model: 'ok' if it has the hormone-domain motif and covers >= --min-cover of the
query; otherwise the failing checks (no_hormone_motif, low_cover_X, length_N aa,
internal_stop, not_collinear).

Usage
  python3 4e_build_models.py --species-table species.tsv --genomedir ${GENOMEDIR} \
      --curated reciprocal/curated_loci.tsv --loci-dir loci \
      --queries queries/forward_queries.faa --query-table queries/query_table.tsv --outdir models

Writes to --outdir:
  curated_aa.fa, curated_cds.fa   all models, named as in curated_loci.tsv
  models.tsv                      name query method queries_tried hsps q_cover aa_len domain_start status
                                  (query = the one whose model was kept)
"""
import argparse, csv, os, re, subprocess, sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import lib_exendin_names as N  # noqa: E402

RESCUE = os.path.join(HERE, "lib_rescue_locus.py")
MOTIF = re.compile(r"[HY][ASG][DEQ][AGP][ATVIGSM]F[TSI]")      # as in 2f / 2d
_B, _AA = "TCAG", "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODE = {a + b + c: _AA[16 * i + 4 * j + k]
        for i, a in enumerate(_B) for j, b in enumerate(_B) for k, c in enumerate(_B)}


def translate(nt):
    return "".join(CODE.get(nt[i:i + 3], "X") for i in range(0, len(nt) - len(nt) % 3, 3))


def safe(ctg):
    return re.sub(r"[^A-Za-z0-9._-]", "_", ctg)          # as in 4c


def read_hits(path):
    out = []
    if not os.path.exists(path):
        return out
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if len(f) < 12:
            continue
        ss, se = int(f[3]), int(f[4])
        out.append(dict(q=f[0], qs=int(f[1]), qe=int(f[2]), ss=ss, se=se, bits=float(f[8]),
                        qseq=f[10], sseq=f[11], strand="+" if se > ss else "-",
                        lo=min(ss, se), hi=max(ss, se)))
    return out


def columns(h):
    """Aligned subject residues of an HSP: [(query position, residue, subject residue index)]."""
    qpos, k, cols = h["qs"] - 1, 0, []
    for qc, sc in zip(h["qseq"], h["sseq"]):
        if qc != "-":
            qpos += 1
        if sc != "-":
            cols.append((qpos, sc, k)); k += 1
    return cols


def faidx(fasta, regions, minus):
    cmd = ["samtools", "faidx"] + (["-i"] if minus else []) + [fasta] + regions
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode:
        raise RuntimeError(r.stderr.strip()[-200:])
    seqs, cur = [], None
    for line in r.stdout.splitlines():
        if line.startswith(">"):
            cur = []; seqs.append(cur)
        else:
            cur.append(line.strip().upper())
    return ["".join(s) for s in seqs]


def hsp_model(hits, fasta, contig, qlen):
    """(protein, cds, n_hsps, cover, collinear, note) from one query's HSPs at a locus."""
    chosen = []
    for h in sorted(hits, key=lambda h: -h["bits"]):
        used = set(range(h["qs"], h["qe"] + 1))
        if any(used <= set(range(c["qs"], c["qe"] + 1)) for c in chosen):
            continue                                   # wholly inside a better HSP
        chosen.append(h)
    chosen.sort(key=lambda h: h["qs"])
    minus = chosen[0]["strand"] == "-"
    # collinear: later query part further along the gene (downstream on its strand)
    pos = [h["lo"] for h in chosen]
    collinear = pos == sorted(pos, reverse=minus)
    regs = [f"{contig}:{h['lo']}-{h['hi']}" for h in chosen]
    nts = faidx(fasta, regs, minus)
    prot, cds, last_q, covered = [], [], 0, set()
    for h, nt in zip(chosen, nts):
        for qpos, res, k in columns(h):
            if qpos <= last_q:
                continue                               # trimmed: already used by an earlier HSP
            prot.append(res); cds.append(nt[3 * k:3 * k + 3]); covered.add(qpos)
        last_q = max(last_q, h["qe"])
    prot, cds = "".join(prot), "".join(cds)
    note = ""
    if translate(cds) != prot.replace("X", "") and translate(cds).replace("*", "X") != prot.replace("*", "X"):
        note = "cds_translation_mismatch"
    return prot, cds, len(chosen), len(covered) / qlen if qlen else 0, collinear, note


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--species-table", required=True)
    ap.add_argument("--genomedir", required=True)
    ap.add_argument("--curated", required=True, help="curated_loci.tsv from 4d")
    ap.add_argument("--loci-dir", help="4c output (its scan/ holds the forward HSPs); needed for --method hsp")
    ap.add_argument("--queries", required=True, help="forward_queries.faa from 4a")
    ap.add_argument("--query-table", help="query_table.tsv from 4a (fallback query of the same group)")
    ap.add_argument("--method", choices=["hsp", "rescue"], default="hsp")
    ap.add_argument("--min-cover", type=float, default=0.7, help="query coverage for status ok")
    ap.add_argument("--flank", type=int, default=2000, help="rescue method only")
    ap.add_argument("--min-identity", type=float, default=30.0, help="rescue method only")
    ap.add_argument("--evalue", default="1e-2", help="rescue method only")
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()
    if a.method == "hsp" and not a.loci_dir:
        sys.exit("[error] --method hsp needs --loci-dir (the 4c output)")
    os.makedirs(a.outdir, exist_ok=True)

    table = {r["code"]: r for r in N.read_species_table(a.species_table)}
    qs = N.read_fasta(a.queries)
    qgroup = {}
    if a.query_table:
        with open(a.query_table) as fh:
            qgroup = {r["id"]: r["group"] for r in csv.DictReader(fh, delimiter="\t")}
    with open(a.curated) as fh:
        cur = list(csv.DictReader(fh, delimiter="\t"))
    rows, aa_out, cds_out, cache = [], [], [], {}
    for r in cur:
        code, name, q = r["species"], r["name"], r["best_forward_query"]
        files = N.find_genome_files(a.genomedir, code, table[code]["accession"])
        lo, hi, st = int(r["g_start"]), int(r["g_end"]), r["strand"]
        prot = cds = ""
        nh, cover, collinear, note, tried = 0, 0.0, True, "", 1
        try:
            if a.method == "hsp":
                key = (code, r["contig"])
                if key not in cache:
                    cache[key] = read_hits(os.path.join(a.loci_dir, "scan", code, safe(r["contig"]) + ".blast.tsv"))
                here = [h for h in cache[key] if h["strand"] == st and h["hi"] >= lo and h["lo"] <= hi]
                byq = defaultdict(list)
                for h in here:
                    byq[h["q"]].append(h)
                if not byq:
                    raise RuntimeError("no forward HSPs at this locus (re-run 4c?)")
                # try every query of the locus's own group (then any query) and keep the best
                # model: a long but faulty query (e.g. a mis-annotated first exon) can win on
                # bitscore and still give a worse model than a correct shorter one
                same = [x for x in byq if qgroup.get(x) == r["group"]] or list(byq)
                best = None
                for x in same:
                    try:
                        mod = hsp_model(byq[x], files["fasta"], r["contig"], len(qs[x][1]) if x in qs else 0)
                    except RuntimeError as e:
                        last_err = e; continue
                    p_, c_, n_, cov_, col_, note_ = mod
                    key = (MOTIF.search(p_) is not None, cov_ >= a.min_cover, 60 <= len(p_) <= 150,
                           col_, cov_, len(p_))
                    if best is None or key > best[0]:
                        best = (key, x, mod)
                if best is None:
                    raise last_err
                _, q, (prot, cds, nh, cover, collinear, note) = best
                tried = len(same)
            else:
                if files["fasta"] not in cache:
                    cache[files["fasta"]] = {l.split("\t")[0]: int(l.split("\t")[1])
                                             for l in open(files["fasta"] + ".fai")}
                clen = cache[files["fasta"]].get(r["contig"], 10 ** 12)
                qf = os.path.join(a.outdir, f"{name}.query.faa")
                open(qf, "w").write(f">{q}\n{qs[q][1]}\n")
                p = subprocess.run([sys.executable, RESCUE, "--genome", files["fasta"], "--contig", r["contig"],
                                    "--region", f"{max(1, lo - a.flank)}-{min(clen, hi + a.flank)}",
                                    "--query", qf, "--protein", "--name", name, "--outdir", a.outdir,
                                    "--min-identity", str(a.min_identity), "--evalue", a.evalue],
                                   capture_output=True, text=True)
                if p.returncode:
                    raise RuntimeError((p.stderr.strip().splitlines() or ["failed"])[-1])
                prot = N.read_fasta(os.path.join(a.outdir, f"{name}.aa.fa"))[name][1]
                cds = N.read_fasta(os.path.join(a.outdir, f"{name}.cds.fa"))[name][1]
                nh = sum(1 for _ in open(os.path.join(a.outdir, f"{name}.exons.tsv"))) - 1
                cover = len(prot) / len(qs[q][1]) if q in qs else 0
        except Exception as e:                          # noqa: BLE001
            rows.append([name, q, a.method, 0, 0, "0.00", 0, "", f"failed: {str(e)[:80]}"])
            N.warn(f"{name}: no model: {e}")
            continue
        m = MOTIF.search(prot)
        dom = f"{prot[m.start()]}@{m.start() + 1}" if m else "none"
        st_ = []
        if not m:
            st_.append("no_hormone_motif")
        if cover < a.min_cover:
            st_.append(f"low_cover_{cover:.2f}")
        if not (60 <= len(prot) <= 150):
            st_.append(f"length_{len(prot)}aa")
        if "*" in prot:
            st_.append("internal_stop")
        if not collinear:
            st_.append("not_collinear")
        if note:
            st_.append(note)
        rows.append([name, q, a.method, tried, nh, f"{cover:.2f}", len(prot), dom, ";".join(st_) or "ok"])
        aa_out.append((name, prot)); cds_out.append((name, cds))

    with open(os.path.join(a.outdir, "models.tsv"), "w") as fh:
        fh.write("name\tquery\tmethod\tqueries_tried\thsps\tq_cover\taa_len\tdomain_start\tstatus\n")
        for row in rows:
            fh.write("\t".join(map(str, row)) + "\n")
    for fn, recs in (("curated_aa.fa", aa_out), ("curated_cds.fa", cds_out)):
        with open(os.path.join(a.outdir, fn), "w") as fh:
            for n, s in recs:
                fh.write(f">{n}\n{s}\n")
    ok = sum(1 for r in rows if r[-1] == "ok")
    N.info(f"{len(aa_out)}/{len(cur)} loci modelled ({a.method}); {ok} ok, {len(rows) - ok} need a look "
           f"(models.tsv status != ok)")


if __name__ == "__main__":
    main()
