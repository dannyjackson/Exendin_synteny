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

Missing exons. A short exon cannot pass the forward search's e-value against a whole
genome: the Heloderma exendins' ~10-residue, Pro-rich last exon (...GPSSGAPPPS) was missing
from every model, which cut each mature peptide ~30 residues in and pushed exendin-1/-2/-3
below 5a's 80% identity although every residue present matched. So for each model, every
stretch of >= --fill-min-gap query residues that no HSP covers (a missing first or last
exon, or a gap between HSPs) is searched with tblastn in its own small window: up to
--fill-extend-bp beyond the first/last HSP, or between the two HSPs around an internal gap.
The best hit on the gene's strand (>= --fill-min-bits) is added and the model rebuilt, and
the filled model is kept only if it ranks at least as well as the unfilled one. The
'filled' column of models.tsv counts the HSPs added. --no-fill turns this off.

Method 'rescue' runs lib_rescue_locus.py (the old repo's 2e) in protein mode, as before.
Use it only to compare.

Hormone domain: the old repo's motif [HY][ASG][DEQ][AGP][ATVIGSM]F[TSI] (fitted to the
Heloderma exendins), else a relaxed one, [HY][ASGE][DEQ][AGP].[FYLI], which squamate GCGL
and diverged copies need; domain_start says '(relaxed)' when only that matched. Residues
before the first Met (within --trim-to-met) are cut: they come from a query whose model
starts upstream of the start codon.

Status per model: 'ok' if it has a hormone domain and covers >= --min-cover of the
query; otherwise the failing checks (no_hormone_motif, low_cover_X, length_N aa,
internal_stop, not_collinear).

Usage
  python3 4e_build_models.py --species-table species.tsv --genomedir ${GENOMEDIR} \
      --curated reciprocal/curated_loci.tsv --loci-dir loci \
      --queries queries/forward_queries.faa --query-table queries/query_table.tsv --outdir models

Writes to --outdir:
  curated_aa.fa, curated_cds.fa   all models, named as in curated_loci.tsv
  models.tsv                      name query method queries_tried hsps filled q_cover aa_len domain_start status
                                  (query = the one whose model was kept)
"""
import argparse, csv, os, re, subprocess, sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import lib_exendin_names as N  # noqa: E402

RESCUE = os.path.join(HERE, "lib_rescue_locus.py")
MOTIF = re.compile(r"[HY][ASG][DEQ][AGP][ATVIGSM]F[TSI]")      # as in 2f / 2d (fitted to Heloderma)
# across squamates positions 5-7 vary (GCGL YSEGTLAS, PhrGui Ex4 YSEADLGS, ShiCro HSDGTYTS);
# His/Tyr, positions 2-4 and a hydrophobic residue at 6 stay conserved
MOTIF_RELAXED = re.compile(r"[HY][ASGE][DEQ][AGP].[FYLI]")


def find_domain(prot):
    """(match, 'strict'|'relaxed') for the hormone-domain start, or (None, '')."""
    m = MOTIF.search(prot)
    if m:
        return m, "strict"
    m = MOTIF_RELAXED.search(prot)
    return (m, "relaxed") if m else (None, "")
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


def select(hits):
    """HSPs best-first, skipping any wholly inside a better one; returned in query order."""
    chosen = []
    for h in sorted(hits, key=lambda h: -h["bits"]):
        used = set(range(h["qs"], h["qe"] + 1))
        if any(used <= set(range(c["qs"], c["qe"] + 1)) for c in chosen):
            continue
        chosen.append(h)
    return sorted(chosen, key=lambda h: h["qs"])


def fill_gaps(chosen, fasta, contig, qseq, clen, a, tmpdir):
    """Search the parts of the query no HSP covers, each in its own small genomic window.

    A short exon (the ~10-residue, Pro-rich last exon of the exendin precursors) cannot
    reach a useful e-value against a whole genome, so the genome-wide forward search
    misses it; against the few kb where it must be, it can. Windows, in gene orientation:
    leading gap = up to --fill-extend-bp upstream of the first HSP; trailing gap = the same
    downstream of the last; internal gap = between the two HSPs around it.
    Returns HSP dicts in the forward-hit format (contig coordinates).
    """
    qlen, minus = len(qseq), chosen[0]["strand"] == "-"
    first, last, ext = chosen[0], chosen[-1], a.fill_extend_bp
    gaps = []
    if first["qs"] - 1 >= a.fill_min_gap:
        win = (first["hi"] + 1, first["hi"] + ext) if minus else (first["lo"] - ext, first["lo"] - 1)
        gaps.append((1, first["qs"] - 1, win))
    if qlen - last["qe"] >= a.fill_min_gap:
        win = (last["lo"] - ext, last["lo"] - 1) if minus else (last["hi"] + 1, last["hi"] + ext)
        gaps.append((last["qe"] + 1, qlen, win))
    for p, n in zip(chosen, chosen[1:]):
        if n["qs"] - p["qe"] - 1 >= a.fill_min_gap:
            win = (n["hi"] + 1, p["lo"] - 1) if minus else (p["hi"] + 1, n["lo"] - 1)
            gaps.append((p["qe"] + 1, n["qs"] - 1, win))
    out = []
    for k, (gs, ge, (wlo, whi)) in enumerate(gaps):
        wlo, whi = max(1, wlo), min(clen, whi)
        if whi - wlo < 30:
            continue
        q0 = max(1, gs - a.fill_pad)
        seg = qseq[q0 - 1:min(qlen, ge + a.fill_pad)]
        qf, sf = os.path.join(tmpdir, "gap.faa"), os.path.join(tmpdir, "win.fa")
        open(qf, "w").write(f">gap\n{seg}\n")
        open(sf, "w").write(f">win\n{faidx(fasta, [f'{contig}:{wlo}-{whi}'], False)[0]}\n")
        r = subprocess.run(["tblastn", "-query", qf, "-subject", sf, "-seg", "no", "-evalue", str(a.fill_evalue),
                            "-word_size", "2", "-comp_based_stats", "0", "-outfmt",
                            "6 qseqid qstart qend sstart send pident length evalue bitscore sframe qseq sseq"],
                           capture_output=True, text=True)
        if r.returncode:
            raise RuntimeError("tblastn (gap fill) failed: " + r.stderr.strip()[-150:])
        best = None
        for line in r.stdout.splitlines():
            f = line.split("\t")
            if len(f) < 12:
                continue
            ss, se, bits = int(f[3]), int(f[4]), float(f[8])
            if ((se < ss) != minus) or bits < a.fill_min_bits:
                continue
            if best is None or bits > best[0]:
                best = (bits, f)
        if best:
            bits, f = best
            ss, se = wlo + int(f[3]) - 1, wlo + int(f[4]) - 1
            out.append(dict(q="fill", qs=q0 + int(f[1]) - 1, qe=q0 + int(f[2]) - 1, ss=ss, se=se, bits=bits,
                            qseq=f[10], sseq=f[11], strand="-" if minus else "+", lo=min(ss, se), hi=max(ss, se)))
    return out


def hsp_model(hits, fasta, contig, qlen):
    """(protein, cds, n_hsps, cover, collinear, note) from one query's HSPs at a locus."""
    chosen = select(hits)
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
    ap.add_argument("--no-fill", dest="fill", action="store_false",
                    help="do not search for missing exons (default: search, needs tblastn)")
    ap.add_argument("--fill-evalue", type=float, default=10.0, help="e-value cut-off of the gap search")
    ap.add_argument("--fill-min-bits", type=float, default=15.0, help="weakest gap hit accepted")
    ap.add_argument("--fill-extend-bp", type=int, default=10000,
                    help="how far beyond the first/last HSP to look for a missing first/last exon")
    ap.add_argument("--fill-min-gap", type=int, default=5, help="unaligned query residues worth a search")
    ap.add_argument("--fill-pad", type=int, default=3, help="query residues of overlap added to each gap search")
    ap.add_argument("--trim-to-met", type=int, default=20,
                    help="cut residues before the first Met if it lies within this many residues (0 = never)")
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
    rows, aa_out, cds_out, cache, clen_cache = [], [], [], {}, {}
    for r in cur:
        code, name, q = r["species"], r["name"], r["best_forward_query"]
        files = N.find_genome_files(a.genomedir, code, table[code]["accession"])
        lo, hi, st = int(r["g_start"]), int(r["g_end"]), r["strand"]
        prot = cds = ""
        nh, cover, collinear, note, tried, filled = 0, 0.0, True, "", 1, 0
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
                def rank(mod):
                    p_, c_, n_, cov_, col_, note_ = mod
                    return (find_domain(p_)[0] is not None, cov_ >= a.min_cover, 60 <= len(p_) <= 150,
                            col_, cov_, len(p_))
                if a.fill and (code, r["contig"]) not in clen_cache:
                    fai = files["fasta"] + ".fai"
                    clen_cache[(code, r["contig"])] = next((int(l.split("\t")[1]) for l in open(fai)
                                                    if l.split("\t")[0] == r["contig"]), 10 ** 12)
                best = None
                for x in same:
                    ql = len(qs[x][1]) if x in qs else 0
                    try:
                        mod = hsp_model(byq[x], files["fasta"], r["contig"], ql)
                    except RuntimeError as e:
                        last_err = e; continue
                    nfill = 0
                    if a.fill and x in qs:
                        extra = fill_gaps(select(byq[x]), files["fasta"], r["contig"], qs[x][1],
                                          clen_cache[(code, r["contig"])], a, a.outdir)
                        if extra:
                            try:
                                mod2 = hsp_model(byq[x] + extra, files["fasta"], r["contig"], ql)
                                if rank(mod2) >= rank(mod):
                                    mod, nfill = mod2, len(extra)
                            except RuntimeError:
                                pass
                    if best is None or (rank(mod), nfill) > best[0]:
                        best = ((rank(mod), nfill), x, mod, nfill)
                if best is None:
                    raise last_err
                _, q, (prot, cds, nh, cover, collinear, note), filled = best
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
            rows.append([name, q, a.method, 0, 0, 0, "0.00", 0, "", f"failed: {str(e)[:80]}"])
            N.warn(f"{name}: no model: {e}")
            continue
        if a.trim_to_met and prot and not prot.startswith("M") and "M" in prot[:a.trim_to_met]:
            k = prot.index("M")                # residues before the start Met come from the query's own
            prot, cds = prot[k:], cds[3 * k:]  # upstream N-terminus, not from this gene
        m, how = find_domain(prot)
        dom = f"{prot[m.start()]}@{m.start() + 1}{'' if how == 'strict' else '(relaxed)'}" if m else "none"
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
        rows.append([name, q, a.method, tried, nh, filled, f"{cover:.2f}", len(prot), dom, ";".join(st_) or "ok"])
        aa_out.append((name, prot)); cds_out.append((name, cds))

    with open(os.path.join(a.outdir, "models.tsv"), "w") as fh:
        fh.write("name\tquery\tmethod\tqueries_tried\thsps\tfilled\tq_cover\taa_len\tdomain_start\tstatus\n")
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
