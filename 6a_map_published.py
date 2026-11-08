#!/usr/bin/env python3
"""
Map each published Heloderma venom peptide (6_published_peptides.fa: helospectin I and II =
exendin-1, helodermin = exendin-2, exendin-4) onto ONE genome (default HelSus, the Gila
monster), codon by codon, for the figure of 6b.

Why a separate step. Steps 4-5 name a gene by comparing its predicted mature peptide with
the published ones. This step goes the other way: it starts from the purified peptide and
asks where in the genome each residue is encoded, so the gene behind every published name
can be seen directly, with no gene model in between.

Method, per peptide
  1. Genome-wide tblastn (-seg no, -comp_based_stats 0, e <= --evalue). Short queries
     need these: with the defaults the Pro-rich residues are masked. Only the main exon
     has to be found here; the short exons are found in step 3.
  2. Hits on the same contig and strand within --cluster-kb are grouped into candidate
     loci; each is scored by its best collinear chain of HSPs around the locus's strongest
     HSP (the anchor): query order matches genomic order on the strand, introns
     <= --max-intron-bp, at most --max-exons HSPs, each other HSP >= --min-exon-bits.
  3. The best --top-loci candidates are searched again, each in its own window (chain
     +/- --refine-bp) with a sensitive tblastn (word size 2, e <= 1000). This finds the
     ~10-residue last exon that a genome-wide search cannot (as in 4e's gap fill).
  4. The candidate with the most identical residues (then the highest chain score) is
     the mapping. Each query residue gets one codon: where two HSPs overlap at an exon
     junction, the HSP whose residue matches the query wins, then the higher-scoring one.
  4b. tblastn stops at mismatched ends, so up to 6 unaligned residues at either end of
     the peptide get the codons that continue the aligned exon (read from the genome).
  5. The codons just before residue 1 and just after the last residue are translated
     from the genome (same exon): the processing site before the peptide and, after it,
     whether the genome encodes more residues (e.g. helospectin I's Ser38 after
     helospectin II's end) or a Gly amidation signal.
Genes from the NCBI GFF that overlap each mapping are reported (which exon), so an
'uncharacterized protein LOC...' carrying a published peptide is named here.

Usage
  python3 6a_map_published.py --species-table ../species.tsv --genomedir ${GENOMEDIR} \
      --code HelSus --peptides ${SCRIPTDIR}/6_published_peptides.fa --threads 8 --outdir map

Writes to --outdir
  peptides.faa         the queries used
  genome_hits.tsv      genome-wide tblastn HSPs
  candidates.tsv       peptide rank contig strand start end genome_bits refined_bits identical mapped chosen
  map.tsv              one row per peptide: the chosen locus, identity, exons, introns,
                       flanking residues, overlapping NCBI genes, the runner-up locus
  hsps.tsv             the HSPs (exons) of each mapping
  residues.tsv         one row per peptide residue: codon coordinates, genome residue, status
                       (identical / mismatch / stop / gap / unmapped), genome residues
                       inserted after it, ext=yes if the codon was read by extending the
                       alignment over an unaligned end (step 4b)
Needs tblastn, makeblastdb and samtools (exendin_env).
"""
import argparse, os, subprocess, sys, tempfile
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import lib_exendin_names as N  # noqa: E402
import lib_gff as G            # noqa: E402

OUTFMT = "6 qseqid saccver qstart qend sstart send pident length evalue bitscore qseq sseq"
CODE = {a + b + c: aa for (a, b, c), aa in zip(
    [(x, y, z) for x in "TCAG" for y in "TCAG" for z in "TCAG"],
    "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG")}
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


# ---------------- small helpers ----------------
def run(cmd):
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode:
        why = r.stderr.strip()[-600:] or r.stdout.strip()[-300:]
        if r.returncode < 0 or r.returncode in (137, 139):
            why += (f"\n(killed: exit {r.returncode}. Usually out of memory or the job's time limit: "
                    "ask for more, e.g. sinteractive --mem=32g --cpus-per-task=8)")
        sys.exit(f"[error] {' '.join(cmd)}\nexit {r.returncode}\n{why}")
    return r.stdout


def translate(nt):
    return "".join(CODE.get(nt[i:i + 3].upper(), "X") for i in range(0, len(nt) - len(nt) % 3, 3))


def faidx(fasta, contig, lo, hi):
    """Plus-strand sequence of contig:lo-hi (1-based, inclusive)."""
    out = run(["samtools", "faidx", fasta, f"{contig}:{lo}-{hi}"])
    return "".join(l.strip() for l in out.splitlines() if not l.startswith(">")).upper()


def read_peptides(path):
    peps = []
    for pid, (head, seq) in N.read_fasta(path).items():
        f = [x.strip() for x in head.split("|")]
        peps.append(dict(id=pid, group=f[1] if len(f) > 1 else "published",
                         label=f[2] if len(f) > 2 else pid, source=f[3] if len(f) > 3 else "",
                         seq=seq.upper().replace("*", "")))
    return peps


def acc_of(sid):
    """'ref|NC_139459.1|' (sseqid of a -parse_seqids database) -> 'NC_139459.1'."""
    parts = [x for x in sid.split("|") if x]
    return parts[1] if len(parts) >= 2 and parts[0] in ("ref", "gb", "emb", "dbj", "lcl", "gnl") else sid


def parse(out, offset=None):
    """tblastn -outfmt OUTFMT -> HSP dicts; offset=(contig, window_start) for -subject searches."""
    hits = []
    for line in out.splitlines():
        f = line.rstrip("\n").split("\t")
        if len(f) < 12:
            continue
        ss, se = int(f[4]), int(f[5])
        contig = acc_of(f[1])
        if offset:
            contig, w0 = offset
            ss, se = w0 + ss - 1, w0 + se - 1
        hits.append(dict(q=f[0], contig=contig, qs=int(f[2]), qe=int(f[3]), ss=ss, se=se,
                         pident=float(f[6]), evalue=float(f[8]), bits=float(f[9]), qseq=f[10], sseq=f[11],
                         strand="+" if se >= ss else "-", lo=min(ss, se), hi=max(ss, se)))
    return hits


# ---------------- chaining and residue mapping ----------------
def chain(hits, strand, max_intron, anchor=None, max_exons=4, min_bits=14.0, min_aa=4, max_qoverlap=8):
    """Best collinear chain (sum of bits) of at most max_exons HSPs: query order = genomic
    order along the strand, consecutive HSPs <= max_intron apart. If anchor is given the
    chain must contain it, and every other HSP needs >= min_bits and >= min_aa residues
    (so chance 3-5-residue matches in the window cannot be strung into extra 'exons')."""
    hs = [h for h in hits if h["strand"] == strand and
          (h is anchor or (h["bits"] >= min_bits and h["qe"] - h["qs"] + 1 >= min_aa))]
    hs.sort(key=lambda h: (h["qs"], h["qe"]))
    if not hs:
        return []

    def ok(a, b):                                       # can b follow a?
        if not (a["qs"] < b["qs"] and a["qe"] < b["qe"] and a["qe"] - b["qs"] + 1 <= max_qoverlap):
            return False
        gap = (b["lo"] - a["hi"]) if strand == "+" else (a["lo"] - b["hi"])
        return -30 <= gap <= max_intron
    n = len(hs)
    NEG = float("-inf")
    # fw[i][k]: best chain of k HSPs ending at i; bw[i][k]: best chain of k HSPs starting at i
    fw = [[NEG] * (max_exons + 1) for _ in range(n)]
    fp = [[None] * (max_exons + 1) for _ in range(n)]
    for i in range(n):
        fw[i][1] = hs[i]["bits"]
        for j in range(i):
            if ok(hs[j], hs[i]):
                for k in range(2, max_exons + 1):
                    if fw[j][k - 1] + hs[i]["bits"] > fw[i][k]:
                        fw[i][k], fp[i][k] = fw[j][k - 1] + hs[i]["bits"], j
    bw = [[NEG] * (max_exons + 1) for _ in range(n)]
    bp = [[None] * (max_exons + 1) for _ in range(n)]
    for i in range(n - 1, -1, -1):
        bw[i][1] = hs[i]["bits"]
        for j in range(i + 1, n):
            if ok(hs[i], hs[j]):
                for k in range(2, max_exons + 1):
                    if bw[j][k - 1] + hs[i]["bits"] > bw[i][k]:
                        bw[i][k], bp[i][k] = bw[j][k - 1] + hs[i]["bits"], j
    if anchor is not None and any(h is anchor for h in hs):
        a = next(i for i, h in enumerate(hs) if h is anchor)
        best = max(((fw[a][k1] + bw[a][k2] - hs[a]["bits"], k1, k2)
                    for k1 in range(1, max_exons + 1) for k2 in range(1, max_exons + 2 - k1)
                    if fw[a][k1] > NEG and bw[a][k2] > NEG), key=lambda t: t[0])
        _, k1, k2 = best
        left, i, k = [], a, k1
        while i is not None and k >= 1:
            left.append(hs[i]); i, k = fp[i][k], k - 1
        right, i, k = [], bp[a][k2], k2 - 1
        while i is not None and k >= 1:
            right.append(hs[i]); i, k = bp[i][k], k - 1
        return left[::-1] + right
    _, i, k = max((fw[i][k], i, k) for i in range(n) for k in range(1, max_exons + 1) if fw[i][k] > NEG)
    out = []
    while i is not None and k >= 1:
        out.append(hs[i]); i, k = fp[i][k], k - 1
    return out[::-1]


def columns(h):
    """[(qpos, subject residue or '-', codon_lo, codon_hi)] and {qpos: genome residues inserted after it}."""
    qpos, k, cols, ins = h["qs"] - 1, 0, [], defaultdict(int)
    for qc, sc in zip(h["qseq"], h["sseq"]):
        if qc != "-" and sc != "-":
            qpos += 1
            if h["strand"] == "+":
                c0 = h["ss"] + 3 * k
            else:
                c0 = h["ss"] - 3 * k - 2
            cols.append((qpos, sc.upper(), c0, c0 + 2)); k += 1
        elif qc != "-":
            qpos += 1
            cols.append((qpos, "-", None, None))
        elif sc != "-":
            ins[qpos] += 1; k += 1
    return cols, ins


def residue_map(ch, seq):
    """One codon per query residue. Returns (rows, identical, mapped).

    Consecutive HSPs that overlap in the query (an exon junction, or an HSP running on into
    the intron) are split at ONE cut: positions up to the cut come from the earlier HSP,
    the rest from the later one, the cut placed where the two together match the query
    best. So every exon is a contiguous run of residues."""
    ch = sorted(ch, key=lambda h: h["qs"])
    cols = [columns(h) for h in ch]
    bounds = [[0, 10 ** 9] for _ in ch]
    for i in range(len(ch) - 1):
        A, B = ch[i], ch[i + 1]
        if A["qe"] < B["qs"]:
            cut = A["qe"]
        else:
            best = None
            for c in range(max(B["qs"] - 1, bounds[i][0]), A["qe"] + 1):
                sc = sum(1 for q, r, _, _ in cols[i][0] if B["qs"] <= q <= c and r == seq[q - 1]) + \
                     sum(1 for q, r, _, _ in cols[i + 1][0] if c < q <= A["qe"] and r == seq[q - 1])
                if best is None or sc > best[0]:
                    best = (sc, c)
            cut = best[1]
        bounds[i][1] = cut
        bounds[i + 1][0] = cut
    pick, ins_total = {}, defaultdict(int)
    for ei, ((lo_b, hi_b), (cl, ins)) in enumerate(zip(bounds, cols), 1):
        for qpos, res, c0, c1 in cl:
            if lo_b < qpos <= hi_b:
                pick[qpos] = (ei, res, c0, c1)
        for q, n in ins.items():
            if lo_b < q <= hi_b:
                ins_total[q] += n
    rows, ident, mapped = [], 0, 0
    for q in range(1, len(seq) + 1):
        qr = seq[q - 1]
        if q not in pick:
            rows.append(dict(qpos=q, qres=qr, exon="", gres="", codon_start="", codon_end="",
                             status="unmapped", ins_after=ins_total.get(q, 0)))
            continue
        ei, res, c0, c1 = pick[q]
        st = ("identical" if res == qr else "gap" if res == "-" else "stop" if res == "*" else "mismatch")
        if res != "-":
            mapped += 1
        ident += st == "identical"
        rows.append(dict(qpos=q, qres=qr, exon=ei, gres=res, codon_start=c0 or "", codon_end=c1 or "",
                         status=st, ins_after=ins_total.get(q, 0)))
    return rows, ident, mapped


def extend_ends(rows, fasta, contig, strand, clen, max_ext=6):
    """tblastn stops an alignment at mismatched ends (exendin-3's HSD against an exendin-4
    gene's HGE). Unaligned residues at either end of the peptide (up to max_ext) are given
    the codons that continue the first / last aligned exon, read from the genome, so they
    show as identical or mismatch instead of unmapped. Marked ext=yes."""
    m = [i for i, r in enumerate(rows) if r["codon_start"] != ""]
    if not m:
        return rows, 0
    added = 0
    for side in ("lead", "trail"):
        if side == "lead":
            k = m[0]
            idx = list(range(k - 1, -1, -1))
            ref = rows[k]
        else:
            k = m[-1]
            idx = list(range(k + 1, len(rows)))
            ref = rows[k]
        if not idx or len(idx) > max_ext:
            continue
        for n, i in enumerate(idx, 1):
            step = -3 * n if side == "lead" else 3 * n
            if strand == "-":
                step = -step
            c0 = ref["codon_start"] + step
            c1 = c0 + 2
            if c0 < 1 or c1 > clen:
                break
            nt = faidx(fasta, contig, c0, c1)
            if strand == "-":
                nt = nt.translate(COMP)[::-1]
            res = translate(nt)
            r = rows[i]
            r.update(codon_start=c0, codon_end=c1, exon=ref["exon"], gres=res,
                     status="identical" if res == r["qres"] else "stop" if res == "*" else "mismatch", ext="yes")
            added += 1
    return rows, added


def flanks(fasta, contig, strand, rows, clen, n_before=2, n_after=3):
    """Genome residues just before residue 1 and just after the last mapped residue (same exon)."""
    m = [r for r in rows if r["codon_start"] != ""]
    if not m:
        return "", ""
    first, last = m[0], m[-1]
    if strand == "+":
        b0, b1 = first["codon_start"] - 3 * n_before, first["codon_start"] - 1
        a0, a1 = last["codon_end"] + 1, last["codon_end"] + 3 * n_after
    else:
        b0, b1 = first["codon_end"] + 1, first["codon_end"] + 3 * n_before
        a0, a1 = last["codon_start"] - 3 * n_after, last["codon_start"] - 1
    out = []
    for lo, hi in ((b0, b1), (a0, a1)):
        lo, hi = max(1, lo), min(clen, hi)
        if hi - lo < 2:
            out.append(""); continue
        nt = faidx(fasta, contig, lo, hi)
        if strand == "-":
            nt = nt.translate(COMP)[::-1]
        out.append(translate(nt))
    return out[0], out[1]


def cluster(hits, gap_bp):
    groups = []
    for key in sorted({(h["contig"], h["strand"]) for h in hits}):
        hs = sorted((h for h in hits if (h["contig"], h["strand"]) == key), key=lambda h: h["lo"])
        cur = [hs[0]]
        for h in hs[1:]:
            if h["lo"] - max(x["hi"] for x in cur) <= gap_bp:
                cur.append(h)
            else:
                groups.append(cur); cur = [h]
        groups.append(cur)
    return groups


# ---------------- main ----------------
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--species-table", required=True)
    ap.add_argument("--genomedir", required=True)
    ap.add_argument("--code", default="HelSus", help="genome to map onto [HelSus]")
    ap.add_argument("--peptides", default=os.path.join(HERE, "6_published_peptides.fa"))
    ap.add_argument("--blastdb", help="nucleotide BLAST db [${GENOMEDIR}/blastdb/<code>_<accession>, as 4b]")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--evalue", type=float, default=1.0, help="genome-wide e-value cut-off [1]")
    ap.add_argument("--max-targets", type=int, default=500, help="genome-wide -max_target_seqs [500]")
    ap.add_argument("--min-bits", type=float, default=20.0, help="genome-wide HSPs below this are ignored [20]")
    ap.add_argument("--cluster-kb", type=float, default=10.0, help="hits this close form one candidate locus [10]")
    ap.add_argument("--max-intron-bp", type=int, default=10000)
    ap.add_argument("--max-exons", type=int, default=4, help="exons a peptide may span [4]")
    ap.add_argument("--min-exon-bits", type=float, default=14.0,
                    help="bits an HSP other than the strongest needs to join a mapping [14]")
    ap.add_argument("--top-loci", type=int, default=5, help="candidate loci refined per peptide [5]")
    ap.add_argument("--refine-bp", type=int, default=5000)
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()

    os.makedirs(a.outdir, exist_ok=True)
    table = {r["code"]: r for r in N.read_species_table(a.species_table)}
    if a.code not in table:
        sys.exit(f"[error] {a.code} not in {a.species_table}")
    acc = table[a.code]["accession"]
    files = N.find_genome_files(a.genomedir, a.code, acc)
    if not files["fasta"]:
        sys.exit(f"[error] no genome FASTA for {a.code} under {a.genomedir} (run 1b)")
    if not os.path.exists(files["fasta"] + ".fai"):
        run(["samtools", "faidx", files["fasta"]])
    clen = {l.split("\t")[0]: int(l.split("\t")[1]) for l in open(files["fasta"] + ".fai")}

    peps = read_peptides(a.peptides)
    qfa = os.path.join(a.outdir, "peptides.faa")
    with open(qfa, "w") as fh:
        for p in peps:
            fh.write(f">{p['id']}\n{p['seq']}\n")
    N.info(f"{len(peps)} peptides: " + ", ".join(f"{p['id']} ({len(p['seq'])} aa, {p['group']})" for p in peps))

    db = a.blastdb or os.path.join(a.genomedir, "blastdb", f"{a.code}_{acc}")
    if not any(os.path.exists(db + x) for x in (".ndb", ".nal", ".nin", ".00.nin")):
        os.makedirs(os.path.dirname(db), exist_ok=True)
        N.info(f"makeblastdb {db}")
        run(["makeblastdb", "-in", files["fasta"], "-dbtype", "nucl", "-parse_seqids", "-out", db,
             "-title", f"{a.code}_{acc}"])

    raw = os.path.join(a.outdir, "genome_hits.tsv")
    if not os.path.exists(raw):
        N.info(f"genome-wide tblastn against {db} ...")
        run(["tblastn", "-query", qfa, "-db", db, "-seg", "no", "-comp_based_stats", "0",
             "-evalue", str(a.evalue), "-max_target_seqs", str(a.max_targets), "-num_threads", str(a.threads),
             "-outfmt", OUTFMT, "-out", raw + ".tmp"])
        os.replace(raw + ".tmp", raw)
    else:
        N.info(f"re-using {raw} (delete it to search again)")
    allhits = [h for h in parse(open(raw).read()) if h["bits"] >= a.min_bits]

    gff_cache = {}

    def genes_on(contig):
        if contig not in gff_cache:
            gff_cache[contig] = G.load_region(files["gff"], contig, 1, 10 ** 12) if files["gff"] else []
        return gff_cache[contig]

    cand_rows, map_rows, hsp_rows, res_rows = [], [], [], []
    tmp = tempfile.mkdtemp(prefix="6a_")
    for p in peps:
        hits = [h for h in allhits if h["q"] == p["id"]]
        if not hits:
            N.warn(f"{p['id']}: no genome-wide hit >= {a.min_bits} bits")
            map_rows.append(dict(peptide=p["id"], group=p["group"], label=p["label"], length=len(p["seq"]),
                                 contig="NA", note="no genome-wide hit"))
            continue
        loci = []
        for grp in cluster(hits, a.cluster_kb * 1000):
            anc = max(grp, key=lambda h: h["bits"])
            ch = chain(grp, grp[0]["strand"], a.max_intron_bp, anchor=anc, max_exons=a.max_exons,
                       min_bits=a.min_exon_bits)
            loci.append(dict(contig=grp[0]["contig"], strand=grp[0]["strand"], lo=min(h["lo"] for h in ch),
                             hi=max(h["hi"] for h in ch), anchor=anc, gbits=sum(h["bits"] for h in ch)))
        loci.sort(key=lambda L: -L["gbits"])
        results = []
        for rank, L in enumerate(loci[:a.top_loci], 1):
            w0 = max(1, L["lo"] - a.refine_bp)
            w1 = min(clen.get(L["contig"], L["hi"] + a.refine_bp), L["hi"] + a.refine_bp)
            sf = os.path.join(tmp, "win.fa")
            open(sf, "w").write(f">win\n{faidx(files['fasta'], L['contig'], w0, w1)}\n")
            pf = os.path.join(tmp, "pep.faa")
            open(pf, "w").write(f">{p['id']}\n{p['seq']}\n")
            loc_hits = parse(run(["tblastn", "-query", pf, "-subject", sf, "-seg", "no", "-comp_based_stats", "0",
                                  "-word_size", "2", "-evalue", "1000", "-outfmt", OUTFMT]),
                             offset=(L["contig"], w0))
            # local anchor: the window HSP at the genome-wide anchor's position (else the best)
            anc = L["anchor"]
            same = [h for h in loc_hits if h["strand"] == L["strand"] and h["hi"] >= anc["lo"] and h["lo"] <= anc["hi"]]
            lanc = max(same or [h for h in loc_hits if h["strand"] == L["strand"]] or [None],
                       key=lambda h: h["bits"] if h else 0)
            ch = chain(loc_hits, L["strand"], a.max_intron_bp, anchor=lanc, max_exons=a.max_exons,
                       min_bits=a.min_exon_bits) if lanc else []
            rows, ident, mapped = residue_map(ch, p["seq"]) if ch else ([], 0, 0)
            results.append(dict(rank=rank, L=L, chain=ch, rows=rows, ident=ident, mapped=mapped,
                                rbits=sum(h["bits"] for h in ch)))
        results.sort(key=lambda r: (-r["ident"], -r["rbits"]))
        best = results[0]
        for r in results:
            L = r["L"]
            cand_rows.append(dict(peptide=p["id"], rank=r["rank"], contig=L["contig"], strand=L["strand"],
                                  start=L["lo"], end=L["hi"], genome_bits=f"{L['gbits']:.1f}",
                                  refined_bits=f"{r['rbits']:.1f}", identical=r["ident"], mapped=r["mapped"],
                                  length=len(p["seq"]), chosen="yes" if r is best else "no"))
        if not best["chain"]:
            N.warn(f"{p['id']}: no chain survived refinement")
            continue
        ch, L = best["chain"], best["L"]
        lo, hi = min(h["lo"] for h in ch), max(h["hi"] for h in ch)
        order = ch                                         # query order = gene order
        introns = []
        for x, y in zip(order, order[1:]):
            introns.append((y["lo"] - x["hi"] - 1) if L["strand"] == "+" else (x["lo"] - y["hi"] - 1))
        best["rows"], n_ext = extend_ends(best["rows"], files["fasta"], L["contig"], L["strand"], clen[L["contig"]])
        if n_ext:
            best["ident"] = sum(r["status"] == "identical" for r in best["rows"])
            best["mapped"] = sum(r["codon_start"] != "" and r["status"] != "gap" for r in best["rows"])
            lo = min(lo, min(r["codon_start"] for r in best["rows"] if r["codon_start"] != ""))
            hi = max(hi, max(r["codon_end"] for r in best["rows"] if r["codon_end"] != ""))
        fb, fa_ = flanks(files["fasta"], L["contig"], L["strand"], best["rows"], clen[L["contig"]])
        ov = G.overlaps(genes_on(L["contig"]), lo, hi)
        ru = results[1] if len(results) > 1 else None
        map_rows.append(dict(
            peptide=p["id"], group=p["group"], label=p["label"], length=len(p["seq"]), contig=L["contig"],
            strand=L["strand"], start=lo, end=hi, exons=len(ch), introns=",".join(map(str, introns)) or "-",
            identical=best["ident"], mapped=best["mapped"], identity=f"{best['ident'] / len(p['seq']):.3f}",
            bits=f"{best['rbits']:.1f}", before=fb or "-", after=fa_ or "-",
            genes="; ".join(f"{G.gene_label(g)}: {w}" for g, w in ov) or "none annotated",
            runner_up=(f"{ru['L']['contig']}:{ru['L']['lo']}-{ru['L']['hi']}({ru['L']['strand']}) "
                       f"{ru['ident']}/{len(p['seq'])}") if ru else "-",
            source=p["source"], note=""))
        for ei, h in enumerate(ch, 1):
            hsp_rows.append(dict(peptide=p["id"], exon=ei, qstart=h["qs"], qend=h["qe"], contig=h["contig"],
                                 gstart=h["lo"], gend=h["hi"], strand=h["strand"], bits=f"{h['bits']:.1f}",
                                 qseq=h["qseq"], sseq=h["sseq"]))
        for r in best["rows"]:
            res_rows.append(dict(peptide=p["id"], contig=L["contig"], strand=L["strand"], **r))
        N.info(f"{p['id']}: {L['contig']}:{lo}-{hi} ({L['strand']}) {best['ident']}/{len(p['seq'])} identical, "
               f"{len(ch)} exon(s); before {fb or '-'} | after {fa_ or '-'}; "
               + ("; ".join(f"{g['name']} {w}" for g, w in ov) or "no annotated gene")
               + (f"; runner-up {ru['ident']}/{len(p['seq'])}" if ru else ""))

    def write(fn, cols, rows):
        with open(os.path.join(a.outdir, fn), "w") as fh:
            fh.write("\t".join(cols) + "\n")
            for r in rows:
                fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")
    write("candidates.tsv", ["peptide", "rank", "contig", "strand", "start", "end", "genome_bits", "refined_bits",
                             "identical", "mapped", "length", "chosen"], cand_rows)
    write("map.tsv", ["peptide", "group", "label", "length", "contig", "strand", "start", "end", "exons", "introns",
                      "identical", "mapped", "identity", "bits", "before", "after", "genes", "runner_up", "source",
                      "note"], map_rows)
    write("hsps.tsv", ["peptide", "exon", "qstart", "qend", "contig", "gstart", "gend", "strand", "bits", "qseq",
                       "sseq"], hsp_rows)
    write("residues.tsv", ["peptide", "qpos", "qres", "contig", "strand", "exon", "codon_start", "codon_end", "gres",
                           "status", "ins_after", "ext"], res_rows)
    with open(os.path.join(a.outdir, "genome.txt"), "w") as fh:
        fh.write(f"code\t{a.code}\naccession\t{acc}\nfasta\t{files['fasta']}\ngff\t{files['gff'] or ''}\n")
    N.info(f"wrote {a.outdir}/map.tsv, residues.tsv, hsps.tsv, candidates.tsv")


if __name__ == "__main__":
    main()
