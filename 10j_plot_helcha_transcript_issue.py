#!/usr/bin/env python3
"""
HelCha's SLC4A8-AQP5 interval, and the TOGA2 annotation error behind its exendin-3 model.

TOGA2 copies the reference transcript's exon structure onto the query. For the HelSus
exendin-4 gene (LOC_00012325) the EviAnn reference has two isoforms. mRNA-2 starts with a
419-bp first exon whose CDS has no protein support; mRNA-1 has the canonical first exon.
TOGA2 projected only mRNA-2 onto the HelCha exendin-3 locus: the HelCha model took the
wrong first exon (~110 extra residues, 198 aa) and missed the real one, 5.9 kb away. Step 10
excludes mRNA-2 and rebuilds the locus from the published H. horridum exendin-3 cDNA
(AJ580309, lib_rescue_locus.py).

  A  HelCha interval: TOGA2 genes (HelSus reference onto HelCha, step 3's genes table),
     above the curated names (names.tsv), as in 10i. The projection of --bad-tx is
     labelled red; the exendin-3 locus of panel C is boxed. One pair of tracks per HelCha
     contig (the assembly may split the interval).
  B  The cause, in the HelSus reference: every transcript of --ref-gene in the EviAnn GFF,
     exon by exon, 5' -> 3'. Exons of --bad-tx that no other isoform has are outlined red.
     Thin box = UTR, tall box = CDS.
  C  The evidence that exendin-3 is at the curated locus: the HelCha genome itself is read
     at the curated exons (--rescued, from lib_rescue_locus.py), spliced, translated, and
     aligned residue by residue to the published precursor (--published, P20394). Each
     residue sits on the colour of the exon its codon starts in; codons split by an intron
     are marked; each exon is labelled with its genomic coordinates. The mature peptide
     (--mature-id) is bracketed. Differences are red.

The genome is read through its samtools index (<genome>.fai), in Python: no samtools needed.
Standard library only; read-only on its inputs. SVG, and PDF / PNG through rsvg-convert.

Usage
  python3 10j_plot_helcha_transcript_issue.py --species HelCha \
      --genes ../synteny_toga2/interval/interval_HelSus_genes.tsv \
      --toga2 ${TOGADIR}/target_HelSus/query_HelCha --ref-gff ${HELSUS_EVIANN_GFF} \
      --names names/names.tsv --ref-gene LOC_00012325 --bad-tx LOC_00012325-mRNA-2 \
      --rescued rescued/HelCha_Ex3.exons.tsv --rescued-cds rescued/HelCha_Ex3.cds.fa \
      --genome ${HELCHA_FASTA} --genome-contig JANEZZ010000386.1 \
      --published ${SCRIPTDIR}/4_reference_peptides.fa \
      --svg helcha_transcript_issue.svg --pdf helcha_transcript_issue.pdf \
      --png helcha_transcript_issue.png 2> helcha_transcript_issue.log
"""
import argparse, csv, html, importlib.util, math, os, re, sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import lib_exendin_names as N  # noqa: E402
import lib_gff as G  # noqa: E402
import lib_pub_synteny as P  # noqa: E402


def load(fn, name):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fn))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


T3C = load("3c_plot_toga2_synteny.py", "toga2_interval")
F10I = load("10i_plot_helsus_nomenclature.py", "helsus_nomenclature")   # spread, colours, short_name
CLASS_COL, BOUNDARY, OTHER, UNNAMED = F10I.CLASS_COL, F10I.BOUNDARY, F10I.OTHER, F10I.UNNAMED
BAD = "#D7263D"
ANGLE = F10I.ANGLE
esc, overlap, spread = F10I.esc, F10I.overlap, F10I.spread


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t"))


def read_fa(path):
    seqs, k = {}, None
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            k = line[1:].split()[0]; seqs[k] = []
        elif k:
            seqs[k].append(line.upper())
    return {k: "".join(v) for k, v in seqs.items()}


def canon(c):
    return re.sub(r"[._]+", "_", (c or "").lower())


def bp(iv):
    return sum(e - s + 1 for s, e in iv)


# ------------------------------------------------------------------ data
def rescued_model(path):
    """curated exons in transcript order: [(g_start, g_end, strand)], and their span."""
    rows = read_tsv(path)
    if not rows:
        sys.exit(f"[error] {path}: no exons")
    ex = [(min(int(r["g_start"]), int(r["g_end"])), max(int(r["g_start"]), int(r["g_end"])), r["strand"])
          for r in sorted(rows, key=lambda r: int(r["exon"]))]
    return ex, min(e[0] for e in ex), max(e[1] for e in ex), ex[0][2]


def read_fai(fasta):
    p = fasta + ".fai"
    if not os.path.exists(p):
        sys.exit(f"[error] no index {p}: samtools faidx {fasta}")
    out = {}
    for line in open(p):
        t = line.rstrip("\n").split("\t")
        out[t[0]] = tuple(int(x) for x in t[1:5])        # length, offset, bases/line, bytes/line
    return out


def fetch(fasta, fai, contig, s, e):
    """bases s..e (1-based, inclusive) of contig, upper case."""
    ln, off, lb, lw = fai[contig]
    if s < 1 or e > ln:
        sys.exit(f"[error] {contig}:{s}-{e} is outside the sequence (length {ln})")
    with open(fasta, "rb") as fh:
        p0 = off + (s - 1) // lb * lw + (s - 1) % lb
        p1 = off + (e - 1) // lb * lw + (e - 1) % lb
        fh.seek(p0)
        raw = fh.read(p1 - p0 + 1).decode()
    return re.sub(r"\s", "", raw).upper()


COMP = str.maketrans("ACGTN", "TGCAN")
_B, _A = "TCAG", "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODE = {x + y + z: _A[16 * i + 4 * j + k] for i, x in enumerate(_B) for j, y in enumerate(_B) for k, z in enumerate(_B)}


def spliced(fasta, fai, contig, exons):
    """spliced sequence 5'->3', the genomic position of every base, and its exon number."""
    seq, pos, exn = [], [], []
    for n, (s, e, st) in enumerate(exons, 1):
        x = fetch(fasta, fai, contig, s, e)
        ps = list(range(s, e + 1))
        if st == "-":
            x, ps = x.translate(COMP)[::-1], ps[::-1]
        seq.append(x); pos += ps; exn += [n] * len(ps)
    return "".join(seq), pos, exn


def align(a, b):
    """global alignment (match 2, mismatch -1, gap -2) -> two gapped strings."""
    n, m = len(a), len(b)
    S = [[0] * (m + 1) for _ in range(n + 1)]
    for i in range(1, n + 1):
        S[i][0] = -2 * i
    for j in range(1, m + 1):
        S[0][j] = -2 * j
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            S[i][j] = max(S[i - 1][j - 1] + (2 if a[i - 1] == b[j - 1] else -1), S[i - 1][j] - 2, S[i][j - 1] - 2)
    i, j, x, y = n, m, [], []
    while i or j:
        if i and j and S[i][j] == S[i - 1][j - 1] + (2 if a[i - 1] == b[j - 1] else -1):
            x.append(a[i - 1]); y.append(b[j - 1]); i -= 1; j -= 1
        elif i and S[i][j] == S[i - 1][j] - 2:
            x.append(a[i - 1]); y.append("-"); i -= 1
        else:
            x.append("-"); y.append(b[j - 1]); j -= 1
    return "".join(x[::-1]), "".join(y[::-1])


def ref_transcripts(gff, gene):
    """transcripts of one reference gene: [(id, exons, cds)], contig, strand."""
    hit = G.find_genes(gff, [f"g={gene}"]).get("g")
    if not hit:
        sys.exit(f"[error] --ref-gene {gene} not found in {gff}")
    contig, s, e, strand, gid = hit
    genes = [g for g in G.load_region(gff, contig, s, e) if g["id"].split(":")[-1] == gene or g["id"] == gid]
    if not genes:
        sys.exit(f"[error] {gene}: no transcripts loaded")
    txs = sorted(genes[0]["tx"], key=lambda t: t["id"])
    return [(t["id"], t["exons"], t["cds"]) for t in txs], contig, strand, (s, e)


# ------------------------------------------------------------------ drawing
def transcript(body, y, exons, cds, X, h, col, red=(), dashed=False, title=""):
    """one transcript: intron line, UTR boxes (h/2), CDS boxes (h); exons in red outlined."""
    xs = [X(v) for s, e in exons for v in (s, e + 1)]
    dash = ' stroke-dasharray="1.5,1"' if dashed else ""
    body.append(f'<line x1="{min(xs):.2f}" y1="{y:.2f}" x2="{max(xs):.2f}" y2="{y:.2f}" stroke="#555" '
                f'stroke-width="0.5"{dash}/>')
    for s, e in exons:
        a, b = sorted((X(s), X(e + 1)))
        b = max(b, a + 0.6)
        isred = (s, e) in red
        body.append(f'<rect x="{a:.2f}" y="{y - h / 4:.2f}" width="{b - a:.2f}" height="{h / 2:.2f}" '
                    f'fill="{col}" fill-opacity="0.45" stroke="{BAD if isred else col}" '
                    f'stroke-width="{1.0 if isred else 0.4}"><title>{esc(title)} exon {s}-{e}</title></rect>')
    for s, e in cds:
        a, b = sorted((X(s), X(e + 1)))
        b = max(b, a + 0.6)
        isred = any(overlap(s, e, r0, r1) for r0, r1 in red)
        body.append(f'<rect x="{a:.2f}" y="{y - h / 2:.2f}" width="{b - a:.2f}" height="{h:.2f}" '
                    f'fill="{col}" stroke="{BAD if isred else col}" stroke-width="{1.0 if isred else 0.4}">'
                    f'<title>{esc(title)} CDS {s}-{e}</title></rect>')


def gene_track(body, feats, y, X, lo, hi, fs, gene_h, x0, x1, fill_of, flip, bold_kinds=("named",)):
    """10i-style gene track with rotated, spread labels; returns the track's bottom y."""
    mm = P.PT_PER_MM
    sin, cos = math.sin(math.radians(ANGLE)), math.cos(math.radians(ANGLE))
    labelled = sorted([f for f in feats if f["label"]], key=lambda f: X((max(f["start"], lo) + min(f["end"], hi)) / 2))
    tw = max([P.text_w(f["label"], fs) for f in labelled] + [0])
    ygene = y + tw * sin + fs * cos + 3 * mm + gene_h / 2
    body.append(f'<line x1="{x0:.2f}" y1="{ygene:.2f}" x2="{x1:.2f}" y2="{ygene:.2f}" stroke="#000" stroke-width="0.5"/>')
    for f in sorted(feats, key=lambda f: f["kind"] not in ("other",)):
        xa, xb = sorted((X(max(f["start"], lo)), X(min(f["end"], hi) + 1)))
        st = P.flipstrand(f["strand"]) if flip else f["strand"]
        pts, _ = P.arrow(xa, xb, ygene, gene_h, st, 0.8)
        col = fill_of(f)
        dash = ' stroke-dasharray="1.2,0.8"' if f.get("dashed") else ""
        body.append(f'<polygon points="{pts}" fill="{col if f["filled"] else "#FFFFFF"}" '
                    f'stroke="{col if col != OTHER else "#9A9A9A"}" stroke-width="{0.5 if f["filled"] else 0.8}"{dash}>'
                    f'<title>{esc(f["id"])} | {esc(f.get("desc", ""))} | {f["start"]}-{f["end"]} {f["strand"]}</title></polygon>')
    anchors = [X((max(f["start"], lo) + min(f["end"], hi)) / 2) for f in labelled]
    ybase = ygene - gene_h / 2 - 3 * mm
    for f, ax, px in zip(labelled, anchors, spread(anchors, fs * 1.75, x0, x1 - fs * cos)):
        body.append(f'<polyline points="{ax:.2f},{ygene - gene_h / 2 - 0.3:.2f} {ax:.2f},{ygene - gene_h / 2 - 1.2 * mm:.2f} '
                    f'{px:.2f},{ybase:.2f}" fill="none" stroke="#7F7F7F" stroke-width="0.3"/>')
        w = (' font-weight="bold"' if f["kind"] in bold_kinds else "") + (f' fill="{BAD}"' if f.get("bad") else "")
        body.append(f'<text transform="translate({px + 0.6:.2f},{ybase - 0.6:.2f}) rotate(-{ANGLE})" '
                    f'font-size="{fs:.1f}"{w}>{esc(f["label"])}</text>')
    return ygene, ygene + gene_h / 2


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--species", default="HelCha")
    ap.add_argument("--genes", required=True, help="3c --genes-out of the HelSus run (interval_HelSus_genes.tsv)")
    ap.add_argument("--toga2", required=True, help="target_HelSus/query_<species>")
    ap.add_argument("--ref-gff", required=True, help="the EviAnn GFF TOGA2 projected (HelSus)")
    ap.add_argument("--names", required=True, help="10f names/names.tsv")
    ap.add_argument("--ref-gene", default="LOC_00012325")
    ap.add_argument("--bad-tx", default="LOC_00012325-mRNA-2")
    ap.add_argument("--rescued", required=True, help="lib_rescue_locus.py <name>.exons.tsv")
    ap.add_argument("--rescued-cds", help="<name>.cds.fa: picks the CDS in the spliced genome "
                                          "(default: the longest ORF, as lib_rescue_locus.py)")
    ap.add_argument("--genome", required=True, help="the HelCha genome FASTA (indexed: <genome>.fai)")
    ap.add_argument("--genome-contig", help="the contig's name in the genome FASTA "
                                            "(default: TOGA2's name with __1 read as .1)")
    ap.add_argument("--published", required=True, help="FASTA with the published precursor (4_reference_peptides.fa)")
    ap.add_argument("--published-id", default="EXENDIN3_precursor", help="its record (first field of the header)")
    ap.add_argument("--mature-id", default="EXENDIN3_mature", help="mature peptide record, bracketed in C ('' for none)")
    ap.add_argument("--flank-kb", type=float, default=1.0, help="panel B: shown beyond the transcripts")
    ap.add_argument("--pad-kb", type=float, default=5.0, help="panel A: shown beyond the genes on each contig")
    ap.add_argument("--width-mm", type=float, default=180.0)
    ap.add_argument("--font-pt", type=float, default=7.0)
    ap.add_argument("--svg"); ap.add_argument("--pdf"); ap.add_argument("--png")
    ap.add_argument("--dpi", type=int, default=600)
    a = ap.parse_args()
    sp = a.species
    mm, fs = P.PT_PER_MM, a.font_pt
    W = a.width_mm * mm

    # ---- data: panel C first (it defines the locus boxed in A) ----
    r_ex, r_lo, r_hi, r_st = rescued_model(a.rescued)
    r_ctg = None
    notes = G.ref_gene_notes(a.ref_gff)
    proj = G.load_toga2_query(a.toga2, notes, merge_bp=0)
    # contig of the rebuilt model: the TOGA2 contig with projections over it, else from --genes
    grows = [r for r in read_tsv(a.genes) if r["species"] == sp]
    for c, gs in proj.items():
        if any(overlap(r_lo, r_hi, g["start"], g["end"]) for g in gs) and any(canon(c) == canon(r["contig"]) for r in grows):
            r_ctg = c
    if r_ctg is None:
        cands = {r["contig"] for r in grows if overlap(r_lo, r_hi, int(r["start"]), int(r["end"]))}
        r_ctg = next(iter(cands), None)
    if r_ctg is None:
        sys.exit(f"[error] no {sp} TOGA2 gene overlaps the rebuilt model ({r_lo}-{r_hi}); check --rescued and --genes")
    N.info(f"C: curated model on {r_ctg}:{r_lo}-{r_hi} {r_st}, {len(r_ex)} exon(s)")

    loc_tx = []                                   # TOGA2 projections over the locus
    for g in proj.get(r_ctg, []):
        for t in g["tx"]:
            if overlap(r_lo, r_hi, t["start"], t["end"]):
                loc_tx.append(t)
    loc_tx.sort(key=lambda t: (not t["id"].startswith(a.bad_tx), t["id"]))
    other_iso = sorted({t["id"] for gs in proj.values() for g in gs for t in g["tx"]
                        if t["id"].startswith(a.ref_gene + "-") and not t["id"].startswith(a.bad_tx)})
    iso_here = [t["id"] for t in loc_tx if t["id"].startswith(a.ref_gene + "-") and not t["id"].startswith(a.bad_tx)]
    for t in loc_tx:
        aa = bp(t["cds"]) // 3 - 1 if t["cds"] else 0
        miss = [e for e in t["exons"] if not any(overlap(*e, s_, e_) for s_, e_, _ in r_ex)]
        N.info(f"A: TOGA2 {t['id']} {t['status']} {t['start']}-{t['end']}, {len(t['exons'])} exon(s), {aa} aa"
               + (f"; exon(s) not in the curated model: {', '.join(f'{x}-{y} ({y - x + 1} bp)' for x, y in miss)}"
                  if miss else ""))
    if not loc_tx:
        N.warn(f"A: no TOGA2 projection overlaps the curated model on {r_ctg}")
    N.info(f"A: other {a.ref_gene} isoforms projected onto {sp}: "
           + (", ".join(other_iso) if other_iso else "none") + f"; at this locus: {', '.join(iso_here) or 'none'}")

    ref_tx, ref_ctg, ref_st, ref_span = ref_transcripts(a.ref_gff, a.ref_gene)
    if not any(t[0].startswith(a.bad_tx) for t in ref_tx):
        N.warn(f"B: {a.bad_tx} is not a transcript of {a.ref_gene} in {a.ref_gff}")
    good = {e for tid, ex, _ in ref_tx if not tid.startswith(a.bad_tx) for e in ex}
    bad_only = {e for tid, ex, _ in ref_tx if tid.startswith(a.bad_tx) for e in ex
                if not any(overlap(*e, *g) for g in good)}
    for tid, ex, cds in ref_tx:
        N.info(f"B: {tid}: {len(ex)} exon(s), CDS {bp(cds)} bp" + (f" = {bp(cds) // 3 - 1} aa" if cds else "")
               + (f"; exon(s) no other isoform has: {', '.join(f'{s}-{e} ({e - s + 1} bp)' for s, e in sorted(bad_only))}"
                  if tid.startswith(a.bad_tx) else ""))

    # ---- panel A data ----
    tq, how = T3C.read_genes(a.toga2)
    names = [r for r in read_tsv(a.names) if r.get("species") == sp and r.get("class") not in ("GCG", "GIP")
             and (r.get("start") or "").isdigit()]
    by_ctg = defaultdict(list)
    for r in grows:
        qg = r["query_gene"]
        tx = sorted({T3C.parse_projection(p)[0] for p in (tq or {}).get(qg, {}).get("proj", [])})
        txs = tx[0] + (f" +{len(tx) - 1}" if len(tx) > 1 else "") if tx else ""
        isb = r.get("boundary") == "boundary"
        grp = r.get("group", "")
        fam = N.is_family(grp) or grp == "proglucagon"
        by_ctg[r["contig"]].append(dict(
            start=int(r["start"]), end=int(r["end"]), strand=r["strand"], id=qg, group=grp,
            kind="boundary" if isb else ("family" if fam else "other"),
            label=(r["symbol"].split("/")[0] if isb else (f"{qg} · {txs}" if txs else qg) if fam else ""),
            filled=r.get("status", "") in ("FI", "I", ""), dashed=qg.startswith("retro_"),
            status=r.get("status", ""), desc=f"{r['symbol']} | {r.get('status', '')} | {', '.join(tx)}",
            bad=any(t.startswith(a.bad_tx) for t in tx)))
    if not by_ctg:
        sys.exit(f"[error] no {sp} rows in {a.genes}")
    fam_groups = sorted({f["group"] for fs_ in by_ctg.values() for f in fs_ if f["kind"] == "family"})
    fcol = N.colour_map(fam_groups)

    def fill_of(f):
        return {"boundary": BOUNDARY, "other": OTHER, "unnamed": UNNAMED}.get(
            f["kind"], fcol.get(f["group"]) if f["kind"] == "family" else CLASS_COL.get(f["group"], CLASS_COL["copy"]))

    # ---- geometry ----
    x0 = max(P.text_w(s, fs) for s in ("TOGA2", "curated", "HelSus reference")) + 4 * mm
    x1 = W - 4 * mm
    # panel A: room on the right for the rightmost rotated label of each contig, and the 'B, C' tag
    cos_ = math.cos(math.radians(ANGLE))
    rA = 4 * mm
    for c, fs_ in by_ctg.items():
        labs = [f for f in fs_ if f["label"]]
        if labs:
            mid = (min(f["start"] for f in fs_) + max(f["end"] for f in fs_)) / 2
            b = {f["label"]: (f["start"] + f["end"]) / 2 for f in fs_ if f["kind"] == "boundary"}
            fl_ = (b["SLC4A8"] > mid) if "SLC4A8" in b else (b["AQP5"] < mid) if "AQP5" in b else False
            last = min(labs, key=lambda f: f["start"]) if fl_ else max(labs, key=lambda f: f["end"])
            rA = max(rA, P.text_w(last["label"], fs) * cos_ + 2 * mm)
    rA = max(rA, P.text_w("C", fs) + 3 * mm)
    x1A = W - rA
    body, y = [], 0.0

    def heading(letter, text):
        nonlocal y
        y += fs * 1.6
        body.append(f'<text x="0" y="{y:.2f}" font-size="{fs + 1:.1f}" font-weight="bold">{letter}</text>'
                    f'<text x="{fs * 1.4:.2f}" y="{y:.2f}" font-size="{fs:.1f}" font-weight="bold">{esc(text)}</text>')
        y += fs * 0.6

    # panel A: one block per contig, shared scale
    heading("A", f"{sp} SLC4A8–AQP5 interval: TOGA2 (HelSus reference) and curated names")
    wins = {}
    for c, fs_ in by_ctg.items():
        lo = min(f["start"] for f in fs_) - int(a.pad_kb * 1000)
        hi = max(f["end"] for f in fs_) + int(a.pad_kb * 1000)
        b = {f["label"]: (f["start"] + f["end"]) / 2 for f in fs_ if f["kind"] == "boundary"}
        mid = (lo + hi) / 2
        flip = (b["SLC4A8"] > mid) if "SLC4A8" in b else (b["AQP5"] < mid) if "AQP5" in b else False
        wins[c] = (max(1, lo), hi, flip)
    scaleA = (x1A - x0) / max(hi - lo for lo, hi, _ in wins.values())
    order = sorted(wins, key=lambda c: (not any(f["kind"] == "boundary" and f["label"] == "SLC4A8" for f in by_ctg[c]), c))
    for c in order:
        lo, hi, flip = wins[c]
        span = (hi - lo) * scaleA
        labs = {f["label"] for f in by_ctg[c] if f["kind"] == "boundary"}
        # as step 2: a contig with only AQP5 is right-aligned, everything else left-aligned
        xa0, xa1 = (x1A - span, x1A) if labs == {"AQP5"} else (x0, x0 + span)

        def X(p, lo=lo, hi=hi, flip=flip, xa0=xa0, xa1=xa1):
            r = (p - lo) / (hi - lo)
            return xa1 - r * (xa1 - xa0) if flip else xa0 + r * (xa1 - xa0)
        feats = by_ctg[c]
        nm = []
        for r in names:
            if canon(r["contig"]) != canon(c):
                continue
            s, e = int(r["start"]), int(r["end"])
            hit = sorted((f for f in feats if overlap(s, e, f["start"], f["end"])), key=lambda f: -overlap(s, e, f["start"], f["end"]))
            nm.append(dict(start=s, end=e, strand=hit[0]["strand"] if hit else "", kind="named", group=r["class"],
                           label=F10I.short_name(r["new_name"], sp), id=r["new_name"], filled=True, desc=r.get("rule", "")))
        for f in feats:
            if f["kind"] == "family" and not any(overlap(f["start"], f["end"], n["start"], n["end"]) for n in nm):
                nm.append(dict(f, kind="unnamed", label="", filled=False))
        nm += [dict(f) for f in feats if f["kind"] == "boundary"]
        ytop = y
        body.append(f'<text x="{xa0:.2f}" y="{y + fs:.2f}" font-size="{fs:.1f}" fill="#555">{esc(c)}'
                    f'{" (reversed)" if flip else ""}</text>')
        y += fs * 1.2
        yg1, yb1 = gene_track(body, feats, y, X, lo, hi, fs, 2.4 * mm, xa0, xa1, fill_of, flip, bold_kinds=())
        body.append(f'<text x="{x0 - 2 * mm:.2f}" y="{yg1 + fs * 0.35:.2f}" font-size="{fs:.1f}" text-anchor="end">TOGA2</text>')
        yg2, yb2 = gene_track(body, nm, yb1 + fs * 0.4, X, lo, hi, fs, 2.4 * mm, xa0, xa1, fill_of, flip)
        body.append(f'<text x="{x0 - 2 * mm:.2f}" y="{yg2 + fs * 0.35:.2f}" font-size="{fs:.1f}" text-anchor="end">curated</text>')
        for n in nm:                                   # bands
            if n["kind"] == "named":
                p, q = sorted((X(n["start"]), X(n["end"] + 1)))
                w = max(q - p, 0.8)
                body.insert(0, f'<rect x="{(p + q) / 2 - w / 2:.2f}" y="{ytop + fs * 1.4:.2f}" width="{w:.2f}" '
                               f'height="{yb2 - ytop - fs * 1.2:.2f}" fill="{fill_of(n)}" fill-opacity="0.12"/>')
        if canon(c) == canon(r_ctg):                   # box the panel-C locus
            p, q = sorted((X(r_lo), X(r_hi + 1)))
            p, q = p - 1.2 * mm, q + 1.2 * mm
            body.append(f'<rect x="{p:.2f}" y="{yg1 - 2.2 * mm:.2f}" width="{q - p:.2f}" height="{yg2 - yg1 + 4.4 * mm:.2f}" '
                        f'fill="none" stroke="#000" stroke-width="0.6" stroke-dasharray="2,1"/>'
                        f'<text x="{q + 0.8:.2f}" y="{yg2 + 2.2 * mm:.2f}" font-size="{fs:.1f}" font-weight="bold">C</text>')
        y = yb2 + fs * 0.8
    sbA = P.nice_scale((x1A - x0) / scaleA)
    body.append(f'<line x1="{x0:.2f}" y1="{y:.2f}" x2="{x0 + sbA * scaleA:.2f}" y2="{y:.2f}" stroke="#000" stroke-width="0.7"/>'
                f'<text x="{x0 + sbA * scaleA + 1.5 * mm:.2f}" y="{y + fs * 0.35:.2f}" font-size="{fs:.1f}">{P.fmt_kb(sbA)}</text>')
    y += fs * 0.8

    # panel B: exon level, 5' -> 3'
    fl = int(a.flank_kb * 1000)
    b_lo = min(s for _, ex, _ in ref_tx for s, _ in ex) - fl
    b_hi = max(e for _, ex, _ in ref_tx for _, e in ex) + fl
    lab_w = max(P.text_w(f"{t} · 999 aa", fs) for t, _, _ in ref_tx) + 2 * mm
    xs0, xs1 = x0, x1 - lab_w
    scaleB = (xs1 - xs0) / (b_hi - b_lo)
    h = 2.6 * mm
    rowh = fs * 2.3

    def XB(p):
        return xs0 + ((b_hi - p) if ref_st == "-" else (p - b_lo)) * scaleB

    heading("B", f"Cause: HelSus reference {a.ref_gene} (EviAnn), every isoform; red = exon only {a.bad_tx} has")
    y += fs * 0.6
    for tid, ex, cds in ref_tx:
        y += rowh / 2
        isbad = tid.startswith(a.bad_tx)
        transcript(body, y, ex, cds, XB, h, BAD if isbad else "#4D4D4D", red=bad_only if isbad else (), title=tid)
        aa = bp(cds) // 3 - 1 if cds else 0
        xe = max(XB(v) for s_, e_ in ex for v in (s_, e_ + 1))
        body.append(f'<text x="{xe + 1.5 * mm:.2f}" y="{y + fs * 0.35:.2f}" font-size="{fs:.1f}">{esc(tid)}'
                    f'{f" · {aa} aa" if aa else ""}</text>')
        for s_, e_ in (bad_only if isbad else ()):
            body.append(f'<text x="{(XB(s_) + XB(e_ + 1)) / 2:.2f}" y="{y - h / 2 - 1:.2f}" font-size="{fs - 1:.1f}" '
                        f'fill="{BAD}" text-anchor="middle">{e_ - s_ + 1} bp</text>')
        y += rowh / 2
    body.append(f'<text x="{x0 - 2 * mm:.2f}" y="{y - rowh * len(ref_tx) / 2 + fs * 0.35:.2f}" font-size="{fs:.1f}" '
                f'text-anchor="end">HelSus reference</text>')
    y += fs * 0.2
    sb = P.nice_scale((xs1 - xs0) / scaleB)
    body.append(f'<line x1="{xs0:.2f}" y1="{y:.2f}" x2="{xs0 + sb * scaleB:.2f}" y2="{y:.2f}" stroke="#000" stroke-width="0.7"/>'
                f'<text x="{xs0 + sb * scaleB + 1.5 * mm:.2f}" y="{y + fs * 0.35:.2f}" font-size="{fs:.1f}">'
                f'{P.fmt_kb(sb)} · 5′ → 3′ · thin = UTR, thick = CDS</text>')
    y += fs * 0.8

    # panel C: the genome, read at the curated exons, translated and aligned to the published precursor
    fai = read_fai(a.genome)
    gctg = a.genome_contig or re.sub(r"__(\d+)$", r".\1", r_ctg)
    if gctg not in fai:
        alt = [k for k in fai if canon(k) == canon(r_ctg)]
        if not alt:
            sys.exit(f"[error] {gctg} is not in {a.genome}.fai; give its name with --genome-contig")
        gctg = alt[0]
    seq, gpos, gexn = spliced(a.genome, fai, gctg, r_ex)
    if a.rescued_cds:
        cds = next(iter(read_fa(a.rescued_cds).values()), "")
        i0 = seq.find(cds) if cds else -1
        if i0 < 0:
            sys.exit(f"[error] the CDS of {a.rescued_cds} is not in the genome spliced at {a.rescued}: "
                     f"wrong --genome / --genome-contig, or the exon table and CDS are from different runs")
        how = f"the CDS of {os.path.basename(a.rescued_cds)}"
    else:
        orfs = [(f + m.start(), m.group(0)) for f in range(3)
                for m in re.finditer(r"ATG(?:[ACGT]{3})*?(?:TAA|TAG|TGA)", seq[f:])]
        if not orfs:
            sys.exit(f"[error] the genome spliced at {a.rescued} has no complete ORF")
        i0, cds = max(orfs, key=lambda x: len(x[1]))
        cds = cds[:-3]
        how = "the longest ORF"
    prot = "".join(CODE.get(cds[k:k + 3], "X") for k in range(0, len(cds) - len(cds) % 3, 3)).rstrip("*")
    res_exon = [gexn[i0 + 3 * k] for k in range(len(prot))]
    res_split = [len(set(gexn[i0 + 3 * k:i0 + 3 * k + 3])) > 1 for k in range(len(prot))]
    ex_cds = {}
    for k in range(len(cds)):
        n_ = gexn[i0 + k]
        lo_, hi_ = ex_cds.get(n_, (gpos[i0 + k], gpos[i0 + k]))
        ex_cds[n_] = (min(lo_, gpos[i0 + k]), max(hi_, gpos[i0 + k]))
    ex_res = {}
    for k, n_ in enumerate(res_exon):
        ex_res.setdefault(n_, [k + 1, k + 1])[1] = k + 1
    pubs = {k.split("|")[0]: (k, v) for k, v in read_fa(a.published).items()}
    if a.published_id not in pubs:
        sys.exit(f"[error] no {a.published_id} in {a.published} (records: {', '.join(pubs)})")
    phead, pub = pubs[a.published_id]
    hf = phead.split("|")
    plab = f"{hf[1]} {hf[2].replace('_', ' ')}" if len(hf) > 2 else hf[0]
    mature = pubs[a.mature_id][1] if a.mature_id and a.mature_id in pubs else ""
    A_, B_ = align(pub, prot)
    ident = sum(1 for x, y_ in zip(A_, B_) if x == y_ and x != "-")
    ngap = sum(1 for x, y_ in zip(A_, B_) if "-" in (x, y_))
    N.info(f"C: {gctg} spliced at {len(r_ex)} exon(s) = {len(seq)} bp; {how} at base {i0 + 1}: "
           f"{len(cds)} bp = {len(prot)} aa")
    N.info(f"C: {plab}: {len(pub)} aa; identical at {ident}/{len(pub)} positions ({100 * ident / len(pub):.1f}%), "
           f"{ngap} gap column(s)")
    for n_ in sorted(ex_cds):
        N.info(f"C: exon {n_}: residues {ex_res.get(n_, ['-', '-'])[0]}-{ex_res.get(n_, ['-', '-'])[1]}, "
               f"coding bases {gctg}:{ex_cds[n_][0]}-{ex_cds[n_][1]} ({r_st})")
    msplit = sum(res_split)

    gl = f"{sp} genome, curated exons"
    x0c = max(x0, max(P.text_w(plab, fs), P.text_w(gl, fs)) + 4 * mm)
    ncol = len(A_)
    step = max(min(fs * 0.95, (x1 - x0c) / ncol), fs * 0.62)
    per = max(10, int((x1 - x0c) / step))
    EXCOL = ["#FCE3CF", "#D6E6F5", "#DFF0D8", "#F3E0EE", "#FFF2B3", "#E5E5E5"]
    mcols = set()
    if mature:
        mi = pub.find(mature)
        if mi < 0:
            N.warn(f"C: {a.mature_id} is not inside {a.published_id}: no bracket")
        else:
            k_, pn = 0, 0
            for c_, ch in enumerate(A_):
                if ch != "-":
                    if mi <= pn < mi + len(mature):
                        mcols.add(c_)
                    pn += 1
    heading("C", f"Evidence: the {sp} genome at the curated exons, translated, against the published "
                 f"exendin-3 precursor")
    pn, gn = 0, 0                                 # residue counters (published, genome)
    for l0 in range(0, ncol, per):
        cols = range(l0, min(ncol, l0 + per))
        X = lambda c_: x0c + (c_ - l0 + 0.5) * step
        y += fs * 0.4
        mc = [c_ for c_ in cols if c_ in mcols]
        if mc:
            ya = y + fs * 0.9
            body.append(f'<path d="M{X(mc[0]) - step / 2:.2f},{ya + 1.6:.2f} V{ya:.2f} H{X(mc[-1]) + step / 2:.2f} V{ya + 1.6:.2f}" '
                        f'fill="none" stroke="#000" stroke-width="0.5"/>'
                        f'<text x="{(X(mc[0]) + X(mc[-1])) / 2:.2f}" y="{ya - 1:.2f}" font-size="{fs - 0.5:.1f}" '
                        f'text-anchor="middle" font-style="italic">mature exendin-3</text>')
            y += fs * 1.2
        yn, yp, ym, yg, ye = y + fs * 0.9, y + fs * 2.0, y + fs * 2.9, y + fs * 3.9, y + fs * 4.9
        body.append(f'<text x="{x0c - 2 * mm:.2f}" y="{yp:.2f}" font-size="{fs:.1f}" text-anchor="end">{esc(plab)}</text>'
                    f'<text x="{x0c - 2 * mm:.2f}" y="{yg:.2f}" font-size="{fs:.1f}" text-anchor="end">{esc(gl)}</text>')
        runs = []
        for c_ in cols:
            pc, gc = A_[c_], B_[c_]
            if pc != "-":
                pn += 1
                if pn % 10 == 0:
                    body.append(f'<text x="{X(c_):.2f}" y="{yn:.2f}" font-size="{fs - 1.5:.1f}" fill="#777" '
                                f'text-anchor="middle">{pn}</text>')
            same = pc == gc and pc != "-"
            col = "#000" if same else BAD
            body.append(f'<text x="{X(c_):.2f}" y="{yp:.2f}" font-size="{fs:.1f}" text-anchor="middle" '
                        f'fill="{col}">{esc(pc)}</text>')
            body.append(f'<text x="{X(c_):.2f}" y="{ym:.2f}" font-size="{fs:.1f}" text-anchor="middle" fill="#777">'
                        f'{"|" if same else ("·" if "-" not in (pc, gc) else " ")}</text>')
            if gc != "-":
                en = res_exon[gn]
                body.append(f'<rect x="{X(c_) - step / 2:.2f}" y="{yg - fs * 0.85:.2f}" width="{step:.2f}" '
                            f'height="{fs * 1.15:.2f}" fill="{EXCOL[(en - 1) % len(EXCOL)]}"><title>residue {gn + 1}, '
                            f'exon {en}, {gctg}:{gpos[i0 + 3 * gn]}</title></rect>')
                if res_split[gn]:
                    body.append(f'<line x1="{X(c_) - step / 2:.2f}" y1="{yg + fs * 0.4:.2f}" x2="{X(c_) + step / 2:.2f}" '
                                f'y2="{yg + fs * 0.4:.2f}" stroke="#000" stroke-width="0.9"/>')
                if runs and runs[-1][0] == en:
                    runs[-1][2] = c_
                else:
                    runs.append([en, c_, c_])
                gn += 1
            body.append(f'<text x="{X(c_):.2f}" y="{yg:.2f}" font-size="{fs:.1f}" text-anchor="middle" '
                        f'font-weight="bold" fill="{col}">{esc(gc)}</text>')
        for en, c0, c1 in runs:
            xa_, xb_ = X(c0) - step / 2 + 0.4, X(c1) + step / 2 - 0.4
            body.append(f'<path d="M{xa_:.2f},{ye - fs * 0.9:.2f} V{ye - fs * 0.6:.2f} H{xb_:.2f} V{ye - fs * 0.9:.2f}" '
                        f'fill="none" stroke="#555" stroke-width="0.4"/>'
                        f'<text x="{(xa_ + xb_) / 2:.2f}" y="{ye + fs * 0.2:.2f}" font-size="{fs - 0.5:.1f}" '
                        f'text-anchor="middle">exon {en}</text>')
        y = ye + fs * 0.6
    y += fs * 0.6
    lines = [f"{gctg}, curated exons spliced 5′→3′ and translated ({how}): {len(prot)} aa; identical to "
             f"{plab} at {ident} of {len(pub)} positions ({100 * ident / len(pub):.0f}%)"
             + (f", {ngap} gap(s)" if ngap else "") + "."]
    lines += [f"exon {n_}: residues {ex_res[n_][0]}–{ex_res[n_][1]}, coding bases {gctg}:{ex_cds[n_][0]:,}–"
              f"{ex_cds[n_][1]:,} ({'minus' if r_st == '-' else 'plus'} strand)" for n_ in sorted(ex_cds) if n_ in ex_res]
    lines.append("Background = exon of the codon's first base; underlined = codon split by an intron"
                 + (f" ({msplit})" if msplit else "") + "; red = differs from the published precursor.")
    for t in lines:                                # word-wrapped to the figure width
        words, cur = t.split(" "), ""
        for w_ in words + [None]:
            if w_ is not None and P.text_w((cur + " " + w_).strip(), fs - 0.5) <= x1 - x0c:
                cur = (cur + " " + w_).strip()
                continue
            y += fs * 1.3
            body.append(f'<text x="{x0c:.2f}" y="{y:.2f}" font-size="{fs - 0.5:.1f}">{esc(cur)}</text>')
            cur = w_ or ""
    H = y + fs * 1.2
    svg = (f'<svg xmlns="http://www.w3.org/2000/svg" width="{W:.2f}pt" height="{H:.2f}pt" viewBox="0 0 {W:.2f} {H:.2f}" '
           f'font-family="{P.FONT}">\n<rect width="100%" height="100%" fill="#FFFFFF"/>\n' + "\n".join(body) + "\n</svg>\n")
    P.write(svg, a)
    if a.svg:
        N.info(f"wrote {a.svg}")
    N.info(f"figure {W / mm:.0f} x {H / mm:.0f} mm")


if __name__ == "__main__":
    main()
