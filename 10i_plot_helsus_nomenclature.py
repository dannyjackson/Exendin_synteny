#!/usr/bin/env python3
"""
One figure of the HelSus SLC4A8-AQP5 interval, drawn three times in three nomenclatures,
on one shared coordinate axis:

  A  the genome's own annotation, original names: the EviAnn annotation TOGA2 projected
     from (one track per --gff; give more than one to compare annotations). Exendin-family
     genes are labelled with the annotation's own ID and description.
  B  TOGA2 with HelSus as both reference and query (target_HelSus/query_HelSus, the
     self-projection): every TOGA2 query gene, labelled in TOGA2's nomenclature, i.e. its
     query-gene ID and the reference transcript(s) it was projected from.
  C  the step-10 names (names/names.tsv: GCGL, Exendin_anc, A, B ... Z, Y ...), in
     10f's class colours. TOGA2 exendin-family genes that have no name (not in the tree)
     are drawn hollow, as in 10f panel A.

A faint vertical band runs from each C locus up through A and B, so the three names of
one locus line up. --crosswalk writes the same correspondence as a table: for each named
locus, the overlapping gene IDs of every annotation.

Coordinates are TOGA2's (= EviAnn's: chr2). An extra GFF whose sequences are named
differently (e.g. NCBI: NC_139459.1) is matched through the NCBI sequence report next to it
(sequence_report.jsonl, as 10g), or --chrom-map chr2=NC_139459.1.

The view runs from SLC4A8 to AQP5 (SLC4A8 on the left; the axis is flipped if SLC4A8 lies
to the right in the assembly), with the boundary genes clipped to --clip-kb on their
outer side and --pad-kb added at both ends.

Colours: exendin-family groups (panels A, B) use the project's fixed colours
(lib_exendin_names); SLC4A8/AQP5 dark grey; other genes light grey, unlabelled unless
--label-all. Hollow = pseudogene (A) / TOGA2 not intact, status other than FI or I (B).
Dashed outline = TOGA2 paralogous or retrogene projection only (B).

Style as the other figures: --width-mm (180), every text --font-pt (7), Arial/Helvetica,
SVG; PDF and PNG through rsvg-convert. Standard library only; read-only on its inputs.

Usage
  python3 10i_plot_helsus_nomenclature.py \
      --gff EviAnn=${HELSUS_EVIANN_GFF} \
      --toga2 ${TOGADIR}/target_HelSus/query_HelSus --names names/names.tsv \
      --boundaries SLC4A8,AQP5=LOC_00012326 \
      --svg helsus_nomenclature.svg --pdf helsus_nomenclature.pdf --png helsus_nomenclature.png \
      --crosswalk names/HelSus_nomenclature_crosswalk.tsv 2> helsus_nomenclature.log
"""
import argparse, csv, html, importlib.util, math, os, re, sys
from collections import Counter, defaultdict

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


T3C = load("3c_plot_toga2_synteny.py", "toga2_interval")     # read_genes, read_status, names, describe
T10G = load("10g_gene_table.py", "gene_table")                # seq_report_map

CLASS_COL = {"GCGL": "#009E73", "Exendin_anc": "#0072B2", "copy": "#D55E00"}   # as 10f
CLASS_LEGEND = {"GCGL": "GCGL", "Exendin_anc": "Exendin_anc", "copy": "lettered copy (A, B … / Z, Y …)"}
BOUNDARY = "#4D4D4D"
OTHER = "#D3D3D3"
UNNAMED = "#9A9A9A"
ANGLE = 40                       # label rotation, degrees above horizontal


def esc(x):
    return html.escape(str(x), quote=True)


def overlap(a0, a1, b0, b1):
    return max(0, min(a1, b1) - max(a0, b0) + 1)


def short_name(new, sp):
    n = re.sub(rf"^{re.escape(sp)}_", "", new)
    return re.sub(r"^Exendin_", "", n)


# ------------------------------------------------------------------ inputs
def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t"))


def gff_track(label, path, contig, lo, hi, manual, bnd):
    sr = T10G.seq_report_map(path)
    c2 = re.sub(r"__(\d+)$", r".\1", contig)
    alias = manual.get(contig) or (contig if contig in sr.values() else None) or sr.get(contig) or sr.get(c2) or c2
    genes = G.load_region(path, alias, lo, hi)
    if not genes and alias == contig:
        # sequence names may still differ: try the report's chromosome number
        alt = sr.get(re.sub(r"^chr", "", contig))
        if alt:
            alias, genes = alt, G.load_region(path, alt, lo, hi)
    if not genes:
        N.warn(f"A/{label}: no genes on {alias}:{lo}-{hi} in {path}; "
               f"give the sequence name with --chrom-map {contig}=<name in that GFF>")
    else:
        N.info(f"A/{label}: {len(genes)} gene(s) on {alias}:{lo}-{hi} ({path})")
    feats = []
    for g in genes:
        text = f"{g['name']} {g['description']} {g['tx'][0]['product']}"
        grp = N.name_group(text)
        fam = N.is_family(grp) or grp == "proglucagon"
        kind = boundary_of(g, bnd) or ("family" if fam else "other")
        desc = re.sub(r"\s+", " ", g["description"] or g["tx"][0]["product"] or "").strip()
        short = re.sub(r"^similar to\s+\S+\s*:\s*", "", desc, flags=re.I)      # EviAnn: 'Similar to EXE4_HELSU: Exendin-4'
        short = short if len(short) <= 26 else short[:25] + "…"
        lab = g["name"] + (f" · {short}" if short and short.lower() != g["name"].lower() and kind not in bnd else "")
        feats.append(dict(start=g["start"], end=g["end"], strand=g["strand"], kind=kind, group=grp,
                          label=lab, id=g["id"] or g["name"], filled="pseudo" not in g["biotype"].lower(),
                          dashed=False, status=g["biotype"], desc=desc))
    return feats, alias


def toga2_track(qdir, ref_gffs, contig, lo, hi, bnd, label_mode):
    genes, how = T3C.read_genes(qdir)
    if genes is None:
        sys.exit(f"[error] {qdir}: {how}")
    status = T3C.read_status(qdir)
    names = T3C.read_ref_names(ref_gffs, None)
    unresolved = Counter()
    feats = []
    for qg, g in genes.items():
        if g["contig"] != contig or g["end"] < lo or g["start"] > hi:
            continue
        d = T3C.describe(qg, g, names, status, unresolved)
        tx = sorted({T3C.parse_projection(p)[0] for p in g["proj"]})
        txs = tx[0] + (f" +{len(tx) - 1}" if len(tx) > 1 else "") if tx else ""
        kind = boundary_of(d, bnd) or ("family" if d["family"] else "other")
        lab = {"query_gene": qg, "projection": txs or qg, "both": f"{qg} · {txs}" if txs else qg}[label_mode]
        if kind in bnd:
            lab = qg                                   # boundaries: the TOGA2 ID only
        feats.append(dict(start=g["start"], end=g["end"], strand=g["strand"], kind=kind, group=d["group"],
                          label=lab, id=qg, filled=d["filled"], dashed=d["dashed"], status=d["status"],
                          desc=f"{d['symbol']} | {d['status']} | {', '.join(tx)}"))
    if unresolved:
        N.warn(f"B: {sum(unresolved.values())} projection(s) whose reference gene has no name "
               f"({', '.join(list(unresolved)[:4])} ...): check --ref-gff")
    N.info(f"B: {len(feats)} TOGA2 query gene(s) in the view ({how}), "
           f"{sum(1 for f in feats if f['kind'] == 'family')} exendin-family")
    return feats


def boundary_of(f, bnd):
    """'SLC4A8' / 'AQP5' if f covers most of that boundary gene (same locus), else ''."""
    for n, b in bnd.items():
        ov = overlap(f["start"], f["end"], b[1], b[2])
        if ov >= 0.5 * min(f["end"] - f["start"] + 1, b[2] - b[1] + 1):
            return n
    return ""


def names_track(path, sp, contig, lo, hi, bfeats, toga):
    rows = [r for r in read_tsv(path) if r.get("species") == sp and r.get("class") not in ("GCG", "GIP")]
    feats, skipped = [], []
    for r in rows:
        if not (r.get("start") or "").isdigit():
            skipped.append(f"{r['new_name']} (no coordinates)"); continue
        if re.sub(r"[._]+", "_", r["contig"].lower()) != re.sub(r"[._]+", "_", contig.lower()):
            skipped.append(f"{r['new_name']} ({r['contig']})"); continue
        s, e = int(r["start"]), int(r["end"])
        if e < lo or s > hi:
            skipped.append(f"{r['new_name']} (outside the view)"); continue
        # strand from the TOGA2 gene at the locus (names.tsv has none)
        hits = sorted((f for f in toga if overlap(s, e, f["start"], f["end"])),
                      key=lambda f: -overlap(s, e, f["start"], f["end"]))
        feats.append(dict(start=s, end=e, strand=hits[0]["strand"] if hits else "", kind="named",
                          group=r["class"], label=short_name(r["new_name"], sp), id=r["new_name"],
                          filled=True, dashed=False, status=r.get("rule", ""), desc=r.get("rule", "")))
    for f in toga:                                 # TOGA2 family loci left out of the tree
        if f["kind"] == "family" and not any(overlap(f["start"], f["end"], n["start"], n["end"]) for n in feats):
            feats.append(dict(f, kind="unnamed", label="", filled=False, group="unnamed"))
    feats += bfeats
    if skipped:
        N.info(f"C: not drawn: {', '.join(skipped)}")
    nn = sum(1 for f in feats if f["kind"] == "named")
    N.info(f"C: {nn} named {sp} locus/loci in the view; "
           f"{sum(1 for f in feats if f['kind'] == 'unnamed')} TOGA2 exendin-family gene(s) without a name (hollow)")
    if not nn:
        N.warn(f"C: no {sp} rows of {path} fall in the view")
    return feats


# ------------------------------------------------------------------ layout helpers
def spread(anchors, d, x0, x1):
    """label x positions for sorted anchors: >= d apart, each crowded cluster centred on its anchors."""
    def place(tot, n):
        return min(max(tot / n - (n - 1) * d / 2, x0), x1 - (n - 1) * d)
    cl = []                                          # [sum of anchors, n]
    for x in anchors:
        cl.append([x, 1])
        while len(cl) > 1 and place(*cl[-2]) + cl[-2][1] * d > place(*cl[-1]):
            t, n = cl.pop(); cl[-1][0] += t; cl[-1][1] += n
    out = []
    for t, n in cl:
        f0 = place(t, n)
        out += [f0 + i * d for i in range(n)]
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--gff", action="append", default=[], metavar="LABEL=PATH",
                    help="panel A track; the first one also gives SLC4A8/AQP5 (use EviAnn, TOGA2's reference)")
    ap.add_argument("--toga2", required=True, help="TOGA2 HelSus self-projection: target_HelSus/query_HelSus")
    ap.add_argument("--ref-gff", action="append", default=[],
                    help="reference annotation(s) naming TOGA2's reference genes (default: the first --gff)")
    ap.add_argument("--names", required=True, help="10f names/names.tsv")
    ap.add_argument("--species", default="HelSus")
    ap.add_argument("--boundaries", default="SLC4A8,AQP5=LOC_00012326")
    ap.add_argument("--chrom-map", action="append", default=[], metavar="CHR=NAME",
                    help="sequence name in the other GFFs, e.g. chr2=NC_139459.1")
    ap.add_argument("--clip-kb", type=float, default=15.0, help="boundary genes are drawn up to this far out")
    ap.add_argument("--pad-kb", type=float, default=0.0, help="extend the view this far beyond the boundaries")
    ap.add_argument("--toga2-label", choices=("both", "query_gene", "projection"), default="both",
                    help="panel B labels: TOGA2 query gene, reference transcript, or both (default)")
    ap.add_argument("--label-all", action="store_true", help="label every gene in A and B, not only the family")
    ap.add_argument("--no-bands", action="store_true", help="no vertical bands linking the panels")
    ap.add_argument("--width-mm", type=float, default=180.0)
    ap.add_argument("--font-pt", type=float, default=7.0)
    ap.add_argument("--svg"); ap.add_argument("--pdf"); ap.add_argument("--png")
    ap.add_argument("--dpi", type=int, default=600)
    ap.add_argument("--crosswalk", help="TSV: every named locus with the IDs of every annotation at it")
    a = ap.parse_args()
    if not a.gff:
        sys.exit("[error] give at least one --gff (EviAnn first)")
    gffs = []
    for s in a.gff:
        lab, _, p = s.partition("=")
        if not p or not os.path.exists(p):
            sys.exit(f"[error] --gff {s}: expected LABEL=existing file")
        gffs.append((lab, p))
    manual = dict(m.split("=", 1) for m in a.chrom_map)

    # ---- view: SLC4A8 .. AQP5 on the first GFF ----
    bspec = [x.strip() for x in a.boundaries.split(",") if x.strip()]
    if len(bspec) != 2:
        sys.exit("[error] --boundaries needs two genes, e.g. SLC4A8,AQP5=LOC_00012326")
    found = G.find_genes(gffs[0][1], bspec)
    bn = [x.partition("=")[0] for x in bspec]
    miss = [n for n in bn if n not in found]
    if miss:
        sys.exit(f"[error] {', '.join(miss)} not found in {gffs[0][1]}; give its ID as NAME=ID")
    (cl, l0, l1, lst, lid), (cr, r0, r1, rst, rid) = found[bn[0]], found[bn[1]]
    if cl != cr:
        sys.exit(f"[error] {bn[0]} ({cl}) and {bn[1]} ({cr}) are on different sequences")
    contig, flip = cl, l0 > r0
    clip, pad = a.clip_kb * 1000, a.pad_kb * 1000
    if not flip:
        lo, hi = max(l0, l1 - clip) - pad, min(r1, r0 + clip) + pad
    else:
        lo, hi = max(r0, r1 - clip) - pad, min(l1, l0 + clip) + pad
    lo = max(1, int(lo)); hi = int(hi)
    bnd = {bn[0]: (contig, l0, l1, lst, lid), bn[1]: (contig, r0, r1, rst, rid)}
    N.info(f"view {contig}:{lo}-{hi} ({(hi - lo + 1) / 1000:.0f} kb), {bn[0]} {lid} {l0}-{l1} {lst}, "
           f"{bn[1]} {rid} {r0}-{r1} {rst}; " + ("axis flipped so that " + bn[0] + " is on the left" if flip
                                                   else bn[0] + " on the left"))

    # ---- the three panels ----
    tracksA, aliases = [], {}
    for lab, p in gffs:
        f, al = gff_track(lab, p, contig, lo, hi, manual, bnd)
        tracksA.append((lab, f)); aliases[lab] = al
    toga = toga2_track(a.toga2, a.ref_gff or [gffs[0][1]], contig, lo, hi, bnd, a.toga2_label)
    bfeats = [dict(start=b[1], end=b[2], strand=b[3], kind=n, group=n, label=n, id=b[4], filled=True,
                   dashed=False, status="", desc="") for n, b in bnd.items()]
    named = names_track(a.names, a.species, contig, lo, hi, bfeats, toga)
    panels = [("A", "Genome annotation, original names" if len(tracksA) == 1 else "Genome annotations, original names",
               tracksA),
              ("B", f"TOGA2 projection, {a.species} reference onto {a.species} (target_{a.species}/query_{a.species})",
               [("TOGA2", toga)]),
              ("C", "Names from the gene tree (step 10)", [("step 10", named)])]

    # ---- colours ----
    fam_groups = sorted({f["group"] for _, tr in tracksA + [("", toga)] for f in tr if f["kind"] == "family"})
    fcol = N.colour_map(fam_groups)

    def fill_of(f):
        if f["kind"] in bn:
            return BOUNDARY
        if f["kind"] == "family":
            return fcol.get(f["group"], UNNAMED)
        if f["kind"] == "named":
            return CLASS_COL.get(f["group"], CLASS_COL["copy"])
        if f["kind"] == "unnamed":
            return UNNAMED
        return OTHER

    # ---- geometry (pt) ----
    mm = P.PT_PER_MM
    fs = a.font_pt
    W = a.width_mm * mm
    left = max(P.text_w(lab, fs) for _, _, trs in panels for lab, _ in trs) + 4 * mm
    cos_ = math.cos(math.radians(ANGLE))

    def is_lab(f):
        return f["label"] and (f["kind"] in ("family", "named") or f["kind"] in bn or a.label_all)
    edge = []
    for _, _, trs in panels:
        for _, feats in trs:
            labs = [f for f in feats if is_lab(f)]
            if labs:
                last = max(labs, key=lambda f: -(f["start"] + f["end"]) if flip else (f["start"] + f["end"]))
                edge.append(P.text_w(last["label"], fs) * cos_)
    right = max([6 * mm] + [e + 2 * mm for e in edge])
    x0, x1 = left, W - right
    span = hi - lo + 1

    def X(pos):
        r = (pos - lo) / span
        return x1 - r * (x1 - x0) if flip else x0 + r * (x1 - x0)

    gene_h = 2.4 * mm
    dlab = fs * 1.75                                  # min horizontal gap between rotated labels
    sin, cos = math.sin(math.radians(ANGLE)), math.cos(math.radians(ANGLE))
    out, y = [], 0.0
    band_spans = [(min(X(f["start"]), X(f["end"])), max(X(f["start"]), X(f["end"])), fill_of(f))
                  for f in named if f["kind"] == "named"]
    panel_tops, body = [], []
    for letter, title, trs in panels:
        y += fs * 1.6
        panel_tops.append(y - fs * 1.2)
        body.append(f'<text x="0" y="{y:.2f}" font-size="{fs + 1:.1f}" font-weight="bold">{letter}</text>'
                    f'<text x="{fs * 1.4:.2f}" y="{y:.2f}" font-size="{fs:.1f}" font-weight="bold">{esc(title)}</text>')
        y += fs * 0.6
        for lab, feats in trs:
            labelled = [f for f in feats if is_lab(f)]
            labelled.sort(key=lambda f: X((max(f["start"], lo) + min(f["end"], hi)) / 2))
            tw = max([P.text_w(f["label"], fs) for f in labelled] + [0])
            lab_h = tw * sin + fs * cos + 3 * mm      # rotated text height + leader rise
            ygene = y + lab_h + gene_h / 2
            # track label and axis line
            body.append(f'<text x="{x0 - 2 * mm:.2f}" y="{ygene + fs * 0.35:.2f}" font-size="{fs:.1f}" '
                        f'text-anchor="end">{esc(lab)}</text>')
            body.append(f'<line x1="{x0:.2f}" y1="{ygene:.2f}" x2="{x1:.2f}" y2="{ygene:.2f}" '
                        f'stroke="#000" stroke-width="0.5"/>')
            # genes: others first, then family on top
            order = sorted(feats, key=lambda f: (f["kind"] in ("family", "named", "unnamed"), f["kind"] in bn))
            for f in order:
                xa, xb = sorted((X(max(f["start"], lo)), X(min(f["end"], hi) + 1)))
                st = P.flipstrand(f["strand"]) if flip else f["strand"]
                pts, _ = P.arrow(xa, xb, ygene, gene_h, st, 0.8)
                col = fill_of(f)
                fill = col if f["filled"] else "#FFFFFF"
                dash = ' stroke-dasharray="1.2,0.8"' if f["dashed"] else ""
                body.append(f'<polygon points="{pts}" fill="{fill}" stroke="{col if col != OTHER else "#9A9A9A"}" '
                            f'stroke-width="{0.5 if f["filled"] else 0.8}"{dash}><title>'
                            f'{esc(f["id"])} | {esc(f["desc"])} | {contig}:{f["start"]}-{f["end"]} {f["strand"]}'
                            f'</title></polygon>')
            # labels: spread, leader, rotated text
            anchors = [X((max(f["start"], lo) + min(f["end"], hi)) / 2) for f in labelled]
            pos = spread(anchors, dlab, x0, x1 - fs * cos)
            ybase = ygene - gene_h / 2 - 3 * mm
            for f, ax, px in zip(labelled, anchors, pos):
                body.append(f'<polyline points="{ax:.2f},{ygene - gene_h / 2 - 0.3:.2f} {ax:.2f},'
                            f'{ygene - gene_h / 2 - 1.2 * mm:.2f} {px:.2f},{ybase:.2f}" fill="none" '
                            f'stroke="#7F7F7F" stroke-width="0.3"/>')
                weight = ' font-weight="bold"' if f["kind"] == "named" else ""
                body.append(f'<text transform="translate({px + 0.6:.2f},{ybase - 0.6:.2f}) rotate(-{ANGLE})" '
                            f'font-size="{fs:.1f}"{weight}>{esc(f["label"])}</text>')
            y = ygene + gene_h / 2 + fs * 0.8
        y += fs * 0.4
    y_bottom = y

    # ---- bands linking the panels ----
    bands = []
    if not a.no_bands:
        for xa, xb, col in band_spans:
            w = max(xb - xa, 0.8)
            bands.append(f'<rect x="{(xa + xb) / 2 - w / 2:.2f}" y="{panel_tops[0]:.2f}" width="{w:.2f}" '
                         f'height="{y_bottom - panel_tops[0]:.2f}" fill="{col}" fill-opacity="0.12"/>')

    # ---- scale bar and legend ----
    y += fs * 0.6
    sb = P.nice_scale(span)
    sw = sb / span * (x1 - x0)
    foot = [f'<line x1="{x0:.2f}" y1="{y:.2f}" x2="{x0 + sw:.2f}" y2="{y:.2f}" stroke="#000" stroke-width="0.7"/>',
            f'<text x="{x0 + sw + 1.5 * mm:.2f}" y="{y + fs * 0.35:.2f}" font-size="{fs:.1f}">{P.fmt_kb(sb)}</text>']
    rowAB = [(fcol[g], P.legend_name(g), True) for g in fam_groups]
    rowAB += [(BOUNDARY, f"{bn[0]} / {bn[1]}", True), (OTHER, "other gene", True),
              ("#555555", "pseudogene / TOGA2 not intact", False)]
    rowC = [(CLASS_COL[k], v, True) for k, v in CLASS_LEGEND.items()
            if any(f["kind"] == "named" and f["group"] == k for f in named)]
    if any(f["kind"] == "unnamed" for f in named):
        rowC.append((UNNAMED, "exendin-family, not named (not in tree)", False))
    xs = x0 + sw + 12 * mm
    xl, yl = xs, y
    for head, items in (("A, B: annotation name", rowAB), ("C: step-10 class", rowC)):
        foot.append(f'<text x="{xl:.2f}" y="{yl + fs * 0.35:.2f}" font-size="{fs:.1f}" font-weight="bold">'
                    f'{esc(head)}</text>')
        xl += P.text_w(head, fs) + 3 * mm
        for col, text, filled in items:
            w = P.text_w(text, fs) + 6 * mm
            if xl + w > W - 2 * mm:
                xl, yl = xs + P.text_w(head, fs) + 3 * mm, yl + fs * 1.5
            foot.append(f'<rect x="{xl:.2f}" y="{yl - fs * 0.45:.2f}" width="{2.6 * mm:.2f}" height="{fs * 0.9:.2f}" '
                        f'fill="{col if filled else "#FFFFFF"}" stroke="{col}" stroke-width="0.6"/>'
                        f'<text x="{xl + 3.4 * mm:.2f}" y="{yl + fs * 0.35:.2f}" font-size="{fs:.1f}">{esc(text)}</text>')
            xl += w
        xl, yl = xs, yl + fs * 1.6
    yl -= fs * 1.6
    H = yl + fs * 1.2

    svg = (f'<svg xmlns="http://www.w3.org/2000/svg" width="{W:.2f}pt" height="{H:.2f}pt" '
           f'viewBox="0 0 {W:.2f} {H:.2f}" font-family="{P.FONT}">\n'
           f'<rect width="100%" height="100%" fill="#FFFFFF"/>\n'
           + "\n".join(bands + body + foot) + "\n</svg>\n")
    P.write(svg, a)
    if a.svg:
        N.info(f"wrote {a.svg}")
    N.info(f"figure {W / mm:.0f} x {H / mm:.0f} mm")

    # ---- crosswalk ----
    if a.crosswalk:
        loci = [f for f in named if f["kind"] in ("named", "unnamed")]
        tracks = [(lab, f) for lab, f in tracksA] + [("TOGA2", toga)]
        with open(a.crosswalk, "w") as fh:
            fh.write("\t".join(["name", "class", "contig", "start", "end", "strand"]
                               + [f"{lab}_ids" for lab, _ in tracks] + [f"{lab}_labels" for lab, _ in tracks]) + "\n")
            for f in sorted(loci, key=lambda f: f["start"]):
                ids, labs = [], []
                for lab, tr in tracks:
                    hit = [g for g in tr if overlap(f["start"], f["end"], g["start"], g["end"])]
                    ids.append(",".join(g["id"] for g in hit) or "none")
                    labs.append(" | ".join(g["label"] for g in hit) or "none")
                fh.write("\t".join(map(str, [f["id"] if f["kind"] == "named" else "(not named)",
                                             f["group"] if f["kind"] == "named" else "unnamed",
                                             contig, f["start"], f["end"], f["strand"]] + ids + labs)) + "\n")
        N.info(f"wrote {a.crosswalk}")
        for f in sorted(loci, key=lambda f: f["start"]):
            parts = []
            for lab, tr in tracks:
                hit = [g["id"] for g in tr if overlap(f["start"], f["end"], g["start"], g["end"])]
                parts.append(f"{lab} {','.join(hit) or 'none'}")
            N.info(f"  {f['label'] or '(not named)':12s} {contig}:{f['start']}-{f['end']}  " + "; ".join(parts))


if __name__ == "__main__":
    main()
