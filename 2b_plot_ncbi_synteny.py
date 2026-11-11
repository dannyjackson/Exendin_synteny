#!/usr/bin/env python3
"""
Step-2 synteny figure, publication layout: the exendin-family genes that each species'
NCBI annotation puts between SLC4A8 and AQP5, species in Cactus-tree order.

What is drawn (only the interval)
  - SLC4A8 and AQP5 on one contig: one line from SLC4A8 to AQP5. Family genes outside the
    pair (2a's --pad-kb margin) are NOT drawn; they are listed in the log.
  - SLC4A8 and AQP5 on different contigs (assembly break): two lines. The SLC4A8 contig is
    drawn left-aligned with the family genes on its interval side, the AQP5 contig
    right-aligned with the family genes on its interval side, each within --near-kb of the
    boundary. Family genes on the far side of a boundary are left out and logged.
  - Only one boundary annotated: one line, as for a split genome.
  - Family genes on a contig with no boundary gene (2a 'unplaced'): left out unless
    --unplaced, then drawn on their own line. Logged either way.
  Every line is oriented with SLC4A8 on the left. A contig with one boundary is oriented by
  that gene's strand in --reference-query (HelSus).

Style
  Width --width-mm (180 = double column), all text --font-pt (7) at that size, Arial/
  Helvetica, no title, no contig names or notes. Species are labelled by name (italic,
  'H. suspectum'; --full-names for 'Heloderma suspectum'), from --species-table. A scale bar
  and one legend replace the axis and the per-gene labels (--labels puts labels back).
  Boundary genes longer than --clip-kb are drawn clipped on their outer side. Stretches
  longer than --max-gap-kb with nothing drawn are shortened and marked //.
  Filled = protein-coding; hollow = pseudogene / other biotype.

Outputs: --svg, --pdf and --png (--dpi, default 600) via rsvg-convert; --layout TSV of
every drawn feature with true coordinates.

Usage
  python3 2b_plot_ncbi_synteny.py --loci ncbi_loci/ncbi_loci.tsv \
      --species-table ../species_toga2.tsv --species-tree ${SPTREE} --top GalGal \
      --svg ncbi_synteny.svg --pdf ncbi_synteny.pdf --png ncbi_synteny.png \
      --layout ncbi_synteny.tsv 2> ncbi_synteny.log
"""
import argparse, csv, html, os, shutil, subprocess, sys, tempfile
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import lib_exendin_names as N  # noqa: E402

S = N.load_engine()            # parse_newick / match_tip / prune / species_order only

PT_PER_MM = 72 / 25.4
BOUNDARY_FILL = {0: "#4D4D4D", 1: "#B3B3B3"}       # left (SLC4A8), right (AQP5)
FONT = "Arial, Helvetica, sans-serif"


def legend_name(group):
    return {"GCGL (glucagon-like)": "GCGL"}.get(group, group)


def species_label(code, sp, full):
    name = (sp or "").replace("_", " ").split()
    if len(name) < 2:
        return code, False
    return (" ".join(name[:2]) if full else f"{name[0][0]}. {name[1]}"), True


# ---------------- lines: which features, oriented, clipped ----------------
def oriented(f, flip):
    a, b = (-f["end"], -f["start"]) if flip else (f["start"], f["end"])
    s = f["strand"] if not flip else {"+": "-", "-": "+"}.get(f["strand"], f["strand"])
    return dict(f, a=a, b=b, ostrand=s)


def make_lines(code, rows, bnames, ref_bstrand, near, clip, unplaced):
    """-> list of lines: dict(align='left'|'right'|'full', feats=[...]); and dropped [(row, why)]."""
    left, right = bnames
    bnd = {r["group"]: r for r in rows if r["kind"] == "boundary"}
    fam = [r for r in rows if r["kind"] == "exendin" and r["placement"] != "elsewhere"]
    for r in rows:
        if r["kind"] in ("boundary", "exendin"):
            r["start"], r["end"] = int(r["start"]), int(r["end"])
    lines, dropped, used = [], [], set()

    def clip_boundary(f, outer_left):
        if clip and f["b"] - f["a"] > clip:
            f["true_len"] = f["b"] - f["a"]
            if outer_left:
                f["a"] = f["b"] - clip
            else:
                f["b"] = f["a"] + clip
        return f

    # both boundaries on one contig (a pair, or two singles far apart): one line between them
    if left in bnd and right in bnd and bnd[left]["contig"] == bnd[right]["contig"]:
        L, R = bnd[left], bnd[right]
        flip = L["start"] > R["start"]
        lo, hi = min(L["start"], R["start"]), max(L["end"], R["end"])
        feats = [clip_boundary(oriented(L, flip), True), clip_boundary(oriented(R, flip), False)]
        for g in fam:
            if g["contig"] == L["contig"] and g["end"] >= lo and g["start"] <= hi:
                feats.append(oriented(g, flip)); used.add(id(g))
        lines.append(dict(align="full", contig=L["contig"], feats=feats))
    else:
        for side, name in enumerate(bnames):
            if name not in bnd:
                continue
            B = bnd[name]
            o = oriented(B, False)
            flip = (ref_bstrand.get(name) is not None and o["ostrand"] != ref_bstrand[name])
            ob = oriented(B, flip)
            feats = [clip_boundary(ob, side == 0)]
            for g in fam:
                if g["contig"] != B["contig"]:
                    continue
                og = oriented(g, flip)
                inside = og["a"] >= ob["b"] if side == 0 else og["b"] <= ob["a"]
                dist = (og["a"] - ob["b"]) if side == 0 else (ob["a"] - og["b"])
                if inside and dist <= near:
                    feats.append(og); used.add(id(g))
                elif not inside:
                    dropped.append((g, f"on the far side of {name}"))
                    used.add(id(g))
                else:
                    dropped.append((g, f"{dist / 1000:.0f} kb from {name} (> --near-kb)"))
                    used.add(id(g))
            lines.append(dict(align="left" if side == 0 else "right", contig=B["contig"], feats=feats))
    rest = [g for g in fam if id(g) not in used]
    by_ctg = defaultdict(list)
    for g in rest:
        on_line = any(g["contig"] == ln["contig"] for ln in lines)
        if on_line:
            dropped.append((g, "outside the SLC4A8-AQP5 interval"))
        elif unplaced:
            by_ctg[g["contig"]].append(g)
        else:
            dropped.append((g, "contig with no boundary gene (use --unplaced to draw)"))
    for c, gs in sorted(by_ctg.items()):
        lines.append(dict(align="left", contig=c, feats=[oriented(g, False) for g in gs], unplaced=True))
    return lines, dropped


def compress(line, max_gap, brk):
    """Shorten empty stretches > max_gap to brk (axis units); x = compressed start, w = width."""
    fs = sorted(line["feats"], key=lambda f: f["a"])
    shift, end, breaks = 0, None, []
    for f in fs:
        if end is not None and max_gap and f["a"] - shift - end > max_gap:
            g = f["a"] - shift - end
            breaks.append(end + brk / 2)
            shift += g - brk
        f["x0"], f["x1"] = f["a"] - shift, f["b"] - shift
        end = f["x1"] if end is None else max(end, f["x1"])
    base = min(f["x0"] for f in fs)
    for f in fs:
        f["x0"] -= base; f["x1"] -= base
    line["feats"], line["breaks"] = fs, [b - base for b in breaks]
    line["span"] = max(f["x1"] for f in fs)


def nice_scale(span):
    best = 1e3
    for v in (1e3, 2e3, 5e3, 1e4, 2e4, 5e4, 1e5, 2e5, 5e5, 1e6, 2e6, 5e6, 1e7):
        if v <= span / 4:
            best = v
    return best


def fmt_kb(v):
    return f"{v / 1e6:g} Mb" if v >= 1e6 else f"{v / 1e3:g} kb"


# ---------------- drawing ----------------
def arrow(x0, x1, y, h, strand, min_w):
    if x1 - x0 < min_w:
        c = (x0 + x1) / 2; x0, x1 = c - min_w / 2, c + min_w / 2
    head = min(h * 0.45, (x1 - x0) * 0.6)
    t, b = y - h / 2, y + h / 2
    if strand == "+":
        pts = [(x0, t), (x1 - head, t), (x1, y), (x1 - head, b), (x0, b)]
    elif strand == "-":
        pts = [(x1, t), (x0 + head, t), (x0, y), (x0 + head, b), (x1, b)]
    else:
        pts = [(x0, t), (x1, t), (x1, b), (x0, b)]
    return " ".join(f"{px:.2f},{py:.2f}" for px, py in pts), (x0, x1)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--loci", required=True, help="ncbi_loci.tsv from 2a")
    ap.add_argument("--species-table", help="species.tsv / species_toga2.tsv, for species names")
    ap.add_argument("--species-tree", help="Newick species tree (order + cladogram)")
    ap.add_argument("--tree-map", help="TSV tree_tip<TAB>code for tips that do not match automatically")
    ap.add_argument("--top", default="", help="comma-separated codes drawn first, above the tree (GalGal)")
    ap.add_argument("--reference-query", default="HelSus", help="sets the strand of single boundaries")
    ap.add_argument("--boundary-names", default="SLC4A8,AQP5", help="as in 2a --boundaries, left one first")
    ap.add_argument("--near-kb", type=float, default=1000.0,
                    help="split genomes: family genes this far from a boundary, on its interval side, are drawn")
    ap.add_argument("--unplaced", action="store_true", help="draw family genes on contigs with no boundary")
    ap.add_argument("--clip-kb", type=float, default=15.0, help="boundary genes drawn at most this long")
    ap.add_argument("--max-gap-kb", type=float, default=1000.0, help="empty stretches longer than this are cut (//)")
    ap.add_argument("--break-kb", type=float, default=20.0, help="drawn width of a cut stretch (kb)")
    ap.add_argument("--width-mm", type=float, default=180.0, help="figure width (180 double, 88 single column)")
    ap.add_argument("--font-pt", type=float, default=7.0, help="size of all text")
    ap.add_argument("--full-names", action="store_true", help="'Heloderma suspectum' instead of 'H. suspectum'")
    ap.add_argument("--labels", action="store_true", help="label each family gene (Ex2, GCGL ...)")
    ap.add_argument("--no-tree", action="store_true")
    ap.add_argument("--svg"); ap.add_argument("--pdf"); ap.add_argument("--png")
    ap.add_argument("--dpi", type=int, default=600)
    ap.add_argument("--layout", help="TSV of every drawn feature (true coordinates)")
    a = ap.parse_args()

    bnames = a.boundary_names.split(",")
    near, clip = a.near_kb * 1000, a.clip_kb * 1000
    names = {}
    if a.species_table:
        for r in N.read_species_table(a.species_table):
            names[r["code"]] = r["species"]

    with open(a.loci) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    by_sp, queries, empty = defaultdict(list), [], {}
    for r in rows:
        if r["species"] not in queries:
            queries.append(r["species"])
        if r["kind"] == "none":
            empty[r["species"]] = {"no_annotation": "no NCBI annotation", "no_genome": "no genome"}.get(
                r["placement"], r["placement"])
        else:
            by_sp[r["species"]].append(r)

    # strand of each boundary in the reference, oriented SLC4A8-left
    ref_bstrand = {}
    for q in [a.reference_query] + queries:
        b = {r["group"]: r for r in by_sp.get(q, []) if r["kind"] == "boundary"}
        if len(b) == 2 and b[bnames[0]]["contig"] == b[bnames[1]]["contig"]:
            flip = int(b[bnames[0]]["start"]) > int(b[bnames[1]]["start"])
            for n in bnames:
                s = b[n]["strand"]
                ref_bstrand[n] = {"+": "-", "-": "+"}.get(s, s) if flip else s
            N.info(f"orientation of single boundaries from {q}: " +
                   ", ".join(f"{n} {s}" for n, s in ref_bstrand.items()))
            break

    lines_by, n_drop = {}, 0
    for q in queries:
        if q in empty:
            continue
        lines, dropped = make_lines(q, by_sp[q], bnames, ref_bstrand, near, clip, a.unplaced)
        for g, why in dropped:
            n_drop += 1
            N.warn(f"{q}: {g['symbol']} ({legend_name(g['group'])}, {g['contig']}:{g['start']}-{g['end']}) "
                   f"not drawn: {why}")
        for ln in lines:
            compress(ln, a.max_gap_kb * 1000, a.break_kb * 1000)
        lines_by[q] = lines
        kinds = "+".join(ln["align"] for ln in lines)
        nfam = sum(f["kind"] == "exendin" for ln in lines for f in ln["feats"])
        N.info(f"{q}: {len(lines)} line(s) [{kinds}], {nfam} family gene(s) drawn")
        if not lines:
            empty[q] = "no SLC4A8 or AQP5 annotated"

    # order and tree
    tree = None
    top = [q for q in a.top.split(",") if q in queries]
    if a.species_tree:
        order, tree = S.species_order(queries, a.reference_query, a.species_tree, a.tree_map, None)
    else:
        order = list(queries)
    order = top + [q for q in order if q not in top]
    if a.no_tree:
        tree = None
    if tree is not None:
        intree = {lf.name for lf in S.leaves(tree)}
        tree = S.prune(tree, intree - set(top))

    # ---------- geometry (points) ----------
    fs = a.font_pt
    W = a.width_mm * PT_PER_MM
    margin = 2.0
    tree_w = 0 if tree is None else 55.0
    labels = {q: species_label(q, names.get(q), a.full_names) for q in order}
    name_w = max(len(t) for t, _ in labels.values()) * fs * 0.46 + 8
    x0 = margin + tree_w + name_w
    x1 = W - margin - 2
    line_h, line_gap, row_gap = fs * 0.95, fs * 1.25, fs * 0.75
    gene_h = fs * 0.95
    header_h = fs * 1.6
    span = max([ln["span"] for ls in lines_by.values() for ln in ls] or [1])
    scale = (x1 - x0) / span
    min_w = 1.6

    y = margin + header_h
    rows_y, centre = [], {}
    for i, q in enumerate(order):
        nl = max(1, len(lines_by.get(q, [])))
        if i and q not in top and order[i - 1] in top:
            y += row_gap * 1.5                       # gap between --top rows and the tree
        ys = [y + line_h / 2 + k * line_gap for k in range(nl)]
        rows_y.append((q, ys))
        centre[q] = sum(ys) / len(ys)
        y = ys[-1] + line_h / 2 + row_gap
    groups = sorted({f["group"] for ls in lines_by.values() for ln in ls for f in ln["feats"]
                     if f["kind"] == "exendin"},
                    key=lambda g: (not g.startswith("GCGL"), not g.startswith("exendin-"), g))
    cmap = N.colour_map(groups)

    out, lay = [], []
    T = lambda x, yy, s, **kw: out.append(
        f'<text x="{x:.2f}" y="{yy:.2f}" font-size="{fs}"'
        + "".join(f' {k.replace("_", "-")}="{v}"' for k, v in kw.items()) + f">{html.escape(s)}</text>")

    # column headers: the boundary genes
    T(x0, margin + fs, bnames[0], font_style="italic")
    T(x1, margin + fs, bnames[1], font_style="italic", text_anchor="end")

    for q, ys in rows_y:
        lab, ital = labels[q]
        T(margin + tree_w + 3, centre[q] + fs * 0.35, lab, **({"font_style": "italic"} if ital else {}))
        lines = lines_by.get(q, [])
        if not lines:
            T(x0, centre[q] + fs * 0.35, empty.get(q, "nothing annotated"), fill="#808080")
            continue
        for ln, yy in zip(lines, ys):
            off = x0 if ln["align"] in ("left", "full") else x1 - ln["span"] * scale
            X = lambda v: off + v * scale
            # backbone, interrupted at breaks
            segs, cur = [], 0.0
            for bx in ln["breaks"]:
                segs.append((cur, bx - a.break_kb * 250)); cur = bx + a.break_kb * 250
            segs.append((cur, ln["span"]))
            for s0, s1 in segs:
                out.append(f'<line x1="{X(s0):.2f}" x2="{X(s1):.2f}" y1="{yy:.2f}" y2="{yy:.2f}" '
                           f'stroke="#000" stroke-width="0.5"/>')
            for bx in ln["breaks"]:
                for d in (-1.0, 1.0):
                    cx = X(bx) + d
                    out.append(f'<line x1="{cx - 1:.2f}" x2="{cx + 1:.2f}" y1="{yy + 2.2:.2f}" '
                               f'y2="{yy - 2.2:.2f}" stroke="#000" stroke-width="0.5"/>')
            for f in ln["feats"]:
                pts, (gx0, gx1) = arrow(X(f["x0"]), X(f["x1"]), yy, gene_h, f["ostrand"], min_w)
                tip = f'{f["symbol"]} | {f["description"]} | {f["contig"]}:{f["start"]}-{f["end"]}'
                if f["kind"] == "boundary":
                    col = BOUNDARY_FILL[bnames.index(f["group"]) if f["group"] in bnames else 0]
                    out.append(f'<polygon points="{pts}" fill="{col}" stroke="none">'
                               f'<title>{html.escape(tip)}</title></polygon>')
                else:
                    col = cmap[f["group"]]
                    filled = f["biotype"] == "protein_coding"
                    out.append(f'<polygon points="{pts}" fill="{col if filled else "#fff"}" stroke="{col}" '
                               f'stroke-width="{0.4 if filled else 0.7}"><title>{html.escape(tip)}</title></polygon>')
                    if a.labels:
                        T((gx0 + gx1) / 2, yy - gene_h / 2 - 1.2, f["label"], text_anchor="middle")
                lay.append([q, ln["contig"], ln["align"], f["kind"], f["symbol"], f["group"], f["start"],
                            f["end"], f["strand"], f["ostrand"], f"{gx0:.2f}", f"{gx1:.2f}"])

    # cladogram (topology), leaves at row centres
    if tree is not None:
        def depth(n):
            return 0 if not n.children else 1 + max(depth(c) for c in n.children)
        D = depth(tree) or 1
        TX = lambda d: margin + d / D * (tree_w - 4)

        def place(n, d):
            if not n.children:
                n.x, n.y = D, centre[n.name]
            else:
                for c in n.children:
                    place(c, d + 1)
                n.x = d
                n.y = (min(c.y for c in n.children) + max(c.y for c in n.children)) / 2

        place(tree, 0)

        def draw(n):
            for c in n.children:
                out.append(f'<path d="M{TX(n.x):.2f},{c.y:.2f} H{TX(c.x):.2f}" stroke="#000" '
                           f'stroke-width="0.5" fill="none"/>')
                draw(c)
            if n.children:
                ys_ = [c.y for c in n.children]
                out.append(f'<line x1="{TX(n.x):.2f}" x2="{TX(n.x):.2f}" y1="{min(ys_):.2f}" '
                           f'y2="{max(ys_):.2f}" stroke="#000" stroke-width="0.5"/>')
        draw(tree)

    # scale bar and legend, one line under the panel
    y += fs * 0.6
    sb = nice_scale(span)
    out.append(f'<line x1="{x0:.2f}" x2="{x0 + sb * scale:.2f}" y1="{y:.2f}" y2="{y:.2f}" '
               f'stroke="#000" stroke-width="0.8"/>')
    T(x0 + sb * scale + 3, y + fs * 0.35, fmt_kb(sb))
    items = [(legend_name(g), cmap[g], True) for g in groups]
    if any(f["kind"] == "exendin" and f["biotype"] != "protein_coding"
           for ls in lines_by.values() for ln in ls for f in ln["feats"]):
        items.append(("pseudogene", "#000", False))
    items += [(bnames[0], BOUNDARY_FILL[0], True), (bnames[1], BOUNDARY_FILL[1], True)]
    lx = x0 + sb * scale + 3 + len(fmt_kb(sb)) * fs * 0.55 + 14
    ly = y
    for text, col, filled in items:
        w = 9 + 3 + len(text) * fs * 0.52 + 10
        if lx + w > x1 + 2:
            lx, ly = margin + tree_w + 3, ly + fs * 1.5
        out.append(f'<rect x="{lx:.2f}" y="{ly - gene_h / 2:.2f}" width="9" height="{gene_h:.2f}" '
                   f'fill="{col if filled else "#fff"}" stroke="{col}" stroke-width="{0.4 if filled else 0.7}"/>')
        T(lx + 12, ly + fs * 0.35, text)
        lx += w
    H = ly + fs * 0.9 + margin

    svg = (f'<svg xmlns="http://www.w3.org/2000/svg" width="{W:.1f}pt" height="{H:.1f}pt" '
           f'viewBox="0 0 {W:.2f} {H:.2f}" font-family="{FONT}">\n'
           f'<rect width="100%" height="100%" fill="#fff"/>\n' + "\n".join(out) + "\n</svg>\n")
    if a.svg:
        open(a.svg, "w").write(svg)
    else:
        sys.stdout.write(svg)
    for path, fmt in ((a.pdf, "pdf"), (a.png, "png")):
        if not path:
            continue
        if not shutil.which("rsvg-convert"):
            N.warn(f"rsvg-convert not found: {path} not written (mamba install -c conda-forge librsvg)")
            continue
        with tempfile.NamedTemporaryFile("w", suffix=".svg", delete=False) as fh:
            fh.write(svg); tmp = fh.name
        cmd = ["rsvg-convert", "-f", fmt, "-o", path, tmp]
        if fmt == "png":
            cmd[3:3] = ["-d", str(a.dpi), "-p", str(a.dpi)]
        r = subprocess.run(cmd, capture_output=True, text=True)
        os.unlink(tmp)
        (N.info if r.returncode == 0 else N.warn)(f"{path}: " + ("written" if r.returncode == 0 else r.stderr.strip()))
    if a.layout:
        with open(a.layout, "w") as fh:
            fh.write("species\tcontig\tline\tkind\tsymbol\tgroup\tstart\tend\tstrand\tstrand_drawn\tx0_pt\tx1_pt\n")
            for r in lay:
                fh.write("\t".join(map(str, r)) + "\n")
    N.info(f"{len(order)} species, {sum(len(v) for v in lines_by.values())} lines, axis {fmt_kb(span)}, "
           f"{W / PT_PER_MM:.0f} x {H / PT_PER_MM:.0f} mm, text {fs:g} pt; {n_drop} family gene(s) not drawn")


if __name__ == "__main__":
    main()
