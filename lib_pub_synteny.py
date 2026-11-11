#!/usr/bin/env python3
"""
Publication-layout synteny figure of the SLC4A8-AQP5 interval. Not run on its own; used
by 2b_plot_ncbi_synteny.py (NCBI annotations) and 3d_plot_toga2_interval.py (TOGA2 runs).

Input to draw(): one entry per species
    {code: dict(bnd={'SLC4A8': feat, 'AQP5': feat}, genes=[feat, ...], empty='reason' or '')}
  feat: dict(contig, start, end, strand, symbol, group, family (bool), filled (bool),
             dashed (bool), label, hover)
  group is the colour key: a family name group (lib_exendin_names.name_group), a gene
  symbol, or 'uncharacterised'.

What is drawn (only the interval)
  - SLC4A8 and AQP5 on one contig: one line from SLC4A8 to AQP5, every gene in between.
  - On different contigs (assembly break), or only one of them annotated: one line per
    boundary contig. SLC4A8 line left-aligned, AQP5 line right-aligned, each with the genes
    on its interval side within --near-kb of the boundary (--split-genes family|all).
  - Family genes on a contig with no boundary: left out unless --unplaced.
  Every line is oriented with SLC4A8 on the left; a contig with only one boundary is
  oriented by that gene's strand in --reference-query. Family genes left out are logged.

Style: --width-mm wide, all text --font-pt, Arial/Helvetica, no title, no contig names.
Species names in italics from --species-table. A scale bar and one legend; boundary
genes clipped to --clip-kb on their outer side; empty stretches over --max-gap-kb cut
and marked //. Outputs SVG, and PDF / PNG through rsvg-convert.
"""
import html, os, re, shutil, subprocess, sys, tempfile
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import lib_exendin_names as N  # noqa: E402

S = N.load_engine()            # parse_newick / leaves / prune / species_order only

PT_PER_MM = 72 / 25.4
FONT = "Arial, Helvetica, sans-serif"
BOUNDARY_FILL = ("#4D4D4D", "#B3B3B3")          # left (SLC4A8), right (AQP5)
UNCHAR = "uncharacterised"
UNCHAR_FILL = "#E6E6E6"
# light tab20 halves: other (non-family) genes, kept paler than the Okabe-Ito family colours
OTHER_FILLS = ["#AEC7E8", "#FFBB78", "#98DF8A", "#FF9896", "#C5B0D5", "#C49C94", "#F7B6D2",
               "#DBDB8D", "#9EDAE5", "#BCBD22", "#17BECF", "#8C564B"]


def add_args(ap, near_default=1000.0, split_default="family"):
    ap.add_argument("--species-table", help="species.tsv / species_toga2.tsv, for species names")
    ap.add_argument("--species-tree", help="Newick species tree (order + cladogram)")
    ap.add_argument("--tree-map", help="TSV tree_tip<TAB>code for tips that do not match automatically")
    ap.add_argument("--no-tree", action="store_true")
    ap.add_argument("--top", default="", help="comma-separated codes drawn first, above the tree")
    ap.add_argument("--reference-query", default="HelSus", help="sets the strand of single boundaries")
    ap.add_argument("--near-kb", type=float, default=near_default,
                    help="split genomes: genes this far from a boundary, on its interval side, are drawn "
                         "(0 = the longest complete interval among the species)")
    ap.add_argument("--split-genes", choices=("family", "all"), default=split_default,
                    help="split genomes: draw only family genes beside a boundary, or every gene")
    ap.add_argument("--unplaced", action="store_true", help="draw family genes on contigs with no boundary")
    ap.add_argument("--clip-kb", type=float, default=15.0, help="boundary genes drawn at most this long")
    ap.add_argument("--max-gap-kb", type=float, default=1000.0, help="empty stretches longer than this are cut (//)")
    ap.add_argument("--break-kb", type=float, default=20.0, help="drawn width of a cut stretch (kb)")
    ap.add_argument("--width-mm", type=float, default=180.0, help="figure width (180 double, 88 single column)")
    ap.add_argument("--font-pt", type=float, default=7.0, help="size of all text")
    ap.add_argument("--full-names", action="store_true", help="'Heloderma suspectum' instead of 'H. suspectum'")
    ap.add_argument("--labels", action="store_true", help="label every non-boundary gene")
    ap.add_argument("--svg"); ap.add_argument("--pdf"); ap.add_argument("--png")
    ap.add_argument("--dpi", type=int, default=600)
    ap.add_argument("--layout", help="TSV of every drawn feature (true coordinates)")


def text_w(s, fs):
    """Approximate width of s in Arial at fs pt."""
    return fs * sum(0.67 if c.isupper() else 0.56 if c.isdigit() else 0.28 if c in " .,-/()il" else 0.52
                    for c in s)


def species_label(code, sp, full):
    name = (sp or "").replace("_", " ").split()
    if len(name) < 2:
        return code, False
    return (" ".join(name[:2]) if full else f"{name[0][0]}. {name[1]}"), True


def legend_name(group):
    return {"GCGL (glucagon-like)": "GCGL"}.get(group, group)


def flipstrand(s):
    return {"+": "-", "-": "+"}.get(s, s)


def oriented(f, flip):
    a, b = (-f["end"], -f["start"]) if flip else (f["start"], f["end"])
    return dict(f, a=a, b=b, ostrand=flipstrand(f["strand"]) if flip else f["strand"])


def ref_strands(species, bnames, ref_query):
    """Strand of each boundary, oriented SLC4A8-left, in the first species with both on one contig."""
    for q in [ref_query] + list(species):
        b = species.get(q, {}).get("bnd", {})
        if all(n in b for n in bnames) and b[bnames[0]]["contig"] == b[bnames[1]]["contig"]:
            flip = b[bnames[0]]["start"] > b[bnames[1]]["start"]
            out = {n: (flipstrand(b[n]["strand"]) if flip else b[n]["strand"]) for n in bnames}
            N.info(f"single boundaries are oriented by {q}: " + ", ".join(f"{n} {s}" for n, s in out.items()))
            return out
    N.warn("no species has both boundaries on one contig: single boundaries are not reoriented")
    return {}


def make_lines(code, sp, bnames, ref_bstrand, near, clip, split_genes, unplaced):
    """-> lines [dict(align, contig, feats)], dropped family genes [(feat, why)]."""
    left, right = bnames
    bnd, genes = sp["bnd"], sp["genes"]
    lines, dropped, used = [], [], set()

    def clip_b(f, outer_left):
        f = dict(f, kind="boundary", bname=bnames[0 if outer_left else 1])
        if clip and f["b"] - f["a"] > clip:
            if outer_left:
                f["a"] = f["b"] - clip
            else:
                f["b"] = f["a"] + clip
        return f

    if left in bnd and right in bnd and bnd[left]["contig"] == bnd[right]["contig"]:
        L, R = bnd[left], bnd[right]
        flip = L["start"] > R["start"]
        lo, hi = min(L["end"], R["end"]), max(L["start"], R["start"])
        feats = [clip_b(oriented(L, flip), True), clip_b(oriented(R, flip), False)]
        for g in genes:
            if g["contig"] == L["contig"] and g["end"] > lo and g["start"] < hi:
                feats.append(dict(oriented(g, flip), kind="gene")); used.add(id(g))
        lines.append(dict(align="full", contig=L["contig"], feats=feats))
    else:
        for side, name in enumerate(bnames):
            if name not in bnd:
                continue
            B = bnd[name]
            flip = name in ref_bstrand and B["strand"] != ref_bstrand[name]
            ob = oriented(B, flip)
            feats = [clip_b(ob, side == 0)]
            for g in genes:
                if g["contig"] != B["contig"]:
                    continue
                og = oriented(g, flip)
                inside = og["a"] >= ob["b"] if side == 0 else og["b"] <= ob["a"]
                dist = (og["a"] - ob["b"]) if side == 0 else (ob["a"] - og["b"])
                used.add(id(g))
                if inside and dist <= near and (split_genes == "all" or g["family"]):
                    feats.append(dict(og, kind="gene"))
                elif g["family"]:
                    dropped.append((g, f"on the far side of {name}" if not inside
                                    else f"{dist / 1000:.0f} kb from {name} (> --near-kb)"))
            lines.append(dict(align="left" if side == 0 else "right", contig=B["contig"], feats=feats))
    by_ctg = defaultdict(list)
    for g in genes:
        if id(g) in used or not g["family"]:
            continue
        if any(g["contig"] == ln["contig"] for ln in lines):
            dropped.append((g, "outside the SLC4A8-AQP5 interval"))
        elif unplaced:
            by_ctg[g["contig"]].append(g)
        else:
            dropped.append((g, "contig with no boundary gene (use --unplaced to draw)"))
    for c, gs in sorted(by_ctg.items()):
        lines.append(dict(align="left", contig=c, feats=[dict(oriented(g, False), kind="gene") for g in gs]))
    return lines, dropped


def compress(line, max_gap, brk):
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


def colours(lines_by, order):
    """family groups: fixed colours; other genes: OTHER_FILLS in order of first appearance."""
    fam, other = set(), []
    for q in order:
        for ln in lines_by.get(q, []):
            for f in sorted(ln["feats"], key=lambda f: f["x0"]):
                if f["kind"] != "gene":
                    continue
                if f["family"]:
                    fam.add(f["group"])
                elif f["group"] != UNCHAR and f["group"] not in other:
                    other.append(f["group"])
    famg = sorted(fam, key=lambda g: (not g.startswith("GCGL"), not g.startswith("exendin-"), g))
    cm = N.colour_map(famg)
    if len(other) > len(OTHER_FILLS):
        N.warn(f"{len(other)} distinct non-family genes; colours repeat after {len(OTHER_FILLS)} "
               f"(add --labels to name them on the figure)")
    for i, g in enumerate(other):
        cm[g] = OTHER_FILLS[i % len(OTHER_FILLS)]
    cm[UNCHAR] = UNCHAR_FILL
    return famg, other, cm


def draw(species, a, bnames, legend_status=(), empty_default="nothing annotated"):
    """Lay out and write the figure. legend_status: [(text, filled, dashed)] after the gene colours."""
    names = {}
    if a.species_table:
        for r in N.read_species_table(a.species_table):
            names[r["code"]] = r["species"]
    queries = list(species)
    ref_bstrand = ref_strands(species, bnames, a.reference_query)

    near = a.near_kb * 1000
    if not near:
        spans = [abs(sp["bnd"][bnames[1]]["start"] - sp["bnd"][bnames[0]]["end"]) for sp in species.values()
                 if all(n in sp["bnd"] for n in bnames)
                 and sp["bnd"][bnames[0]]["contig"] == sp["bnd"][bnames[1]]["contig"]]
        near = max(spans) if spans else 1e6
        N.info(f"--near-kb 0: split genomes drawn up to {near / 1000:.0f} kb from each boundary "
               f"(longest complete interval)")

    lines_by, n_drop, empty = {}, 0, {}
    for q in queries:
        sp = species[q]
        if sp.get("empty"):
            empty[q] = sp["empty"]
            continue
        lines, dropped = make_lines(q, sp, bnames, ref_bstrand, near, a.clip_kb * 1000,
                                    a.split_genes, a.unplaced)
        for g, why in dropped:
            n_drop += 1
            N.warn(f"{q}: {g['symbol']} ({legend_name(g['group'])}, {g['contig']}:{g['start']}-{g['end']}) "
                   f"not drawn: {why}")
        for ln in lines:
            compress(ln, a.max_gap_kb * 1000, a.break_kb * 1000)
        lines_by[q] = lines
        nfam = sum(f["kind"] == "gene" and f["family"] for ln in lines for f in ln["feats"])
        ngen = sum(f["kind"] == "gene" for ln in lines for f in ln["feats"])
        N.info(f"{q}: {len(lines)} line(s) [{'+'.join(ln['align'] for ln in lines) or '-'}], "
               f"{ngen} gene(s) drawn, {nfam} of them family")
        if not lines:
            empty[q] = f"no {bnames[0]} or {bnames[1]} found"

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
    if tree is not None and top:
        tree = S.prune(tree, {lf.name for lf in S.leaves(tree)} - set(top))

    # geometry, in points
    fs = a.font_pt
    W = a.width_mm * PT_PER_MM
    margin = 2.0
    tree_w = 0 if tree is None else 55.0
    labels = {q: species_label(q, names.get(q), a.full_names) for q in order}
    name_w = max([text_w(t, fs) for t, _ in labels.values()] or [20]) + 10
    x0, x1 = margin + tree_w + name_w, W - margin - 2
    line_h, row_gap = fs * 0.95, fs * 0.75
    line_gap = fs * (2.35 if a.labels else 1.25)
    gene_h = fs * 0.95
    header_h = fs * (2.6 if a.labels else 1.6)
    span = max([ln["span"] for ls in lines_by.values() for ln in ls] or [1])
    scale = (x1 - x0) / span
    min_w = 1.6
    famg, other, cmap = colours(lines_by, order)

    y = margin + header_h
    rows_y, centre = [], {}
    for i, q in enumerate(order):
        nl = max(1, len(lines_by.get(q, [])))
        if i and q not in top and order[i - 1] in top:
            y += row_gap * 1.5
        if a.labels and i:
            y += fs * 1.1
        ys = [y + line_h / 2 + k * line_gap for k in range(nl)]
        rows_y.append((q, ys))
        centre[q] = sum(ys) / len(ys)
        y = ys[-1] + line_h / 2 + row_gap

    out, lay = [], []

    def T(x, yy, s, **kw):
        out.append(f'<text x="{x:.2f}" y="{yy:.2f}" font-size="{fs}"'
                   + "".join(f' {k.replace("_", "-")}="{v}"' for k, v in kw.items())
                   + f">{html.escape(s)}</text>")

    T(x0, margin + fs, bnames[0], font_style="italic")
    T(x1, margin + fs, bnames[1], font_style="italic", text_anchor="end")

    for q, ys in rows_y:
        lab, ital = labels[q]
        T(margin + tree_w + 3, centre[q] + fs * 0.35, lab, **({"font_style": "italic"} if ital else {}))
        lines = lines_by.get(q, [])
        if not lines:
            T(x0, centre[q] + fs * 0.35, empty.get(q, empty_default), fill="#808080")
            continue
        for ln, yy in zip(lines, ys):
            off = x0 if ln["align"] in ("left", "full") else x1 - ln["span"] * scale
            X = lambda v, off=off: off + v * scale
            segs, cur, half = [], 0.0, a.break_kb * 250
            for bx in ln["breaks"]:
                segs.append((cur, bx - half)); cur = bx + half
            segs.append((cur, ln["span"]))
            for s0, s1 in segs:
                out.append(f'<line x1="{X(s0):.2f}" x2="{X(s1):.2f}" y1="{yy:.2f}" y2="{yy:.2f}" '
                           f'stroke="#000" stroke-width="0.5"/>')
            for bx in ln["breaks"]:
                for d in (-1.0, 1.0):
                    cx = X(bx) + d
                    out.append(f'<line x1="{cx - 1:.2f}" x2="{cx + 1:.2f}" y1="{yy + 2.2:.2f}" '
                               f'y2="{yy - 2.2:.2f}" stroke="#000" stroke-width="0.5"/>')
            label_slots = []
            for f in sorted(ln["feats"], key=lambda f: f["x0"]):
                pts, (gx0, gx1) = arrow(X(f["x0"]), X(f["x1"]), yy, gene_h, f["ostrand"], min_w)
                tip = html.escape(f.get("hover") or f["symbol"])
                if f["kind"] == "boundary":
                    col = BOUNDARY_FILL[bnames.index(f["bname"])]
                    out.append(f'<polygon points="{pts}" fill="{col}"><title>{tip}</title></polygon>')
                else:
                    col = cmap.get(f["group"], UNCHAR_FILL)
                    edge = col if f["family"] else "#555"
                    dash = ' stroke-dasharray="1.2,0.8"' if f.get("dashed") else ""
                    fill = col if f["filled"] else "#fff"
                    sw = 0.7 if (not f["filled"] or f.get("dashed")) else 0.35
                    out.append(f'<polygon points="{pts}" fill="{fill}" stroke="{edge}" stroke-width="{sw}"{dash}>'
                               f'<title>{tip}</title></polygon>')
                    if a.labels:
                        label_slots.append(((gx0 + gx1) / 2, f["label"]))
                lay.append([q, ln["contig"], ln["align"], f["kind"], f["symbol"], f["group"], f["start"],
                            f["end"], f["strand"], f["ostrand"], f"{gx0:.2f}", f"{gx1:.2f}"])
            if label_slots:
                tiers = S.place_labels([(cx, text_w(t, fs)) for cx, t in label_slots], max_tiers=2)
                for (cx, t), k in zip(label_slots, tiers):
                    sym = bool(re.match(r"^[A-Z][A-Z0-9-]*\+?$", t)) and t != "GCGL"
                    T(cx, yy - gene_h / 2 - 1.5 - k * fs * 1.05, t, text_anchor="middle",
                      **({"font_style": "italic"} if sym else {}))

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

        def draw_n(n):
            for c in n.children:
                out.append(f'<path d="M{TX(n.x):.2f},{c.y:.2f} H{TX(c.x):.2f}" stroke="#000" '
                           f'stroke-width="0.5" fill="none"/>')
                draw_n(c)
            if n.children:
                ys_ = [c.y for c in n.children]
                out.append(f'<line x1="{TX(n.x):.2f}" x2="{TX(n.x):.2f}" y1="{min(ys_):.2f}" '
                           f'y2="{max(ys_):.2f}" stroke="#000" stroke-width="0.5"/>')
        place(tree, 0)
        draw_n(tree)

    # scale bar + legend
    y += fs * 0.6
    sb = nice_scale(span)
    out.append(f'<line x1="{x0:.2f}" x2="{x0 + sb * scale:.2f}" y1="{y:.2f}" y2="{y:.2f}" '
               f'stroke="#000" stroke-width="0.8"/>')
    T(x0 + sb * scale + 3, y + fs * 0.35, fmt_kb(sb))
    # (text, fill, edge, dashed, italic)
    items = [(legend_name(g), cmap[g], cmap[g], False, False) for g in famg]
    items += [(g, cmap[g], "#555", False, bool(re.match(r"^[A-Z][A-Z0-9-]*\+?$", g))) for g in other]
    if any(f["kind"] == "gene" and f["group"] == UNCHAR for ls in lines_by.values() for ln in ls for f in ln["feats"]):
        items.append((UNCHAR, UNCHAR_FILL, "#555", False, False))
    drawn = [f for ls in lines_by.values() for ln in ls for f in ln["feats"] if f["kind"] == "gene"]
    for t, filled, dashed in legend_status:
        if (dashed and any(f.get("dashed") for f in drawn)) or \
           (not dashed and not filled and any(not f["filled"] for f in drawn)) or (filled and not dashed):
            items.append((t, "#888" if filled else "#fff", "#000", dashed, False))
    items += [(bnames[0], BOUNDARY_FILL[0], BOUNDARY_FILL[0], False, True),
              (bnames[1], BOUNDARY_FILL[1], BOUNDARY_FILL[1], False, True)]
    lx = x0 + sb * scale + 3 + text_w(fmt_kb(sb), fs) + 12
    ly = y
    for text, fill, edge, dashed, italic in items:
        w = 12 + text_w(text, fs) + 9
        if lx + w > x1 + 2:
            lx, ly = margin + tree_w + 3, ly + fs * 1.5
        dash = ' stroke-dasharray="1.2,0.8"' if dashed else ""
        out.append(f'<rect x="{lx:.2f}" y="{ly - gene_h / 2:.2f}" width="9" height="{gene_h:.2f}" '
                   f'fill="{fill}" stroke="{edge}" stroke-width="0.6"{dash}/>')
        T(lx + 12, ly + fs * 0.35, text, **({"font_style": "italic"} if italic else {}))
        lx += w
    H = ly + fs * 0.9 + margin

    svg = (f'<svg xmlns="http://www.w3.org/2000/svg" width="{W:.1f}pt" height="{H:.1f}pt" '
           f'viewBox="0 0 {W:.2f} {H:.2f}" font-family="{FONT}">\n'
           f'<rect width="100%" height="100%" fill="#fff"/>\n' + "\n".join(out) + "\n</svg>\n")
    write(svg, a)
    if a.layout:
        with open(a.layout, "w") as fh:
            fh.write("species\tcontig\tline\tkind\tsymbol\tgroup\tstart\tend\tstrand\tstrand_drawn\tx0_pt\tx1_pt\n")
            for r in lay:
                fh.write("\t".join(map(str, r)) + "\n")
    N.info(f"{len(order)} species, {sum(len(v) for v in lines_by.values())} lines, axis {fmt_kb(span)}, "
           f"{W / PT_PER_MM:.0f} x {H / PT_PER_MM:.0f} mm, text {fs:g} pt; {n_drop} family gene(s) not drawn")


def write(svg, a):
    if a.svg:
        open(a.svg, "w").write(svg)
    elif not (a.pdf or a.png):
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
        if r.returncode == 0:
            N.info(f"wrote {path}")
        else:
            N.warn(f"{path}: {r.stderr.strip()}")
