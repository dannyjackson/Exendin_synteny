#!/usr/bin/env python3
"""
Alignment-based synteny of the HelSus SLC4A8-AQP5 interval across every genome of the
Cactus alignment, in the publication layout of steps 2 and 3. Input: the --indir of
9a_hal_lift_region.py.

Each genome is one or more lines (one per contig that carries at least --min-frac of its
aligned sequence, up to --max-lines). On each line:
  band    the syntenic blocks: sequence aligned to the HelSus interval, coloured by its
          position in HelSus (SLC4A8 end -> AQP5 end, a viridis gradient). A block whose
          colours run backwards is inverted relative to HelSus. Blocks drawn faded below the
          band are secondary alignments (duplicated / paralogous sequence).
  arrows  HelSus genes lifted through the alignment: exendins (and GCGL) in their fixed
          colours, every other gene grey, SLC4A8 and AQP5 dark grey. Hollow = less than
          half of the gene lifted; dashed = a secondary lift.
HelSus itself is drawn with its own genes and the full gradient, as the key.

Lines are oriented to match HelSus (majority strand of their blocks) and shifted so their
blocks sit under the matching part of HelSus. Stretches longer than --max-gap-kb in a
target with nothing aligned are cut and marked //.

Usage
  python3 9b_plot_hal_synteny.py --indir hal_lift --species-table ../species.tsv \
      --species-tree ${SPTREE} --svg hal_synteny.svg --pdf hal_synteny.pdf --png hal_synteny.png \
      --layout hal_synteny_layout.tsv 2> hal_synteny.log
"""
import argparse, csv, html, os, sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import lib_exendin_names as N  # noqa: E402
import lib_pub_synteny as P  # noqa: E402

S = P.S
VIRIDIS = [(0.0, (68, 1, 84)), (0.25, (59, 82, 139)), (0.5, (33, 145, 140)), (0.75, (94, 201, 98)), (1.0, (253, 231, 37))]


def viridis(x):
    x = min(1.0, max(0.0, x))
    for (x0, c0), (x1, c1) in zip(VIRIDIS, VIRIDIS[1:]):
        if x <= x1:
            f = (x - x0) / (x1 - x0) if x1 > x0 else 0
            return "#%02X%02X%02X" % tuple(round(a + (b - a) * f) for a, b in zip(c0, c1))
    return "#FDE725"


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def weighted_median(pairs):
    pairs = sorted(pairs)
    tot, acc = sum(w for _, w in pairs), 0
    for v, w in pairs:
        acc += w
        if acc >= tot / 2:
            return v
    return pairs[-1][0] if pairs else 0


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--indir", required=True, help="--outdir of 9a")
    ap.add_argument("--species-table", help="species.tsv, for species names")
    ap.add_argument("--species-tree", help="Newick species tree (order + cladogram)")
    ap.add_argument("--tree-map", help="TSV tree_tip<TAB>code")
    ap.add_argument("--no-tree", action="store_true")
    ap.add_argument("--min-frac", type=float, default=0.05, help="draw a contig carrying at least this fraction")
    ap.add_argument("--max-lines", type=int, default=3, help="at most this many contigs per genome")
    ap.add_argument("--no-secondary", action="store_true", help="leave secondary alignments out")
    ap.add_argument("--colour", choices=("family", "exendin"), default="family",
                    help="family: exendins and GCGL coloured (default); exendin: GCGL grey too")
    ap.add_argument("--pad-kb", type=float, default=5.0, help="drawn beyond the outermost block of a line")
    ap.add_argument("--max-gap-kb", type=float, default=100.0, help="unaligned target stretches longer than this are cut")
    ap.add_argument("--break-kb", type=float, default=10.0)
    ap.add_argument("--overhang", type=float, default=0.3,
                    help="lines may extend this fraction of the HelSus region beyond it; the rest is clipped")
    ap.add_argument("--width-mm", type=float, default=180.0)
    ap.add_argument("--font-pt", type=float, default=7.0)
    ap.add_argument("--full-names", action="store_true")
    ap.add_argument("--svg"); ap.add_argument("--pdf"); ap.add_argument("--png")
    ap.add_argument("--dpi", type=int, default=600)
    ap.add_argument("--layout", help="TSV of every drawn element")
    a = ap.parse_args()

    reg = read_tsv(os.path.join(a.indir, "region.tsv"))[0]
    lo, hi = int(reg["start"]), int(reg["end"])
    L = hi - lo
    left, right = reg["left"], reg["right"]
    refcode = reg["ref"] if len(reg["ref"]) == 6 else "".join(p[:3].capitalize() for p in reg["ref"].split("_")[:2])
    rgenes = read_tsv(os.path.join(a.indir, "ref_genes.tsv"))
    gene_info = {g["gene_id"]: g for g in rgenes}
    blocks = read_tsv(os.path.join(a.indir, "blocks.tsv"))
    glift = read_tsv(os.path.join(a.indir, "genes_lifted.tsv"))
    targets = read_tsv(os.path.join(a.indir, "targets.tsv"))
    names = {r["code"]: r["species"] for r in N.read_species_table(a.species_table)} if a.species_table else {}

    def colour_of(g):
        fam = g["family"] == "1" and (a.colour == "family" or g["group"].startswith("exendin"))
        if g["kind"] == "boundary":
            return P.BOUNDARY_FILL[0], True
        return (N.colour_map([g["group"]])[g["group"]], True) if fam else (P.OTHER_FILL, False)

    # ---- lines per genome, in HelSus-aligned coordinates ----
    lines_by = {}
    lines_by[refcode] = [dict(contig=reg["chrom"], flip=False, offset=0.0, breaks=[], L=lo, R=hi,
                              blocks=[dict(x0=lo, x1=hi, c0=0.0, c1=1.0, sec=False)],
                              genes=[dict(x0=int(g["start"]), x1=int(g["end"]), strand=g["strand"], g=g,
                                          frac=1.0, sec=False) for g in rgenes])]
    for t in targets:
        code = t["code"]
        tb = [b for b in blocks if b["code"] == code and (b["type"] == "primary" or not a.no_secondary)]
        prim = [b for b in tb if b["type"] == "primary"]
        tot = sum(int(b["aligned_bp"]) for b in prim)
        if not prim:
            lines_by[code] = []
            continue
        bp_ctg = defaultdict(int)
        for b in prim:
            bp_ctg[b["contig"]] += int(b["aligned_bp"])
        ctgs = [c for c in sorted(bp_ctg, key=bp_ctg.get, reverse=True) if bp_ctg[c] >= a.min_frac * tot][:a.max_lines]
        dropped = [c for c in bp_ctg if c not in ctgs]
        if dropped:
            N.info(f"{code}: {len(dropped)} contig(s) with {sum(bp_ctg[c] for c in dropped) / 1000:.1f} kb aligned not drawn")
        out = []
        for c in ctgs:
            cb = [b for b in tb if b["contig"] == c]
            cp = [b for b in cb if b["type"] == "primary"]
            minus = sum(int(b["aligned_bp"]) for b in cp if b["strand"] == "-")
            flip = minus * 2 > sum(int(b["aligned_bp"]) for b in cp)
            o = (lambda v: -v) if flip else (lambda v: v)
            cg = [g for g in glift if g["code"] == code and g["contig"] == c]
            items = []
            for b in cb:
                x0, x1 = sorted((o(int(b["tstart"])), o(int(b["tend"]))))
                same = (b["strand"] == "+") != flip
                rs, re_ = (int(b["rstart"]) - lo) / L, (int(b["rend"]) - lo) / L
                items.append(dict(kind="block", x0=x0, x1=x1, c0=rs if same else re_, c1=re_ if same else rs,
                                  sec=b["type"] == "secondary", rmid=(int(b["rstart"]) + int(b["rend"])) / 2,
                                  w=int(b["aligned_bp"]), raw=b))
            wlo = min(i["x0"] for i in items) - a.pad_kb * 1000
            whi = max(i["x1"] for i in items) + a.pad_kb * 1000
            for g in cg:
                x0, x1 = sorted((o(int(g["start"])), o(int(g["end"]))))
                if x1 < wlo or x0 > whi or g["gene_id"] not in gene_info:
                    continue
                st = g["strand"]
                gs = gene_info[g["gene_id"]]["strand"]
                ostrand = gs if (st == "+") != flip else P.flipstrand(gs)
                items.append(dict(kind="gene", x0=x0, x1=x1, strand=ostrand, g=gene_info[g["gene_id"]],
                                  frac=float(g["frac_lifted"]), sec=g["hit"] == "secondary"))
            # cut long unaligned stretches
            items.sort(key=lambda i: i["x0"])
            shift, end, breaks = 0.0, None, []
            gap, brk = a.max_gap_kb * 1000, a.break_kb * 1000
            for i in items:
                if end is not None and i["x0"] - shift - end > gap:
                    g_ = i["x0"] - shift - end
                    breaks.append(end + brk / 2); shift += g_ - brk
                i["x0"] -= shift; i["x1"] -= shift
                end = i["x1"] if end is None else max(end, i["x1"])
            wlo2, whi2 = min(i["x0"] for i in items) - a.pad_kb * 1000, end + a.pad_kb * 1000
            pb = [i for i in items if i["kind"] == "block" and not i["sec"]]
            off = weighted_median([(i["rmid"] - (i["x0"] + i["x1"]) / 2, i["w"]) for i in pb])
            for i in items:
                i["x0"] += off; i["x1"] += off
            out.append(dict(contig=c, flip=flip, offset=off, breaks=[b + off for b in breaks],
                            L=wlo2 + off, R=whi2 + off,
                            blocks=[i for i in items if i["kind"] == "block"],
                            genes=[i for i in items if i["kind"] == "gene"]))
        lines_by[code] = out

    # ---- order ----
    queries = [refcode] + [t["code"] for t in targets]
    tree = None
    if a.species_tree:
        order, tree = S.species_order(queries, refcode, a.species_tree, a.tree_map, None)
    else:
        order = queries
    if a.no_tree:
        tree = None

    # ---- geometry ----
    fs = a.font_pt
    W = a.width_mm * P.PT_PER_MM
    margin = 2.0
    tree_w = 0 if tree is None else 55.0
    labels = {q: P.species_label(q, names.get(q), a.full_names) for q in order}
    name_w = max(P.text_w(t, fs) for t, _ in labels.values()) + 10
    x0p, x1p = margin + tree_w + name_w, W - margin - 2
    allL = [ln["L"] for ls in lines_by.values() for ln in ls] + [lo]
    allR = [ln["R"] for ls in lines_by.values() for ln in ls] + [hi]
    XL = max(min(allL), lo - a.overhang * L)
    XR = min(max(allR), hi + a.overhang * L)
    sc = (x1p - x0p) / (XR - XL)
    X = lambda v: x0p + (min(max(v, XL), XR) - XL) * sc
    gene_h, band_h = fs * 0.85, fs * 0.5
    line_gap, row_gap = fs * 2.0, fs * 0.6
    y = margin + fs * 1.6
    rows_y, centre = [], {}
    for q in order:
        n = max(1, len(lines_by.get(q, [])))
        ys = [y + gene_h / 2 + k * line_gap for k in range(n)]
        rows_y.append((q, ys)); centre[q] = sum(ys) / len(ys) + band_h / 2
        y = ys[-1] + gene_h / 2 + band_h + 2 + row_gap

    out, defs, lay = [], [], []
    gid = [0]

    def T(x, yy, s, **kw):
        out.append(f'<text x="{x:.2f}" y="{yy:.2f}" font-size="{fs}"'
                   + "".join(f' {k.replace("_", "-")}="{v}"' for k, v in kw.items()) + f">{html.escape(s)}</text>")

    def grad(c0, c1):
        gid[0] += 1
        n = f"g{gid[0]}"
        defs.append(f'<linearGradient id="{n}" x1="0" x2="1" y1="0" y2="0">'
                    + "".join(f'<stop offset="{k / 4:.2f}" stop-color="{viridis(c0 + (c1 - c0) * k / 4)}"/>' for k in range(5))
                    + "</linearGradient>")
        return f"url(#{n})"

    T(X(lo), margin + fs, left, font_style="italic")
    T(X(hi), margin + fs, right, font_style="italic", text_anchor="end")

    for q, ys in rows_y:
        lab, ital = labels[q]
        T(margin + tree_w + 3, centre[q] + fs * 0.35, lab,
          **({"font_style": "italic"} if ital else {}), **({"font_weight": "bold"} if q == refcode else {}))
        lines = lines_by.get(q, [])
        if not lines:
            T(x0p, centre[q] + fs * 0.35, "not aligned" if q in lines_by else "not in the alignment", fill="#808080")
            continue
        for ln, yy in zip(lines, ys):
            by = yy + gene_h / 2 + 1.0
            out.append(f'<line x1="{X(ln["L"]):.2f}" x2="{X(ln["R"]):.2f}" y1="{yy:.2f}" y2="{yy:.2f}" '
                       f'stroke="#000" stroke-width="0.4"/>')
            for bx in ln["breaks"]:
                for d in (-1.0, 1.0):
                    cx = X(bx) + d
                    out.append(f'<line x1="{cx - 1:.2f}" x2="{cx + 1:.2f}" y1="{yy + 2:.2f}" y2="{yy - 2:.2f}" '
                               f'stroke="#000" stroke-width="0.5"/>')
            for b in sorted(ln["blocks"], key=lambda b: b["sec"]):
                bx0, bx1 = X(b["x0"]), X(b["x1"])
                if bx1 - bx0 < 0.3:
                    bx1 = bx0 + 0.3
                yb = by + (band_h + 0.6 if b["sec"] else 0)
                op = ' fill-opacity="0.45"' if b["sec"] else ""
                out.append(f'<rect x="{bx0:.2f}" y="{yb:.2f}" width="{bx1 - bx0:.2f}" height="{band_h:.2f}" '
                           f'fill="{grad(b["c0"], b["c1"])}"{op}/>')
                rb = b.get("raw", {})
                lay.append([q, ln["contig"], "block_secondary" if b["sec"] else "block", rb.get("tstart", lo),
                            rb.get("tend", hi), rb.get("strand", "+"), rb.get("rstart", lo), rb.get("rend", hi),
                            f"{bx0:.2f}", f"{bx1:.2f}"])
            for g in sorted(ln["genes"], key=lambda g: g["g"]["family"]):
                fill, fam = colour_of(g["g"])
                pts, (gx0, gx1) = P.arrow(X(g["x0"]), X(g["x1"]), yy, gene_h, g["strand"], 1.6)
                hollow = g["frac"] < 0.5
                edge = fill if (fam or g["g"]["kind"] == "boundary") else "#777"
                dash = ' stroke-dasharray="1.2,0.8"' if g["sec"] else ""
                tip = html.escape(f'{g["g"]["symbol"] or g["g"]["gene_id"]} | {g["g"]["description"]} | '
                                  f'{100 * g["frac"]:.0f}% lifted' + (" | secondary" if g["sec"] else ""))
                out.append(f'<polygon points="{pts}" fill="{"#fff" if hollow else fill}" stroke="{edge}" '
                           f'stroke-width="{0.7 if hollow or g["sec"] else 0.35}"{dash}><title>{tip}</title></polygon>')
                lay.append([q, ln["contig"], "gene", g["g"]["gene_id"], g["g"]["symbol"], g["strand"],
                            f'{g["frac"]:.2f}', "secondary" if g["sec"] else "primary", f"{gx0:.2f}", f"{gx1:.2f}"])

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
                out.append(f'<path d="M{TX(n.x):.2f},{c.y:.2f} H{TX(c.x):.2f}" stroke="#000" stroke-width="0.5" fill="none"/>')
                draw_n(c)
            if n.children:
                ys_ = [c.y for c in n.children]
                out.append(f'<line x1="{TX(n.x):.2f}" x2="{TX(n.x):.2f}" y1="{min(ys_):.2f}" y2="{max(ys_):.2f}" '
                           f'stroke="#000" stroke-width="0.5"/>')
        place(tree, 0)
        draw_n(tree)

    # scale bar, gradient key, legend
    y += fs * 0.6
    sb = P.nice_scale(XR - XL)
    out.append(f'<line x1="{x0p:.2f}" x2="{x0p + sb * sc:.2f}" y1="{y:.2f}" y2="{y:.2f}" stroke="#000" stroke-width="0.8"/>')
    T(x0p + sb * sc + 3, y + fs * 0.35, P.fmt_kb(sb))
    kx = x0p + sb * sc + 3 + P.text_w(P.fmt_kb(sb), fs) + 14
    T(kx, y + fs * 0.35, left, font_style="italic")
    kx += P.text_w(left, fs) + 3
    out.append(f'<rect x="{kx:.2f}" y="{y - band_h / 2 - 0.5:.2f}" width="60" height="{band_h + 1:.2f}" fill="{grad(0, 1)}"/>')
    kx += 63
    T(kx, y + fs * 0.35, right, font_style="italic")
    kx += P.text_w(right, fs) + 4
    T(kx, y + fs * 0.35, "position in")
    kx += P.text_w("position in ", fs)
    T(kx, y + fs * 0.35, labels[refcode][0], **({"font_style": "italic"} if labels[refcode][1] else {}))
    lx, ly = kx + P.text_w(labels[refcode][0], fs) + 14, y
    groups = sorted({g["group"] for g in rgenes if g["family"] == "1" and g["kind"] != "boundary"
                     and (a.colour == "family" or g["group"].startswith("exendin"))},
                    key=lambda g: (not g.startswith("GCGL"), not g.startswith("exendin-"), g))
    cm = N.colour_map(groups)
    items = [(P.legend_name(g), cm[g], cm[g], "", False) for g in groups]
    items += [(P.OTHER, P.OTHER_FILL, "#777", "", False), (f"{left} / {right}", P.BOUNDARY_FILL[0], P.BOUNDARY_FILL[0], "", True),
              ("<50% of gene aligned", "#fff", "#000", "", False)]
    if not a.no_secondary:
        items.append(("secondary alignment", "#fff", "#000", ' stroke-dasharray="1.2,0.8"', False))
    for text, fill, edge, dash, ital in items:
        w = 12 + P.text_w(text, fs) + 9
        if lx + w > x1p + 2:
            lx, ly = margin + tree_w + 3, ly + fs * 1.5
        out.append(f'<rect x="{lx:.2f}" y="{ly - gene_h / 2:.2f}" width="9" height="{gene_h:.2f}" fill="{fill}" '
                   f'stroke="{edge}" stroke-width="0.6"{dash}/>')
        T(lx + 12, ly + fs * 0.35, text, **({"font_style": "italic"} if ital else {}))
        lx += w
    H = ly + fs * 0.9 + margin

    svg = (f'<svg xmlns="http://www.w3.org/2000/svg" width="{W:.1f}pt" height="{H:.1f}pt" viewBox="0 0 {W:.2f} {H:.2f}" '
           f'font-family="{P.FONT}">\n<defs>{"".join(defs)}</defs>\n<rect width="100%" height="100%" fill="#fff"/>\n'
           + "\n".join(out) + "\n</svg>\n")
    P.write(svg, a)
    if a.layout:
        with open(a.layout, "w") as fh:
            fh.write("species\tcontig\telement\tf1\tf2\tf3\tf4\tf5\tx0_pt\tx1_pt\n")
            for r in lay:
                fh.write("\t".join(map(str, r)) + "\n")
    clipped = [q for q, ls in lines_by.items() for ln in ls if ln["L"] < XL or ln["R"] > XR]
    if clipped:
        N.warn(f"clipped at the figure edge (more than --overhang beyond the HelSus region): {', '.join(sorted(set(clipped)))}")
    N.info(f"{len(order)} genomes, {sum(len(v) for v in lines_by.values())} lines, "
           f"{W / P.PT_PER_MM:.0f} x {H / P.PT_PER_MM:.0f} mm, text {fs:g} pt")


if __name__ == "__main__":
    main()
