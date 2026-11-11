#!/usr/bin/env python3
"""
River (riparian) plot of the HelSus SLC4A8-AQP5 interval across the genomes of the Cactus
alignment, in the style of GENESPACE / SVbyEye. Input: the --indir of 9a_hal_lift_region.py.

Rows. One horizontal bar per genome, in species-tree order: the stretch of each contig that
the HelSus interval aligns to. A genome whose interval is on several contigs has several
bars side by side (assembly break or rearrangement). Each genome is oriented to match
HelSus. The HelSus genes are drawn on every bar where they lift: exendins and GCGL in their
fixed colours, every other gene grey, SLC4A8 and AQP5 dark grey.

Ribbons. Between each pair of neighbouring rows, every run of HelSus tiles that is
collinear in both genomes is one ribbon from its position in the upper genome to its
position in the lower one. Positions come from each tile's primary alignment (9a).
  --ribbon-colour direction (default): grey = same orientation in both genomes; brown =
      inverted between them (the ribbon twists); ribbons over exendin-family genes in that
      gene's colour.
  --ribbon-colour position: the GENESPACE look: every ribbon coloured by its position along
      HelSus (viridis, SLC4A8 end purple -> AQP5 end yellow); exendins still highlighted.
--secondary adds faded ribbons from each genome's primary position to the secondary
(duplicated) positions in the next genome.

Genomes with no alignment to the interval are left out of the river and named in the log.

Usage
  python3 9b_plot_hal_river.py --indir hal_lift --species-table ../species.tsv \
      --species-tree ${SPTREE} --svg hal_river.svg --pdf hal_river.pdf --png hal_river.png \
      --layout hal_river_layout.tsv 2> hal_river.log
"""
import argparse, csv, html, os, sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import lib_exendin_names as N  # noqa: E402
import lib_pub_synteny as P  # noqa: E402

S = P.S
VIRIDIS = [(0.0, (68, 1, 84)), (0.25, (59, 82, 139)), (0.5, (33, 145, 140)), (0.75, (94, 201, 98)), (1.0, (253, 231, 37))]
FWD, INV = "#BDBDBD", "#A6611A"


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


class Bar:
    """One contig's stretch in one genome: maps target coordinates to figure bp-coordinates."""

    def __init__(self, contig, flip, ivals, pad, max_gap, brk):
        self.contig, self.flip = contig, flip
        o = sorted((min(self.o(s), self.o(e)), max(self.o(s), self.o(e))) for s, e in ivals)
        self.segs, shift, end = [], 0.0, None             # (oriented start, shift)
        for s, e in o:
            if end is not None and s - end > max_gap:
                shift += (s - end) - brk
                self.segs.append((s, shift))
            end = e if end is None else max(end, e)
        self.lo, self.hi = o[0][0] - pad, end + pad
        self.breaks = []
        self.offset = 0.0
        for s, sh in self.segs:
            self.breaks.append(s - sh - brk / 2)

    def o(self, v):
        return -v if self.flip else v

    def x(self, v):
        ov = self.o(v)
        sh = 0.0
        for s, shift in self.segs:
            if ov >= s:
                sh = shift
        return ov - sh + self.offset

    @property
    def L(self):
        return self.lo + self.offset

    @property
    def R(self):
        sh = self.segs[-1][1] if self.segs else 0.0
        return self.hi - sh + self.offset


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--indir", required=True, help="--outdir of 9a")
    ap.add_argument("--species-table", help="species.tsv, for species names")
    ap.add_argument("--species-tree", help="Newick species tree (row order + cladogram)")
    ap.add_argument("--tree-map", help="TSV tree_tip<TAB>code")
    ap.add_argument("--no-tree", action="store_true")
    ap.add_argument("--ribbon-colour", choices=("direction", "position"), default="direction")
    ap.add_argument("--colour", choices=("family", "exendin"), default="family",
                    help="family: exendins and GCGL highlighted (default); exendin: only exendins")
    ap.add_argument("--secondary", action="store_true", help="also draw ribbons to secondary (duplicated) hits")
    ap.add_argument("--min-frac", type=float, default=0.05, help="draw a contig carrying at least this fraction")
    ap.add_argument("--max-bars", type=int, default=3, help="at most this many contigs per genome")
    ap.add_argument("--max-jump-kb", type=float, default=20.0, help="a ribbon breaks where either genome jumps more")
    ap.add_argument("--pad-kb", type=float, default=5.0)
    ap.add_argument("--max-gap-kb", type=float, default=100.0, help="unaligned stretches longer than this are cut (//)")
    ap.add_argument("--break-kb", type=float, default=10.0)
    ap.add_argument("--bar-gap-kb", type=float, default=20.0, help="space between two contigs of one genome")
    ap.add_argument("--width-mm", type=float, default=180.0)
    ap.add_argument("--font-pt", type=float, default=7.0)
    ap.add_argument("--full-names", action="store_true")
    ap.add_argument("--svg"); ap.add_argument("--pdf"); ap.add_argument("--png")
    ap.add_argument("--dpi", type=int, default=600)
    ap.add_argument("--layout", help="TSV of every drawn bar and ribbon")
    a = ap.parse_args()

    reg = read_tsv(os.path.join(a.indir, "region.tsv"))[0]
    lo, hi = int(reg["start"]), int(reg["end"])
    Lreg = hi - lo
    left, right = reg["left"], reg["right"]
    ref = reg["ref"]
    refcode = ref if len(ref) == 6 else "".join(p[:3].capitalize() for p in ref.split("_")[:2])
    rgenes = read_tsv(os.path.join(a.indir, "ref_genes.tsv"))
    ginfo = {g["gene_id"]: g for g in rgenes}
    hits = read_tsv(os.path.join(a.indir, "tile_hits.tsv"))
    glift = read_tsv(os.path.join(a.indir, "genes_lifted.tsv"))
    targets = read_tsv(os.path.join(a.indir, "targets.tsv"))
    names = {r["code"]: r["species"] for r in N.read_species_table(a.species_table)} if a.species_table else {}
    pad, maxgap, brk = a.pad_kb * 1000, a.max_gap_kb * 1000, a.break_kb * 1000

    def fam_group(g):
        if g["kind"] == "boundary" or g["family"] != "1":
            return None
        return g["group"] if (a.colour == "family" or g["group"].startswith("exendin")) else None

    famgenes = [(int(g["start"]), int(g["end"]), fam_group(g)) for g in rgenes if fam_group(g)]
    groups = sorted({f[2] for f in famgenes}, key=lambda g: (not g.startswith("GCGL"), not g.startswith("exendin-"), g))
    cmap = N.colour_map(groups)

    def tile_group(rs, re_):
        best = None
        for s, e, g in famgenes:
            ov = min(e, re_) - max(s, rs)
            if ov > 0 and (best is None or ov > best[0]):
                best = (ov, g)
        return best[1] if best else None

    # ---- tile positions per genome: {code: {rstart: (contig, tstart, tend, strand)}} ----
    tiles = sorted({int(h["rstart"]): int(h["rend"]) for h in hits}.items())
    tile_end = dict(tiles)
    if not tiles:
        sys.exit("[error] no tile hits in tile_hits.tsv: re-run 9a (it writes tile_hits.tsv since this version)")
    tile_len = tiles[0][1] - tiles[0][0] + 1
    pos = defaultdict(dict)
    sec = defaultdict(lambda: defaultdict(list))
    for s in range(lo, hi + 1, tile_len):
        pos[refcode][s] = (reg["chrom"], s, min(s + tile_len - 1, hi), "+")
    for h in hits:
        rec = (h["contig"], int(h["tstart"]), int(h["tend"]), h["strand"])
        if h["type"] == "primary":
            pos[h["code"]][int(h["rstart"])] = rec
        else:
            sec[h["code"]][int(h["rstart"])].append(rec)

    # ---- bars per genome ----
    bars_by, unaligned = {}, []
    codes = [refcode] + [t["code"] for t in targets]
    for code in codes:
        p = pos.get(code, {})
        if not p:
            unaligned.append(code)
            continue
        bp = defaultdict(int)
        for c, s, e, st in p.values():
            bp[c] += e - s + 1
        tot = sum(bp.values())
        ctgs = [c for c in sorted(bp, key=bp.get, reverse=True) if bp[c] >= a.min_frac * tot][:a.max_bars]
        left_out = [c for c in bp if c not in ctgs]
        if left_out:
            N.info(f"{code}: {len(left_out)} contig(s) with {sum(bp[c] for c in left_out) / 1000:.1f} kb aligned not drawn")
        bars = []
        for c in ctgs:
            recs = [(r, v) for r, v in p.items() if v[0] == c]
            minus = sum(v[2] - v[1] for _, v in recs if v[3] == "-")
            flip = minus * 2 > sum(v[2] - v[1] for _, v in recs)
            ivals = [(v[1], v[2]) for _, v in recs]
            ivals += [(int(g["start"]), int(g["end"])) for g in glift if g["code"] == code and g["contig"] == c]
            b = Bar(c, flip, ivals, pad, maxgap, brk)
            # keep only lifted genes inside the drawn window (others would stretch it)
            b.recs = recs
            b.target = weighted_median([(r + tile_len / 2 - (b.x(v[1]) + b.x(v[2])) / 2, v[2] - v[1]) for r, v in recs])
            b.rmed = weighted_median([(r, v[2] - v[1]) for r, v in recs])
            bars.append(b)
        bars.sort(key=lambda b: b.rmed)
        prev_R = None
        for b in bars:
            b.offset = b.target
            if prev_R is not None and b.L < prev_R + a.bar_gap_kb * 1000:
                b.offset += prev_R + a.bar_gap_kb * 1000 - b.L
            prev_R = b.R
        bars_by[code] = bars
    if unaligned:
        N.warn(f"not aligned to the interval, left out of the river: {', '.join(unaligned)}")

    # ---- order ----
    drawn = [c for c in codes if c in bars_by]
    tree = None
    if a.species_tree:
        order, tree = S.species_order(drawn, refcode, a.species_tree, a.tree_map, None)
    else:
        order = drawn
    if a.no_tree:
        tree = None

    # ---- geometry ----
    fs = a.font_pt
    W = a.width_mm * P.PT_PER_MM
    margin = 2.0
    tree_w = 0 if tree is None else 45.0
    labels = {q: P.species_label(q, names.get(q), a.full_names) for q in order}
    name_w = max(P.text_w(t, fs) for t, _ in labels.values()) + 8
    x0p, x1p = margin + tree_w + name_w, W - margin - 2
    XL = min(b.L for bs in bars_by.values() for b in bs)
    XR = max(b.R for bs in bars_by.values() for b in bs)
    sc = (x1p - x0p) / (XR - XL)
    X = lambda v: x0p + (v - XL) * sc
    bar_h, gene_h, row_step = fs * 0.45, fs * 0.9, fs * 3.4
    top = margin + fs * 1.8
    ys = {q: top + k * row_step for k, q in enumerate(order)}

    out_rib, out_bar, out_txt, lay = [], [], [], []

    def T(x, yy, s, **kw):
        out_txt.append(f'<text x="{x:.2f}" y="{yy:.2f}" font-size="{fs}"'
                       + "".join(f' {k.replace("_", "-")}="{v}"' for k, v in kw.items()) + f">{html.escape(s)}</text>")

    def bar_of(code, contig):
        for b in bars_by[code]:
            if b.contig == contig:
                return b
        return None

    # ribbons between neighbouring rows
    def ribbon(xa0, xa1, ya, xb0, xb1, yb, fill, opacity, twist):
        if twist:
            xb0, xb1 = xb1, xb0
        ym = (ya + yb) / 2
        d = (f"M{xa0:.2f},{ya:.2f} L{xa1:.2f},{ya:.2f} C{xa1:.2f},{ym:.2f} {xb1:.2f},{ym:.2f} {xb1:.2f},{yb:.2f} "
             f"L{xb0:.2f},{yb:.2f} C{xb0:.2f},{ym:.2f} {xa0:.2f},{ym:.2f} {xa0:.2f},{ya:.2f} Z")
        out_rib.append(f'<path d="{d}" fill="{fill}" fill-opacity="{opacity}" stroke="none"/>')

    jump = a.max_jump_kb * 1000
    n_rib = 0
    for qa, qb in zip(order, order[1:]):
        ya, yb = ys[qa] + gene_h / 2, ys[qb] - gene_h / 2
        runs, cur = [], None
        for r, _ in tiles:
            va, vb = pos[qa].get(r), pos[qb].get(r)
            ba = bar_of(qa, va[0]) if va else None
            bb = bar_of(qb, vb[0]) if vb else None
            if not (ba and bb):
                cur = None
                continue
            xa = sorted((ba.x(va[1]), ba.x(va[2])))
            xb = sorted((bb.x(vb[1]), bb.x(vb[2])))
            # orientation of this tile in each drawn bar, relative to HelSus order
            da = 1 if (va[3] == "+") != ba.flip else -1
            db = 1 if (vb[3] == "+") != bb.flip else -1
            grp = tile_group(r, tile_end[r])
            if cur and cur["ba"] is ba and cur["bb"] is bb and cur["da"] == da and cur["db"] == db \
                    and cur["grp"] == grp and abs(xa[0] - cur["xa"][-1][0]) <= jump + tile_len * 2 \
                    and abs(xb[0] - cur["xb"][-1][0]) <= jump + tile_len * 2 and r - cur["r1"] <= jump:
                cur["xa"].append(xa); cur["xb"].append(xb); cur["r1"] = r
            else:
                cur = dict(ba=ba, bb=bb, da=da, db=db, grp=grp, xa=[xa], xb=[xb], r0=r, r1=r)
                runs.append(cur)
        for run in sorted(runs, key=lambda r: r["grp"] is not None):
            a0, a1 = min(x[0] for x in run["xa"]), max(x[1] for x in run["xa"])
            b0, b1 = min(x[0] for x in run["xb"]), max(x[1] for x in run["xb"])
            twist = run["da"] != run["db"]
            if run["grp"]:
                fill, op = cmap[run["grp"]], 0.85
            elif a.ribbon_colour == "position":
                fill, op = viridis(((run["r0"] + run["r1"]) / 2 - lo) / Lreg), 0.7
            else:
                fill, op = (INV if twist else FWD), 0.6
            ribbon(X(a0), X(a1), ya, X(b0), X(b1), yb, fill, op, twist)
            n_rib += 1
            lay.append([f"{qa}->{qb}", run["ba"].contig, run["bb"].contig, run["r0"], run["r1"] + tile_len - 1,
                        "inverted" if twist else "collinear", run["grp"] or ""])
        if a.secondary:
            for r, recs in sec[qb].items():
                va = pos[qa].get(r)
                ba = bar_of(qa, va[0]) if va else None
                for vb in recs:
                    bb = bar_of(qb, vb[0])
                    if ba and bb:
                        xa = sorted((ba.x(va[1]), ba.x(va[2]))); xb = sorted((bb.x(vb[1]), bb.x(vb[2])))
                        grp = tile_group(r, tile_end[r])
                        ribbon(X(xa[0]), X(xa[1]), ya, X(xb[0]), X(xb[1]), yb, cmap.get(grp, "#999999"), 0.35, False)

    # bars, genes, labels
    for q in order:
        y = ys[q]
        lab, ital = labels[q]
        T(margin + tree_w + 3, y + fs * 0.35, lab, **({"font_style": "italic"} if ital else {}),
          **({"font_weight": "bold"} if q == refcode else {}))
        for b in bars_by[q]:
            out_bar.append(f'<rect x="{X(b.L):.2f}" y="{y - bar_h / 2:.2f}" width="{(b.R - b.L) * sc:.2f}" '
                           f'height="{bar_h:.2f}" rx="{bar_h / 2:.2f}" fill="#EFEFEF" stroke="#666" stroke-width="0.4"/>')
            for bx in b.breaks:
                cx = X(bx + b.offset)
                out_bar.append(f'<rect x="{cx - 1.6:.2f}" y="{y - bar_h:.2f}" width="3.2" height="{2 * bar_h:.2f}" fill="#fff"/>')
                for d in (-1.0, 1.0):
                    out_bar.append(f'<line x1="{cx + d - 1:.2f}" x2="{cx + d + 1:.2f}" y1="{y + 2.2:.2f}" '
                                   f'y2="{y - 2.2:.2f}" stroke="#000" stroke-width="0.5"/>')
            lay.append([q, b.contig, "", "", "", "flipped" if b.flip else "as is", f"{X(b.L):.1f}-{X(b.R):.1f}"])
            if q == refcode:
                gl = [(int(g["start"]), int(g["end"]), g, 1.0, "primary") for g in rgenes]
            else:
                gl = [(int(g["start"]), int(g["end"]), ginfo[g["gene_id"]], float(g["frac_lifted"]), g["hit"])
                      for g in glift if g["code"] == q and g["contig"] == b.contig and g["gene_id"] in ginfo]
            for s, e, g, frac, hit in sorted(gl, key=lambda t: fam_group(t[2]) is not None):
                gx0, gx1 = sorted((X(b.x(s)), X(b.x(e))))
                st = g["strand"] if (q == refcode) else g["strand"]
                if q != refcode:
                    v = [h for h in glift if h["code"] == q and h["gene_id"] == g["gene_id"] and h["contig"] == b.contig]
                    if v and v[0]["strand"] == "-":
                        st = P.flipstrand(st)
                if b.flip:
                    st = P.flipstrand(st)
                grp = fam_group(g)
                if g["kind"] == "boundary":
                    fill, edge = P.BOUNDARY_FILL[0], P.BOUNDARY_FILL[0]
                elif grp:
                    fill, edge = cmap[grp], cmap[grp]
                else:
                    fill, edge = P.OTHER_FILL, "#777"
                pts, _ = P.arrow(gx0, gx1, y, gene_h, st, 1.6)
                hollow = frac < 0.5
                dash = ' stroke-dasharray="1.2,0.8"' if hit == "secondary" else ""
                tip = html.escape(f'{g["symbol"] or g["gene_id"]} | {g["description"]} | {100 * frac:.0f}% aligned')
                out_bar.append(f'<polygon points="{pts}" fill="{"#fff" if hollow else fill}" stroke="{edge}" '
                               f'stroke-width="{0.7 if hollow or dash else 0.35}"{dash}><title>{tip}</title></polygon>')

    if tree is not None:
        def depth(n):
            return 0 if not n.children else 1 + max(depth(c) for c in n.children)
        D = depth(tree) or 1
        TX = lambda d: margin + d / D * (tree_w - 4)

        def place(n, d):
            if not n.children:
                n.x, n.y = D, ys[n.name]
            else:
                for c in n.children:
                    place(c, d + 1)
                n.x = d
                n.y = (min(c.y for c in n.children) + max(c.y for c in n.children)) / 2

        def draw_n(n):
            for c in n.children:
                out_txt.append(f'<path d="M{TX(n.x):.2f},{c.y:.2f} H{TX(c.x):.2f}" stroke="#000" stroke-width="0.5" fill="none"/>')
                draw_n(c)
            if n.children:
                yy = [c.y for c in n.children]
                out_txt.append(f'<line x1="{TX(n.x):.2f}" x2="{TX(n.x):.2f}" y1="{min(yy):.2f}" y2="{max(yy):.2f}" '
                               f'stroke="#000" stroke-width="0.5"/>')
        place(tree, 0)
        draw_n(tree)

    # headers from HelSus's boundary genes
    rb = bars_by[refcode][0]
    T(X(rb.x(lo)), margin + fs, left, font_style="italic")
    T(X(rb.x(hi)), margin + fs, right, font_style="italic", text_anchor="end")

    # scale bar + legend
    y = ys[order[-1]] + gene_h + fs * 1.4
    sb = P.nice_scale(XR - XL)
    out_txt.append(f'<line x1="{x0p:.2f}" x2="{x0p + sb * sc:.2f}" y1="{y:.2f}" y2="{y:.2f}" stroke="#000" stroke-width="0.8"/>')
    T(x0p + sb * sc + 3, y + fs * 0.35, P.fmt_kb(sb))
    lx, ly = x0p + sb * sc + 3 + P.text_w(P.fmt_kb(sb), fs) + 14, y
    items = [(P.legend_name(g), cmap[g], cmap[g], "", False) for g in groups]
    items += [(P.OTHER, P.OTHER_FILL, "#777", "", False), (f"{left} / {right}", P.BOUNDARY_FILL[0], P.BOUNDARY_FILL[0], "", True)]
    if a.ribbon_colour == "direction":
        items += [("collinear", FWD, FWD, "", False), ("inverted", INV, INV, "", False)]
    items += [("<50% of gene aligned", "#fff", "#000", "", False)]
    if a.secondary:
        items.append(("duplicated (secondary)", "#999999", "#999999", ' fill-opacity="0.35"', False))
    defs = ""
    if a.ribbon_colour == "position":
        defs = ('<defs><linearGradient id="vir" x1="0" x2="1" y1="0" y2="0">'
                + "".join(f'<stop offset="{k / 4:.2f}" stop-color="{viridis(k / 4)}"/>' for k in range(5))
                + "</linearGradient></defs>")
        T(lx, ly + fs * 0.35, left, font_style="italic"); lx += P.text_w(left, fs) + 3
        out_txt.append(f'<rect x="{lx:.2f}" y="{ly - gene_h / 2:.2f}" width="50" height="{gene_h:.2f}" fill="url(#vir)"/>')
        lx += 53
        T(lx, ly + fs * 0.35, right, font_style="italic"); lx += P.text_w(right, fs) + 14
    for text, fill, edge, extra, ital in items:
        w = 12 + P.text_w(text, fs) + 9
        if lx + w > x1p + 2:
            lx, ly = margin + tree_w + 3, ly + fs * 1.5
        out_txt.append(f'<rect x="{lx:.2f}" y="{ly - gene_h / 2:.2f}" width="9" height="{gene_h:.2f}" fill="{fill}" '
                       f'stroke="{edge}" stroke-width="0.6"{extra}/>')
        T(lx + 12, ly + fs * 0.35, text, **({"font_style": "italic"} if ital else {}))
        lx += w
    H = ly + fs * 0.9 + margin

    svg = (f'<svg xmlns="http://www.w3.org/2000/svg" width="{W:.1f}pt" height="{H:.1f}pt" viewBox="0 0 {W:.2f} {H:.2f}" '
           f'font-family="{P.FONT}">\n{defs}<rect width="100%" height="100%" fill="#fff"/>\n'
           + "\n".join(out_rib + out_bar + out_txt) + "\n</svg>\n")
    P.write(svg, a)
    if a.layout:
        with open(a.layout, "w") as fh:
            fh.write("rows\tcontig_a\tcontig_b\thelsus_start\thelsus_end\ttype\tgroup\n")
            for r in lay:
                fh.write("\t".join(map(str, r)) + "\n")
    N.info(f"{len(order)} genomes, {sum(len(v) for v in bars_by.values())} bars, {n_rib} ribbons, "
           f"{W / P.PT_PER_MM:.0f} x {H / P.PT_PER_MM:.0f} mm, text {fs:g} pt")


if __name__ == "__main__":
    main()
