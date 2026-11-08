#!/usr/bin/env python3
"""
Figure of where the published Heloderma venom peptides are encoded in one genome (from 6a).

  A   Genomic context. Top: the whole contig (chromosome) with the region marked. Below:
      the region, with every NCBI gene in it (exendin-family genes in their name colour,
      genes carrying a mapped peptide outlined) and SLC4A8 / AQP5 if they are within
      --context-max-kb. A marker under each mapped peptide is joined to its pop-out.
  B-  One pop-out per published exendin (exendin-1, -2, -4 ... in the order of the peptide
      file). Each has
        - the locus: NCBI gene models (thin = exon, thick = CDS) and, for each peptide, the
          genomic blocks its residues map to. Minus-strand loci are drawn reversed so the
          peptide reads left to right in every panel;
        - the alignment: peptide (top) over the genome's translation of the codons it maps
          to (bottom), one block per exon, the intron length between blocks; in grey, the
          genome residues just before residue 1 (processing site) and just after the last
          residue;
        - identity, introns, the overlapping NCBI gene and exon, and the runner-up locus.
      When a pop-out has more than one peptide (helospectin I and II), the first (longest)
      is the reference and the others are aligned to it: residues encoded by the SAME
      codons in both are drawn in the exendin colour, residues that map to different codons,
      or exist in only one peptide, in magenta, in the locus blocks and in the alignment.

Usage
  python3 6b_plot_published_map.py --map-dir map --out published_map
Writes <out>.svg and <out>.png.
"""
import argparse, csv, os, sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import lib_exendin_names as N  # noqa: E402
import lib_gff as G            # noqa: E402

DIFF = "#CC79A7"          # maps differently / only in one peptide
GREY = "#9A9A9A"
FONT = 'font-family="Helvetica, Arial, sans-serif"'
MONO = 'font-family="DejaVu Sans Mono, Menlo, Consolas, monospace"'
E = N.html_escape


def read_tsv(p):
    with open(p) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def tint(hexc, f):
    """Mix a colour with white: f=0 white, f=1 the colour."""
    h = hexc.lstrip("#")
    r, g, b = (int(h[i:i + 2], 16) for i in (0, 2, 4))
    return "#%02X%02X%02X" % tuple(round(255 - (255 - c) * f) for c in (r, g, b))


def colour(group):
    return N.FIXED_COLOURS.get(group, "#555555")


def windows(ms, gap=50000):
    """Mappings grouped into loci: same contig, within gap bp."""
    wins = []
    for m in sorted(ms, key=lambda m: (m["contig"], int(m["start"]))):
        if wins and wins[-1]["contig"] == m["contig"] and int(m["start"]) - wins[-1]["hi"] <= gap:
            wins[-1]["hi"] = max(wins[-1]["hi"], int(m["end"])); wins[-1]["ms"].append(m)
        else:
            wins.append(dict(contig=m["contig"], lo=int(m["start"]), hi=int(m["end"]), ms=[m]))
    return wins


def ranges(qs, seq):
    """[3,4,5,9] -> 'S3–A5, P9' (residue letters from seq)."""
    out, qs = [], sorted(qs)
    i = 0
    while i < len(qs):
        j = i
        while j + 1 < len(qs) and qs[j + 1] == qs[j] + 1:
            j += 1
        a, b = qs[i], qs[j]
        out.append(f"{seq[a - 1]}{a}" if a == b else f"{seq[a - 1]}{a}–{seq[b - 1]}{b}")
        i = j + 1
    return ", ".join(out)


def fmt(n):
    return f"{int(n):,}"


class Svg:
    def __init__(self):
        self.el = []

    def add(self, s):
        self.el.append(s)

    def text(self, x, y, s, size=11, anchor="start", weight="normal", fill="#222", mono=False, italic=False):
        ital = ' font-style="italic"' if italic else ""
        self.add(f'<text x="{x:.1f}" y="{y:.1f}" {MONO if mono else FONT} font-size="{size}" text-anchor="{anchor}" '
                 f'font-weight="{weight}" fill="{fill}"{ital}>{E(s)}</text>')

    def rect(self, x, y, w, h, fill, stroke="none", sw=1, rx=0, opacity=1.0, dash=None):
        d = f' stroke-dasharray="{dash}"' if dash else ""
        self.add(f'<rect x="{x:.1f}" y="{y:.1f}" width="{max(w, 0):.1f}" height="{h:.1f}" rx="{rx}" fill="{fill}" '
                 f'stroke="{stroke}" stroke-width="{sw}" opacity="{opacity}"{d}/>')

    def line(self, x0, y0, x1, y1, stroke="#444", sw=1, dash=None):
        d = f' stroke-dasharray="{dash}"' if dash else ""
        self.add(f'<line x1="{x0:.1f}" y1="{y0:.1f}" x2="{x1:.1f}" y2="{y1:.1f}" stroke="{stroke}" '
                 f'stroke-width="{sw}"{d}/>')

    def poly(self, pts, fill, stroke="none", sw=1, opacity=1.0, dash=None):
        d = f' stroke-dasharray="{dash}"' if dash else ""
        self.add(f'<polygon points="{" ".join(f"{x:.1f},{y:.1f}" for x, y in pts)}" fill="{fill}" stroke="{stroke}" '
                 f'stroke-width="{sw}" opacity="{opacity}"{d}/>')


def text_w(s, size):
    return 0.56 * size * len(s)


_CLIP = [0]


def clip_open(svg, x0, x1, y0, y1):
    """Start a group clipped to x0-x1 / y0-y1 (gene models longer than the window stop at its edge)."""
    _CLIP[0] += 1
    cid = f"clip{_CLIP[0]}"
    svg.add(f'<defs><clipPath id="{cid}"><rect x="{x0:.1f}" y="{y0:.1f}" width="{x1 - x0:.1f}" '
            f'height="{y1 - y0:.1f}"/></clipPath></defs><g clip-path="url(#{cid})">')


def clip_close(svg):
    svg.add("</g>")


def draw_gene(svg, g, X, y, flip, colour_fill, outline=None, label=True, label_y=None, max_label=None, xlim=None):
    """One gene model: intron line, exons (thin), CDS (thick), strand chevrons. xlim=(x0, x1):
    the drawing is clipped there, and the label is centred on the visible part."""
    t = g["tx"][0]
    if xlim:
        clip_open(svg, xlim[0], xlim[1], y - 12, y + 12)
    xs = sorted([X(g["start"]), X(g["end"])])
    svg.line(xs[0], y, xs[1], y, stroke=colour_fill, sw=1.2)
    for s, e in t["exons"]:
        a, b = sorted([X(s), X(e + 1)])
        svg.rect(a, y - 4, max(b - a, 1.2), 8, colour_fill)
    for s, e in t["cds"]:
        a, b = sorted([X(s), X(e + 1)])
        svg.rect(a, y - 7, max(b - a, 1.2), 14, colour_fill)
    if outline:
        svg.rect(xs[0] - 3, y - 10, xs[1] - xs[0] + 6, 20, "none", stroke=outline, sw=1.6, rx=3)
    # chevrons pointing in the drawn direction of transcription
    right = (g["strand"] == "+") != flip
    span = xs[1] - xs[0]
    if span > 30:
        step = max(span / 6, 14)
        x = xs[0] + step / 2
        while x < xs[1] - 4:
            d = 3 if right else -3
            svg.add(f'<polyline points="{x - d:.1f},{y - 3} {x + d:.1f},{y} {x - d:.1f},{y + 3}" fill="none" '
                    f'stroke="white" stroke-width="1"/>')
            x += step
    if xlim:
        clip_close(svg)
        xs = [max(xs[0], xlim[0]), min(xs[1], xlim[1])]
    if label:
        s = G.gene_label(g, 40) if max_label is None else G.gene_label(g, max_label)
        if xlim and (X(g["start"]) < xlim[0] - 1 or X(g["end"]) > xlim[1] + 1 or
                     X(g["end"]) < xlim[0] - 1 or X(g["start"]) > xlim[1] + 1):
            s += "  (continues beyond the window)"
        svg.text((xs[0] + xs[1]) / 2, label_y if label_y is not None else y + 22, s, size=10, anchor="middle",
                 fill="#333")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--map-dir", required=True, help="6a output")
    ap.add_argument("--out", required=True, help="output prefix (.svg and .png)")
    ap.add_argument("--gff", help="override the GFF recorded by 6a")
    ap.add_argument("--boundary-genes", default="SLC4A8,AQP5")
    ap.add_argument("--context-max-kb", type=float, default=3000.0,
                    help="panel A reaches out to the boundary genes if they are within this distance [3000]")
    ap.add_argument("--context-pad-kb", type=float, default=15.0)
    ap.add_argument("--locus-pad-bp", type=int, default=600)
    ap.add_argument("--width", type=int, default=1900)
    ap.add_argument("--png-width", type=int, default=3800)
    ap.add_argument("--title", default=None)
    a = ap.parse_args()

    md = a.map_dir
    for f in ("genome.txt", "map.tsv", "residues.tsv"):
        if not os.path.exists(os.path.join(md, f)):
            sys.exit(f"[error] {md}/{f} missing: run 6a_map_published.py to completion first")
    info = dict(l.rstrip("\n").split("\t", 1) for l in open(os.path.join(md, "genome.txt")) if "\t" in l)
    gff = a.gff or info.get("gff", "")
    fasta = info.get("fasta", "")
    clen = {}
    if fasta and os.path.exists(fasta + ".fai"):
        clen = {l.split("\t")[0]: int(l.split("\t")[1]) for l in open(fasta + ".fai")}
    maps = [r for r in read_tsv(os.path.join(md, "map.tsv")) if r["contig"] not in ("NA", "")]
    if not maps:
        sys.exit("[error] map.tsv has no mapped peptide")
    res = defaultdict(list)
    for r in read_tsv(os.path.join(md, "residues.tsv")):
        res[r["peptide"]].append(r)

    # ---------------- panels: one per exendin group, in file order ----------------
    groups = []
    for m in maps:
        if m["group"] not in groups:
            groups.append(m["group"])
    panel_of = {g: [m for m in maps if m["group"] == g] for g in groups}
    for g in groups:                                   # reference = longest peptide first
        panel_of[g].sort(key=lambda m: -int(m["length"]))
    letters = {g: chr(ord("B") + i) for i, g in enumerate(groups)}

    # comparison classes for multi-peptide panels, by codon: a residue is 'same' if the
    # other peptide(s) of the panel have a residue encoded by exactly the same codon,
    # 'diff' if not (it maps elsewhere, is unmapped, or exists in one peptide only)
    cls = {}
    for g in groups:
        ms = panel_of[g]
        if len(ms) == 1:
            cls[ms[0]["peptide"]] = {}
            continue
        cod = {m["peptide"]: {(r["contig"], r["codon_start"]) for r in res[m["peptide"]] if r["codon_start"]}
               for m in ms}
        for m in ms:
            pid = m["peptide"]
            others = set.intersection(*[cod[x["peptide"]] for x in ms if x["peptide"] != pid])
            cls[pid] = {int(r["qpos"]): ("same" if r["codon_start"] and (r["contig"], r["codon_start"]) in others
                                         else "diff") for r in res[pid]}

    # ---------------- genes ----------------
    contigs = sorted({m["contig"] for m in maps}, key=lambda c: min(int(m["start"]) for m in maps if m["contig"] == c))
    genes = {}
    for c in contigs:
        genes[c] = G.load_region(gff, c, 1, 10 ** 12) if gff and os.path.exists(gff) else []
        if not genes[c]:
            N.warn(f"no GFF genes for {c} ({gff or 'no GFF'})")
    bnames = [b.strip().upper() for b in a.boundary_genes.split(",") if b.strip()]

    W, M = a.width, 40
    svg = Svg()
    code = info.get("code", "")
    title = a.title or f"Published Heloderma venom peptides mapped to the {code} genome ({info.get('accession', '')})"
    svg.text(M, 34, title, size=20, weight="bold")

    # ================= panel A =================
    y = 60
    svg.text(M, y + 22, "A", size=24, weight="bold")
    svg.text(M + 30, y + 21, "Genomic context", size=15, weight="bold")
    y += 40
    X0, X1 = M + 150, W - M - 30
    markers = {}                                  # group -> list of (x, y) marker positions in A
    for c in contigs:
        ms = [m for m in maps if m["contig"] == c]
        lo = min(int(m["start"]) for m in ms)
        hi = max(int(m["end"]) for m in ms)
        bs = [g for g in genes[c] if g["name"].upper() in bnames
              and min(abs(g["start"] - hi), abs(g["end"] - lo)) <= a.context_max_kb * 1000]
        # genes that carry a mapped peptide are shown whole (if not absurdly long)
        carry = [g for g in genes[c] if any(g["end"] >= int(m["start"]) and g["start"] <= int(m["end"]) for m in ms)
                 and g["end"] - g["start"] <= a.context_max_kb * 1000]
        rlo = min([lo] + [g["start"] for g in bs + carry]) - int(a.context_pad_kb * 1000)
        rhi = max([hi] + [g["end"] for g in bs + carry]) + int(a.context_pad_kb * 1000)
        rlo = max(1, rlo)
        L = clen.get(c, max(rhi, 1))
        rhi = min(rhi, L)
        # whole contig
        yc = y + 12
        svg.text(X0 - 12, yc + 4, c, size=12, anchor="end", weight="bold")
        svg.text(X0 - 12, yc + 18, f"{L / 1e6:.1f} Mb", size=10, anchor="end", fill="#666")
        svg.rect(X0, yc - 5, X1 - X0, 10, "#E4E0E2", stroke="#9A9399", rx=5)
        step = 10 ** max(0, len(str(int(L / 1e6))) - 1) * 1e6
        step = step if L / step >= 4 else step / 2
        t = step
        while t < L:
            xt = X0 + t / L * (X1 - X0)
            svg.line(xt, yc + 5, xt, yc + 9, stroke="#888")
            svg.text(xt, yc + 20, f"{t / 1e6:.0f}", size=9, anchor="middle", fill="#777")
            t += step
        hx0, hx1 = X0 + rlo / L * (X1 - X0), X0 + rhi / L * (X1 - X0)
        if hx1 - hx0 < 4:
            hx0, hx1 = (hx0 + hx1) / 2 - 2, (hx0 + hx1) / 2 + 2
        svg.rect(hx0, yc - 9, hx1 - hx0, 18, "none", stroke="#D9643A", sw=2, rx=2)
        ya = yc + 120                             # region axis
        svg.poly([(hx0, yc + 9), (hx1, yc + 9), (X1, ya - 62), (X0, ya - 62)], fill="#F6F2F1", stroke="#C9C2C6",
                 dash="3,3")
        svg.text(X0 - 12, ya + 4, f"{rlo / 1e6:.3f}–{rhi / 1e6:.3f} Mb", size=10, anchor="end", fill="#666")

        def XA(p, rlo=rlo, rhi=rhi):
            return X0 + (p - rlo) / (rhi - rlo) * (X1 - X0)
        svg.line(X0, ya, X1, ya, stroke="#BBB", sw=1)
        # scale bar
        span = rhi - rlo
        sb = 10 ** int(len(str(int(span / 5))) - 1)
        sb = sb * 5 if span / (sb * 5) >= 3 else sb
        svg.line(X1 - (X1 - X0) * sb / span, ya + 62, X1, ya + 62, stroke="#444", sw=2)
        svg.text(X1, ya + 76, f"{sb / 1000:g} kb", size=10, anchor="end")
        # genes in region; labels in rows above the axis, no overlaps
        rows_end = [-1e9] * 4
        mapped_genes = set()
        for m in ms:
            for g in genes[c]:
                if g["end"] >= int(m["start"]) and g["start"] <= int(m["end"]):
                    mapped_genes.add(g["id"])
        for g in genes[c]:
            if g["end"] < rlo or g["start"] > rhi:
                continue
            grp = N.name_group(g["description"] + " " + g["name"])
            fam = N.is_family(grp)
            col = colour(grp) if fam else ("#3A3238" if g["name"].upper() in bnames else "#8E878C")
            outl = "#1E1A1D" if g["id"] in mapped_genes else None
            draw_gene(svg, g, XA, ya, False, col, outline=outl, label=False, xlim=(X0, X1))
            xm = (max(XA(g["start"]), X0) + min(XA(g["end"]), X1)) / 2
            lab = g["name"] if not fam else f"{g['name']} ({N.short_label(g['description']) or grp})"
            if "uncharacterized" in g["description"].lower():
                lab += " (uncharacterized)"
            w = text_w(lab, 9.5)
            for k in range(len(rows_end)):
                if xm - w / 2 > rows_end[k] + 4:
                    rows_end[k] = xm + w / 2
                    ly = ya - 16 - 13 * k
                    svg.line(xm, ya - 9, xm, ly + 3, stroke="#CCC", sw=0.6)
                    svg.text(xm, ly, lab, size=9.5, anchor="middle", fill=col if fam else "#444",
                             weight="bold" if (fam or g["id"] in mapped_genes) else "normal")
                    break
        # peptide markers below the axis, one per panel (group) and locus; markers that would
        # overprint (several exendins on one gene) are set side by side
        mk = []
        for g in groups:
            gm = [m for m in ms if m["group"] == g]
            for w in windows(gm):
                mk.append([(XA(w["lo"]) + XA(w["hi"])) / 2, g])
        mk.sort()
        clusters = []
        for x, g in mk:
            if clusters and x - clusters[-1][-1][0] < 18:
                clusters[-1].append((x, g))
            else:
                clusters.append([(x, g)])
        for cl in clusters:
            x0 = sum(x for x, _ in cl) / len(cl)
            for k, (_, g) in enumerate(cl):
                xm = x0 + (k - (len(cl) - 1) / 2) * 18
                col = colour(g)
                svg.poly([(xm, ya + 12), (xm - 7, ya + 24), (xm + 7, ya + 24)], fill=col)
                svg.text(xm, ya + 38, letters[g], size=13, anchor="middle", weight="bold", fill=col)
                markers.setdefault(g, []).append((xm, ya + 42))
        y = ya + 90
    A_bottom = y

    # ================= pop-outs =================
    n = len(groups)
    gap = 28
    pw = (W - 2 * M - (n - 1) * gap) / n
    PY = A_bottom + 70
    panel_bottoms = []
    wedges = Svg()                                     # drawn first, under panel A and the panels
    body = Svg()                                       # panels drawn after the wedges (on top)
    for i, g in enumerate(groups):
        px = M + i * (pw + gap)
        col = colour(g)
        ms = panel_of[g]
        s = body
        yy = PY + 26
        s.text(px + 12, yy, letters[g], size=22, weight="bold", fill=col)
        names = " / ".join(m["label"] for m in ms)
        s.text(px + 38, yy - 1, f"{g[0].upper() + g[1:]}: {names}", size=14, weight="bold")
        yy += 20
        # windows: peptides of this panel grouped by contig and proximity
        for w in windows(ms):
            strand = w["ms"][0]["strand"]
            flip = strand == "-"
            ov = [gg for gg in genes.get(w["contig"], []) if gg["end"] >= w["lo"] and gg["start"] <= w["hi"]]
            lo, hi = w["lo"] - a.locus_pad_bp, w["hi"] + a.locus_pad_bp
            for gg in ov:                               # show overlapping genes whole if not huge
                if gg["end"] - gg["start"] <= 40000:
                    lo, hi = min(lo, gg["start"] - 200), max(hi, gg["end"] + 200)
            lo = max(1, lo)
            near = [gg for gg in genes.get(w["contig"], []) if gg["end"] >= lo and gg["start"] <= hi]
            LX0, LX1 = px + 92, px + pw - 16

            def X(p, lo=lo, hi=hi, flip=flip, LX0=LX0, LX1=LX1):
                f = (p - lo) / (hi - lo)
                return LX0 + (1 - f if flip else f) * (LX1 - LX0)
            s.text(px + 12, yy + 4, f"{w['contig']}:{fmt(lo)}–{fmt(hi)}  ({'−' if flip else '+'} strand"
                   + (", drawn reversed)" if flip else ")"), size=10.5, fill="#555")
            yy += 22
            # axis with end coordinates
            s.line(LX0, yy, LX1, yy, stroke="#CCC")
            s.text(LX0, yy - 4, fmt(hi if flip else lo), size=8.5, fill="#888")
            s.text(LX1, yy - 4, fmt(lo if flip else hi), size=8.5, anchor="end", fill="#888")
            span = hi - lo
            sb = 10 ** int(len(str(max(1, int(span / 4)))) - 1)
            sb = sb * 5 if span / (sb * 5) >= 2.5 else sb
            yy += 22
            # gene models, one row each
            s.text(px + 12, yy + 4, "NCBI genes", size=10, fill="#555")
            for gg in near:
                grp = N.name_group(gg["description"] + " " + gg["name"])
                gcol = colour(grp) if N.is_family(grp) else "#8E878C"
                draw_gene(s, gg, X, yy, flip, gcol, label=True, label_y=yy + 20, max_label=int((LX1 - LX0) / 8),
                          xlim=(LX0, LX1))
                yy += 34
            if not near:
                s.text(LX0, yy + 4, "no NCBI gene annotated here", size=10, fill="#999", italic=True)
                yy += 22
            # peptide rows: residues coloured as runs (same / differs)
            for m in w["ms"]:
                pid = m["peptide"]
                s.text(px + 12, yy + 4, m["label"], size=10.5, weight="bold", fill=col)
                s.line(LX0, yy, LX1, yy, stroke="#EEE")
                runs = []
                for r in res[pid]:
                    if r["codon_start"] == "":
                        continue
                    c = cls.get(pid, {}).get(int(r["qpos"]), "same")
                    c0, c1 = int(r["codon_start"]), int(r["codon_end"])
                    if runs and runs[-1][2] == c and runs[-1][3] == r["exon"]:
                        runs[-1][0] = min(runs[-1][0], c0); runs[-1][1] = max(runs[-1][1], c1)
                    else:
                        runs.append([c0, c1, c, r["exon"]])
                for c0, c1, c, ex in runs:
                    a0, a1 = sorted([X(c0), X(c1 + 1)])
                    if c == "diff":
                        a0, a1 = (a0 + a1) / 2 - max(1.5, (a1 - a0) / 2), (a0 + a1) / 2 + max(1.5, (a1 - a0) / 2)
                    s.rect(a0, yy - 6, max(a1 - a0, 2), 12, col if c == "same" else DIFF)
                # exon labels: query range
                ex_rng = defaultdict(list)
                for r in res[pid]:
                    if r["exon"]:
                        ex_rng[r["exon"]].append(int(r["qpos"]))
                for ex, qs in ex_rng.items():
                    cs = [int(r["codon_start"]) for r in res[pid] if r["exon"] == ex and r["codon_start"]]
                    xm = (X(min(cs)) + X(max(cs))) / 2
                    s.text(xm, yy + 17, f"{min(qs)}–{max(qs)}", size=8.5, anchor="middle", fill="#555")
                yy += 30
            s.line(LX1 - (LX1 - LX0) * sb / span, yy - 4, LX1, yy - 4, stroke="#444", sw=1.5)
            s.text(LX1, yy + 8, f"{sb:,} bp" if sb < 1000 else f"{sb / 1000:g} kb", size=9, anchor="end")
            yy += 14
        # ---------- alignment ----------
        yy += 10
        s.text(px + 12, yy, "Peptide over the genome's translation of its codons", size=11, weight="bold")
        yy += 8
        nmax = max(int(m["length"]) for m in ms)
        nex = max(len([x for x in m["introns"].split(",") if x not in ("-", "")]) for m in ms)
        flank_n = 5
        avail = pw - 24 - 74 - nex * 50 - 16
        cw = max(7.0, min(14.0, avail / (nmax + flank_n)))
        fs = min(12.5, cw * 0.95)
        for m in ms:
            pid = m["peptide"]
            rows = res[pid]
            intr = [x for x in m["introns"].split(",") if x not in ("-", "")]
            yq, yg = yy + 16, yy + 16 + cw + 6
            s.text(px + 12, yq + 2, m["label"], size=10, weight="bold", fill=col)
            s.text(px + 12, yg + 2, code or "genome", size=10, fill="#555")
            x = px + 86
            before = "" if m["before"] in ("-", "") else m["before"]
            for ch in before:
                s.text(x + cw / 2, yg + 2, ch, size=fs, anchor="middle", fill=GREY, mono=True)
                x += cw
            if before:
                s.line(x + 1, yq - cw, x + 1, yg + 6, stroke="#BBB", dash="2,2")
                x += 3
            prev_ex, k_int, seg_start, seg = None, 0, x, []
            segs = []
            for r in rows:
                ex = r["exon"] or "u"
                if prev_ex is not None and ex != prev_ex and ex != "u" and prev_ex != "u":
                    segs.append((seg_start, x, seg)); seg = []
                    il = intr[k_int] if k_int < len(intr) else "?"
                    k_int += 1
                    s.text(x + 23, yq + 1, "//", size=12, anchor="middle", fill="#777", weight="bold")
                    s.text(x + 23, yg + 2, f"{int(il):,} bp" if il.lstrip("-").isdigit() else il, size=8,
                           anchor="middle", fill="#777")
                    x += 46; seg_start = x
                prev_ex = ex
                c = cls.get(pid, {}).get(int(r["qpos"]), "same")
                st = r["status"]
                if len(ms) > 1:
                    s.rect(x + 0.5, yq - cw + 2, cw - 1, cw + 1, tint(col, 0.35) if c == "same" else tint(DIFF, 0.6))
                s.text(x + cw / 2, yq, r["qres"], size=fs, anchor="middle", mono=True, weight="bold")
                if st == "identical":
                    s.rect(x + 0.5, yg - cw + 2, cw - 1, cw + 1, tint(col, 0.3))
                    s.text(x + cw / 2, yg, r["gres"], size=fs, anchor="middle", mono=True)
                elif st == "stop":
                    s.rect(x + 0.5, yg - cw + 2, cw - 1, cw + 1, "#C0392B")
                    s.text(x + cw / 2, yg, "*", size=fs, anchor="middle", mono=True, fill="white")
                elif st == "mismatch":
                    s.text(x + cw / 2, yg, r["gres"], size=fs, anchor="middle", mono=True, fill="#C0392B",
                           weight="bold")
                else:
                    s.text(x + cw / 2, yg, "·", size=fs, anchor="middle", mono=True, fill="#AAA")
                if int(r["ins_after"] or 0):
                    s.line(x + cw, yg - cw, x + cw, yg + 3, stroke="#C0392B", sw=1.5)
                if r["codon_start"]:
                    seg.append(int(r["codon_start"]))
                x += cw
            segs.append((seg_start, x, seg))
            after = "" if m["after"] in ("-", "") else m["after"]
            if after:
                s.line(x + 1, yq - cw, x + 1, yg + 6, stroke="#BBB", dash="2,2")
                x += 3
                for ch in after:
                    s.text(x + cw / 2, yg + 2, ch, size=fs, anchor="middle", fill=GREY, mono=True)
                    x += cw
            # coordinates under each exon block
            for x0, x1, cs in segs:
                if not cs:
                    continue
                a0, a1 = (max(cs), min(cs)) if m["strand"] == "-" else (min(cs), max(cs))
                lab = f"{fmt(a0)}–{fmt(a1)}"
                if text_w(lab, 8) <= x1 - x0 + 30:
                    s.text((x0 + x1) / 2, yg + 14, lab, size=8, anchor="middle", fill="#777")
                else:
                    s.text(x0, yg + 14, fmt(a0), size=8, fill="#777")
            yy = yg + 22
        # ---------- legend and facts ----------
        yy += 4
        lx = px + 12
        items = [(tint(col, 0.3), "identical"), ("#FFFFFF", "mismatch (red letter)"), ("#C0392B", "stop"),
                 (GREY, "genome flank (grey letters)")]
        if len(ms) > 1:
            items = [(tint(col, 0.35), f"same codons in {ms[0]['label']} and the other(s)"),
                     (tint(DIFF, 0.6), "different codons, or in one peptide only")] + items
        for fill, lab in items:
            s.rect(lx, yy - 8, 10, 10, fill, stroke="#999", sw=0.6)
            s.text(lx + 14, yy + 1, lab, size=9.5, fill="#444")
            lx += 24 + text_w(lab, 9.5)
            if lx > px + pw - 120:
                lx = px + 12; yy += 15
        yy += 20
        for m in ms:
            intr = [x for x in m["introns"].split(",") if x not in ("-", "")]
            facts = [f"{m['label']}: {m['identical']}/{m['length']} identical "
                     f"({100 * int(m['identical']) / int(m['length']):.0f}%), {m['exons']} exon(s)"
                     + (f", intron{'s' if len(intr) > 1 else ''} " + ", ".join(f"{int(x):,} bp" for x in intr)
                        if intr else ""),
                     f"   before: …{m['before']}   after: {m['after']}…",
                     f"   gene: {m['genes']}",
                     f"   runner-up locus: {m['runner_up']}"]
            for f in facts:
                lim = int((pw - 24) / (0.52 * 10))
                while f:
                    s.text(px + 12, yy, f[:lim], size=10, fill="#333" if not f.startswith("   ") else "#555")
                    f = ("     " + f[lim:]) if len(f) > lim else ""
                    yy += 13
            yy += 4
        if len(ms) > 1:
            nsame = sum(1 for v in cls[ms[0]["peptide"]].values() if v == "same")
            lines = [f"{' and '.join(m['label'] for m in ms)}: {nsame} residues are encoded by the same codons."]
            for m in ms:
                d = [int(q) for q, v in cls[m["peptide"]].items() if v == "diff"]
                if d:
                    sq = "".join(r["qres"] for r in res[m["peptide"]])
                    lines.append(f"Only in {m['label']} (or encoded elsewhere): "
                                 + ("all residues" if len(d) == len(sq) else ranges(d, sq)))
            lim = int((pw - 24) / (0.55 * 10))
            for f in lines:
                while f:
                    cut = f.rfind(" ", 0, lim) if len(f) > lim else len(f)
                    s.text(px + 12, yy, f[:cut], size=10, weight="bold", fill="#333")
                    f = f[cut:].strip(); yy += 14
        panel_bottoms.append(yy + 10)
        # wedge from marker(s) in A
        for mx, my in markers.get(g, []):
            wedges.poly([(mx - 4, my), (mx + 4, my), (px + pw, PY), (px, PY)], fill=tint(col, 0.08), stroke=tint(col, 0.55),
                     dash="4,3")
    for i, g in enumerate(groups):
        px = M + i * (pw + gap)
        svg.rect(px, PY, pw, panel_bottoms[i] - PY, "#FFFFFF", stroke=tint(colour(g), 0.7), sw=1.5, rx=8)
    svg.el = svg.el[:1] + wedges.el + svg.el[1:] + body.el
    H = int(max(panel_bottoms) + 30)
    svg.text(M, H - 10, "Mapping: 6a_map_published.py (genome-wide tblastn, then a sensitive tblastn in a window "
             "around the best loci; one codon per peptide residue). Genes: NCBI RefSeq GFF.", size=9.5, fill="#777")
    out_svg = (f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}">\n'
               f'<rect width="{W}" height="{H}" fill="white"/>\n' + "\n".join(svg.el) + "\n</svg>\n")
    open(a.out + ".svg", "w").write(out_svg)
    N.load_engine().write_png(out_svg, a.out + ".png", a.png_width)
    N.info(f"wrote {a.out}.svg")


if __name__ == "__main__":
    main()
