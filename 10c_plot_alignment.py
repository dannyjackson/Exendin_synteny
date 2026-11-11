#!/usr/bin/env python3
"""
Draw a protein alignment for curation: every residue as a coloured cell (Clustal X classes),
each row labelled with its tip name and TOGA2 facts, so broken models, extensions, frame
problems and misplaced copies can be seen at a glance (10_toga2_exendin_tree.md).

  rows       reference sequences (alignment rows not in --tips, e.g. the published
             precursors added with mafft --add) first, in grey italics; then the tips,
             grouped by TOGA2 label (Ex2, Ex4, ...), in species-tree order within a group
  columns    tip | status | kind (o = ortholog, p = paralogous projection) | aa length |
             the alignment; X (internal stop or unknown) in red; a ruler every 10 columns
  --region   shade alignment columns, e.g. the mature peptide: --region 48-86

Usage
  python3 10c_plot_alignment.py --aln aln/exendin.aa.fa --tips collect/candidates.tsv \
      --species-tree ${SPTREE} --svg aln/exendin.aa.svg --pdf aln/exendin.aa.pdf
"""
import argparse, csv, html, os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import lib_pub_synteny as P  # noqa: E402

S = P.S
CLASSES = [("AILMFWV", "#80A0F0"), ("KR", "#F01505"), ("DE", "#C048C0"), ("NQST", "#15C015"),
           ("C", "#F08080"), ("G", "#F09048"), ("P", "#C0C000"), ("HY", "#15A4A4")]
COL = {aa: c for s, c in CLASSES for aa in s}


def read_fa(path):
    out, name = [], None
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith(">"):
            name = line[1:].split()[0]; out.append([name, []])
        elif out:
            out[-1][1].append(line.strip())
    return [(n, "".join(s).upper()) for n, s in out]


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--aln", required=True, help="aligned protein FASTA")
    ap.add_argument("--tips", help="candidates.tsv from 10a (status, kind, group per tip)")
    ap.add_argument("--curation-log", help="_log.tsv from 10b, to map renamed tips back")
    ap.add_argument("--species-tree", help="order species within each group")
    ap.add_argument("--region", action="append", default=[], help="alignment columns to shade, e.g. 48-86")
    ap.add_argument("--cell-pt", type=float, default=5.2, help="width of one residue")
    ap.add_argument("--font-pt", type=float, default=6.0)
    ap.add_argument("--svg"); ap.add_argument("--pdf"); ap.add_argument("--png")
    ap.add_argument("--dpi", type=int, default=300)
    a = ap.parse_args()

    aln = read_fa(a.aln)
    if len({len(s) for _, s in aln}) != 1:
        sys.exit("[error] sequences differ in length: not an alignment")
    ncol = len(aln[0][1])
    info = {}
    if a.tips:
        with open(a.tips) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                info[r["tip"]] = r
    if a.curation_log:
        with open(a.curation_log) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r["new_name"] not in ("-", r["tip"]) and r["tip"] in info:
                    info[r["new_name"]] = info[r["tip"]]
    rank = {}
    if a.species_tree:
        tips = {r["query"] for r in info.values()}
        order, _ = S.species_order(sorted(tips), "HelSus", a.species_tree, None, None)
        rank = {q: i for i, q in enumerate(order)}

    def key(row):
        n = row[0]
        r = info.get(n)
        if r is None:
            return (0, "", 0, n)
        return (1, r["group"], rank.get(r["query"], 999), n)
    aln.sort(key=key)

    fs, cw = a.font_pt, a.cell_pt
    rh = fs * 1.25
    margin = 4.0
    name_w = max(P.text_w(n, fs) for n, _ in aln) + 6
    col_w = [P.text_w("status", fs) + 5, P.text_w("kind", fs) + 5, P.text_w("aa", fs) + 12]
    x0 = margin + name_w + sum(col_w)
    W = x0 + ncol * cw + margin
    top = margin + fs * 2.6
    H = top + len(aln) * rh + margin + fs
    out = []
    mono = 'font-family="Menlo, Consolas, DejaVu Sans Mono, monospace"'

    def T(x, y, s, **kw):
        out.append(f'<text x="{x:.2f}" y="{y:.2f}" font-size="{fs}"'
                   + "".join(f' {k.replace("_", "-")}="{v}"' for k, v in kw.items()) + f">{html.escape(s)}</text>")

    for reg in a.region:
        s, _, e = reg.partition("-")
        out.append(f'<rect x="{x0 + (int(s) - 1) * cw:.2f}" y="{top - fs * 0.4:.2f}" width="{(int(e) - int(s) + 1) * cw:.2f}" '
                   f'height="{len(aln) * rh + fs * 0.4:.2f}" fill="#FFF3B0"/>')
    for c in range(10, ncol + 1, 10):
        T(x0 + (c - 0.5) * cw, top - fs * 0.6, str(c), text_anchor="middle")
    hx = margin + name_w
    for lab, w in zip(("status", "kind", "aa"), col_w):
        T(hx + 1, top - fs * 0.6, lab, fill="#555"); hx += w

    prev_group = None
    for i, (n, s) in enumerate(aln):
        y = top + i * rh
        r = info.get(n)
        g = r["group"] if r else "reference"
        if prev_group is not None and g != prev_group:
            out.append(f'<line x1="{margin:.2f}" x2="{W - margin:.2f}" y1="{y - 0.4:.2f}" y2="{y - 0.4:.2f}" '
                       f'stroke="#000" stroke-width="0.4"/>')
        prev_group = g
        if r is None:
            T(margin, y + rh * 0.72, n.split("|")[0], fill="#777", font_style="italic")
        else:
            T(margin, y + rh * 0.72, n)
            hx = margin + name_w
            bad = r["status"] not in ("FI", "I")
            T(hx + 1, y + rh * 0.72, r["status"], **({"fill": "#C00000"} if bad else {})); hx += col_w[0]
            T(hx + 1, y + rh * 0.72, "p" if r["kind"] == "paralog" else "o"); hx += col_w[1]
            aal = sum(ch not in "-." for ch in s)          # length as aligned (after any trim)
            T(hx + 1, y + rh * 0.72, str(aal), **({"fill": "#C00000"} if not 80 <= aal <= 120 else {}))
        for j, ch in enumerate(s):
            if ch in "-.":
                continue
            cx = x0 + j * cw
            fill = "#E8312F" if ch in "X*" else COL.get(ch, "#DDDDDD")
            out.append(f'<rect x="{cx:.2f}" y="{y:.2f}" width="{cw:.2f}" height="{rh - 0.6:.2f}" fill="{fill}" '
                       f'fill-opacity="0.55"/>')
            out.append(f'<text x="{cx + cw / 2:.2f}" y="{y + rh * 0.7:.2f}" font-size="{fs * 0.92:.2f}" {mono} '
                       f'text-anchor="middle">{ch}</text>')
    svg = (f'<svg xmlns="http://www.w3.org/2000/svg" width="{W:.1f}pt" height="{H:.1f}pt" viewBox="0 0 {W:.2f} {H:.2f}" '
           f'font-family="{P.FONT}">\n<rect width="100%" height="100%" fill="#fff"/>\n' + "\n".join(out) + "\n</svg>\n")
    P.write(svg, a)
    print(f"[info] {len(aln)} rows x {ncol} columns; {sum(1 for n, _ in aln if n not in info)} reference row(s); "
          f"{W / P.PT_PER_MM:.0f} x {H / P.PT_PER_MM:.0f} mm", file=sys.stderr)


if __name__ == "__main__":
    main()
