#!/usr/bin/env python3
"""
Draw the step-4 gene tree with the same colour per name group as the synteny plots, and
(with --sources) a code after each tip saying where its name came from.

This runs lib_draw_tree.py (the old repo's 2c_draw_tree.py, copied unchanged), after
replacing its fixed colour table with lib_exendin_names.FIXED_COLOURS, so exendin-2 is
the same colour in the tree and in all three synteny figures. All of 2c's options work:
  --tree --groups --outgroup --title --min-support --width --row --png --png-scale
Node labels are SH-aLRT/UFBoot; nodes with UFBoot >= --min-support are dark. Bold =
Heloderma.

--sources tree/name_sources.tsv (from 4g --reverse-table) adds, after each tip name, the
code of the annotation its name was taken from (its best reverse hit in 4d), and a
legend line under the tree:
  R  reference RefSeq proteome (e.g. HelSus), named by NCBI
  N  NCBI annotation of another species (step 2)
  T  TOGA2 projection, named after the HelSus EviAnn reference gene (step 3)
  P  published precursor (4_reference_peptides.fa)
  O  outgroup sequence added with 4g --outgroup-fasta
Hovering over a code in the SVG shows the species and the protein it came from.

Usage
  python3 4h_draw_tree.py --tree tree/exendin_aa.rooted.nwk --groups tree/groups.tsv \
      --sources tree/name_sources.tsv --outgroup "$(cat tree/outgroup.txt)" \
      --title "..." --png tree/exendin_aa.png > tree/exendin_aa.svg
"""
import html, importlib.util, io, os, re, sys
from contextlib import redirect_stdout

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import lib_exendin_names as N  # noqa: E402

spec = importlib.util.spec_from_file_location("draw_tree", os.path.join(HERE, "lib_draw_tree.py"))
T = importlib.util.module_from_spec(spec)
spec.loader.exec_module(T)
T.FIXED = dict(N.FIXED_COLOURS, outgroup="#333333")
T.PALETTE = N.EXTRA_COLOURS

LEGEND = {"R": "reference RefSeq proteome", "N": "NCBI annotation", "T": "TOGA2 projection",
          "P": "published precursor", "O": "outgroup sequence", "?": "unknown"}


def pop_arg(name):
    """Remove '--name VALUE' (or --name=VALUE) from sys.argv and return VALUE (None if absent)."""
    for i, x in enumerate(sys.argv):
        if x == name and i + 1 < len(sys.argv):
            v = sys.argv[i + 1]; del sys.argv[i:i + 2]; return v
        if x.startswith(name + "="):
            del sys.argv[i]; return x.split("=", 1)[1]
    return None


def add_sources(svg, src):
    """Append [code] after each tip label, widen the figure for it, and add a legend line."""
    used = set()
    for tip, (code, detail) in src.items():
        pat = re.compile(r'(<text [^>]*>)' + re.escape(html.escape(tip)) + r'(</text>)')
        new, n = pat.subn(lambda m: f'{m.group(1)}{html.escape(tip)}<tspan dx="4" fill="#777" font-weight="normal" '
                                    f'font-size="9">[{code}]<title>{html.escape(detail)}</title></tspan>'
                                    f'{m.group(2)}', svg, count=1)
        if n:
            svg = new; used.add(code)
    if not used:
        N.warn("--sources: no tip of the tree is in the sources table")
        return svg, 0, 0
    first = svg.split("\n", 1)[0]
    W = float(re.search(r'width="([\d.]+)"', first).group(1))
    H = float(re.search(r'height="([\d.]+)"', first).group(1))
    dw, dh = 40, 22
    head = re.sub(r'width="[\d.]+"', f'width="{W + dw:.0f}"', first, count=1)
    head = re.sub(r'height="[\d.]+"', f'height="{H + dh:.0f}"', head, count=1)
    legend = "Name source:\u00a0\u00a0" + "\u00a0\u00a0·\u00a0\u00a0".join(f"[{c}] {LEGEND[c]}" for c in "RNTPO?" if c in used)
    body = svg.split("\n", 1)[1].rsplit("</svg>", 1)[0]
    svg = (head + "\n" + body + f'<text x="20" y="{H + dh - 8:.0f}" font-size="10" fill="#555">'
           f'{html.escape(legend)}</text>\n</svg>')
    return svg, dw, dh


def main():
    sources = pop_arg("--sources")
    if not sources:
        T.main()
        return
    src = {}
    with open(sources) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 2 and f[0] != "tip":
                detail = f"{f[2] if len(f) > 2 else ''} {f[3] if len(f) > 3 else ''}: " \
                         f"{f[4] if len(f) > 4 else ''} {f[5] if len(f) > 5 else ''}"
                src[f[0]] = (f[1], detail.strip())
    png_calls = []
    orig_png = T.write_png
    T.write_png = lambda svg, path, w, h: png_calls.append((path, w, h))
    buf = io.StringIO()
    try:
        with redirect_stdout(buf):
            T.main()
    finally:
        T.write_png = orig_png
    first = buf.getvalue().split("\n", 1)[0]
    W0 = float(re.search(r'width="([\d.]+)"', first).group(1))
    svg, dw, dh = add_sources(buf.getvalue().rstrip("\n"), src)
    print(svg)
    for path, w, h in png_calls:
        scale = w / W0
        orig_png(svg, path, int(w + dw * scale), int(h + dh * scale))
    counts = {}
    for code, _ in src.values():
        counts[code] = counts.get(code, 0) + 1
    print("[info] name sources: " + ", ".join(f"[{c}] {LEGEND.get(c, c)} {n}" for c, n in sorted(counts.items())),
          file=sys.stderr)


if __name__ == "__main__":
    main()
