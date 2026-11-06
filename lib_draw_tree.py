#!/usr/bin/env python3
"""
Draw a rooted gene tree (IQ-TREE .treefile) as an SVG, with tips coloured by
gene group and node support shown. No dependencies beyond Python 3.

Tip groups follow the 2b naming:
  <Query>_ANC                 -> ANC (ancestral exendin gene)
  <Query>_A038 / _S118 / _U1  -> exendin copy, AQP5 side / SLC4A8 side / unplaced
  tips given with --outgroup  -> outgroup
(--groups TSV, tip<TAB>group, overrides this.) Heloderma tips are drawn in bold.

Support labels "SH-aLRT/UFBoot" (or a single number) are printed at internal
nodes; values with UFBoot >= --min-support are dark, the rest light grey.

Usage
  python3 2c_draw_tree_svg.py --tree exendin_origin.treefile --outgroup TaeGut_ANC \
      --title "..." > exendin_origin.svg
"""
import argparse, html, os, re, sys

PALETTE = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9",
           "#7F3C8D", "#8C564B", "#B8860B", "#2F4F4F", "#E7298A", "#666666"]
FIXED = {"ANC": "#0072B2", "exendin copy, AQP5 side": "#D55E00",
         "exendin copy, SLC4A8 side": "#CC79A7", "exendin copy, unplaced": "#E69F00", "outgroup": "#333333"}


class Node:
    __slots__ = ("name", "length", "children", "x", "y")

    def __init__(self):
        self.name, self.length, self.children, self.x, self.y = "", 0.0, [], 0.0, 0.0


def parse(s):
    s = re.sub(r"\[[^\]]*\]", "", s.strip()).rstrip(";")
    i, root, stack = 0, Node(), []
    cur = root
    while i < len(s):
        c = s[i]
        if c == "(":
            n = Node(); cur.children.append(n); stack.append(cur); cur = n; i += 1
        elif c == ",":
            n = Node(); stack[-1].children.append(n); cur = n; i += 1
        elif c == ")":
            cur = stack.pop(); i += 1
        elif c == ":":
            m = re.match(r":([^,():;]+)", s[i:]); cur.length = float(m.group(1)); i += len(m.group(0))
        else:
            m = re.match(r"('[^']*'|[^,():;]+)", s[i:]); cur.name = m.group(1).strip("'"); i += len(m.group(0))
    return root.children[0] if len(root.children) == 1 and not root.name else root


def leaves(n):
    return [n] if not n.children else [l for c in n.children for l in leaves(c)]


def infer_group(tip):
    parts = tip.split("_")
    t = parts[1] if len(parts) > 1 else ""
    if t.startswith("ANC"):
        return "ANC"
    side = {"A": "AQP5 side", "S": "SLC4A8 side", "U": "unplaced"}
    if re.match(r"^[ASU]\d", t):
        return "exendin copy, " + side[t[0]]
    return "other"


def support(name):
    """'85.3/97' -> (label, ufboot); '97' -> ('97', 97)."""
    if not name:
        return None, None
    nums = re.findall(r"[\d.]+", name)
    if not nums:
        return None, None
    try:
        return name, float(nums[-1])
    except ValueError:
        return None, None

def write_png(svg, path, w, h):
    """Render the SVG with whatever converter is installed."""
    import shutil, subprocess, tempfile
    with tempfile.NamedTemporaryFile("w", suffix=".svg", delete=False) as fh:
        fh.write(svg)
        tmp = fh.name
    cmds = []
    if shutil.which("rsvg-convert"):
        cmds.append(["rsvg-convert", "-w", str(w), "-h", str(h), "-o", path, tmp])
    if shutil.which("inkscape"):
        cmds.append(["inkscape", tmp, "--export-type=png", f"--export-filename={path}",
                     f"--export-width={w}"])
    if shutil.which("convert"):
        cmds.append(["convert", "-density", "150", "-background", "white", tmp, path])
    for c in cmds:
        try:
            if subprocess.run(c, capture_output=True).returncode == 0 and os.path.exists(path):
                print(f"[info] wrote {path} with {c[0]}", file=sys.stderr)
                os.unlink(tmp)
                return
        except Exception:
            pass
    try:
        import cairosvg
        cairosvg.svg2png(bytestring=svg.encode(), write_to=path, output_width=w)
        print(f"[info] wrote {path} with cairosvg", file=sys.stderr)
        os.unlink(tmp)
        return
    except ImportError:
        pass
    except Exception as e:
        print(f"[warn] cairosvg failed: {e}", file=sys.stderr)
    for b in ("chromium", "chromium-browser", "google-chrome"):
        if shutil.which(b):
            r = subprocess.run([b, "--headless", "--no-sandbox", "--disable-gpu",
                                f"--screenshot={path}", f"--window-size={w},{h}",
                                "file://" + tmp], capture_output=True)
            if os.path.exists(path):
                print(f"[info] wrote {path} with {b}", file=sys.stderr)
                os.unlink(tmp)
                return
    os.unlink(tmp)
    print(f"[warn] no SVG converter found; PNG not written. Install one with:\n"
          f"       mamba install -c conda-forge cairosvg\n"
          f"       (or: module load ImageMagick / inkscape)", file=sys.stderr)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tree", required=True)
    ap.add_argument("--groups")
    ap.add_argument("--outgroup", default="", help="comma-separated outgroup tip names")
    ap.add_argument("--title", default="")
    ap.add_argument("--min-support", type=float, default=95)
    ap.add_argument("--width", type=int, default=1000)
    ap.add_argument("--row", type=int, default=15, help="pixels per tip")
    ap.add_argument("--png", help="also write a PNG here (needs rsvg-convert, inkscape, "
                                  "cairosvg, ImageMagick or chromium)")
    ap.add_argument("--png-scale", type=float, default=2.0, help="PNG resolution multiplier")
    a = ap.parse_args()

    root = parse(open(a.tree).read())
    root.length = 0.0
    groups = {}
    if a.groups:
        for line in open(a.groups):
            f = line.rstrip("\n").split("\t")
            if len(f) >= 2 and f[0] != "tip":
                groups[f[0]] = f[1]

    tips = leaves(root)
    for o in filter(None, a.outgroup.split(",")):
        groups.setdefault(o, "outgroup")
    for i, t in enumerate(tips):
        t.y = i
    def setx(n, x):
        n.x = x + (n.length or 0.0)
        for c in n.children:
            setx(c, n.x)
    setx(root, 0.0)
    def sety(n):
        if n.children:
            for c in n.children:
                sety(c)
            n.y = (n.children[0].y + n.children[-1].y) / 2
    sety(root)

    maxx = max(t.x for t in tips) or 1.0
    label_w = max(len(t.name) for t in tips) * 7 + 20
    left, top = 20, 50 + (20 if a.title else 0)
    plot_w = a.width - left - label_w - 20
    X = lambda v: left + v / maxx * plot_w
    Y = lambda v: top + v * a.row

    gnames = sorted({groups.get(t.name, infer_group(t.name)) for t in tips})
    colours, k = {}, 0
    for g in gnames:
        if g in FIXED:
            colours[g] = FIXED[g]
        else:
            while PALETTE[k % len(PALETTE)] in FIXED.values():
                k += 1
            colours[g] = PALETTE[k % len(PALETTE)]; k += 1

    H = top + len(tips) * a.row + 60
    out = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{a.width}" height="{H}" '
           f'font-family="Helvetica, Arial, sans-serif" font-size="11">',
           f'<rect width="100%" height="100%" fill="white"/>']
    if a.title:
        out.append(f'<text x="{left}" y="22" font-size="15" font-weight="bold">{html.escape(a.title)}</text>')

    def draw(n):
        for c in n.children:
            out.append(f'<line x1="{X(n.x):.1f}" y1="{Y(c.y):.1f}" x2="{X(c.x):.1f}" y2="{Y(c.y):.1f}" stroke="#333" stroke-width="1.2"/>')
            draw(c)
        if n.children:
            out.append(f'<line x1="{X(n.x):.1f}" y1="{Y(n.children[0].y):.1f}" x2="{X(n.x):.1f}" '
                       f'y2="{Y(n.children[-1].y):.1f}" stroke="#333" stroke-width="1.2"/>')
            lab, uf = support(n.name)
            if lab is not None and n is not root:
                col = "#222" if uf >= a.min_support else "#AAA"
                out.append(f'<text x="{X(n.x) - 3:.1f}" y="{Y(n.y) - 3:.1f}" text-anchor="end" font-size="8" '
                           f'fill="{col}">{html.escape(lab)}</text>')
        else:
            g = groups.get(n.name, infer_group(n.name))
            bold = ' font-weight="bold"' if re.match(r"^Hel(Sus|Cha)_", n.name) else ""
            out.append(f'<circle cx="{X(n.x):.1f}" cy="{Y(n.y):.1f}" r="3" fill="{colours[g]}"/>')
            out.append(f'<text x="{X(n.x) + 6:.1f}" y="{Y(n.y) + 4:.1f}" fill="{colours[g]}"{bold}>'
                       f'{html.escape(n.name)}</text>')
    draw(root)

    # legend and scale bar
    ly = top - 22
    lx = left
    for g in gnames:
        out.append(f'<circle cx="{lx + 4}" cy="{ly}" r="4" fill="{colours[g]}"/>'
                   f'<text x="{lx + 12}" y="{ly + 4}">{html.escape(g)}</text>')
        lx += 20 + 7 * len(g) + 10
    sb = 10 ** int(f"{maxx / 5:e}".split("e")[1])
    sy = top + len(tips) * a.row + 25
    out.append(f'<line x1="{left}" y1="{sy}" x2="{X(sb):.1f}" y2="{sy}" stroke="#333" stroke-width="1.5"/>'
               f'<text x="{left}" y="{sy + 14}">{sb:g} substitutions/site</text>')
    out.append(f'<text x="{a.width - 20}" y="{sy + 14}" text-anchor="end" fill="#777">'
               f'support: SH-aLRT/UFBoot; dark = UFBoot ≥ {a.min_support:g}. Bold = Heloderma</text>')
    out.append("</svg>")
    svg = "\n".join(out)
    print(svg)
    if a.png:
        write_png(svg, a.png, int(a.width * a.png_scale), int(H * a.png_scale))

if __name__ == "__main__":
    main()
