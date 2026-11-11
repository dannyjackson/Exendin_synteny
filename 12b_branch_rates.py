#!/usr/bin/env python3
"""
Map the FAST partition onto the SLOW topology: compare the branch lengths two partitions
give for the same tree, branch by branch, and test whether a named set of branches is
faster than the rest (12_partitioned_tree.md).

Input: two Newick trees with the SAME topology and different branch lengths -- the slow
partition's ML tree, and the fast partition's branch lengths estimated on that fixed
topology (IQ-TREE -te). Branches are matched by the set of tips below them, so tip order
and rooting do not matter.

Per branch it reports the slow length, the fast length, the ratio, and log2 of the ratio
relative to the tree's median, which is the branch's rate shift. With --highlight it tests
Irwin & Prentice's (2011) prediction directly: that the mature peptide, but not the signal
and pro-regions, evolved rapidly on the Heloderma lineage. A permutation test compares the
median log-ratio of the highlighted branches with the rest.

Outputs (--out-prefix): _branches.tsv, and an SVG/PDF/PNG of the slow tree with each branch
coloured by its rate shift (red = faster in the mature peptide, blue = slower).

Usage
  python3 12b_branch_rates.py --slow-tree partitioned/slow.treefile \
      --fast-tree partitioned/fast_on_slow.treefile \
      --highlight HelSus_A,HelSus_B,HelCha_A,HelCha_B --highlight-name "Heloderma venom copies" \
      --out-prefix partitioned/rates
"""
import argparse, html, os, random, re, statistics, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import lib_pub_synteny as P  # noqa: E402  (write(), text_w(), FONT, PT_PER_MM)


class Node:
    __slots__ = ("name", "length", "children", "parent", "clade", "x", "y")

    def __init__(self):
        self.name, self.length, self.children, self.parent = "", 0.0, [], None


def parse(s):
    s = re.sub(r"\[[^\]]*\]", "", s.strip()).rstrip(";")
    root = Node(); cur = root; i = 0
    while i < len(s):
        c = s[i]
        if c == "(":
            n = Node(); n.parent = cur; cur.children.append(n); cur = n; i += 1
        elif c == ",":
            n = Node(); n.parent = cur.parent; cur.parent.children.append(n); cur = n; i += 1
        elif c == ")":
            cur = cur.parent; i += 1
        elif c == ":":
            m = re.match(r":([-\d.eE+]+)", s[i:]); cur.length = float(m.group(1)); i += len(m.group(0))
        else:
            m = re.match(r"('[^']*'|[^,():;]+)", s[i:]); cur.name = m.group(1).strip("'"); i += len(m.group(0))
    return root


def walk(n):
    yield n
    for c in n.children:
        yield from walk(c)


def set_clades(root):
    def rec(n):
        n.clade = frozenset([n.name]) if not n.children else frozenset().union(*(rec(c) for c in n.children))
        return n.clade
    rec(root)


def branches(root):
    """{clade: length} for every branch except the root's own."""
    set_clades(root)
    return {n.clade: n.length for n in walk(root) if n is not root}


def colour(z, lo, hi):
    """diverging blue -> grey -> red over [lo, hi], z is log2 rate shift."""
    if hi == lo:
        t = 0.5
    else:
        t = (z - lo) / (hi - lo)
    t = min(1.0, max(0.0, t))
    a, b, c = (49, 84, 162), (190, 190, 190), (178, 34, 34)
    if t < 0.5:
        f = t * 2; rgb = [round(x + (y - x) * f) for x, y in zip(a, b)]
    else:
        f = (t - 0.5) * 2; rgb = [round(x + (y - x) * f) for x, y in zip(b, c)]
    return "#%02X%02X%02X" % tuple(rgb)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--slow-tree", required=True)
    ap.add_argument("--fast-tree", required=True)
    ap.add_argument("--highlight", default="", help="comma-separated tips; their terminal branches and any "
                                                    "internal branch inside their clade are the test set")
    ap.add_argument("--highlight-name", default="highlighted")
    ap.add_argument("--min-slow", type=float, default=1e-6,
                    help="branches shorter than this in the slow partition give meaningless ratios and are "
                         "reported but left out of the statistics")
    ap.add_argument("--permutations", type=int, default=10000)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--width-mm", type=float, default=180.0)
    ap.add_argument("--font-pt", type=float, default=7.0)
    ap.add_argument("--out-prefix", required=True)
    a = ap.parse_args()

    slow_t, fast_t = parse(open(a.slow_tree).read()), parse(open(a.fast_tree).read())
    sb, fb = branches(slow_t), branches(fast_t)
    tips_s = {n.name for n in walk(slow_t) if not n.children}
    tips_f = {n.name for n in walk(fast_t) if not n.children}
    if tips_s != tips_f:
        only = ", ".join(sorted(tips_s ^ tips_f)[:8])
        sys.exit(f"[error] the trees hold different tips (e.g. {only}). Both must come from the same alignment")
    if set(sb) != set(fb):
        sys.exit(f"[error] the trees have different topologies ({len(set(sb) ^ set(fb))} branches differ). "
                 f"Estimate the fast partition on the fixed slow topology (IQ-TREE -te slow.treefile)")

    hl = {x for x in a.highlight.split(",") if x}
    missing = hl - tips_s
    if missing:
        sys.exit(f"[error] --highlight names that are not tips: {', '.join(sorted(missing))}")

    rows, zs = [], {}
    usable = []
    for clade in sb:
        s, f = sb[clade], fb[clade]
        ok = s >= a.min_slow and f > 0
        ratio = (f / s) if ok else None
        if ok:
            usable.append((clade, ratio))
        rows.append([clade, s, f, ratio])
    if not usable:
        sys.exit("[error] no branch has usable lengths in both partitions")
    med = statistics.median(r for _, r in usable)
    import math
    for clade, ratio in usable:
        zs[clade] = math.log2(ratio / med)

    def is_hl(clade):
        return bool(hl) and clade <= hl

    out = []
    for clade, s, f, ratio in rows:
        tip = next(iter(clade)) if len(clade) == 1 else ""
        out.append([tip or f"clade of {len(clade)}", "terminal" if len(clade) == 1 else "internal",
                    f"{s:.6f}", f"{f:.6f}", f"{ratio:.3f}" if ratio else "-",
                    f"{zs[clade]:+.2f}" if clade in zs else "-",
                    "yes" if is_hl(clade) else "no",
                    ",".join(sorted(clade)) if len(clade) > 1 else tip])
    out.sort(key=lambda r: (-float(r[5]) if r[5] != "-" else 99))
    with open(a.out_prefix + "_branches.tsv", "w") as fh:
        fh.write("branch\ttype\tslow_length\tfast_length\tratio\tlog2_rate_shift\thighlighted\ttips\n")
        for r in out:
            fh.write("\t".join(r) + "\n")

    print(f"[info] {len(rows)} branches; median fast/slow ratio {med:.2f} "
          f"({len(rows) - len(usable)} branch(es) too short in the slow partition to rate)", file=sys.stderr)
    print(f"{'branch':<34}{'slow':>10}{'fast':>10}{'ratio':>9}{'log2 shift':>12}  set")
    for r in out[:12]:
        print(f"{r[0][:33]:<34}{r[2]:>10}{r[3]:>10}{r[4]:>9}{r[5]:>12}  {'*' if r[6] == 'yes' else ''}")
    if len(out) > 12:
        print(f"... {len(out) - 12} more in {a.out_prefix}_branches.tsv")

    if hl:
        A = [zs[c] for c in zs if is_hl(c)]
        B = [zs[c] for c in zs if not is_hl(c)]
        if not A or not B:
            print("[warn] --highlight selects every branch or none: no test run", file=sys.stderr)
        else:
            obs = statistics.median(A) - statistics.median(B)
            pool = A + B
            rnd = random.Random(a.seed)
            hits = 0
            for _ in range(a.permutations):
                rnd.shuffle(pool)
                d = statistics.median(pool[:len(A)]) - statistics.median(pool[len(A):])
                hits += abs(d) >= abs(obs)
            p = (hits + 1) / (a.permutations + 1)
            print(f"\n{a.highlight_name}: {len(A)} branch(es), median log2 shift {statistics.median(A):+.2f}")
            print(f"rest of the tree:  {len(B)} branch(es), median log2 shift {statistics.median(B):+.2f}")
            print(f"difference {obs:+.2f} log2 units ({2 ** obs:.2f}x), permutation p = {p:.4f} "
                  f"({a.permutations} permutations)")
            if obs > 0 and p < 0.05:
                print(f"=> the mature peptide evolves faster, relative to signal+pro, on the "
                      f"{a.highlight_name} than elsewhere")

    # ---- figure: slow tree, branches coloured by rate shift ----
    fs = a.font_pt
    W = a.width_mm * P.PT_PER_MM
    leaves = [n for n in walk(slow_t) if not n.children]
    row = fs * 1.3
    name_w = max(P.text_w(l.name, fs) for l in leaves) + 6
    margin = 3.0
    plot_w = W - name_w - margin * 2 - 60
    depth = {}

    def setd(n, d=0.0):
        depth[id(n)] = d
        for c in n.children:
            setd(c, d + c.length)
    setd(slow_t)
    maxd = max(depth[id(l)] for l in leaves) or 1.0
    X = lambda d: margin + d / maxd * plot_w
    for i, l in enumerate(leaves):
        l.y = margin + fs * 1.6 + i * row

    def sety(n):
        if n.children:
            for c in n.children:
                sety(c)
            n.y = (n.children[0].y + n.children[-1].y) / 2
    sety(slow_t)
    vals = sorted(zs.values())
    lo, hi = (vals[int(0.05 * len(vals))], vals[int(0.95 * (len(vals) - 1))]) if vals else (-1, 1)
    lim = max(abs(lo), abs(hi)) or 1.0
    svg = []
    for n in walk(slow_t):
        for c in n.children:
            col = colour(zs.get(c.clade, 0.0), -lim, lim)
            wdt = 2.2 if is_hl(c.clade) else 1.4
            svg.append(f'<path d="M{X(depth[id(n)]):.2f},{c.y:.2f} H{X(depth[id(c)]):.2f}" stroke="{col}" '
                       f'stroke-width="{wdt}" fill="none" stroke-linecap="butt"><title>'
                       f'{html.escape(",".join(sorted(c.clade))[:120])}: log2 shift '
                       f'{zs.get(c.clade, float("nan")):+.2f}</title></path>')
        if n.children:
            svg.append(f'<line x1="{X(depth[id(n)]):.2f}" x2="{X(depth[id(n)]):.2f}" '
                       f'y1="{n.children[0].y:.2f}" y2="{n.children[-1].y:.2f}" stroke="#555" stroke-width="0.6"/>')
    for l in leaves:
        bold = ' font-weight="bold"' if l.name in hl else ""
        svg.append(f'<text x="{X(depth[id(l)]) + 3:.2f}" y="{l.y + fs * 0.35:.2f}" font-size="{fs}"{bold}>'
                   f'{html.escape(l.name)}</text>')
    y = margin + fs * 1.6 + len(leaves) * row + fs
    gx, gw = margin, 90.0
    for k in range(30):
        svg.append(f'<rect x="{gx + k * gw / 30:.2f}" y="{y - fs * 0.5:.2f}" width="{gw / 30 + 0.4:.2f}" '
                   f'height="{fs:.2f}" fill="{colour(-lim + 2 * lim * k / 29, -lim, lim)}"/>')
    svg.append(f'<text x="{gx:.2f}" y="{y - fs * 0.9:.2f}" font-size="{fs}">{-lim:+.1f}</text>')
    svg.append(f'<text x="{gx + gw:.2f}" y="{y - fs * 0.9:.2f}" font-size="{fs}" text-anchor="end">{lim:+.1f}</text>')
    svg.append(f'<text x="{gx + gw + 6:.2f}" y="{y + fs * 0.35:.2f}" font-size="{fs}">'
               f'log2 rate shift, mature vs signal+pro (median = 0)</text>')
    H = y + fs * 1.6
    doc = (f'<svg xmlns="http://www.w3.org/2000/svg" width="{W:.1f}pt" height="{H:.1f}pt" '
           f'viewBox="0 0 {W:.2f} {H:.2f}" font-family="{P.FONT}">\n<rect width="100%" height="100%" fill="#fff"/>\n'
           + "\n".join(svg) + "\n</svg>\n")

    class Out:
        svg, pdf, png, dpi = a.out_prefix + ".svg", a.out_prefix + ".pdf", a.out_prefix + ".png", 600
    P.write(doc, Out)
    print(f"[info] branch table {a.out_prefix}_branches.tsv; figure {a.out_prefix}.pdf", file=sys.stderr)


if __name__ == "__main__":
    main()
