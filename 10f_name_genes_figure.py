#!/usr/bin/env python3
"""
Name the exendin-family genes from the tree and their position, and draw one publication
figure: (A) the SLC4A8-AQP5 interval of every genome with each gene labelled by its new
name, (B) the gene tree with the same names, and with --aln (C) the protein alignment of
the exendin genes only (Exendin_anc and the lettered copies; no GCGL, GCG or GIP), in tree
order, columns that are gaps in every exendin removed (10_toga2_exendin_tree.md, step 5).

Naming rule
  GCGL         every tip in the GCGL clade: the smallest clade holding the --gcgl-tips
               (default: tips already named *_GCGL). Name <Sp>_GCGL.
  Exendin_anc  the ancestral anguimorph exendin: the smallest clade holding the --anc-tips
               (e.g. one Heloderma, one anguid and the Shinisaurus copy). One gene per
               genome in that clade becomes <Sp>_Exendin_anc. If a genome has several
               there, the one closest to the clade's base (shortest path from the clade's
               node) is chosen and the choice is logged; --names overrides it.
  copies       every other exendin gene of a genome: a letter by its order from that
               genome's Exendin_anc along the chromosome, counting exendin genes only:
               toward AQP5 A, B, C ...; toward SLC4A8 Z, Y, X ...  (<Sp>_A, <Sp>_Z)
               The AQP5 direction comes from the boundary genes on the same contig.
               A genome with no Exendin_anc is counted from its GCGL instead (logged).
               A copy on another contig than the anchor cannot be ordered: <Sp>_U1, U2 ...
  GCG, GIP     outgroups keep their names.
--names TSV (old_tip<TAB>new_name) overrides any of these. Every decision is written to
names.tsv with the rule that made it.

Inputs
  --tree         IQ-TREE .treefile (support labels SH-aLRT/UFBoot), rooted with -o on GIP
  --candidates   collect/candidates.tsv (10a): coordinates of every TOGA2 tip
  --curation-log curated/exendin_log.tsv (10b): renames and added sequences
  --coords       TSV tip, contig, start, end, strand: coordinates of added tips; or
  --rescued      TIP=rescued/<name>.exons.tsv from lib_rescue_locus.py
  --genes        the 3c gene table of the run the tree's genes come from (e.g.
                 interval_HelSus_genes.tsv): SLC4A8 / AQP5 coordinates per genome, and
                 the other genes drawn grey in panel A

Outputs: --outdir names.tsv, renamed.treefile; --svg/--pdf/--png the figure.

Usage
  python3 10f_name_genes_figure.py --tree tree/exendin_aa.treefile --outgroup "$(cat tree/outgroup.txt)" \
      --anc-tips HelSus_Ex2_Ex4,ShiCro_Ex2_Ex4,ElgMul_Ex2_Ex4 --gcgl-tips GalGal_GCGL,HelSus_GCGL \
      --candidates collect/candidates.tsv --curation-log curated/exendin_log.tsv \
      --rescued HelCha_Ex3=rescued/HelCha_Ex3.exons.tsv \
      --genes ../synteny_toga2/interval/interval_HelSus_genes.tsv \
      --species-table ../species_toga2.tsv --species-tree ${SPTREE} \
      --outdir names --svg exendin_figure.svg --pdf exendin_figure.pdf --png exendin_figure.png
"""
import argparse, csv, html, os, re, sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import lib_exendin_names as N  # noqa: E402
import lib_pub_synteny as P  # noqa: E402

COL = {"GIP": "#7F3C8D", "GCG": "#7F7F7F", "GCGL": "#009E73", "Exendin_anc": "#0072B2", "copy": "#D55E00"}
SLC_RX = re.compile(r"SLC4A8|solute carrier family 4 member 8", re.I)
AQP_RX = re.compile(r"AQP5|aquaporin[- ]?5", re.I)


# ---------------- tree ----------------
class Node:
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
            m = re.match(r":([^,():;]+)", s[i:]); cur.length = float(m.group(1)); i += len(m.group(0))
        else:
            m = re.match(r"('[^']*'|[^,():;]+)", s[i:]); cur.name = m.group(1).strip("'"); i += len(m.group(0))
    return root


def leaves(n):
    return [n] if not n.children else [l for c in n.children for l in leaves(c)]


def mrca(root, names):
    names = set(names)
    best = None

    def walk(n):
        nonlocal best
        here = {l.name for l in leaves(n)}
        if names <= here:
            best = n
            for c in n.children:
                walk(c)
    walk(root)
    return best


def root_on(root, og):
    """IQ-TREE -o leaves the outgroup tips as children of a multifurcating root: gather them."""
    og = set(og)
    if not og:
        return root
    oc = [c for c in root.children if {l.name for l in leaves(c)} <= og]
    ic = [c for c in root.children if c not in oc]
    if {l.name for c in oc for l in leaves(c)} != og or not ic:
        N.warn("outgroup tips are not grouped at the root of this tree: drawn as given (root with IQ-TREE -o)")
        return root
    new = Node()
    for side in (oc, ic):
        if len(side) == 1:
            k = side[0]
        else:
            k = Node(); k.children = side
            for x in side:
                x.parent = k
        k.parent = new; new.children.append(k)
    return new


def prune(n, gone):
    """Remove leaves in gone; a node left with one child is replaced by it (branch lengths added)."""
    n.children = [c for c in (prune(c, gone) for c in n.children) if c is not None]
    if not n.children:
        return None if (n.name in gone or n.name == "") else n
    if len(n.children) == 1:
        c = n.children[0]
        c.length += n.length
        c.parent = n.parent
        return c
    for c in n.children:
        c.parent = n
    return n


def support(name):
    m = re.match(r"^([\d.]+)/([\d.]+)$", name or "")
    if m:
        return float(m.group(1)), float(m.group(2))
    m = re.match(r"^([\d.]+)$", name or "")
    return (None, float(m.group(1))) if m else (None, None)


def ladderize(n):
    for c in n.children:
        ladderize(c)
    n.children.sort(key=lambda c: len(leaves(c)))


def dist_to(node, anc):
    d = 0.0
    while node is not anc and node is not None:
        d += node.length; node = node.parent
    return d


def write_newick(n, rename):
    if not n.children:
        return f"{rename.get(n.name, n.name)}:{n.length:g}"
    return "(" + ",".join(write_newick(c, rename) for c in n.children) + f"){n.name}:{n.length:g}"


# ---------------- coordinates ----------------
def norm_ctg(c):
    return re.sub(r"__(\d+)$", r".\1", c or "")


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t"))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--tree", required=True)
    ap.add_argument("--outgroup", default="", help="comma-separated GIP tips (tree/outgroup.txt)")
    ap.add_argument("--anc-tips", required=True, help="tips whose smallest clade is the ancestral anguimorph exendin")
    ap.add_argument("--gcgl-tips", default="", help="tips whose smallest clade is GCGL (default: tips named *_GCGL)")
    ap.add_argument("--candidates", required=True)
    ap.add_argument("--curation-log")
    ap.add_argument("--coords", help="TSV tip contig start end strand, for added tips")
    ap.add_argument("--rescued", action="append", default=[], help="TIP=exons.tsv from lib_rescue_locus.py")
    ap.add_argument("--genes", required=True, help="3c --genes-out table of the same TOGA2 run")
    ap.add_argument("--names", help="TSV old_tip<TAB>new_name: overrides")
    ap.add_argument("--no-context", action="store_true", help="panel A: leave out the grey non-exendin genes")
    ap.add_argument("--left-mm", type=float, default=105.0, help="width of panel A")
    ap.add_argument("--drop-species", default="",
                    help="comma-separated species codes left out of the figure: their row in A, their tips "
                         "pruned from the tree in B (not re-inferred) and their rows in C, e.g. BarImb")
    ap.add_argument("--aln", help="protein alignment of the tree (curated tip names): adds panel C, the exendin "
                                  "genes only, in tree order, under their new names")
    ap.add_argument("--aln-region", action="append", default=[],
                    help="START-END[:LABEL] in panel C's columns to shade and label, e.g. 52-90:hormone domain")
    ap.add_argument("--outdir", required=True)
    P.add_args(ap, near_default=0.0, split_default="all")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    tree = root_on(parse(open(a.tree).read()), [t for t in a.outgroup.split(",") if t])
    drop = {x for x in a.drop_species.split(",") if x}
    if drop:
        gone = [l.name for l in leaves(tree) if l.name.split("_")[0] in drop]
        tree = prune(tree, set(gone))
        N.info(f"--drop-species {','.join(sorted(drop))}: {len(gone)} tip(s) pruned from the tree "
               f"({', '.join(gone)}); the tree is not re-inferred")
    ladderize(tree)
    tips = [l.name for l in leaves(tree)]
    tipnode = {l.name: l for l in leaves(tree)}

    # ---- classes ----
    cls = {}
    for t in tips:
        base = t.split("_", 1)[1] if "_" in t else t
        if re.sub(r"[b-z]$", "", base) == "GIP":
            cls[t] = "GIP"
        elif re.sub(r"[b-z]$", "", base) == "GCG":
            cls[t] = "GCG"
    gt = [t for t in a.gcgl_tips.split(",") if t] or [t for t in tips if re.search(r"_GCGL[b-z]?$", t)]
    if gt:
        g = mrca(tree, gt)
        for l in leaves(g):
            if l.name in cls:
                N.warn(f"GCGL clade (MRCA of {','.join(gt)}) contains {l.name} ({cls[l.name]}): check --gcgl-tips")
            cls.setdefault(l.name, "GCGL")
    anc_t = [t for t in a.anc_tips.split(",") if t]
    missing = [t for t in anc_t + gt if t not in tipnode]
    if missing:
        sys.exit(f"[error] not tips of the tree: {', '.join(missing)}")
    anc_clade = mrca(tree, anc_t)
    anc_members = [l.name for l in leaves(anc_clade)]
    for t in anc_members:
        if t in cls:
            N.warn(f"Exendin_anc clade (MRCA of {','.join(anc_t)}) contains {t} ({cls[t]}): check --anc-tips")
    for t in tips:
        cls.setdefault(t, "exendin")
    N.info(f"classes: " + ", ".join(f"{k} {sum(v == k for v in cls.values())}" for k in ("GIP", "GCG", "GCGL", "exendin")))
    N.info(f"Exendin_anc clade: {len(anc_members)} tips: {', '.join(anc_members)}")

    # ---- coordinates of tree tips ----
    cand = {r["tip"]: r for r in read_tsv(a.candidates)}
    newname = {}
    if a.curation_log:
        for r in read_tsv(a.curation_log):
            if r["new_name"] != "-":
                newname[r["new_name"]] = r["tip"]
    pos = {}
    for t in tips:
        r = cand.get(newname.get(t, t))
        if r:
            pos[t] = dict(contig=r["contig"], start=int(r["start"]), end=int(r["end"]), strand=r["strand"])
    if a.coords:
        for r in read_tsv(a.coords):
            pos[r["tip"]] = dict(contig=r["contig"], start=int(r["start"]), end=int(r["end"]), strand=r["strand"])
    for spec in a.rescued:
        t, _, path = spec.partition("=")
        ex = read_tsv(path)
        ctg, tfa = None, path.replace(".exons.tsv", ".target.fa")     # lib_rescue_locus writes both
        if os.path.exists(tfa):
            with open(tfa) as fh:
                ctg = fh.readline()[1:].split()[0].split(":")[0]
        pos[t] = dict(contig=ctg or "", start=min(int(e["g_start"]) for e in ex), end=max(int(e["g_end"]) for e in ex),
                      strand=ex[0]["strand"])
        if not ctg:
            N.warn(f"{t}: contig name not found next to {path}; give it with --coords")
    for t in tips:
        if cls[t] not in ("GIP", "GCG") and t not in pos and t.split("_")[0] not in ("GalGal", "HomSap"):
            N.warn(f"{t}: no coordinates (not in candidates.tsv, --coords or --rescued): not drawn, not ordered")

    # ---- boundaries and context genes per genome ----
    grows = read_tsv(a.genes)
    bnd, context = defaultdict(dict), defaultdict(list)
    for r in grows:
        f = dict(contig=r["contig"], start=int(r["start"]), end=int(r["end"]), strand=r["strand"],
                 symbol=r["symbol"], status=r.get("status", ""), desc=r.get("description", ""), group=r["group"])
        if r["boundary"] == "boundary":
            text = " ".join([r["symbol"], r["group"], r.get("description", "")])
            name = "SLC4A8" if SLC_RX.search(text) else "AQP5" if AQP_RX.search(text) else None
            if name:
                bnd[r["species"]][name] = f
        else:
            context[r["species"]].append(f)

    # ---- names ----
    override = {}
    if a.names:
        for line in open(a.names):
            p = line.rstrip("\n").split("\t")
            if len(p) >= 2 and not line.startswith("#") and p[0] != "old_tip":
                override[p[0]] = p[1]
    rows, rename = [], {}
    by_sp = defaultdict(list)
    for t in tips:
        by_sp[t.split("_")[0]].append(t)
    for sp, ts in sorted(by_sp.items()):
        anc = [t for t in ts if t in anc_members]
        chosen = None
        if anc:
            chosen = min(anc, key=lambda t: dist_to(tipnode[t], anc_clade))
            if len(anc) > 1:
                N.warn(f"{sp}: {len(anc)} genes in the Exendin_anc clade ({', '.join(anc)}); {chosen} is closest to the "
                       f"clade's base and is called Exendin_anc. Override with --names if another is the orthologue")
        for t in ts:
            if cls[t] == "GCGL":
                rename[t] = f"{sp}_GCGL"; rows.append([t, rename[t], sp, "GCGL", "GCGL clade"])
            elif cls[t] in ("GIP", "GCG"):
                rename[t] = t; rows.append([t, t, sp, cls[t], "outgroup"])
            elif t == chosen:
                rename[t] = f"{sp}_Exendin_anc"; rows.append([t, rename[t], sp, "Exendin_anc", "Exendin_anc clade"])
        copies = [t for t in ts if cls[t] == "exendin" and t != chosen]
        if not copies:
            continue
        anchor, how = chosen, "Exendin_anc"
        if anchor is None:
            gl = [t for t in ts if cls[t] == "GCGL"]
            anchor, how = (gl[0], "GCGL (no Exendin_anc in this genome)") if gl else (None, "none")
            N.warn(f"{sp}: no Exendin_anc; copies counted from {anchor or 'nothing'}")
        if anchor is None or anchor not in pos:
            for k, t in enumerate(copies, 1):
                rename[t] = f"{sp}_U{k}"; rows.append([t, rename[t], sp, "copy", "no anchor position: unordered"])
            continue
        ap_ = pos[anchor]
        actg = norm_ctg(ap_["contig"])
        amid = (ap_["start"] + ap_["end"]) / 2
        b = bnd.get(sp, {})
        if "AQP5" in b and norm_ctg(b["AQP5"]["contig"]) == actg:
            aqp_dir = 1 if (b["AQP5"]["start"] + b["AQP5"]["end"]) / 2 > amid else -1; dhow = "AQP5 on the same contig"
        elif "SLC4A8" in b and norm_ctg(b["SLC4A8"]["contig"]) == actg:
            aqp_dir = -1 if (b["SLC4A8"]["start"] + b["SLC4A8"]["end"]) / 2 > amid else 1; dhow = "SLC4A8 on the same contig"
        else:
            aqp_dir, dhow = 1, "no boundary on the anchor's contig: AQP5 assumed toward increasing coordinate"
            N.warn(f"{sp}: {dhow}")
        side = {1: [], -1: []}
        unplaced = []
        for t in copies:
            p = pos.get(t)
            if not p or norm_ctg(p["contig"]) != actg:
                unplaced.append(t); continue
            d = (p["start"] + p["end"]) / 2 - amid
            side[1 if d * aqp_dir > 0 else -1].append((abs(d), t))
        for sgn, letters in ((1, "ABCDEFGHIJKLM"), (-1, "ZYXWVUTSRQPON")):
            for k, (d, t) in enumerate(sorted(side[sgn])):
                rename[t] = f"{sp}_{letters[k]}"
                rows.append([t, rename[t], sp, "copy", f"{k + 1} {'toward AQP5' if sgn == 1 else 'toward SLC4A8'} "
                                                       f"from {how} ({d / 1000:.1f} kb; {dhow})"])
        for k, t in enumerate(unplaced, 1):
            rename[t] = f"{sp}_U{k}"; rows.append([t, rename[t], sp, "copy", "not on the anchor's contig: unordered"])
    for old, new in override.items():
        if old in rename:
            rename[old] = new
            for r in rows:
                if r[0] == old:
                    r[1], r[4] = new, r[4] + "; overridden by --names"
    dup = [n for n in set(rename.values()) if list(rename.values()).count(n) > 1]
    if dup:
        sys.exit(f"[error] two tips would get the same name: {', '.join(dup)}")
    with open(os.path.join(a.outdir, "names.tsv"), "w") as fh:
        fh.write("old_tip\tnew_name\tspecies\tclass\trule\tcontig\tstart\tend\n")
        for r in rows:
            p = pos.get(r[0], {})
            fh.write("\t".join(map(str, r + [p.get("contig", ""), p.get("start", ""), p.get("end", "")])) + "\n")
    with open(os.path.join(a.outdir, "renamed.treefile"), "w") as fh:
        fh.write(write_newick(tree, rename) + ";\n")
    for r in rows:
        N.info(f"  {r[0]:24s} -> {r[1]:22s} {r[4]}")

    # ---- panel A: interval with names ----
    klass = {t: ("Exendin_anc" if rename.get(t, "").endswith("_Exendin_anc") else
                 "copy" if cls[t] == "exendin" else cls[t]) for t in tips}
    N.FIXED_COLOURS.update({"Exendin_anc": COL["Exendin_anc"], "exendin copy": COL["copy"],
                            "GCGL (glucagon-like)": COL["GCGL"]})
    species = {}
    queries = sorted((set(bnd) | {t.split("_")[0] for t in tips if t in pos}) - drop)
    for sp in queries:
        genes = []
        drawn = []
        for t in tips:
            if t.split("_")[0] != sp or t not in pos or klass[t] in ("GIP", "GCG"):
                continue
            p = pos[t]
            k = klass[t]
            lab = rename[t].split("_", 1)[1]
            lab = {"Exendin_anc": "anc"}.get(lab, lab)
            genes.append(dict(contig=norm_ctg(p["contig"]), start=p["start"], end=p["end"], strand=p["strand"], symbol=rename[t],
                              group={"GCGL": "GCGL (glucagon-like)", "Exendin_anc": "Exendin_anc", "copy": "exendin copy"}[k],
                              family=True, filled=True, dashed=False, label=lab, label_italic=False,
                              hover=f"{rename[t]} (was {t})"))
            drawn.append((norm_ctg(p["contig"]), p["start"], p["end"]))
        for f in context.get(sp, []):
            ov = any(norm_ctg(f["contig"]) == c and f["start"] <= e and f["end"] >= s for c, s, e in drawn)
            fam = N.is_family(f["group"]) or f["group"].startswith("exendin")
            if ov:
                continue
            if fam:
                genes.append(dict(f, contig=norm_ctg(f["contig"]), family=True, filled=False, dashed=False, group="exendin copy", label="",
                                  hover=f"{f['symbol']} {f['status']}: not in the tree"))
            elif not a.no_context:
                genes.append(dict(f, contig=norm_ctg(f["contig"]), family=False, filled=True, dashed=False, group=P.OTHER, label="", hover=f["symbol"]))
        b = bnd.get(sp, {})
        species[sp] = dict(bnd={k: dict(v, contig=norm_ctg(v["contig"]), family=False, label="", dashed=False, filled=True, hover=k) for k, v in b.items()},
                           genes=genes, empty="" if (b or genes) else "no TOGA2 annotation")
    a.labels, a.colour = True, "family"
    width_total = a.width_mm
    a.width_mm = a.left_mm
    captured = []
    orig_write = P.write
    P.write = lambda svg, args: captured.append(svg)
    try:
        P.draw(species, a, ["SLC4A8", "AQP5"], legend_status=[("not in the tree", False, False)])
    finally:
        P.write = orig_write
    left = captured[0]

    # ---- panel B: tree ----
    fs = a.font_pt
    rw = (width_total - a.left_mm) * P.PT_PER_MM
    # match panel A's height: spread the tips over it (within readable limits)
    lh0 = float(re.search(r'height="([\d.]+)pt"', left.split("\n", 1)[0]).group(1))
    row = min(max((lh0 - fs * 2.2 - fs * 3.2) / max(1, len(tips)), fs * 1.05), fs * 1.8)
    lab_w = max(P.text_w(rename.get(t, t), fs) for t in tips) + 6
    top = fs * 2.2
    xs = {}

    def depth(n, d=0.0):
        xs[id(n)] = d
        for c in n.children:
            depth(c, d + c.length)
    depth(tree)
    maxd = max(xs[id(l)] for l in leaves(tree)) or 1
    plot_w = rw - lab_w - 10
    X = lambda d: 4 + d / maxd * plot_w
    ys = {}
    for i, l in enumerate(leaves(tree)):
        ys[id(l)] = top + i * row

    def sety(n):
        if n.children:
            for c in n.children:
                sety(c)
            ys[id(n)] = (ys[id(n.children[0])] + ys[id(n.children[-1])]) / 2
    sety(tree)
    out = []
    for n in [tree] + [x for x in _walk(tree)]:
        for c in n.children:
            out.append(f'<path d="M{X(xs[id(n)]):.2f},{ys[id(c)]:.2f} H{X(xs[id(c)]):.2f}" stroke="#000" stroke-width="0.5" fill="none"/>')
        if n.children:
            out.append(f'<line x1="{X(xs[id(n)]):.2f}" x2="{X(xs[id(n)]):.2f}" y1="{ys[id(n.children[0])]:.2f}" '
                       f'y2="{ys[id(n.children[-1])]:.2f}" stroke="#000" stroke-width="0.5"/>')
            sh, uf = support(n.name)
            if uf is not None and n is not tree:
                strong = uf >= 95 and (sh is None or sh >= 80)
                if strong or uf >= 70:
                    out.append(f'<circle cx="{X(xs[id(n)]):.2f}" cy="{ys[id(n)]:.2f}" r="1.5" '
                               f'fill="{"#000" if strong else "#fff"}" stroke="#000" stroke-width="0.5">'
                               f'<title>{html.escape(n.name)}</title></circle>')
    for l in leaves(tree):
        t = l.name
        col = COL.get(klass[t], COL["copy"])
        out.append(f'<text x="{X(xs[id(l)]) + 3:.2f}" y="{ys[id(l)] + fs * 0.35:.2f}" font-size="{fs}" fill="{col}">'
                   f'{html.escape(rename.get(t, t))}<title>was {html.escape(t)}</title></text>')
    yb = top + len(tips) * row + fs * 0.3
    sb = _nice(maxd)
    out.append(f'<line x1="4" x2="{4 + sb / maxd * plot_w:.2f}" y1="{yb:.2f}" y2="{yb:.2f}" stroke="#000" stroke-width="0.8"/>')
    out.append(f'<text x="{6 + sb / maxd * plot_w:.2f}" y="{yb + fs * 0.35:.2f}" font-size="{fs}">{sb:g} subst./site</text>')
    yl = yb + fs * 1.5
    lx = 4
    for txt, fill in (("UFBoot ≥ 95, SH-aLRT ≥ 80", "#000"), ("UFBoot ≥ 70", "#fff")):
        out.append(f'<circle cx="{lx + 2:.2f}" cy="{yl:.2f}" r="1.5" fill="{fill}" stroke="#000" stroke-width="0.5"/>')
        out.append(f'<text x="{lx + 6:.2f}" y="{yl + fs * 0.35:.2f}" font-size="{fs}">{html.escape(txt)}</text>')
        lx += 10 + P.text_w(txt, fs) + 8
    rh = yl + fs
    right = (f'<svg xmlns="http://www.w3.org/2000/svg" width="{rw:.1f}pt" height="{rh:.1f}pt" viewBox="0 0 {rw:.2f} {rh:.2f}" '
             f'font-family="{P.FONT}">\n' + "\n".join(out) + "\n</svg>")

    # ---- compose ----
    def dims(svg):
        h = svg.split("\n", 1)[0]
        return (float(re.search(r'width="([\d.]+)pt"', h).group(1)), float(re.search(r'height="([\d.]+)pt"', h).group(1)))
    lw, lh = dims(left)
    rw_, rh_ = dims(right)
    W, H = width_total * P.PT_PER_MM, max(lh, rh_) + fs * 1.4
    def nest(svg, x, y):
        # nested panels in the parent's user units (pt), not 'pt' again, which would scale them twice
        head, rest = svg.split("\n", 1)
        head = re.sub(r'width="([\d.]+)pt"', r'width="\1"', head)
        head = re.sub(r'height="([\d.]+)pt"', r'height="\1"', head)
        return re.sub(r"^<svg ", f'<svg x="{x:.2f}" y="{y:.2f}" ', head, count=1) + "\n" + rest
    panel_c, ch = "", 0.0
    if a.aln:
        keep = [l.name for l in leaves(tree) if klass[l.name] in ("Exendin_anc", "copy")]
        csvg = alignment_panel(a.aln, keep, rename, {t: COL[klass[t]] for t in keep}, W, fs, a.aln_region)
        if csvg:
            _, ch = dims(csvg)
            panel_c = (f'\n<text x="2" y="{H + fs * 1.1:.2f}" font-size="{fs + 1}" font-weight="bold">C</text>\n'
                       + nest(csvg, 0, H + fs * 1.4))
            H += ch + fs * 1.4
    body = (f'<text x="2" y="{fs * 1.1:.2f}" font-size="{fs + 1}" font-weight="bold">A</text>\n'
            f'<text x="{lw + 2:.2f}" y="{fs * 1.1:.2f}" font-size="{fs + 1}" font-weight="bold">B</text>\n'
            + nest(left, 0, fs * 1.4) + "\n" + nest(right, lw, fs * 1.4) + panel_c)
    svg = (f'<svg xmlns="http://www.w3.org/2000/svg" width="{W:.1f}pt" height="{H:.1f}pt" viewBox="0 0 {W:.2f} {H:.2f}" '
           f'font-family="{P.FONT}">\n<rect width="100%" height="100%" fill="#fff"/>\n{body}\n</svg>\n')
    P.write(svg, a)
    N.info(f"figure {W / P.PT_PER_MM:.0f} x {H / P.PT_PER_MM:.0f} mm; names in {a.outdir}/names.tsv, "
           f"tree in {a.outdir}/renamed.treefile")


CLASSES = [("AILMFWV", "#80A0F0"), ("KR", "#F01505"), ("DE", "#C048C0"), ("NQST", "#15C015"),
           ("C", "#F08080"), ("G", "#F09048"), ("P", "#C0C000"), ("HY", "#15A4A4")]
RES_COL = {aa: c for aas, c in CLASSES for aa in aas}


def alignment_panel(path, keep, rename, namecol, W, fs, regions):
    """Panel C: the exendin rows of the alignment, tree order, all-gap columns removed, wrapped to width W."""
    seqs, name = {}, None
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            name = line[1:].split()[0]; seqs[name] = []
        elif name:
            seqs[name].append(line.upper())
    seqs = {k: "".join(v) for k, v in seqs.items()}
    rows = [(t, seqs[t]) for t in keep if t in seqs]
    miss = [t for t in keep if t not in seqs]
    if miss:
        N.warn(f"panel C: {len(miss)} exendin tip(s) not in {path}: {', '.join(miss[:5])}")
    if not rows:
        N.warn(f"panel C: no exendin tip of the tree is in {path}; panel C left out")
        return ""
    L = len(rows[0][1])
    cols = [j for j in range(L) if any(s[j] not in "-." for _, s in rows)]
    rows = [(t, "".join(s[j] for j in cols)) for t, s in rows]
    n = len(cols)
    names = [rename.get(t, t) for t, _ in rows]
    name_w = max(P.text_w(x, fs) for x in names) + 6
    cw, rh = fs * 0.66, fs * 1.2
    per = max(10, int((W - name_w - 6) // cw) // 10 * 10)
    blocks = [(b, min(b + per, n)) for b in range(0, n, per)]
    mono = 'font-family="Menlo, Consolas, DejaVu Sans Mono, monospace"'
    out, y = [], fs * 0.4
    for b0, b1 in blocks:
        top = y + fs * (2.2 if regions else 1.2)
        for reg in regions:
            span, _, lab = reg.partition(":")
            r0, _, r1 = span.partition("-")
            r0, r1 = max(int(r0) - 1, b0), min(int(r1), b1)
            if r0 < r1:
                x = name_w + (r0 - b0) * cw
                out.append(f'<rect x="{x:.2f}" y="{top - fs * 0.3:.2f}" width="{(r1 - r0) * cw:.2f}" '
                           f'height="{len(rows) * rh + fs * 0.3:.2f}" fill="none" stroke="#000" stroke-width="0.6"/>')
                if lab:
                    out.append(f'<text x="{x:.2f}" y="{top - fs * 1.4:.2f}" font-size="{fs}">{html.escape(lab)}</text>')
        for c in range(((b0 // 10) + 1) * 10, b1 + 1, 10):
            out.append(f'<text x="{name_w + (c - b0 - 0.5) * cw:.2f}" y="{top - fs * 0.45:.2f}" font-size="{fs}" '
                       f'text-anchor="middle">{c}</text>')
        for i, ((t, sq), nm) in enumerate(zip(rows, names)):
            yy = top + i * rh
            out.append(f'<text x="0" y="{yy + rh * 0.72:.2f}" font-size="{fs}" fill="{namecol.get(t, "#000")}">'
                       f'{html.escape(nm)}</text>')
            for j, ch in enumerate(sq[b0:b1]):
                if ch in "-.":
                    continue
                x = name_w + j * cw
                fill = "#E8312F" if ch in "X*" else RES_COL.get(ch, "#DDDDDD")
                out.append(f'<rect x="{x:.2f}" y="{yy:.2f}" width="{cw:.2f}" height="{rh - 0.5:.2f}" fill="{fill}" '
                           f'fill-opacity="0.45"/>')
                out.append(f'<text x="{x + cw / 2:.2f}" y="{yy + rh * 0.74:.2f}" font-size="{fs}" {mono} '
                           f'text-anchor="middle">{ch}</text>')
        y = top + len(rows) * rh + fs * 0.8
    N.info(f"panel C: {len(rows)} exendin rows x {n} columns (of {L}; all-gap columns removed), {len(blocks)} block(s)")
    return (f'<svg xmlns="http://www.w3.org/2000/svg" width="{W:.1f}pt" height="{y:.1f}pt" viewBox="0 0 {W:.2f} {y:.2f}" '
            f'font-family="{P.FONT}">\n' + "\n".join(out) + "\n</svg>")


def _walk(n):
    for c in n.children:
        yield c
        yield from _walk(c)


def _nice(m):
    best = 0.01
    for v in (0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0):
        if v <= m / 3:
            best = v
    return best


if __name__ == "__main__":
    main()
