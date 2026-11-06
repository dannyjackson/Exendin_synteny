#!/usr/bin/env python3
"""
Plot, for every species, where each TOGA2-annotated exendin gene sits on the
orthologous chromosome, between the two conserved boundary genes (SLC4A8 and
AQP5), with the species tree alongside. One track per contig, so tandem
duplications show up as several genes between the boundaries, and copies
elsewhere as extra tracks.

Inputs are TOGA2 only, plus the species tree for order:
  - the exendin genes: tips.tsv from 4b_collect_toga2_exendins.py (every row is
    drawn; rows with included = no are drawn hollow, so status filters stay visible);
  - the boundary genes: TOGA2 ortholog projections of the two genes in each query
    (orthology_classification.tsv + query_annotation.bed + loss_summary.tsv);
  - --species-tree (optional): Newick; species are drawn in its tip order with the
    tree (topology only, or --branch-lengths) to the left. Tree tips are matched to
    query codes exactly, or as Genus_species -> GenSpe (Heloderma_suspectum -> HelSus);
    --tree-map TSV (tree_tip<TAB>query) overrides. Unmatched names are reported.

Boundary genes are given as NAME=REFERENCE_GENE_ID. A NAME with no ID is looked up
in --gff: the gene whose attributes contain NAME as a whole word, nearest to the
HelSus exendin genes. The IDs actually used are printed, so check them.

One boundary locus per contig: overlapping projections of a boundary gene (isoforms,
chains) are merged, the locus nearest the exendin genes is kept, and within it the
projection with the best TOGA2 status, then the lowest mRNA number. A boundary gene
longer than --clip-kb is drawn clipped on the side away from the exendins (marked
'//'); its true coordinates are kept in the layout table.

Orientation: every track is drawn with the first --boundaries gene on the left
(default SLC4A8, which keeps HelSus chr2 in its own orientation, as in the
"SLC4A8 - ANC - AQP5" interval). If both boundaries are on a track, their order
decides it; otherwise a track is flipped when its boundary gene (or, failing that,
its exendin genes) lies on the opposite strand from HelSus.

Labels are stacked in tiers above the track so they do not overlap, with a leader
line to the gene when lifted. Labels drop the species prefix of the tip name.

Read-only on the TOGA2 directory. Writes the SVG (stdout or --svg), optional PNG,
and --layout (one row per drawn feature).

Usage
  python3 4d_plot_exendin_synteny.py \
      --run ${TOGADIR}/target_HelSus \
      --tips toga2_exendins/tips.tsv \
      --boundaries SLC4A8,AQP5=LOC_00012326 --gff ${HELSUS_EVIANN_GFF} \
      --species-tree ${CACTUSDIR}/Anguimorpha_species_tree_iteration-01_5-23-26.nwk \
      --svg exendin_synteny.svg --png exendin_synteny.png --layout exendin_synteny.tsv
"""
import argparse, csv, gzip, html, os, re, shutil, subprocess, sys, tempfile
from collections import defaultdict
from urllib.parse import unquote

RANK = {"FI": 0, "I": 1, "PI": 2, "UL": 3, "M": 4, "L": 5, "PG": 6, "N": 7}
PALETTE = ["#D55E00", "#0072B2", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#7F3C8D", "#8C564B"]
BOUNDARY_COL = "#555555"


def opener(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def gene_of(pr):
    return pr.split("-mRNA")[0]


def gff_lookup(gff, name, near):
    """Reference gene ID whose gene/mRNA attributes contain NAME as a word, nearest to `near`."""
    pat = re.compile(rf"(?<![A-Za-z0-9]){re.escape(name)}(?![A-Za-z0-9])", re.I)
    genes, parent, hits = {}, {}, set()
    with opener(gff) as fh:
        for line in fh:
            t = line.rstrip("\n").split("\t")
            if line.startswith("#") or len(t) < 9 or t[2] not in ("gene", "mRNA", "transcript"):
                continue
            at = dict(kv.split("=", 1) for kv in t[8].strip().strip(";").split(";") if "=" in kv)
            fid = at.get("ID", "")
            if t[2] == "gene":
                genes[fid] = (t[0], int(t[3]), int(t[4]))
            else:
                parent[fid] = at.get("Parent", "").split(",")[0]
            if any(pat.search(unquote(v)) for k, v in at.items() if k not in ("ID", "Parent")):
                hits.add(fid if t[2] == "gene" else parent[fid])
    hits = [h for h in hits if h in genes]
    if not hits:
        sys.exit(f"[error] no gene in {gff} matches '{name}'; give it as {name}=GENE_ID")
    chrom, pos = near

    def d(h):
        c, s, e = genes[h]
        return (c != chrom, abs((s + e) / 2 - pos))
    hits.sort(key=d)
    if len(hits) > 1:
        print(f"[info] {name}: {len(hits)} GFF matches; using the one nearest the HelSus exendins "
              f"({hits[0]} at {genes[hits[0]][0]}:{genes[hits[0]][1]}); others: {','.join(hits[1:6])}",
              file=sys.stderr)
    return hits[0]


def boundary_projections(qdir, ref_ids):
    """{ref_gene: [(chrom, start, end, strand, status, projection)]} for TOGA2 ortholog projections."""
    orth = set()
    p = os.path.join(qdir, "orthology_classification.tsv")
    if os.path.exists(p):
        with opener(p) as fh:
            next(fh, None)
            for line in fh:
                t = line.rstrip("\n").split("\t")
                if len(t) >= 5 and t[0] in ref_ids and t[3] != "None":
                    orth.add(t[3])
    status = {}
    p = os.path.join(qdir, "loss_summary.tsv")
    if os.path.exists(p):
        with opener(p) as fh:
            for line in fh:
                t = line.rstrip("\n").split("\t")
                if len(t) >= 3 and t[0] == "PROJECTION" and t[1] in orth:
                    status[t[1]] = t[2]
    out = defaultdict(list)
    p = os.path.join(qdir, "query_annotation.bed")
    if os.path.exists(p):
        with opener(p) as fh:
            for line in fh:
                t = line.rstrip("\n").split("\t")
                if len(t) >= 6 and t[3] in orth:
                    out[gene_of(t[3])].append((t[0], int(t[1]), int(t[2]), t[5], status.get(t[3], "NA"), t[3]))
    return out


def write_png(svg, path, w):
    with tempfile.NamedTemporaryFile("w", suffix=".svg", delete=False) as fh:
        fh.write(svg); tmp = fh.name
    cmds = []
    if shutil.which("rsvg-convert"):
        cmds.append(["rsvg-convert", "-w", str(w), "-o", path, tmp])
    if shutil.which("inkscape"):
        cmds.append(["inkscape", tmp, "--export-type=png", f"--export-filename={path}", f"--export-width={w}"])
    if shutil.which("convert"):
        cmds.append(["convert", "-density", "150", "-background", "white", tmp, path])
    for c in cmds:
        if subprocess.run(c, capture_output=True).returncode == 0 and os.path.exists(path):
            print(f"[info] wrote {path} with {c[0]}", file=sys.stderr); os.unlink(tmp); return
    try:
        import cairosvg
        cairosvg.svg2png(bytestring=svg.encode(), write_to=path, output_width=w)
        print(f"[info] wrote {path} with cairosvg", file=sys.stderr); os.unlink(tmp); return
    except Exception:
        pass
    os.unlink(tmp)
    print("[warn] no SVG converter found; PNG not written (mamba install -c conda-forge cairosvg)", file=sys.stderr)




# ---------------- species tree ----------------
class Node:
    __slots__ = ("name", "length", "children", "x", "y")

    def __init__(self):
        self.name, self.length, self.children, self.x, self.y = "", 0.0, [], 0.0, 0.0


def parse_newick(s):
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
        elif c.isspace():
            i += 1
        else:
            m = re.match(r"('[^']*'|[^,():;]+)", s[i:]); cur.name = m.group(1).strip("'").strip(); i += len(m.group(0))
    return root.children[0] if len(root.children) == 1 and not root.name else root


def leaves(n):
    return [n] if not n.children else [l for c in n.children for l in leaves(c)]


def match_tip(name, queries):
    """Tree tip -> query code: exact, case-insensitive, or Genus_species -> GenSpe."""
    if name in queries:
        return name
    low = {q.lower(): q for q in queries}
    if name.lower() in low:
        return low[name.lower()]
    parts = [p for p in re.split(r"[_\s.]+", name) if p]
    if len(parts) >= 2:
        code = (parts[0][:3] + parts[1][:3]).lower()
        if code in low:
            return low[code]
    return None


def prune(n, keep):
    """Drop leaves not in keep; collapse nodes left with one child (lengths are summed)."""
    if not n.children:
        return n if n.name in keep else None
    kids = [k for k in (prune(c, keep) for c in n.children) if k is not None]
    if not kids:
        return None
    if len(kids) == 1:
        kids[0].length += n.length
        return kids[0]
    n.children = kids
    return n


# ---------------- tracks ----------------
def choose_boundary(ps, ref_point):
    """Merge overlapping projections into loci; keep the locus nearest ref_point (or the best one)."""
    ps = sorted(ps, key=lambda p: p[1])
    loci, cur, end = [], [], -1
    for p in ps:
        if cur and p[1] > end:
            loci.append(cur); cur, end = [], -1
        cur.append(p); end = max(end, p[2])
    if cur:
        loci.append(cur)

    def best(locus):
        return sorted(locus, key=lambda p: (RANK.get(p[4], 99), isoform(p[5]), p[2] - p[1]))[0]
    if ref_point is None:
        return best(sorted(loci, key=lambda L: min(RANK.get(p[4], 99) for p in L))[0])
    return best(min(loci, key=lambda L: min(min(abs(p[1] - ref_point), abs(p[2] - ref_point)) for p in L)))


def isoform(projection):
    m = re.search(r"-mRNA-(\d+)", projection)
    return int(m.group(1)) if m else 99


def place_labels(items, max_tiers=8, gap=4):
    """items: [(centre_x, width)] -> tier per item (0 = lowest), greedy left to right."""
    tiers, out = [], []
    for cx, w in items:
        l, r = cx - w / 2, cx + w / 2
        for k, right in enumerate(tiers):
            if l >= right + gap:
                tiers[k] = r; out.append(k); break
        else:
            if len(tiers) < max_tiers:
                tiers.append(r); out.append(len(tiers) - 1)
            else:
                k = min(range(len(tiers)), key=lambda k: tiers[k])
                tiers[k] = max(tiers[k], r); out.append(k)
    return out


# ---------------- shared by 4d (TOGA2) and 3b (step-3 scan) ----------------
# A feature is a dict with: kind ('exendin' | 'boundary'), name, label, group (colour),
# chrom, start, end (genomic), strand, filled (bool), dashed (bool), hover (str).
# contigs: {query: {chrom: [features]}}.

def build_tracks(contigs, queries, bnames, ref_query, clip_kb):
    """Clip long boundary genes, orient every contig (first boundary name on the left)."""
    left_name = bnames[0]
    ref_strand, votes_by_group = {}, defaultdict(lambda: defaultdict(int))
    for fs in contigs.get(ref_query, {}).values():
        for f in fs:
            if f["kind"] == "boundary":
                ref_strand.setdefault(f["name"], f["strand"])
            else:
                votes_by_group[f["group"]][f["strand"]] += 1
    # majority strand of each exendin group in the reference (ties: group not used)
    ex_ref_strand = {g: max(v, key=v.get) for g, v in votes_by_group.items()
                     if sorted(v.values())[-1] > (sorted(v.values())[-2] if len(v) > 1 else 0)}
    clip = int(clip_kb * 1000)
    tracks = []
    for q in queries:
        cq = contigs.get(q, {})
        if not cq:
            tracks.append(dict(query=q, chrom="", feats=[], flip=False, why="nothing found", main=True))
            continue
        for fs in cq.values():
            for f in fs:
                f.setdefault("clip", "")
                if f["kind"] != "boundary" or not clip or f["end"] - f["start"] <= clip:
                    continue
                others = [g for g in fs if g is not f]
                if not others:
                    continue                      # handled after orientation
                ref = sum((g["start"] + g["end"]) / 2 for g in others) / len(others)
                f["true"] = (f["start"], f["end"])
                if (f["start"] + f["end"]) / 2 < ref:
                    f["start"], f["clip"] = f["end"] - clip, "start"
                else:
                    f["end"], f["clip"] = f["start"] + clip, "end"
        main = max(cq, key=lambda c: (sum(f["kind"] == "boundary" for f in cq[c]),
                                      sum(f["kind"] == "exendin" for f in cq[c])))
        for c in sorted(cq, key=lambda c: (c != main, c)):
            fs = cq[c]
            b = {f["name"]: f for f in fs if f["kind"] == "boundary"}
            flip, why = False, "as is"
            if len(b) == 2:
                flip = b[left_name]["start"] > b[bnames[1]]["start"]; why = "boundary order"
            elif b:
                n, f = next(iter(b.items()))
                if n in ref_strand:
                    flip = f["strand"] != ref_strand[n]; why = f"{n} strand"
            else:
                votes = [f["strand"] != ex_ref_strand[f["group"]] for f in fs if f["group"] in ex_ref_strand]
                if votes:
                    flip = sum(votes) > len(votes) / 2; why = "exendin strand"
            if len(fs) == 1 and fs[0]["kind"] == "boundary" and clip and fs[0]["end"] - fs[0]["start"] > clip:
                f = fs[0]
                f["true"] = (f["start"], f["end"])
                if (f["name"] == left_name) != flip:          # left boundary faces right
                    f["start"], f["clip"] = f["end"] - clip, "start"
                else:
                    f["end"], f["clip"] = f["start"] + clip, "end"
            tracks.append(dict(query=q, chrom=c, feats=fs, flip=flip, why=why, main=c == main))
    return tracks


def species_order(queries, ref_query, species_tree=None, tree_map=None, order=None):
    """(order, pruned tree or None)."""
    if species_tree:
        tree = parse_newick(open(species_tree).read())
        manual = {}
        if tree_map:
            for line in open(tree_map):
                f = line.rstrip("\n").split("\t")
                if len(f) >= 2 and not line.startswith("#"):
                    manual[f[0]] = f[1]
        matched, unmatched = set(), []
        for lf in leaves(tree):
            q = manual.get(lf.name) or match_tip(lf.name, queries)
            if q:
                matched.add(q); lf.name = q
            else:
                unmatched.append(lf.name)
        if unmatched:
            print(f"[info] tree tips with no matching species (pruned): {', '.join(unmatched)}", file=sys.stderr)
        tree = prune(tree, matched)
        out = [lf.name for lf in leaves(tree)] if tree else []
        missing = [q for q in queries if q not in out]
        if missing:
            print(f"[warn] species not in the tree (drawn below it): {', '.join(missing)}; "
                  f"add them with --tree-map", file=sys.stderr)
        return out + missing, tree
    if order:
        return ([l.strip() for l in open(order) if l.strip()] if os.path.exists(order)
                else order.split(",")), None
    return [ref_query] + [q for q in queries if q != ref_query], None


def render(tracks, tree, *, title, note, legend_extra=(), width=1400, pad_kb=10.0, min_kb=0.0,
           branch_lengths=False, bold_re=r"^Hel"):
    """SVG text and layout rows. `tracks` must already be in display order."""
    pad = pad_kb * 1000
    for t in tracks:
        if not t["feats"]:
            t["span"] = 0; continue
        lo = min(f["start"] for f in t["feats"]) - pad
        hi = max(f["end"] for f in t["feats"]) + pad
        t["lo"], t["hi"], t["span"] = lo, hi, hi - lo
        for f in t["feats"]:
            if t["flip"]:
                f["x0"], f["x1"] = hi - f["end"], hi - f["start"]
                f["dir"] = {"+": "-", "-": "+"}.get(f["strand"], ".")
                f["clipside"] = {"start": "right", "end": "left"}.get(f.get("clip", ""), "")
            else:
                f["x0"], f["x1"] = f["start"] - lo, f["end"] - lo
                f["dir"] = f["strand"]
                f["clipside"] = {"start": "left", "end": "right"}.get(f.get("clip", ""), "")
    maxspan = max([t["span"] for t in tracks] + [min_kb * 1000, 1])

    tree_w = 170 if tree else 0
    left, sp_w, ctg_w = 15, 70, 175
    plot_x0 = left + tree_w + sp_w + ctg_w
    plot_w = width - plot_x0 - 25
    X = lambda bp: plot_x0 + bp / maxspan * plot_w
    CH = 4.7

    y = 70
    for t in tracks:
        t["feats"] = sorted(t["feats"], key=lambda f: f["x0"])
        ex = [f for f in t["feats"] if f["kind"] == "exendin"]
        tiers = place_labels([((X(f["x0"]) + X(f["x1"])) / 2, len(f["label"]) * CH) for f in ex])
        for f, k in zip(ex, tiers):
            f["tier"] = k
        n_t = (max(tiers) + 1) if tiers else 0
        t["y"] = y + 12 + n_t * 10
        t["h"] = 12 + n_t * 10 + 24
        y += t["h"]
    H = y + 70 + (20 if legend_extra else 0)

    groups = sorted({f["group"] for t in tracks for f in t["feats"] if f["kind"] == "exendin"})
    col = {g: PALETTE[i % len(PALETTE)] for i, g in enumerate(groups)}
    out = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{H}" '
           f'font-family="Helvetica, Arial, sans-serif" font-size="11">',
           '<rect width="100%" height="100%" fill="white"/>',
           f'<text x="{left}" y="24" font-size="15" font-weight="bold">{html.escape(title)}</text>',
           f'<text x="{left}" y="42" font-size="10" fill="#555">{html.escape(note)}</text>']

    first = {}
    for i, t in enumerate(tracks):
        top = t["y"] - (t["h"] - 24)
        if t["query"] not in first:
            first[t["query"]] = t
            if i:
                out.append(f'<line x1="{left + tree_w}" x2="{width - 15}" y1="{top - 2:.1f}" y2="{top - 2:.1f}" '
                           f'stroke="#e3e3e3"/>')
            bold = ' font-weight="bold"' if re.match(bold_re, t["query"]) else ""
            out.append(f'<text x="{left + tree_w + 6}" y="{t["y"] + 4:.1f}"{bold}>{html.escape(t["query"])}</text>')
        if not t["feats"]:
            out.append(f'<text x="{plot_x0}" y="{t["y"] + 4:.1f}" fill="#999" font-style="italic">{t["why"]}</text>')
            continue
        a_mb, b_mb = (t["hi"], t["lo"]) if t["flip"] else (t["lo"], t["hi"])
        out.append(f'<text x="{left + tree_w + sp_w}" y="{t["y"] + 1:.1f}" font-size="9" fill="#444">'
                   f'{html.escape(t["chrom"][:26])}{" (rev)" if t["flip"] else ""}</text>')
        out.append(f'<text x="{left + tree_w + sp_w}" y="{t["y"] + 11:.1f}" font-size="8" fill="#888">'
                   f'{max(a_mb, 0) / 1e6:.3f}–{max(b_mb, 0) / 1e6:.3f} Mb</text>')

    if tree:
        for lf in leaves(tree):
            lf.y = first[lf.name]["y"] if lf.name in first else 0
        def height(n):
            return 0 if not n.children else 1 + max(height(c) for c in n.children)
        if branch_lengths:
            def setx(n, d):
                n.x = d
                for c in n.children:
                    setx(c, d + c.length)
            setx(tree, 0.0)
        else:
            H0 = height(tree)
            def setx(n, d=None):
                n.x = H0 - height(n)
                for c in n.children:
                    setx(c)
            setx(tree)
        maxx = max(l.x for l in leaves(tree)) or 1.0
        TX = lambda v: left + v / maxx * (tree_w - 10)
        def sety(n):
            if n.children:
                for c in n.children:
                    sety(c)
                n.y = (n.children[0].y + n.children[-1].y) / 2
        sety(tree)
        def draw(n):
            for c in n.children:
                out.append(f'<line x1="{TX(n.x):.1f}" y1="{c.y:.1f}" x2="{TX(c.x):.1f}" y2="{c.y:.1f}" '
                           f'stroke="#333" stroke-width="1.2"/>')
                draw(c)
            if n.children:
                out.append(f'<line x1="{TX(n.x):.1f}" y1="{n.children[0].y:.1f}" x2="{TX(n.x):.1f}" '
                           f'y2="{n.children[-1].y:.1f}" stroke="#333" stroke-width="1.2"/>')
            else:
                out.append(f'<line x1="{TX(n.x):.1f}" y1="{n.y:.1f}" x2="{left + tree_w + 3}" y2="{n.y:.1f}" '
                           f'stroke="#ccc" stroke-dasharray="2,2"/>')
        draw(tree)

    lay = []
    for t in tracks:
        if not t["feats"]:
            continue
        y = t["y"]
        out.append(f'<line x1="{X(0):.1f}" x2="{X(t["span"]):.1f}" y1="{y}" y2="{y}" stroke="#999" stroke-width="2"/>')
        for f in t["feats"]:
            x0, x1 = X(f["x0"]), X(f["x1"])
            if x1 - x0 < 4:
                m = (x0 + x1) / 2; x0, x1 = m - 2, m + 2
            h, tip = 12, min(5, (x1 - x0) * 0.6)
            c = BOUNDARY_COL if f["kind"] == "boundary" else col[f["group"]]
            if f["dir"] == "-":
                pts = [(x0, y), (x0 + tip, y - h / 2), (x1, y - h / 2), (x1, y + h / 2), (x0 + tip, y + h / 2)]
            elif f["dir"] == "+":
                pts = [(x0, y - h / 2), (x1 - tip, y - h / 2), (x1, y), (x1 - tip, y + h / 2), (x0, y + h / 2)]
            else:
                pts = [(x0, y - h / 2), (x1, y - h / 2), (x1, y + h / 2), (x0, y + h / 2)]
            fill = c if f["filled"] else "white"
            dash = ' stroke-dasharray="2,1.5"' if f["dashed"] else ""
            s0, e0 = f.get("true", (f["start"], f["end"]))
            title_t = f'{f["name"]} {f["chrom"]}:{s0}-{e0} ({f["strand"]}) {f.get("hover", "")}'
            out.append(f'<polygon points="{" ".join(f"{px:.1f},{py:.1f}" for px, py in pts)}" fill="{fill}" '
                       f'stroke="{c}" stroke-width="1.2"{dash}><title>{html.escape(title_t)}</title></polygon>')
            if f["clipside"]:
                cx = x0 + 4 if f["clipside"] == "left" else x1 - 10
                for dx in (0, 4):
                    out.append(f'<line x1="{cx + dx:.1f}" y1="{y + h / 2 + 2:.1f}" x2="{cx + dx + 3:.1f}" '
                               f'y2="{y - h / 2 - 2:.1f}" stroke="white" stroke-width="2"/>')
            if f["kind"] == "boundary":
                out.append(f'<text x="{(x0 + x1) / 2:.1f}" y="{y + h / 2 + 10:.1f}" font-size="8" '
                           f'text-anchor="middle" fill="#444">{html.escape(f["name"])}</text>')
            else:
                ly = y - h / 2 - 3 - f["tier"] * 10
                if f["tier"]:
                    out.append(f'<line x1="{(x0 + x1) / 2:.1f}" y1="{y - h / 2:.1f}" x2="{(x0 + x1) / 2:.1f}" '
                               f'y2="{ly + 2:.1f}" stroke="{c}" stroke-width="0.5" opacity="0.6"/>')
                out.append(f'<text x="{(x0 + x1) / 2:.1f}" y="{ly:.1f}" font-size="8" text-anchor="middle" '
                           f'fill="{c}">{html.escape(f["label"])}</text>')
            lay.append([t["query"], t["chrom"], "rev" if t["flip"] else "fwd", t["why"], f["kind"], f["name"],
                        f["group"], s0, e0, f["strand"], f.get("status", ""),
                        "yes" if f["filled"] else "no", f.get("clip") or "no",
                        f"{f['x0'] / 1000:.1f}", f"{f['x1'] / 1000:.1f}"])

    yb = H - 40 - (20 if legend_extra else 0)
    step = next(s for s in (1e3, 2e3, 5e3, 1e4, 2e4, 5e4, 1e5, 2e5, 5e5, 1e6, 2e6, 5e6, 1e7)
                if s >= maxspan / 6)
    out.append(f'<line x1="{X(0):.1f}" x2="{X(step):.1f}" y1="{yb}" y2="{yb}" stroke="#333" stroke-width="2"/>')
    out.append(f'<text x="{X(0):.1f}" y="{yb + 13}" font-size="9">{step / 1000:g} kb</text>')
    lx = X(step) + 40
    for g, c in list(col.items()) + [("boundary gene", BOUNDARY_COL)]:
        out.append(f'<rect x="{lx:.1f}" y="{yb - 6}" width="12" height="10" fill="{c}"/>')
        out.append(f'<text x="{lx + 16:.1f}" y="{yb + 3}" font-size="10">{html.escape(g)}</text>')
        lx += 30 + 6.5 * len(g)
    lx2 = X(0)
    for text, filled, dashed in legend_extra:
        yy = yb + 34
        fill = "#777" if filled else "white"
        dash = ' stroke-dasharray="2,1.5"' if dashed else ""
        out.append(f'<rect x="{lx2:.1f}" y="{yy - 8}" width="12" height="10" fill="{fill}" stroke="#777"{dash}/>')
        out.append(f'<text x="{lx2 + 16:.1f}" y="{yy + 1}" font-size="10">{html.escape(text)}</text>')
        lx2 += 40 + 5.8 * len(text)
    out.append("</svg>")
    return "\n".join(out), lay, maxspan


LAYOUT_COLS = ["query", "contig", "orientation", "oriented_by", "kind", "name", "group", "start", "end",
               "strand", "status", "filled", "clipped", "x0_kb", "x1_kb"]


def write_outputs(svg, lay, svg_path=None, png_path=None, layout_path=None, width=1400):
    if svg_path:
        open(svg_path, "w").write(svg)
    else:
        sys.stdout.write(svg)
    if png_path:
        write_png(svg, png_path, width * 2)
    if layout_path:
        with open(layout_path, "w") as fh:
            fh.write("\t".join(LAYOUT_COLS) + "\n")
            for r in lay:
                fh.write("\t".join(map(str, r)) + "\n")


def add_common_args(ap):
    ap.add_argument("--species-tree", help="Newick species tree; sets the order and is drawn on the left")
    ap.add_argument("--tree-map", help="TSV tree_tip<TAB>species, for tips that do not match automatically")
    ap.add_argument("--branch-lengths", action="store_true", help="draw tree branch lengths (default: topology)")
    ap.add_argument("--order", help="comma-separated species order, or a file with one per line "
                                    "(ignored with --species-tree)")
    ap.add_argument("--reference-query", default="HelSus", help="species that sets the orientation")
    ap.add_argument("--pad-kb", type=float, default=10.0, help="margin around each track")
    ap.add_argument("--clip-kb", type=float, default=50.0,
                    help="draw boundary genes longer than this clipped (0 = never)")
    ap.add_argument("--min-kb", type=float, default=0.0, help="minimum axis length (kb)")
    ap.add_argument("--width", type=int, default=1400)
    ap.add_argument("--svg", help="output SVG (default stdout)")
    ap.add_argument("--png")
    ap.add_argument("--layout", help="TSV of every drawn feature")


def finish(a, contigs, queries, bnames, title, note, legend_extra=()):
    tracks = build_tracks(contigs, queries, bnames, a.reference_query, a.clip_kb)
    order, tree = species_order(queries, a.reference_query, a.species_tree, a.tree_map, a.order)
    rank = {q: i for i, q in enumerate(order)}
    tracks.sort(key=lambda t: (rank.get(t["query"], len(rank)), not t["main"], t["chrom"]))
    svg, lay, span = render(tracks, tree, title=title, note=note, legend_extra=legend_extra, width=a.width,
                            pad_kb=a.pad_kb, min_kb=a.min_kb, branch_lengths=a.branch_lengths)
    write_outputs(svg, lay, a.svg, a.png, a.layout, a.width)
    print(f"[info] {len(tracks)} tracks from {len({t['query'] for t in tracks})} species; "
          f"axis {span / 1000:.0f} kb", file=sys.stderr)


# ---------------- TOGA2 input ----------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", required=True, help="${TOGADIR}/target_HelSus")
    ap.add_argument("--tips", required=True, help="tips.tsv from 4b")
    ap.add_argument("--boundaries", default="SLC4A8,AQP5=LOC_00012326",
                    help="two NAME[=REF_GENE_ID], left one first")
    ap.add_argument("--gff", help="reference GFF, to look up boundary names given without an ID")
    ap.add_argument("--title", default="Exendin genes between the conserved boundary genes (TOGA2, HelSus reference)")
    add_common_args(ap)
    a = ap.parse_args()

    with open(a.tips) as fh:
        tips = [r for r in csv.DictReader(fh, delimiter="\t") if r["q_chrom"] != "NA"]
    if not tips:
        sys.exit("[error] no tips with coordinates in " + a.tips)

    ref_tips = [t for t in tips if t["query"] == a.reference_query]
    near = ((ref_tips[0]["q_chrom"], sum(int(t["q_start"]) for t in ref_tips) / len(ref_tips))
            if ref_tips else ("", 0))
    bnames, bids = [], {}
    for item in a.boundaries.split(","):
        name, _, gid = item.partition("=")
        if not gid:
            if not a.gff:
                sys.exit(f"[error] boundary {name} has no gene ID; give {name}=GENE_ID or --gff")
            gid = gff_lookup(a.gff, name, near)
        bnames.append(name); bids[gid] = name
        print(f"[info] boundary {name} = {gid}", file=sys.stderr)
    if len(bnames) != 2:
        sys.exit("[error] --boundaries needs exactly two genes")

    contigs = defaultdict(lambda: defaultdict(list))
    for t in tips:
        contigs[t["query"]][t["q_chrom"]].append(dict(
            kind="exendin", name=t["tip"], label=t["tip"].split("_", 1)[1] if "_" in t["tip"] else t["tip"],
            group=t["group"], chrom=t["q_chrom"], start=int(t["q_start"]), end=int(t["q_end"]),
            strand=t["q_strand"], status=t["status"], filled=t["included"] == "yes",
            dashed=t.get("kind") == "paralog",
            hover=f'status {t["status"]} {t.get("orthology_class", "")} {t.get("q_gene", "")}'))
    queries = sorted(d[6:] for d in os.listdir(a.run) if d.startswith("query_"))
    for q in queries:
        bnd = boundary_projections(os.path.join(a.run, "query_" + q), set(bids))
        ex_mid = {c: sum((f["start"] + f["end"]) / 2 for f in fs) / len(fs) for c, fs in contigs[q].items()}
        for gid, name in bids.items():
            ps = bnd.get(gid, [])
            on = defaultdict(list)
            for p in ps:
                if p[0] in ex_mid:
                    on[p[0]].append(p)
            picks = ([choose_boundary(v, ex_mid[c]) for c, v in on.items()] if on
                     else ([choose_boundary(ps, None)] if ps else []))
            for p in picks:
                contigs[q][p[0]].append(dict(kind="boundary", name=name, label=name, group=name, chrom=p[0],
                                             start=p[1], end=p[2], strand=p[3], status=p[4], filled=True,
                                             dashed=False, hover=f"status {p[4]} {p[5]}"))

    finish(a, contigs, queries, bnames, a.title,
           f"Tracks oriented with {bnames[0]} on the left; arrows give strand after orientation. "
           f"Hollow = not in the tree (status/length filter); dashed = paralogous projection; "
           f"// = boundary gene clipped to {a.clip_kb:g} kb.")


if __name__ == "__main__":
    main()
