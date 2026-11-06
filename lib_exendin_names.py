#!/usr/bin/env python3
"""
Shared helpers for the synteny plots of steps 2, 3 and 4. Not run on its own.

  name_group(text)     annotation text -> colour group. 'exendin-2' and 'exendin-2-like'
                       give the same group, so a gene is coloured by its name whatever
                       the annotation pipeline appended.
  short_label(text)    compact label for the plot ('exendin-2-like' -> 'Ex2-like').
  is_family(group)     True for exendin / GCGL groups (the ligands of this gene family);
                       False for proglucagon and anything else.
  load_engine()        the drawing code of the old repo's 4d_plot_exendin_synteny.py,
                       copied here verbatim as lib_synteny_engine.py.
  finish(...)          engine.finish(), with two additions: every group always gets the
                       same colour in every figure (FIXED_COLOURS), and species with
                       nothing to draw can say why ('no NCBI annotation').
  read_species_table() species.tsv from step 1.
  find_genome_files()  FASTA / GFF / protein.faa of one species under ${GENOMEDIR}.
"""
import csv, glob, importlib.util, os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))

# Okabe-Ito based. Groups not listed here get EXTRA colours in sorted order.
FIXED_COLOURS = {
    "exendin-1": "#E69F00",
    "exendin-2": "#D55E00",
    "exendin-3": "#56B4E9",
    "exendin-4": "#0072B2",
    "exendin (unnumbered)": "#CC79A7",
    "GCGL (glucagon-like)": "#009E73",
    "proglucagon": "#7F7F7F",
}
EXTRA_COLOURS = ["#7F3C8D", "#8C564B", "#B8860B", "#2F4F4F", "#E7298A", "#1B9E77", "#A6761D", "#666666"]

_EX_NUM = re.compile(r"exendin[\s_-]*(\d+)", re.I)
_EX_ANY = re.compile(r"exendin", re.I)
_GCGL = re.compile(r"glucagon[\s_-]*like(?![\s_-]*peptide)|\bGCGL\b", re.I)
_PROGCG = re.compile(r"pro[\s_-]*glucagon|\bglucagon\b|\bGCG\b", re.I)
_RECEPTOR = re.compile(r"receptor", re.I)


def name_group(text):
    """Annotation text (gene symbol, description, product, label) -> colour group.

    exendin-N ... (any suffix: -like, transcript variant ...)   -> 'exendin-N'
    several numbers ('exendin-3/4', 'EX2_EX4')                    -> 'exendin-2/4'
    exendin without a number                                       -> 'exendin (unnumbered)'
    glucagon-like (not '-like peptide', not a receptor), GCGL     -> 'GCGL (glucagon-like)'
    glucagon / pro-glucagon / GCG (not a receptor)                -> 'proglucagon'
    TOGA2 labels EX2, EX4, EXN (from 3a_list_reference_exendins)  -> as exendin-N / unnumbered
    anything else                                                  -> 'other: <text>'
    """
    t = text or ""
    if _EX_ANY.search(t):
        nums = sorted(set(_EX_NUM.findall(t)), key=int)
        return "exendin-" + "/".join(nums) if nums else "exendin (unnumbered)"
    labs = re.findall(r"(?<![A-Za-z])EX(\d+|N)(?![A-Za-z])", t)
    if labs:
        nums = sorted({x for x in labs if x != "N"}, key=int)
        return "exendin-" + "/".join(nums) if nums else "exendin (unnumbered)"
    if _RECEPTOR.search(t):
        return "other: " + t.strip()[:40]
    if _GCGL.search(t):
        return "GCGL (glucagon-like)"
    if _PROGCG.search(t):
        return "proglucagon"
    return "other: " + t.strip()[:40]


def is_family(group):
    return group.startswith("exendin") or group.startswith("GCGL")


def short_label(text):
    """'exendin-2-like' -> 'Ex2-like'; 'exendin-4' -> 'Ex4'; 'glucagon-like' -> 'GCGL'."""
    t = (text or "").strip()
    m = _EX_NUM.search(t)
    if m:
        like = "-like" if re.search(r"like", t[m.end():m.end() + 8], re.I) else ""
        return f"Ex{m.group(1)}{like}"
    if _EX_ANY.search(t):
        return "Ex-like" if re.search(r"like", t, re.I) else "Ex"
    g = name_group(t)
    if g.startswith("GCGL"):
        return "GCGL"
    if g == "proglucagon":
        return "GCG"
    return t[:12]


def group_abbrev(group):
    """'exendin-2' -> 'Ex2', 'exendin-2/4' -> 'Ex2_4', 'exendin (unnumbered)' -> 'ExN', GCGL -> 'GCGL'."""
    if group.startswith("exendin-"):
        return "Ex" + group[8:].replace("/", "_")
    if group.startswith("exendin"):
        return "ExN"
    if group.startswith("GCGL"):
        return "GCGL"
    if group == "proglucagon":
        return "GCG"
    return re.sub(r"[^A-Za-z0-9]+", "", group.replace("other:", ""))[:10] or "other"


def colour_map(groups):
    col, k = {}, 0
    for g in sorted(groups):
        if g in FIXED_COLOURS:
            col[g] = FIXED_COLOURS[g]
        else:
            col[g] = EXTRA_COLOURS[k % len(EXTRA_COLOURS)]; k += 1
    return col


def load_engine():
    path = os.path.join(HERE, "lib_synteny_engine.py")
    spec = importlib.util.spec_from_file_location("synteny_engine", path)
    S = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(S)
    return S


def add_split_arg(ap):
    ap.add_argument("--max-gap-kb", type=float, default=1000.0,
                    help="within one contig, a stretch longer than this with no drawn feature is "
                         "shortened and drawn as a // break (e.g. boundary genes tens of Mb apart "
                         "after a rearrangement); 0 = never")
    ap.add_argument("--break-kb", type=float, default=20.0,
                    help="width a shortened stretch is drawn at, in kb of the axis")
    ap.add_argument("--inset", default="",
                    help="comma-separated species codes drawn in a box at the top right instead of in the "
                         "tree-ordered panel (outgroups not on the tree, e.g. GalGal)")
    ap.add_argument("--inset-title", default="Not on the tree", help="heading of the --inset box")
    ap.add_argument("--inset-width", type=int, default=640, help="width of the --inset box (px)")


def compress_gaps(tracks, max_gap_kb, break_kb):
    """Shorten long empty stretches inside each track; record where, for the // marks.

    One track spans from its first to its last feature, and the axis of the whole figure is
    the longest track, so two features far apart on one contig (AnoSag: SLC4A8 and AQP5
    57 Mb apart, nothing annotated between) squash every other species. Features after a
    gap > max_gap_kb are shifted left so the gap is drawn break_kb wide. True coordinates
    are kept in f['true'], which the engine already uses for hover text and the layout
    table. Features on DIFFERENT contigs (an assembly break, as in AnoCar) are untouched:
    they stay separate tracks.
    """
    if not max_gap_kb:
        return
    gap, width = max_gap_kb * 1000, break_kb * 1000
    for t in tracks:
        fs = sorted(t["feats"], key=lambda f: f["start"])
        t["breaks"], shift, end = [], 0, None
        for f in fs:
            f.setdefault("true", (f["start"], f["end"]))
            if end is not None and f["start"] - shift - end > gap:
                g = f["start"] - shift - end
                t["breaks"].append((end + width / 2, g))
                shift += g - width
            f["start"] -= shift; f["end"] -= shift
            end = f["end"] if end is None else max(end, f["end"])
        t["shift"] = shift
        if t["breaks"]:
            info(f"{t['query']}: {t['chrom']}: " + ", ".join(f"{g / 1e6:.1f} Mb" for _, g in t["breaks"])
                 + f" with nothing drawn shortened to {break_kb:g} kb (// on the track)")


def _draw_breaks(svg, tracks, tree, width, maxspan):
    """Add // marks and the true length at each shortened stretch (positions as in the engine's render)."""
    tree_w = 170 if tree else 0
    plot_x0 = 15 + tree_w + 70 + 175                      # left + tree_w + sp_w + ctg_w in render()
    plot_w = width - plot_x0 - 25
    X = lambda bp: plot_x0 + bp / maxspan * plot_w
    add = []
    for t in tracks:
        if not t.get("breaks") or not t["feats"]:
            continue
        y = t["y"]
        for pos, g in t["breaks"]:
            x = X(t["hi"] - pos) if t["flip"] else X(pos - t["lo"])
            add.append(f'<rect x="{x - 5:.1f}" y="{y - 3}" width="10" height="6" fill="white"/>')
            for dx in (-3, 2):
                add.append(f'<line x1="{x + dx - 2:.1f}" y1="{y + 5}" x2="{x + dx + 2:.1f}" y2="{y - 5}" '
                           f'stroke="#555" stroke-width="1.3"/>')
            add.append(f'<text x="{x:.1f}" y="{y + 15}" font-size="8" text-anchor="middle" fill="#777">'
                       f'{g / 1e6:.1f} Mb</text><title>{g / 1e6:.2f} Mb on {html_escape(t["chrom"])} '
                       f'with nothing drawn (shortened)</title>')
        # the engine printed the range from the shortened coordinates: put the true end back
        lo, hi = t["lo"], t["hi"]
        a_mb, b_mb = (hi, lo) if t["flip"] else (lo, hi)
        a2, b2 = (hi + t["shift"], lo) if t["flip"] else (lo, hi + t["shift"])
        old = f'{max(a_mb, 0) / 1e6:.3f}–{max(b_mb, 0) / 1e6:.3f} Mb'
        svg = svg.replace(old, f'{max(a2, 0) / 1e6:.3f}–{max(b2, 0) / 1e6:.3f} Mb (//)', 1)
    return svg.replace("</svg>", "\n".join(add) + "\n</svg>") if add else svg


def html_escape(x):
    import html
    return html.escape(str(x))


def _track_boxes(tracks, tree, width, maxspan):
    """(x0, x1, y0, y1) drawn by each track in the engine's render (same geometry as _draw_breaks)."""
    tree_w = 170 if tree else 0
    plot_x0 = 15 + tree_w + 70 + 175
    plot_w = width - plot_x0 - 25
    out = []
    for t in tracks:
        top = t["y"] - (t["h"] - 24)
        if t["feats"]:
            x1 = plot_x0 + t["span"] / maxspan * plot_w
        else:
            x1 = plot_x0 + 6 * len(t.get("why", ""))          # italic 'no ...' text
        out.append((plot_x0 - 5, x1 + 30, top, t["y"] + 22))   # +30 for labels/arrow tips
    return out


def _compose_inset(main_svg, inset_svg, main_boxes, width, box_w, heading):
    """Put inset_svg in a framed box at the top right of main_svg; push the panel down if they would overlap."""
    import re as _re
    ih = int(float(_re.search(r'height="([\d.]+)"', inset_svg.split("\n", 1)[0]).group(1)))
    trim = 52                                             # the nested figure's empty title/note band
    bw, bh = box_w, ih - trim + 22
    bx, by = width - bw - 12, 52                          # below the title and the note line
    clash = any(x1 > bx and x0 < bx + bw and y1 > by and y0 < by + bh for x0, x1, y0, y1 in main_boxes)
    lines = main_svg.split("\n")
    head, body, tail = lines[:4], lines[4:-1], lines[-1:]
    H = float(_re.search(r'height="([\d.]+)"', lines[0]).group(1))
    dy = 0
    if clash:
        dy = bh + 12
        head[0] = _re.sub(r'height="[\d.]+"', f'height="{H + dy:.0f}"', head[0], count=1)
        body = [f'<g transform="translate(0,{dy})">'] + body + ["</g>"]
        info(f"inset box would overlap the panel's top rows: panel moved down {dy} px")
    nested = _re.sub(r"^<svg ", f'<svg x="{bx + 1}" y="{by + 18 - trim}" ', inset_svg, count=1)
    # drop the nested figure's empty title/note band (it starts drawing tracks at y=70)
    nested = nested.replace('<rect width="100%" height="100%" fill="white"/>', "", 1)
    box = [f'<rect x="{bx}" y="{by}" width="{bw}" height="{bh}" fill="white" stroke="#888" stroke-width="1" rx="4"/>',
           f'<text x="{bx + 8}" y="{by + 14}" font-size="11" font-weight="bold" fill="#444">{html_escape(heading)}</text>',
           nested]
    return "\n".join(head + body + box + tail)


def finish(S, a, contigs, queries, bnames, title, note, legend_extra=(), empty_reason=None):
    """engine.finish() with fixed colours per group and a reason for empty species.

    The engine colours groups by PALETTE[i] in sorted order and labels empty species
    'nothing found'. Both are set here by wrapping the engine's own functions for one
    call, so lib_synteny_engine.py stays a verbatim copy of the old 4d script. Long empty
    stretches within a contig are shortened and marked // (compress_gaps, --max-gap-kb).
    Species in --inset are drawn in a box at the top right, not in the tree-ordered panel.
    """
    empty_reason = empty_reason or {}
    inset = [q for q in getattr(a, "inset", "").split(",") if q]
    for q in inset:
        if q not in queries:
            warn(f"--inset {q}: not a species of this figure; ignored")
    inset = [q for q in inset if q in queries]
    main_q = [q for q in queries if q not in inset]
    orig_render, orig_build = S.render, S.build_tracks

    def build(*args, **kw):
        tracks = orig_build(*args, **kw)
        compress_gaps(tracks, getattr(a, "max_gap_kb", 1000.0), getattr(a, "break_kb", 20.0))
        for t in tracks:
            if not t["feats"] and t["query"] in empty_reason:
                t["why"] = empty_reason[t["query"]]
        return tracks

    def palette(tracks):
        groups = sorted({f["group"] for t in tracks for f in t["feats"] if f["kind"] == "exendin"})
        cm = colour_map(groups)
        S.PALETTE = [cm[g] for g in groups] or S.PALETTE

    def render(tracks, tree, **kw):
        palette(tracks)
        svg, lay, span = orig_render(tracks, tree, **kw)
        svg = _draw_breaks(svg, tracks, tree, kw.get("width", 1400), span)
        if inset:
            it = build(contigs, inset, bnames, a.reference_query, a.clip_kb)
            rank = {q: i for i, q in enumerate(inset)}
            it.sort(key=lambda t: (rank[t["query"]], not t["main"], t["chrom"]))
            palette(it)
            ikw = dict(kw, width=a.inset_width, legend_extra=())
            isvg, ilay, ispan = orig_render(it, None, title="", note="", **{k: v for k, v in ikw.items()
                                                                            if k not in ("title", "note")})
            isvg = _draw_breaks(isvg, it, None, a.inset_width, ispan)
            svg = _compose_inset(svg, isvg, _track_boxes(tracks, tree, kw.get("width", 1400), span),
                                 kw.get("width", 1400), a.inset_width, a.inset_title)
            lay = lay + ilay
            info(f"inset: {', '.join(inset)} ({len(it)} track(s), own scale {ispan / 1000:.0f} kb)")
        return svg, lay, span

    S.render, S.build_tracks = render, build
    try:
        S.finish(a, contigs, main_q, bnames, title, note, legend_extra=legend_extra)
    finally:
        S.render, S.build_tracks = orig_render, orig_build


# ---------------- TOGA2 names (steps 3 and 4) ----------------
def toga2_reference_groups(reference_table, name_map=None):
    """reference_exendins.tsv (3a) -> {gene_id: group}, from each gene's annotation text.

    name_map: optional TSV gene_id<TAB>group that overrides the name of a reference gene.
    """
    override = {}
    if name_map:
        for line in open(name_map):
            f = line.rstrip("\n").split("\t")
            if len(f) >= 2 and not line.startswith("#") and f[0] != "gene_id":
                override[f[0]] = f[1]
    out, text = {}, {}
    with open(reference_table) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            text[r["gene_id"]] = r.get("matched_text", "") or r.get("label", "")
            out[r["gene_id"]] = override.get(r["gene_id"]) or name_group(text[r["gene_id"]])
    for g in sorted(out):
        info(f"TOGA2 reference {g}: '{text[g][:70]}' -> {out[g]}")
    return out


def toga2_tip_group(ref_genes, ref_group):
    """Group of a 3b tip from its ref_genes (comma-separated); many2many -> 'exendin-2/4'."""
    gs = sorted({ref_group.get(g, g) for g in ref_genes.split(",") if g})
    if len(gs) == 1:
        return gs[0]
    return name_group(" ".join(gs)) if all(x.startswith("exendin") for x in gs) else " / ".join(gs)


# ---------------- species table and genome files ----------------
TABLE_COLS = ["code", "species", "accession", "toga2_query", "tree_tip", "source", "note"]


def read_species_table(path):
    """species.tsv from 1a_make_species_table.py -> list of dicts (comment lines skipped)."""
    with open(path) as fh:
        rows = [r for r in csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t")]
    for r in rows:
        for c in TABLE_COLS:
            r.setdefault(c, "")
    return [r for r in rows if r["code"]]


def find_genome_files(genomedir, code, accession):
    """{'fasta','gff','protein'} under ${GENOMEDIR}/<code>_<accession>/ncbi_dataset/data/<acc>/ (None if absent)."""
    base = os.path.join(genomedir, f"{code}_{accession}", "ncbi_dataset", "data")

    def first(pattern):
        hits = sorted(glob.glob(os.path.join(base, "*", pattern)))
        return hits[0] if hits else None
    return dict(fasta=first("*_genomic.fna") or first("*.fna"),
                gff=first("genomic.gff") or first("genomic.gff.gz"),
                protein=first("protein.faa") or first("protein.faa.gz"))


def read_fasta(path):
    import gzip
    op = gzip.open(path, "rt") if path.endswith(".gz") else open(path)
    out, name, head, buf = {}, None, "", []
    with op as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name:
                    out[name] = (head, "".join(buf))
                head = line[1:].strip()
                name, buf = head.split()[0] if head else "", []
            elif line:
                buf.append(line.strip())
    if name:
        out[name] = (head, "".join(buf))
    return out


def warn(msg):
    print(f"[warn] {msg}", file=sys.stderr)


def info(msg):
    print(f"[info] {msg}", file=sys.stderr)
