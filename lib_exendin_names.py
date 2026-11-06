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


def finish(S, a, contigs, queries, bnames, title, note, legend_extra=(), empty_reason=None):
    """engine.finish() with fixed colours per group and a reason for empty species.

    The engine colours groups by PALETTE[i] in sorted order and labels empty species
    'nothing found'. Both are set here by wrapping the engine's own functions for one
    call, so lib_synteny_engine.py stays a verbatim copy of the old 4d script.
    """
    empty_reason = empty_reason or {}
    orig_render, orig_build = S.render, S.build_tracks

    def render(tracks, tree, **kw):
        groups = sorted({f["group"] for t in tracks for f in t["feats"] if f["kind"] == "exendin"})
        cm = colour_map(groups)
        S.PALETTE = [cm[g] for g in groups] or S.PALETTE
        return orig_render(tracks, tree, **kw)

    def build_tracks(*args, **kw):
        tracks = orig_build(*args, **kw)
        for t in tracks:
            if not t["feats"] and t["query"] in empty_reason:
                t["why"] = empty_reason[t["query"]]
        return tracks

    S.render, S.build_tracks = render, build_tracks
    try:
        S.finish(a, contigs, queries, bnames, title, note, legend_extra=legend_extra)
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
