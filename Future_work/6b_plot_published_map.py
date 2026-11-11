#!/usr/bin/env python3
"""
Figure of where the published Heloderma venom peptides are encoded in one genome (from 6a).

  A   Genomic context. Top: the whole contig (chromosome) with the region marked. Below:
      the region from SLC4A8 to AQP5 (+/- --context-pad-kb), the two boundary genes taken
      from the TOGA2 reference annotation (--toga2-gff, the EviAnn GFF that TOGA2 projects
      from; --boundaries as in 3c). Tracks:
        - TOGA2 genes: every gene of that annotation in the region (exendin-family genes in
          their name colour, genes carrying a mapped peptide outlined);
        - Reciprocal BLAST: the step-4 curated loci of this species (--curated), labelled
          and coloured by their step-5 names if --names is given (EX1-EX4 in the exendin
          colours, GCGL green, positional letters grey);
        - NCBI genes only with --ncbi (NCBI's exendin models are unreliable).
      A marker under each mapped peptide is joined to its pop-out.
      For a genome other than HelSus (step 7, HelCha), --toga2-run draws TOGA2's
      PROJECTIONS of the reference genes onto that genome (query_<code>/query_annotation.bed),
      each named after its HelSus reference gene with its TOGA2 status; the boundaries are
      the projections of SLC4A8 and AQP5 nearest the mappings.
      TOGA2 contig names (chr2) are translated to the genome's (NC_139459.1) with the NCBI
      sequence report, or --contig-map chr2=NC_139459.1. The two must be the same assembly:
      contig lengths in the TOGA2 GFF's ##sequence-region lines are checked against the
      genome, and the log lists the TOGA2 gene each peptide falls in (a published exendin
      landing in no TOGA2 gene means the coordinates do not match).
  B-  One pop-out per published exendin (exendin-1, -2, -4 ... in the order of the peptide
      file). Each has
        - the locus: TOGA2 (and with --ncbi, NCBI) gene models (thin = exon, thick = CDS),
          the curated reciprocal-BLAST loci, and, for each peptide, the
          genomic blocks its residues map to. Minus-strand loci are drawn reversed so the
          peptide reads left to right in every panel;
        - the alignment: peptide (top) over the genome's translation of the codons it maps
          to (bottom), one block per exon, the intron length between blocks; in grey, the
          genome residues just before residue 1 (processing site) and just after the last
          residue;
        - identity, introns, the TOGA2 gene and exon and the reciprocal-BLAST locus each
          peptide falls in, and the runner-up locus.
      When a pop-out has more than one peptide (helospectin I and II), the first (longest)
      is the reference and the others are aligned to it: residues encoded by the SAME
      codons in both are drawn in the exendin colour, residues that map to different codons,
      or exist in only one peptide, in magenta, in the locus blocks and in the alignment.

Usage
  python3 6b_plot_published_map.py --map-dir map --out published_map \
      --toga2-gff ${HELSUS_EVIANN_GFF} --boundaries SLC4A8,AQP5=LOC_00012326 \
      --curated ../reciprocal_blast/reciprocal/curated_loci.tsv --names ../nomenclature/names_names.tsv
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
    ap.add_argument("--toga2-gff", help="the HelSus annotation TOGA2 projects from (EviAnn GFF3, ${HELSUS_EVIANN_GFF}): "
                    "drawn as the gene track and used for the boundaries")
    ap.add_argument("--toga2-run", help="${TOGADIR}/target_HelSus: draw the TOGA2 PROJECTIONS of a query genome "
                    "(query_<code>/query_annotation.bed) instead of the reference annotation; for genomes other "
                    "than HelSus, e.g. HelCha. --toga2-gff then only names the reference genes")
    ap.add_argument("--toga2-query", help="query name in --toga2-run [the 6a genome code]")
    ap.add_argument("--boundaries", default="SLC4A8,AQP5=LOC_00012326",
                    help="boundary genes in --toga2-gff, NAME or NAME=ID, as in 3c [SLC4A8,AQP5=LOC_00012326]")
    ap.add_argument("--contig-map", action="append", default=[], metavar="TOGA2=NCBI",
                    help="TOGA2 contig name -> genome contig (e.g. chr2=NC_139459.1). Default: from the NCBI "
                         "sequence report (chromosome names), else the names as they are")
    ap.add_argument("--curated", help="step 4 curated loci (reciprocal_blast/reciprocal/curated_loci.tsv): a track of "
                    "the curated genes")
    ap.add_argument("--names", help="step 5 names (nomenclature/names_names.tsv): label curated genes by their new names")
    ap.add_argument("--ncbi", action="store_true", help="also draw the NCBI genes (off by default: NCBI exendin models "
                    "are unreliable)")
    ap.add_argument("--force", action="store_true", help="draw even if the TOGA2 coordinates fail the assembly check")
    ap.add_argument("--gff", help="override the NCBI GFF recorded by 6a")
    ap.add_argument("--boundary-genes", default="SLC4A8,AQP5", help="NCBI boundary names, used only without --toga2-gff")
    ap.add_argument("--context-max-kb", type=float, default=3000.0,
                    help="without --toga2-gff: reach out to NCBI boundary genes within this distance [3000]")
    ap.add_argument("--context-pad-kb", type=float, default=5.0, help="padding outside the boundary genes [5]")
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

    # ---------------- gene tracks ----------------
    import json
    fai_names = set(clen)
    rep_path = os.path.join(os.path.dirname(fasta), "sequence_report.jsonl") if fasta else ""
    to_ncbi = {}
    if rep_path and os.path.exists(rep_path):
        for line in open(rep_path):
            d = json.loads(line)
            acc = d.get("refseqAccession") if d.get("refseqAccession") in fai_names else d.get("genbankAccession")
            if not acc:
                continue
            for k in ("chrName", "sequenceName", "genbankAccession", "refseqAccession", "ucscStyleName"):
                v = d.get(k)
                if v:
                    to_ncbi[v] = acc
                    if k == "chrName":
                        to_ncbi["chr" + v] = acc; to_ncbi["Chr" + v] = acc
    for cm in a.contig_map:
        t, _, n_ = cm.partition("=")
        to_ncbi[t.strip()] = n_.strip()
    to_toga = defaultdict(list)
    for t, n_ in to_ncbi.items():
        to_toga[n_].append(t)

    def ncbi_name(t):
        """TOGA2 contig name -> genome contig: explicit map / sequence report first, then the
        renamings TOGA2 makes ('JANEZZ010002025__1' for 'JANEZZ010002025.1')."""
        if t in to_ncbi:
            return to_ncbi[t]
        if t in fai_names or not fai_names:
            return t
        import re as _re
        for c_ in (t.replace("__", "."), _re.sub(r"_(\d+)$", r".\1", t), t + ".1"):
            if c_ in fai_names:
                return c_
        return t

    contigs = sorted({m["contig"] for m in maps}, key=lambda c: min(int(m["start"]) for m in maps if m["contig"] == c))
    tracks = []                                        # (label, {ncbi contig: [genes]})
    bnd = {}                                           # boundary name -> (ncbi contig, start, end)
    if a.toga2_run:
        # TOGA2 projections onto this (query) genome
        if not a.toga2_gff:
            sys.exit("[error] --toga2-run needs --toga2-gff (the reference GFF, for gene names and boundary IDs)")
        q = a.toga2_query or info.get("code", "")
        qdir = os.path.join(a.toga2_run, "query_" + q)
        if not os.path.isdir(qdir):
            sys.exit(f"[error] no {qdir} (give --toga2-query)")
        notes = G.ref_gene_notes(a.toga2_gff)
        raw = G.load_toga2_query(qdir, notes)
        tg = defaultdict(list)
        for tc, gs in raw.items():
            tg[ncbi_name(tc)].extend(gs)
        miss = [tc for tc in raw if ncbi_name(tc) not in clen]
        if clen and miss and len(miss) == len(raw):
            sys.exit(f"[error] none of the TOGA2 contig names of {q} (e.g. {miss[0]}) is in the genome FASTA; "
                     "give --contig-map TOGA2=GENOME")
        N.info(f"TOGA2 {q}: {sum(len(v) for v in tg.values())} projected genes on {len(tg)} contigs")
        specs = [x.strip() for x in a.boundaries.split(",") if x.strip()]
        found = G.find_genes(a.toga2_gff, [x for x in specs if "=" not in x])
        mids = {c: sum((int(m["start"]) + int(m["end"])) / 2 for m in maps if m["contig"] == c) /
                max(1, sum(1 for m in maps if m["contig"] == c)) for c in contigs}
        for sp in specs:
            n_, _, gid = sp.partition("=")
            gid = gid or (found[n_][4] if n_ in found else "")
            gid = gid.split(":")[-1]
            if not gid:
                N.warn(f"boundary {n_} not found in {a.toga2_gff}"); continue
            cand = [(c, g_) for c, gs in tg.items() for g_ in gs if g_["name"] == gid]
            if not cand:
                N.warn(f"boundary {n_} ({gid}): no TOGA2 projection in {q}"); continue
            c, g_ = min(cand, key=lambda x: (x[0] not in mids, abs((x[1]["start"] + x[1]["end"]) / 2
                                                                    - mids.get(x[0], 0))))
            bnd[n_] = (c, g_["start"], g_["end"])
            N.info(f"TOGA2 boundary {n_}: {gid} projected to {c}:{g_['start']}-{g_['end']} ({g_['status']})"
                   + (f"; {len(cand) - 1} other projection(s)" if len(cand) > 1 else ""))
        tracks.append((f"TOGA2 genes ({q})", tg))
    elif a.toga2_gff:
        specs = [x.strip() for x in a.boundaries.split(",") if x.strip()]
        found = G.find_genes(a.toga2_gff, specs)
        for sp in specs:
            n_ = sp.split("=")[0]
            if n_ in found:
                tc, st, en, strand_, gid = found[n_]
                bnd[n_] = (ncbi_name(tc), st, en)
                N.info(f"TOGA2 boundary {n_}: {gid} {tc}:{st}-{en} ({strand_}) -> {ncbi_name(tc)}")
            else:
                N.warn(f"TOGA2 boundary {sp} not found in {a.toga2_gff}")
        # assembly check: TOGA2 contig lengths (##sequence-region) against the genome
        sr = G.sequence_regions(a.toga2_gff)
        bad = [(t, l, clen.get(ncbi_name(t))) for t, l in sr.items()
               if ncbi_name(t) in clen and clen[ncbi_name(t)] != l]
        if bad and not a.force:
            sys.exit("[error] the TOGA2 annotation and the genome differ in contig length, so their coordinates "
                     "are not interchangeable: " + ", ".join(f"{t} {l:,} vs {ncbi_name(t)} {g_:,}" for t, l, g_ in bad[:3])
                     + " (use --force to draw anyway)")
        tg = {}
        need = set(contigs) | {v[0] for v in bnd.values()}
        for c in need:
            tg[c] = []
            for t in to_toga.get(c, []) + [c]:
                tg[c] = G.load_region(a.toga2_gff, t, 1, 10 ** 12)
                if tg[c]:
                    break
            if not tg[c]:
                N.warn(f"no TOGA2 genes on {c} (TOGA2 names tried: {', '.join(to_toga.get(c, []) + [c])}); "
                       "give --contig-map TOGA2=NCBI")
        tracks.append(("TOGA2 genes", tg))
    if a.ncbi or not (a.toga2_gff or a.toga2_run):
        ng = {}
        for c in contigs + [v[0] for v in bnd.values()]:
            ng[c] = G.load_region(gff, c, 1, 10 ** 12) if gff and os.path.exists(gff) else []
        tracks.append(("NCBI genes", ng))
    if not (a.toga2_gff or a.toga2_run):              # NCBI boundaries, as before
        bnames = [b.strip().upper() for b in a.boundary_genes.split(",") if b.strip()]
        for c in contigs:
            ms_ = [m for m in maps if m["contig"] == c]
            lo_, hi_ = min(int(m["start"]) for m in ms_), max(int(m["end"]) for m in ms_)
            for g_ in tracks[-1][1].get(c, []):
                if g_["name"].upper() in bnames and \
                        min(abs(g_["start"] - hi_), abs(g_["end"] - lo_)) <= a.context_max_kb * 1000:
                    bnd.setdefault(g_["name"], (c, g_["start"], g_["end"]))
    bnames_all = {k.upper() for k in bnd}

    # mappings against the TOGA2 genes (consistency of the two coordinate systems)
    if a.toga2_gff or a.toga2_run:
        hits = 0
        for m in maps:
            ov = G.overlaps(tracks[0][1].get(m["contig"], []), int(m["start"]), int(m["end"]))
            hits += bool(ov)
            N.info(f"{m['peptide']} {m['contig']}:{m['start']}-{m['end']} in TOGA2 genes: "
                   + ("; ".join(f"{G.gene_label(g_)}: {w}" for g_, w in ov) or "none"))
        if not hits:
            N.warn("no mapped peptide falls in a TOGA2 gene: check that the TOGA2 annotation is on this assembly")

    # curated reciprocal-BLAST loci (step 4), named by step 5 if given
    code = info.get("code", "")
    newname = {}
    if a.names:
        for r in read_tsv(a.names):
            newname[r["old_name"]] = r["new_name"]
    curated = defaultdict(list)
    if a.curated:
        for r in read_tsv(a.curated):
            if r["species"] != code:
                continue
            nm = newname.get(r["name"], r["name"])
            suf = nm[len(code) + 1:] if nm.startswith(code + "_") else nm
            if a.names:
                import re as _re
                mm = _re.match(r"EX(\d)", suf)
                grp = f"exendin-{mm.group(1)}" if mm else ("GCGL (glucagon-like)" if suf.startswith("GCGL")
                                                           else "no close published match")
            else:
                grp = r["group"]
            curated[r["contig"]].append(dict(name=nm, short=suf, start=int(r["g_start"]), end=int(r["g_end"]),
                                             strand=r["strand"], group=grp, old=r["name"],
                                             step4=r["group"], placement=r.get("placement", "")))
        N.info(f"{sum(len(v) for v in curated.values())} curated {code} loci"
               + (f", named from {a.names}" if a.names else ""))

    def curated_over(c, lo, hi):
        return [k for k in curated.get(c, []) if k["end"] >= lo and k["start"] <= hi]

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
    regions = []                                  # (contig, lo, hi, bounded)
    bctg = {v[0] for v in bnd.values()}
    if len(bnd) >= 2 and len(bctg) == 1:
        c = bctg.pop()
        rlo = min(v[1] for v in bnd.values()) - int(a.context_pad_kb * 1000)
        rhi = max(v[2] for v in bnd.values()) + int(a.context_pad_kb * 1000)
        inside = [m for m in maps if m["contig"] == c]
        out_ = [m for m in inside if int(m["end"]) < rlo or int(m["start"]) > rhi]
        if out_:
            N.warn("mapped outside the boundary genes (region widened): " + ", ".join(m["peptide"] for m in out_))
            rlo = min([rlo] + [int(m["start"]) - 2000 for m in out_])
            rhi = max([rhi] + [int(m["end"]) + 2000 for m in out_])
        regions.append((c, max(1, rlo), rhi, True))
    else:
        if a.toga2_gff:
            N.warn("boundary genes not both found on one contig: panel A is drawn around the mappings"
                   + (" (" + "; ".join(f"{k} on {v[0]}:{v[1]:,}" for k, v in bnd.items()) + ")" if bnd else ""))
        for c in contigs:
            ms_ = [m for m in maps if m["contig"] == c]
            bs_ = [v for v in bnd.values() if v[0] == c]
            regions.append((c, max(1, min([int(m["start"]) for m in ms_] + [v[1] for v in bs_]) - 15000),
                            max([int(m["end"]) for m in ms_] + [v[2] for v in bs_]) + 15000, False))
    elsewhere = [m for m in maps if not any(m["contig"] == r[0] and int(m["end"]) >= r[1] and int(m["start"]) <= r[2]
                                            for r in regions)]
    for c, rlo, rhi, bounded in regions:
        ms = [m for m in maps if m["contig"] == c and int(m["end"]) >= rlo and int(m["start"]) <= rhi]
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
        ytop = yc + 46                            # top of the zoomed region

        def XA(p, rlo=rlo, rhi=rhi):
            return X0 + (p - rlo) / (rhi - rlo) * (X1 - X0)
        rows_y = []
        yr = ytop + 58                            # first gene track axis (labels above it)
        for lab_, src in tracks:
            rows_y.append((lab_, src, yr)); yr += 76
        ycur = yr - 76 + (40 if curated else 0)   # curated track axis
        svg.poly([(hx0, yc + 9), (hx1, yc + 9), (X1, ytop), (X0, ytop)], fill="#F6F2F1", stroke="#C9C2C6",
                 dash="3,3")
        svg.text(X0 - 12, ytop + 14, f"{rlo / 1e6:.3f}–{rhi / 1e6:.3f} Mb", size=10, anchor="end", fill="#666")
        if bounded:
            svg.text(X0 - 12, ytop + 27, "SLC4A8 to AQP5" if {"SLC4A8", "AQP5"} <= bnames_all else
                     " to ".join(bnd), size=9.5, anchor="end", fill="#666", italic=True)
        mapped_ids = set()
        for lab_, src, ya in rows_y:
            svg.text(X0 - 12, ya + 4, lab_, size=10.5, anchor="end", weight="bold", fill="#444")
            svg.line(X0, ya, X1, ya, stroke="#BBB", sw=1)
            rows_end = [-1e9] * 3
            for g in src.get(c, []):
                if g["end"] < rlo or g["start"] > rhi:
                    continue
                carries = any(g["end"] >= int(m["start"]) and g["start"] <= int(m["end"]) for m in ms)
                if carries:
                    mapped_ids.add(g["id"])
                grp = N.name_group(g["description"] + " " + g["name"])
                fam = N.is_family(grp)
                isb = g["name"].upper() in bnames_all or any(g["start"] == v[1] and g["end"] == v[2] and v[0] == c
                                                            for v in bnd.values())
                col = colour(grp) if fam else ("#3A3238" if isb else "#8E878C")
                draw_gene(svg, g, XA, ya, False, col, outline="#1E1A1D" if carries else None, label=False,
                          xlim=(X0, X1))
                xm = (max(XA(g["start"]), X0) + min(XA(g["end"]), X1)) / 2
                bname = next((k for k, v in bnd.items() if v[0] == c and v[1] == g["start"] and v[2] == g["end"]), None)
                lab = g["name"]
                if bname and bname.upper() != g["name"].upper():
                    lab = f"{bname} ({g['name']})"
                elif fam:
                    lab = f"{g['name']} ({N.short_label(g['description']) or grp})"
                elif "uncharacterized" in g["description"].lower():
                    lab += " (uncharacterized)"
                w = text_w(lab, 9.5)
                for k in range(len(rows_end)):
                    if xm - w / 2 > rows_end[k] + 4:
                        rows_end[k] = xm + w / 2
                        ly = ya - 16 - 13 * k
                        svg.line(xm, ya - 9, xm, ly + 3, stroke="#CCC", sw=0.6)
                        svg.text(xm, ly, lab, size=9.5, anchor="middle", fill=col if (fam or isb) else "#444",
                                 weight="bold" if (fam or isb or carries) else "normal")
                        break
        if curated:
            svg.text(X0 - 12, ycur + 4, "Reciprocal BLAST", size=10.5, anchor="end", weight="bold", fill="#444")
            svg.text(X0 - 12, ycur + 17, "(step 4 curated" + (", step 5 names)" if a.names else ")"), size=9,
                     anchor="end", fill="#666")
            svg.line(X0, ycur, X1, ycur, stroke="#BBB", sw=1)
            rows_end = [-1e9] * 2
            for k_ in sorted(curated.get(c, []), key=lambda k_: k_["start"]):
                if k_["end"] < rlo or k_["start"] > rhi:
                    continue
                x0_, x1_ = XA(k_["start"]), XA(k_["end"])
                kc = colour(k_["group"])
                svg.rect(x0_, ycur - 7, max(x1_ - x0_, 4), 14, kc, rx=2)
                xm = (x0_ + x1_) / 2
                w = text_w(k_["short"], 9.5)
                for k in range(len(rows_end)):
                    if xm - w / 2 > rows_end[k] + 4:
                        rows_end[k] = xm + w / 2
                        svg.text(xm, ycur + 21 + 12 * k, k_["short"], size=9.5, anchor="middle", weight="bold",
                                 fill=kc if kc not in ("#A0A0A0", "#C8C8C8") else "#666")
                        break
            ylast = ycur + 34
        else:
            ylast = rows_y[-1][2] if rows_y else ytop + 58
        # scale bar
        span = rhi - rlo
        sb = 10 ** int(len(str(int(span / 5))) - 1)
        sb = sb * 5 if span / (sb * 5) >= 3 else sb
        ya = ylast
        svg.line(X1 - (X1 - X0) * sb / span, ya + 62, X1, ya + 62, stroke="#444", sw=2)
        svg.text(X1, ya + 76, f"{sb / 1000:g} kb", size=10, anchor="end")
        # peptide markers below the last track, one per panel (group) and locus; markers that
        # would overprint (several exendins on one gene) are set side by side
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
    unbounded = [k for k, v in bnd.items() if not any(v[0] == r[0] and r[3] for r in regions)]
    if unbounded and not any(r[3] for r in regions):
        svg.text(X0, y - 4, "Boundary genes (TOGA2): " + "; ".join(
            f"{k} on {v[0]}:{fmt(v[1])}" for k, v in bnd.items())
            + (" (different contigs: the region is drawn around the mappings)" if len(bnd) > 1 else ""),
            size=10.5, fill="#7A5C00")
        y += 16
    if elsewhere:
        svg.text(X0, y - 4, "Mapped outside this region: " + "; ".join(
            f"{m['label']} {m['contig']}:{fmt(m['start'])}" for m in elsewhere), size=10.5, fill="#C0392B")
        y += 12
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
            lo, hi = w["lo"] - a.locus_pad_bp, w["hi"] + a.locus_pad_bp
            for _, src in tracks:                       # show overlapping genes whole if not huge
                for gg in src.get(w["contig"], []):
                    if gg["end"] >= w["lo"] and gg["start"] <= w["hi"] and gg["end"] - gg["start"] <= 40000:
                        lo, hi = min(lo, gg["start"] - 200), max(hi, gg["end"] + 200)
            for k_ in curated_over(w["contig"], w["lo"], w["hi"]):
                if k_["end"] - k_["start"] <= 40000:
                    lo, hi = min(lo, k_["start"] - 200), max(hi, k_["end"] + 200)
            lo = max(1, lo)
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
            # gene models, one row each, per track
            for tlab, src in tracks:
                near = [gg for gg in src.get(w["contig"], []) if gg["end"] >= lo and gg["start"] <= hi]
                s.text(px + 12, yy + 4, tlab.split(" (")[0], size=10, fill="#555")
                for gg in near:
                    grp = N.name_group(gg["description"] + " " + gg["name"])
                    gcol = colour(grp) if N.is_family(grp) else "#8E878C"
                    draw_gene(s, gg, X, yy, flip, gcol, label=True, label_y=yy + 20,
                              max_label=int((LX1 - LX0) / 8), xlim=(LX0, LX1))
                    yy += 34
                if not near:
                    s.text(LX0, yy + 4, "no gene annotated here", size=10, fill="#999", italic=True)
                    yy += 22
            # curated reciprocal-BLAST loci
            if curated:
                kk = [k_ for k_ in curated.get(w["contig"], []) if k_["end"] >= lo and k_["start"] <= hi]
                s.text(px + 12, yy + 4, "Recip. BLAST", size=10, fill="#555")
                clip_open(s, LX0, LX1, yy - 10, yy + 26)
                for k_ in kk:
                    a0_, a1_ = sorted([X(k_["start"]), X(k_["end"])])
                    kc = colour(k_["group"])
                    s.rect(a0_, yy - 6, max(a1_ - a0_, 3), 12, kc, rx=2)
                    s.text((max(a0_, LX0) + min(a1_, LX1)) / 2, yy + 20, k_["name"], size=10, anchor="middle",
                           weight="bold", fill=kc if kc not in ("#A0A0A0", "#C8C8C8") else "#555")
                clip_close(s)
                if not kk:
                    s.text(LX0, yy + 4, "no curated locus here", size=10, fill="#999", italic=True)
                yy += 34
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
                     *[f"   {tlab.split(' genes')[0]} gene: "
                       + ("; ".join(f"{G.gene_label(g_)}: {w_}" for g_, w_ in
                                    G.overlaps(src.get(m["contig"], []), int(m["start"]), int(m["end"])))
                          or "none") for tlab, src in tracks],
                     *([f"   reciprocal BLAST locus: "
                        + ("; ".join(f"{k_['name']} (step 4: {k_['old']}, {k_['step4']})"
                                     for k_ in curated_over(m["contig"], int(m["start"]), int(m["end"])))
                           or "none")] if curated else []),
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
             "around the best loci; one codon per peptide residue). Genes: "
             + ", ".join(("TOGA2 projections from the HelSus reference" if t_.startswith("TOGA2 genes (")
                          else "TOGA2 reference annotation (EviAnn)" if t_.startswith("TOGA2") else "NCBI RefSeq GFF")
                         for t_, _ in tracks)
             + (". Reciprocal BLAST: step 4 curated loci" + (", step 5 names" if a.names else "") if curated else "")
             + ".", size=9.5, fill="#777")
    out_svg = (f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}">\n'
               f'<rect width="{W}" height="{H}" fill="white"/>\n' + "\n".join(svg.el) + "\n</svg>\n")
    open(a.out + ".svg", "w").write(out_svg)
    N.load_engine().write_png(out_svg, a.out + ".png", a.png_width)
    N.info(f"wrote {a.out}.svg")


if __name__ == "__main__":
    main()
