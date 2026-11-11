#!/usr/bin/env python3
"""
Every TOGA2 gene in the SLC4A8-AQP5 interval, in every query of one TOGA2 run, drawn as
the publication figure of step 2 (lib_pub_synteny.py). One run per call: target_HelSus,
target_HemCap, the AnoCar and CroPor runs, or the integrated annotation.

Where the files are: --dirs is a path with {q} where the query code goes, e.g.
  ${TOGADIR}/target_HelSus/query_{q}                        (target_* runs)
  ${TOGADIR}/query_{q}/target_AnoCar/output                 (per-query runs)
  ${TOGADIR}/integrated_annotations/query_{q}/integrated    (integration)
Queries are every directory the pattern matches (or --queries).

Genes. A gene is a TOGA2 query gene: query_genes.tsv (query_gene -> projections) and
query_genes.bed (its locus). If those are missing, projections in query_annotation.bed
are grouped by reference gene and overlap instead. Projection names are
transcript#gene#chain[#paralog|#retro]; the reference gene of each projection is named
from --ref-gff (the annotation(s) the run projected: symbol from gene=/Name=, description
from description=/Note=/product=), or from --name-map (ref_gene<TAB>symbol).
  colour   exendin-family genes (exendin / GCGL / proglucagon by name): the fixed family
           colours of every figure; other named genes: one pale colour per symbol, the
           same in every species; genes without a usable name: 'uncharacterised'
  filled   best projection status FI or I (loss_summary.tsv); hollow = PI, UL, M, L ...
  dashed   only paralogous or retrogene projections (#paralog, #retro, retro_ genes)

Boundaries. The query genes named SLC4A8 / AQP5 (symbol, or 2a's description regexes;
NAME=REFID|REFID adds reference gene IDs, e.g. AQP5=LOC_00012326 for the HelSus EviAnn
annotation). As in 2a: the pair on one contig within --max-interval-kb with the most
family genes between; otherwise each boundary on its own (split genome).

Read-only on the TOGA2 directory. Writes the figure files and --genes-out.

Usage
  python3 3c_plot_toga2_synteny.py --dirs "${TOGADIR}/target_HelSus/query_{q}" \
      --ref-gff ${HELSUS_EVIANN_GFF} --boundaries SLC4A8,AQP5=LOC_00012326 \
      --species-table ../species_toga2.tsv --species-tree ${SPTREE} \
      --svg interval_HelSus.svg --pdf interval_HelSus.pdf --png interval_HelSus.png \
      --genes-out interval_HelSus_genes.tsv 2> interval_HelSus.log
"""
import argparse, glob, gzip, os, re, sys
from collections import defaultdict, Counter
from urllib.parse import unquote

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import lib_exendin_names as N  # noqa: E402
import lib_pub_synteny as P  # noqa: E402

RANK = {"FI": 0, "I": 1, "PI": 2, "UL": 3, "M": 4, "L": 5, "PG": 6, "N": 7}
DEFAULT_RX = {"SLC4A8": r"solute carrier family 4 member 8", "AQP5": r"aquaporin[- ]?5"}
ID_LIKE = re.compile(r"^(LOC_?\d+|reg_\d+|retro_|gene-|rna-|[XN][MR]_\d+|ENS[A-Z]*\d+|g\d+$)", re.I)


def opener(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)


def first(d, *names):
    for n in names:
        for c in (n, n + ".gz"):
            if os.path.exists(os.path.join(d, c)):
                return os.path.join(d, c)
    return None


# ---------------- reference names ----------------
def read_ref_names(gffs, name_map):
    """key (gene ID, Name, transcript ID, with and without gene-/rna- prefix and version) -> (symbol, description)."""
    names = {}
    for gff in gffs:
        genes, tx = {}, []
        with opener(gff) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                t = line.rstrip("\n").split("\t")
                if len(t) < 9:
                    continue
                at = {k: unquote(v) for k, v in (kv.split("=", 1) for kv in t[8].strip().strip(";").split(";")
                                                  if "=" in kv)}
                if t[2] in ("gene", "pseudogene"):
                    sym = at.get("gene") or at.get("Name") or ""
                    if sym == at.get("ID") or ID_LIKE.match(sym):
                        sym = ""
                    desc = at.get("description") or at.get("Note") or at.get("product") or ""
                    genes[at.get("ID", "")] = [sym, desc]
                    for k in ("ID", "Name", "gene", "locus_tag", "gene_id"):
                        if at.get(k):
                            names.setdefault(at[k], genes[at.get("ID", "")])
                elif t[2] in ("mRNA", "transcript", "lnc_RNA", "ncRNA"):
                    tx.append(at)
        for at in tx:
            g = genes.get(at.get("Parent", "").split(",")[0])
            if g is None:
                continue
            if not g[1] and (at.get("product") or at.get("Note")):
                g[1] = at.get("product") or at.get("Note")
            for k in ("ID", "Name", "transcript_id"):
                if at.get(k):
                    names.setdefault(at[k], g)
        N.info(f"{gff}: {len(genes)} reference genes named")
    for k in list(names):
        for kk in (re.sub(r"^(gene|rna)-", "", k), re.sub(r"\.\d+$", "", re.sub(r"^(gene|rna)-", "", k))):
            names.setdefault(kk, names[k])
    if name_map:
        for line in open(name_map):
            f = line.rstrip("\n").split("\t")
            if len(f) >= 2 and not line.startswith("#"):
                names[f[0]] = [f[1], f[2] if len(f) > 2 else ""]
    return names


def parse_projection(p):
    """-> (transcript, reference gene, flag)."""
    parts = p.split("#")
    flag = parts[-1] if parts[-1] in ("paralog", "retro", "processed_pseudogene") else ""
    if flag:
        parts = parts[:-1]
    tx = parts[0]
    gene = parts[1] if len(parts) >= 3 else re.sub(r"-mRNA.*$|-T\d+$", "", tx)
    return tx, gene, flag


def name_of(tx, gene, names):
    for k in (gene, tx, re.sub(r"\.\d+$", "", tx), re.sub(r"-mRNA.*$", "", tx)):
        if k in names:
            return names[k], True
    if gene and not ID_LIKE.match(gene) and not gene.isdigit():
        return [gene, ""], True                      # the run already uses symbols
    return ["", ""], False


# ---------------- one query ----------------
def read_status(d):
    p = first(d, "loss_summary.tsv")
    out = {}
    if p:
        with opener(p) as fh:
            for line in fh:
                t = line.rstrip("\n").split("\t")
                if len(t) >= 3 and t[0] == "PROJECTION":
                    out[t[1]] = t[2]
    return out


def read_genes(d):
    """{query_gene: dict(contig, start, end, strand, projections)} from query_genes.*, else from the annotation."""
    gbed, gtsv = first(d, "query_genes.bed"), first(d, "query_genes.tsv")
    genes = {}
    if gbed and gtsv:
        with opener(gbed) as fh:
            for line in fh:
                t = line.rstrip("\n").split("\t")
                if len(t) >= 6 and not line.startswith(("#", "track")):
                    genes[t[3]] = dict(contig=t[0], start=int(t[1]) + 1, end=int(t[2]), strand=t[5], proj=[])
        with opener(gtsv) as fh:
            for line in fh:
                t = line.rstrip("\n").split("\t")
                if len(t) >= 2 and t[0] in genes:
                    genes[t[0]]["proj"].append(t[1])
        return genes, "query_genes"
    abed = first(d, "query_annotation.bed")
    if not abed:
        return None, "no query_genes.* or query_annotation.bed"
    by = defaultdict(list)
    with opener(abed) as fh:
        for line in fh:
            t = line.rstrip("\n").split("\t")
            if len(t) >= 6 and not line.startswith(("#", "track")):
                by[(t[0], t[5], parse_projection(t[3])[1])].append((int(t[1]) + 1, int(t[2]), t[3]))
    for (c, s, g), ps in by.items():
        ps.sort()
        cur = None
        for st, en, p in ps:
            if cur and st <= cur["end"]:
                cur["end"] = max(cur["end"], en); cur["proj"].append(p)
            else:
                cur = dict(contig=c, start=st, end=en, strand=s, proj=[p])
                genes[f"{g}@{c}:{st}"] = cur
    return genes, "query_annotation (no query_genes.*)"


def describe(qg, g, names, status, unresolved):
    syms, descs, flags, best = [], [], [], 99
    for p in g["proj"]:
        tx, gene, flag = parse_projection(p)
        (sym, desc), ok = name_of(tx, gene, names)
        if not ok:
            unresolved[gene] += 1
        syms.append(sym); descs.append(desc); flags.append(flag)
        best = min(best, RANK.get(status.get(p, ""), 99))
    if not g["proj"]:
        flags = [""]
    sym = "/".join(sorted({s for s in syms if s})) or ""
    text = " | ".join([sym] + sorted({d for d in descs if d}))
    grp = N.name_group(text) if text.strip(" |") else "other"
    family = N.is_family(grp) or grp == "proglucagon"
    if family:
        group, label = grp, N.group_abbrev(grp)
    elif sym:
        group = label = sym.split("/")[0].upper() + ("+" if "/" in sym else "")
    elif descs and any(descs) and not re.search(r"uncharacteri[sz]ed|hypothetical|^protein$", descs[0] or "", re.I):
        group = label = re.sub(r"\s*\[.*?\]|,? ?transcript variant.*$", "", sorted(d for d in descs if d)[0])[:28]
    else:
        group, label = P.UNCHAR, ""
    dashed = qg.startswith("retro_") or (bool(flags) and all(f in ("paralog", "retro") for f in flags))
    st = {v: k for k, v in RANK.items()}.get(best, "NA")
    return dict(contig=g["contig"], start=g["start"], end=g["end"], strand=g["strand"], symbol=sym or qg,
                description="; ".join(sorted({d for d in descs if d}))[:120], group=group, label=label,
                family=family, filled=best <= 1 or best == 99, dashed=dashed, status=st, query_gene=qg,
                ref_genes=",".join(sorted({parse_projection(p)[1] for p in g["proj"]})),
                hover=f"{qg} | {sym or 'unnamed'} | {st} | {g['contig']}:{g['start']}-{g['end']}")


def choose_boundaries(q, feats, bspec, max_interval, near):
    cands = {}
    for name, ids, rx in bspec:
        cands[name] = []
        for f in feats:
            refs = set(f["ref_genes"].split(","))
            exact = any(s.upper() == name.upper() for s in f["symbol"].split("/")) or bool(refs & ids)
            if exact or rx.search(f["description"]):
                cands[name].append(dict(f, exact=exact))
    fam = [f for f in feats if f["family"]]
    (n0, *_), (n1, *_) = bspec

    def between(c, lo, hi):
        return sum(1 for g in fam if g["contig"] == c and g["end"] >= lo and g["start"] <= hi)

    pairs = []
    for s in cands[n0]:
        for t in cands[n1]:
            if s["contig"] == t["contig"]:
                lo, hi = min(s["start"], t["start"]), max(s["end"], t["end"])
                if hi - lo <= max_interval:
                    pairs.append(((-between(s["contig"], lo, hi), not (s["exact"] and t["exact"]), hi - lo), s, t))
    if pairs:
        _, s, t = min(pairs, key=lambda p: p[0])
        N.info(f"{q}: {n0} {s['query_gene']} and {n1} {t['query_gene']} on {s['contig']}, "
               f"{abs(t['start'] - s['start']) / 1000:.0f} kb apart")
        return {n0: s, n1: t}
    out = {}
    for name, *_ in bspec:
        if cands[name]:
            b = min(cands[name], key=lambda b: (-sum(1 for g in fam if g["contig"] == b["contig"] and
                                                      min(abs(g["start"] - b["end"]), abs(b["start"] - g["end"])) <= near),
                                                not b["exact"]))
            out[name] = b
            N.info(f"{q}: {name} = {b['query_gene']} on {b['contig']} (no pair on one contig)")
        else:
            N.info(f"{q}: no {name} gene")
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--dirs", required=True, help="query directory with {q} for the query code")
    ap.add_argument("--queries", help="comma-separated codes (default: every match of --dirs)")
    ap.add_argument("--ref-gff", action="append", default=[], help="reference annotation(s) the run projected")
    ap.add_argument("--name-map", help="TSV ref_gene<TAB>symbol[<TAB>description], overrides --ref-gff")
    ap.add_argument("--boundaries", default="SLC4A8,AQP5",
                    help="left,right; NAME=REFID|REFID adds reference gene IDs (AQP5=LOC_00012326)")
    ap.add_argument("--max-interval-kb", type=float, default=5000.0, help="pairs further apart are not a pair")
    ap.add_argument("--genes-out", help="TSV of the boundary genes, every gene between or near them, and every family gene")
    P.add_args(ap, near_default=0.0, split_default="all")
    a = ap.parse_args()

    if "{q}" not in a.dirs:
        sys.exit("[error] --dirs needs {q} where the query code goes")
    if a.queries:
        queries = a.queries.split(",")
    else:
        rx = re.compile("^" + re.escape(a.dirs).replace(re.escape("{q}"), "([^/]+)") + "/?$")
        queries = sorted(m.group(1) for m in (rx.match(p.rstrip("/")) for p in glob.glob(a.dirs.replace("{q}", "*")))
                         if m and os.path.isdir(m.group(0)))
    if not queries:
        sys.exit(f"[error] no directory matches {a.dirs}")
    N.info(f"{len(queries)} queries: {', '.join(queries)}")

    bspec, bnames = [], []
    for item in a.boundaries.split(","):
        name, _, ids = item.partition("=")
        bnames.append(name)
        bspec.append((name, set(i for i in ids.split("|") if i),
                      re.compile(rf"(?<![A-Za-z0-9])(?:{DEFAULT_RX.get(name, re.escape(name))})(?![A-Za-z0-9])", re.I)))
    names = read_ref_names(a.ref_gff, a.name_map)
    near_sel = (a.near_kb or 1000.0) * 1000

    species, rows = {}, []
    for q in queries:
        d = a.dirs.replace("{q}", q)
        genes, how = read_genes(d)
        if genes is None:
            N.warn(f"{q}: {how} in {d}")
            species[q] = dict(bnd={}, genes=[], empty="no TOGA2 annotation found")
            continue
        status = read_status(d)
        if not status:
            N.warn(f"{q}: no loss_summary.tsv; every gene drawn filled")
        unresolved = Counter()
        feats = [describe(qg, g, names, status, unresolved) for qg, g in genes.items()]
        if unresolved:
            N.warn(f"{q}: {sum(unresolved.values())} projections of {len(unresolved)} reference genes have no name "
                   f"(e.g. {', '.join(list(unresolved)[:3])}); give --ref-gff or --name-map")
        bnd = choose_boundaries(q, feats, bspec, a.max_interval_kb * 1000, near_sel)
        bq = {b["query_gene"] for b in bnd.values()}
        win = defaultdict(list)                       # contig -> [(lo, hi)] around the boundaries
        bl = list(bnd.values())
        if len(bl) == 2 and bl[0]["contig"] == bl[1]["contig"]:
            win[bl[0]["contig"]].append((min(b["start"] for b in bl), max(b["end"] for b in bl)))
        else:
            for b in bl:
                win[b["contig"]].append((b["start"] - near_sel, b["end"] + near_sel))
        rest = [f for f in feats if f["query_gene"] not in bq and
                (f["family"] or any(f["end"] >= lo and f["start"] <= hi for lo, hi in win.get(f["contig"], [])))]
        species[q] = dict(bnd=bnd, genes=rest, empty="")
        N.info(f"{q}: {len(genes)} genes ({how}); {sum(f['family'] for f in feats)} family")
        for f in rest + list(bnd.values()):
            rows.append([q, f["query_gene"], f["symbol"], f["group"], f["contig"], f["start"], f["end"], f["strand"],
                         f["status"], "boundary" if f["query_gene"] in bq else "", f["ref_genes"]])

    if a.genes_out:
        with open(a.genes_out, "w") as fh:
            fh.write("species\tquery_gene\tsymbol\tgroup\tcontig\tstart\tend\tstrand\tstatus\tboundary\tref_genes\n")
            for r in sorted(rows, key=lambda r: (r[0], r[4], r[5])):
                fh.write("\t".join(map(str, r)) + "\n")
    P.draw(species, a, bnames,
           legend_status=[("not intact", False, False), ("paralog / retrogene", False, True)],
           empty_default="no TOGA2 annotation")


if __name__ == "__main__":
    main()
