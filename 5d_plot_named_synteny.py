#!/usr/bin/env python3
"""
Synteny plot of the curated step-4 genes under their step-5 names, in the layout and
colours of the step-2/3/4 figures.

  label    the new name without the species prefix: GCGL, A, B ..., EX1-EX4
  colour   5b's class: the published peptide the gene encodes at >= 80% identity (5a);
           grey = no close published match (most copies outside Heloderma)
  filled   intact model
  hollow   pseudogene: a premature stop (in 2f's hits or in the 4e model)
  dashed   partial or flagged model (4e status not ok, no stop)
  boundary SLC4A8 / AQP5 as located by tblastn in 4c/4d

--top GalGal draws chicken as the first row, above the tree-ordered panel, as the
reference row. A curated locus that 5b left out of the names (5b --drop, e.g. chicken's
partial locus GalGal_Ex3, the same gene as the outgroup protein) is drawn under the name
given with --rename OLD=NEW (e.g. GalGal_Ex3=GalGal_GCGL); without one it is drawn with
its old name and a warning.

Usage
  python3 5d_plot_named_synteny.py --species-table species.tsv \
      --loci reciprocal/reciprocal_loci.tsv --boundaries-used reciprocal/boundaries_used.tsv \
      --models models/models.tsv --names names_names.tsv --groups names_groups.tsv \
      --rename GalGal_Ex3=GalGal_GCGL --top GalGal \
      --species-tree ${CACTUSDIR}/Anguimorpha_species_tree_iteration-01_5-23-26.nwk \
      --svg named_synteny.svg --png named_synteny.png --layout named_synteny.tsv
"""
import argparse, csv, os, sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import lib_exendin_names as N  # noqa: E402

S = N.load_engine()


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--species-table", required=True)
    ap.add_argument("--loci", required=True, help="4d reciprocal_loci.tsv")
    ap.add_argument("--boundaries-used", required=True, help="4d boundaries_used.tsv")
    ap.add_argument("--models", help="4e models.tsv (for intact / pseudogene / partial)")
    ap.add_argument("--names", required=True, help="5b _names.tsv")
    ap.add_argument("--groups", required=True, help="5b _groups.tsv (colour class per new name)")
    ap.add_argument("--rename", default="", help="OLD=NEW,... for curated loci 5b did not name")
    ap.add_argument("--confidence", default="high", help="4d confidence levels to draw (high = the curated set)")
    ap.add_argument("--placements", default="interval,near_boundary,unplaced",
                    help="4d placements to draw (add 'elsewhere' to see everything)")
    ap.add_argument("--boundary-names", default="SLC4A8,AQP5", help="left one first")
    ap.add_argument("--title", default="Exendin-family genes between SLC4A8 and AQP5, standard names")
    S.add_common_args(ap)
    N.add_split_arg(ap)
    a = ap.parse_args()

    names = {r["old_name"]: r for r in read_tsv(a.names)}
    groups = {r["tip"]: r["group"] for r in read_tsv(a.groups)}
    models = {r["name"]: r for r in read_tsv(a.models)} if a.models else {}
    extra = dict(x.split("=", 1) for x in a.rename.split(",") if "=" in x)
    conf, keep = set(a.confidence.split(",")), set(a.placements.split(","))
    bnames = a.boundary_names.split(",")
    queries = [r["code"] for r in N.read_species_table(a.species_table)]
    contigs = defaultdict(lambda: defaultdict(list))
    n_drawn, unnamed = 0, []
    for r in read_tsv(a.loci):
        if r["class"] not in ("ortholog", "paralog") or r["confidence"] not in conf:
            continue
        if r["placement"] not in keep:
            continue
        old = r["name"]
        if old in names:
            new = names[old]["new_name"]
        elif old in extra:
            new = extra[old]
        else:
            new = old; unnamed.append(old)
        m = models.get(old, {})
        stop = "premature_stop" in r["state"] or "internal_stop" in m.get("status", "")
        flagged = bool(m) and m.get("status", "ok") != "ok" and not stop
        lab = new.split("_", 1)[1] if "_" in new else new
        nr = names.get(old, {})
        contigs[r["species"]][r["contig"]].append(dict(
            kind="exendin", name=new, label=lab, group=groups.get(new, "not compared"), chrom=r["contig"],
            start=int(r["g_start"]), end=int(r["g_end"]), strand=r["strand"], status=r["class"],
            filled=not stop, dashed=flagged,
            hover=f'{new} (was {old}) | rule {nr.get("rule", "-")} {nr.get("order_basis", "")} | '
                  f'{nr.get("best_match", "-")} {nr.get("identity", "")}% | model {m.get("status", "-")} '
                  f'{m.get("aa_len", "")} aa'))
        n_drawn += 1
    if unnamed:
        N.warn(f"not in {a.names} and no --rename, drawn with their old names: {', '.join(unnamed)}")
    for b in read_tsv(a.boundaries_used):
        contigs[b["species"]][b["contig"]].append(dict(
            kind="boundary", name=b["boundary"], label=b["boundary"], group=b["boundary"],
            chrom=b["contig"], start=int(b["start"]), end=int(b["end"]), strand=b["strand"],
            status="", filled=True, dashed=False,
            hover=f'tblastn {b["query"]} {b["bits"]} bits, query cover {b["query_cov"]}'))
    for q in list(contigs):
        if q not in queries:
            queries.append(q)
    empty = {q: "no curated locus" for q in queries if not contigs.get(q)}
    N.info(f"{n_drawn} genes drawn")
    N.finish(S, a, contigs, queries, bnames, a.title,
             f"Names from step 5 (single copy = GCGL; letters = order from {bnames[0]} to {bnames[1]}, not "
             f"orthology; Heloderma = published peptide matched). Colour = published peptide encoded at >= 80% "
             f"(grey: none). Boundaries located by tblastn; tracks oriented with {bnames[0]} on the left.",
             legend_extra=[("intact", True, False), ("pseudogene (stop)", False, False),
                           ("partial / flagged model", True, True)],
             empty_reason=empty)


if __name__ == "__main__":
    main()
