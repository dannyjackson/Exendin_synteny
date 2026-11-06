#!/usr/bin/env python3
"""
Step-4 synteny plot: the curated reciprocal-BLAST set (4d), drawn exactly like the
step-2 (NCBI) and step-3 (TOGA2) figures, with the same colours per name group.

  colour   name group of the locus (its best reverse hit): exendin-2, exendin-4, GCGL ...
  label    group + 'p' and a number for paralogs (Ex4, Ex2p1 ...), as in curated_loci.tsv
  filled   ortholog (reciprocal best hit for its group)
  hollow   paralog (named after its group, not the best forward hit for it)
  dashed   disabled copy (premature stop, or hormone domain not starting His/Tyr)
  boundary SLC4A8 / AQP5 as found by tblastn (4c/4d), so unannotated genomes have them too

By default only high-confidence loci are drawn; --confidence high,medium adds the
medium ones (their label gets a '?').

Usage
  python3 4f_plot_curated_synteny.py --species-table species.tsv \
      --loci reciprocal/reciprocal_loci.tsv --boundaries-used reciprocal/boundaries_used.tsv \
      --species-tree ${CACTUSDIR}/Anguimorpha_species_tree_iteration-01_5-23-26.nwk \
      --svg curated_synteny.svg --png curated_synteny.png --layout curated_synteny.tsv
"""
import argparse, csv, os, sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import lib_exendin_names as N  # noqa: E402

S = N.load_engine()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--species-table", required=True)
    ap.add_argument("--loci", required=True, help="reciprocal_loci.tsv from 4d")
    ap.add_argument("--boundaries-used", required=True, help="boundaries_used.tsv from 4d")
    ap.add_argument("--confidence", default="high", help="comma-separated: high,medium")
    ap.add_argument("--placements", default="interval,near_boundary,unplaced",
                    help="4d placements to draw (add 'elsewhere' to see everything)")
    ap.add_argument("--boundary-names", default="SLC4A8,AQP5", help="left one first")
    ap.add_argument("--title", default="Exendin-family genes between SLC4A8 and AQP5 "
                                       "(curated reciprocal-BLAST orthologs and paralogs)")
    S.add_common_args(ap)
    a = ap.parse_args()

    conf, keep = set(a.confidence.split(",")), set(a.placements.split(","))
    bnames = a.boundary_names.split(",")
    queries = [r["code"] for r in N.read_species_table(a.species_table)]
    contigs = defaultdict(lambda: defaultdict(list))
    with open(a.loci) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    skipped = defaultdict(int)
    for r in rows:
        if r["class"] not in ("ortholog", "paralog") or r["confidence"] not in conf:
            continue
        if r["placement"] not in keep:
            skipped[r["placement"]] += 1
            continue
        lab = r["name"].split("_", 1)[1] + ("?" if r["confidence"] == "medium" else "")
        contigs[r["species"]][r["contig"]].append(dict(
            kind="exendin", name=r["name"], label=lab, group=r["group"], chrom=r["contig"],
            start=int(r["g_start"]), end=int(r["g_end"]), strand=r["strand"], status=r["class"],
            filled=r["class"] == "ortholog", dashed=r["disabled"] == "yes",
            hover=f'{r["class"]} ({r["confidence"]}) | {r["state"]} | reverse best {r["rev_best_desc"][:40]} '
                  f'{r["rev_best_bits"]} vs non-family {r["rev_best_nonfamily_desc"][:30]} '
                  f'{r["rev_nonfamily_bits"]} (margin {r["margin"]}) | {r["reason"]}'))
    for p, n in skipped.items():
        N.info(f"{n} selected loci not drawn: placement '{p}' (add it to --placements to draw them)")
    with open(a.boundaries_used) as fh:
        for b in csv.DictReader(fh, delimiter="\t"):
            contigs[b["species"]][b["contig"]].append(dict(
                kind="boundary", name=b["boundary"], label=b["boundary"], group=b["boundary"],
                chrom=b["contig"], start=int(b["start"]), end=int(b["end"]), strand=b["strand"],
                status="", filled=True, dashed=False,
                hover=f'tblastn {b["query"]} {b["bits"]} bits, query cover {b["query_cov"]}'))
    for q in contigs:
        if q not in queries:
            queries.append(q)
    empty = {q: "no curated locus" for q in queries if not contigs.get(q)}
    N.finish(S, a, contigs, queries, bnames, a.title,
             f"Reciprocal tblastn/blastx; boundaries located by tblastn. Tracks oriented with {bnames[0]} on the "
             f"left; colour = name group of the best reverse hit; // = boundary gene clipped to {a.clip_kb:g} kb.",
             legend_extra=[("ortholog (reciprocal best hit)", True, False), ("paralog", False, False),
                           ("disabled copy", False, True)],
             empty_reason=empty)


if __name__ == "__main__":
    main()
