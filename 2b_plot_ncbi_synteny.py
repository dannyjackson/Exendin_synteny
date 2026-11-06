#!/usr/bin/env python3
"""
Step-2 synteny plot: every exendin-family gene in each species' NCBI annotation, between
the two boundary genes (SLC4A8 ... AQP5), species in the order of the Cactus tree, which
is drawn alongside.

Input: ncbi_loci.tsv from 2a_ncbi_exendin_loci.py. Drawing, orientation, clipping of long
boundary genes and the tree come from lib_synteny_engine.py (the old repo's
4d_plot_exendin_synteny.py, unchanged), so this figure reads like the old 4d / 2g / 3b
figures.

  colour   the gene's NAME group (lib_exendin_names.name_group): exendin-2 and
           exendin-2-like share a colour; the same colours are used in steps 3 and 4
  label    short form of the NCBI description (Ex2-like, Ex4, GCGL ...)
  filled   protein-coding gene
  hollow   pseudogene, or any other non-coding biotype
  tracks   one per contig; 'unplaced' genes (contig with no boundary gene) get their
           own track below the species' main track

Species with no NCBI annotation are drawn as an empty row saying so.

Usage
  python3 2b_plot_ncbi_synteny.py --loci ncbi_loci/ncbi_loci.tsv \
      --species-tree ${CACTUSDIR}/Anguimorpha_species_tree_iteration-01_5-23-26.nwk \
      --svg ncbi_synteny.svg --png ncbi_synteny.png --layout ncbi_synteny.tsv
"""
import argparse, csv, os, sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import lib_exendin_names as N  # noqa: E402

S = N.load_engine()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--loci", required=True, help="ncbi_loci.tsv from 2a")
    ap.add_argument("--placements", default="interval,near_boundary,unplaced",
                    help="which 2a placements to draw (add 'elsewhere' to see everything)")
    ap.add_argument("--boundary-names", default="SLC4A8,AQP5", help="as in 2a --boundaries, left one first")
    ap.add_argument("--title", default="Exendin-family genes between SLC4A8 and AQP5 (NCBI annotations)")
    S.add_common_args(ap)
    a = ap.parse_args()

    keep = set(a.placements.split(","))
    bnames = a.boundary_names.split(",")
    with open(a.loci) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    queries, empty = [], {}
    contigs = defaultdict(lambda: defaultdict(list))
    for r in rows:
        q = r["species"]
        if q not in queries:
            queries.append(q)
        if r["kind"] == "none":
            empty[q] = {"no_annotation": "no NCBI annotation", "no_genome": "no genome"}.get(
                r["placement"], r["placement"])
            continue
        if r["kind"] == "boundary":
            contigs[q][r["contig"]].append(dict(
                kind="boundary", name=r["group"], label=r["group"], group=r["group"], chrom=r["contig"],
                start=int(r["start"]), end=int(r["end"]), strand=r["strand"], status="", filled=True,
                dashed=False, hover=f'{r["symbol"]} {r["gene_id"]} {r["description"]} ({r["note"]})'))
            continue
        if r["placement"] not in keep:
            continue
        contigs[q][r["contig"]].append(dict(
            kind="exendin", name=f'{q}_{r["symbol"]}', label=r["label"], group=r["group"], chrom=r["contig"],
            start=int(r["start"]), end=int(r["end"]), strand=r["strand"], status=r["biotype"],
            filled=r["biotype"] == "protein_coding", dashed=False,
            hover=f'{r["symbol"]} | {r["description"]} | {r["biotype"]} | {r["placement"]} | '
                  f'{r["protein_id"]} {r["aa_len"]} aa'))
    for q in queries:
        if q not in empty and not any(f["kind"] == "exendin" for fs in contigs[q].values() for f in fs):
            N.info(f"{q}: no family gene drawn")
            if not contigs[q]:
                empty[q] = "no family gene or boundary gene annotated"
    N.finish(S, a, contigs, queries, bnames, a.title,
             f"NCBI RefSeq/GenBank annotations. Tracks oriented with {bnames[0]} on the left; arrows give strand "
             f"after orientation; colour = gene name (-like merged); // = boundary gene clipped to {a.clip_kb:g} kb.",
             legend_extra=[("protein-coding", True, False), ("pseudogene / non-coding", False, False)],
             empty_reason=empty)


if __name__ == "__main__":
    main()
