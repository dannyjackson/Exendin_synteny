#!/usr/bin/env python3
"""
Step-3 synteny plot: the step-2 figure drawn from TOGA2 instead of NCBI. Every TOGA2
query gene orthologous to a reference exendin (tips.tsv from 3b_collect_toga2_exendins.py),
between TOGA2's ortholog projections of the two boundary genes, in Cactus-tree order.

This is the old repo's 4d_plot_exendin_synteny.py main() with one change: tips are
coloured by the NAME of the reference gene(s) they are orthologous to, using the same
name groups and fixed colours as step 2 (lib_exendin_names), instead of by TOGA2 label.
So an exendin-2 ortholog is the same colour in the NCBI and the TOGA2 figure.

Name of a reference gene: name_group() of its annotation text (matched_text column of
reference_exendins.tsv from 3a). A tip orthologous to several reference genes
(many2many) gets the combined name, e.g. 'exendin-2/4'. --name-map (TSV
gene_id<TAB>group) overrides the name of any reference gene.

  filled   tip that passed 3b's status / length filters (included = yes)
  hollow   tip that did not
  dashed   paralogous projection (only with 3b --include-paralogs)

--species-table adds an empty row for every species in species.tsv that is not a query
of this TOGA2 run, so the figure lines up with step 2.

Usage
  python3 3c_plot_toga2_synteny.py --run ${TOGADIR}/target_HelSus \
      --tips toga2_exendins/tips.tsv --reference-table reference_exendins.tsv \
      --boundaries SLC4A8,AQP5=LOC_00012326 --gff ${HELSUS_EVIANN_GFF} \
      --species-table ../species.tsv \
      --species-tree ${CACTUSDIR}/Anguimorpha_species_tree_iteration-01_5-23-26.nwk \
      --svg toga2_synteny.svg --png toga2_synteny.png --layout toga2_synteny.tsv
"""
import argparse, csv, os, sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import lib_exendin_names as N  # noqa: E402

S = N.load_engine()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", required=True, help="${TOGADIR}/target_HelSus")
    ap.add_argument("--tips", required=True, help="tips.tsv from 3b")
    ap.add_argument("--reference-table", required=True, help="reference_exendins.tsv from 3a")
    ap.add_argument("--name-map", help="TSV gene_id<TAB>group, overrides the name of a reference gene")
    ap.add_argument("--boundaries", default="SLC4A8,AQP5=LOC_00012326",
                    help="two NAME[=REF_GENE_ID], left one first")
    ap.add_argument("--gff", help="reference GFF, to look up boundary names given without an ID")
    ap.add_argument("--species-table", help="species.tsv; species not in the TOGA2 run get an empty row")
    ap.add_argument("--title", default="Exendin-family genes between SLC4A8 and AQP5 (TOGA2, HelSus reference)")
    S.add_common_args(ap)
    a = ap.parse_args()

    ref_group = N.toga2_reference_groups(a.reference_table, a.name_map)

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
            gid = S.gff_lookup(a.gff, name, near)
        bnames.append(name); bids[gid] = name
        N.info(f"boundary {name} = {gid}")
    if len(bnames) != 2:
        sys.exit("[error] --boundaries needs exactly two genes")

    contigs = defaultdict(lambda: defaultdict(list))
    for t in tips:
        g = N.toga2_tip_group(t["ref_genes"], ref_group)
        contigs[t["query"]][t["q_chrom"]].append(dict(
            kind="exendin", name=t["tip"], label=t["tip"].split("_", 1)[1] if "_" in t["tip"] else t["tip"],
            group=g, chrom=t["q_chrom"], start=int(t["q_start"]), end=int(t["q_end"]),
            strand=t["q_strand"], status=t["status"], filled=t["included"] == "yes",
            dashed=t.get("kind") == "paralog",
            hover=f'{g} | status {t["status"]} {t.get("orthology_class", "")} {t.get("q_gene", "")} '
                  f'{t.get("aa_len", "")} aa {t.get("reason", "")}'))
    queries = sorted(d[6:] for d in os.listdir(a.run) if d.startswith("query_"))
    for q in queries:
        bnd = S.boundary_projections(os.path.join(a.run, "query_" + q), set(bids))
        ex_mid = {c: sum((f["start"] + f["end"]) / 2 for f in fs) / len(fs) for c, fs in contigs[q].items()}
        for gid, name in bids.items():
            ps = bnd.get(gid, [])
            on = defaultdict(list)
            for p in ps:
                if p[0] in ex_mid:
                    on[p[0]].append(p)
            picks = ([S.choose_boundary(v, ex_mid[c]) for c, v in on.items()] if on
                     else ([S.choose_boundary(ps, None)] if ps else []))
            for p in picks:
                contigs[q][p[0]].append(dict(kind="boundary", name=name, label=name, group=name, chrom=p[0],
                                             start=p[1], end=p[2], strand=p[3], status=p[4], filled=True,
                                             dashed=False, hover=f"status {p[4]} {p[5]}"))
    empty = {}
    if a.species_table:
        for sp in N.read_species_table(a.species_table):
            if sp["code"] not in queries:
                queries.append(sp["code"]); empty[sp["code"]] = "not a query in this TOGA2 run"
    for q in queries:
        if q not in empty and not contigs.get(q):
            empty[q] = "no TOGA2 ortholog of a reference exendin"

    N.finish(S, a, contigs, queries, bnames, a.title,
             f"TOGA2 projections (HelSus reference). Tracks oriented with {bnames[0]} on the left; colour = name "
             f"of the reference exendin(s) the gene is orthologous to (-like merged); "
             f"// = boundary gene clipped to {a.clip_kb:g} kb.",
             legend_extra=[("passed 3b filters", True, False), ("failed 3b filters", False, False),
                           ("paralogous projection", False, True)],
             empty_reason=empty)


if __name__ == "__main__":
    main()
