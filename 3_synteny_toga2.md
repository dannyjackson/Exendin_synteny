# Synteny from TOGA2 annotations
The step-2 figure from TOGA2 instead of NCBI: same rows, layout and colours. TOGA2 projects
the HelSus EviAnn annotation onto each query, so it finds copies NCBI did not annotate, but
a model is only as good as the reference isoform behind it (HelCha exendin-3: TOGA2 used
isoform `mRNA-2`, whose first exon has no protein support, and produced a 198-aa model; see
step 4.1).

Steps 1–2 are the old repo's TOGA2 collection, unchanged; step 3 colours tips by name.

Run the set-up block of `0_Environment_management.md` first.
```
cd ${WORK}/synteny_toga2
```

## 1. Reference genes annotated as exendin (`3a` = old `4a`, unchanged)
```
python3 ${SCRIPTDIR}/3a_list_reference_exendins.py \
    --gff ${HELSUS_EVIANN_GFF} --run ${TOGADIR}/target_HelSus \
    --out reference_exendins.tsv 2> reference_exendins.log
cat reference_exendins.log
cut -f1,2,8,9,10 reference_exendins.tsv
```
Expect `LOC_00012324` (exendin-2) and `LOC_00012325` (exendin-4). `matched_text` names each
gene in the plot. GCGL (`LOC_00012323`) is included only if EviAnn calls it exendin; to add
it, re-run with `--pattern 'exendin|glucagon-like'` and check no real proglucagon got in.

## 2. One TOGA2 gene per exendin ortholog per query (`3b` = old `4b`, unchanged)
```
python3 ${SCRIPTDIR}/3b_collect_toga2_exendins.py \
    --run ${TOGADIR}/target_HelSus --gene-table reference_exendins.tsv \
    --status FI,I,PI --min-cds 150 --outdir toga2_exendins 2> toga2_exendins.log
cat toga2_exendins.log
cut -f1,3,4,7,9,10,11,16,18,19 toga2_exendins/tips.tsv | column -t

# optional: paralogous projections too (drawn dashed)
python3 ${SCRIPTDIR}/3b_collect_toga2_exendins.py \
    --run ${TOGADIR}/target_HelSus --gene-table reference_exendins.tsv \
    --include-paralogs --outdir toga2_exendins_paralogs 2> toga2_exendins_paralogs.log
```
`aa_len` far from 84–113 aa = a model problem (the old `4_toga2_exendin_tree.md` step 2);
many2many tips are copies TOGA2 could not split. `tips.tsv` feeds step 4's query set.

## 3. Plot
```
python3 ${SCRIPTDIR}/3c_plot_toga2_synteny.py --run ${TOGADIR}/target_HelSus \
    --tips toga2_exendins/tips.tsv --reference-table reference_exendins.tsv \
    --boundaries SLC4A8,AQP5=LOC_00012326 --gff ${HELSUS_EVIANN_GFF} \
    --species-table ../species.tsv --species-tree ${SPTREE} --top GalGal \
    --svg toga2_synteny.svg --png toga2_synteny.png --layout toga2_synteny.tsv 2> toga2_synteny.log
cat toga2_synteny.log
```
For the paralog version, use `--tips toga2_exendins_paralogs/tips.tsv` and other output names.
- **Log.** The name each reference gene got (`LOC_00012324 ... -> exendin-2`) and the
  boundary IDs; check the SLC4A8 ID. Override a name with `--name-map`
  (`gene_id<TAB>group`).
- **Drawing.** Colour = name of the reference exendin(s) a tip is orthologous to
  (many2many → combined name, e.g. `exendin-2/4`). Filled = passed 3b's filters, hollow =
  failed, dashed = paralogous projection. Labels are 3b tip names (`EX2`, `EX4b`, `EX2p`).
- **Rows.** `--species-table` gives every species a row in step-2 order ("not a query in
  this TOGA2 run" for chicken and any species TOGA2 did not cover). `//` breaks as in
  step 2.
- **Coordinates** are TOGA2's (HelSus `chr2`), not NCBI's (`NC_139459.1`): the figures line
  up by species and gene order, not by position.
