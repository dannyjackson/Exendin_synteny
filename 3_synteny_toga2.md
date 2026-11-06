# Synteny of exendin-family genes from TOGA2 annotations
The step-2 figure, drawn from TOGA2 instead of NCBI. The rows, layout and colours are the
same, so the two figures can be read side by side. TOGA2 projects the HelSus reference
annotation (EviAnn) onto each query, so it finds copies NCBI did not annotate. A TOGA2 model
is only as good as the reference isoform behind it, though. The old repo's HelCha exendin-3
case is the warning example: TOGA2 placed it 5.9 kb from the real first exon.

Steps 1–2 are the old repo's TOGA2 collection, unchanged. Step 3 is new, so that the tips
are coloured by name like step 2 rather than by TOGA2 label.

Run the set-up block of `0_Environment_management.md` first.
```
HELSUS_EVIANN_GFF=${PROJDIR}/datasets/genomes/Heloderma_suspectum/Heloderma_suspectum.final.reference.gff3
cd ${WORK}/synteny_toga2
```

## 1. Reference genes annotated as exendin
`3a_list_reference_exendins.py` is `4a_list_reference_exendins.py` from
[Exendin_evolution](https://github.com/dannyjackson/Exendin_evolution), copied unchanged
(renamed only).
```
python3 ${SCRIPTDIR}/3a_list_reference_exendins.py \
    --gff ${HELSUS_EVIANN_GFF} --run ${TOGADIR}/target_HelSus \
    --out reference_exendins.tsv 2> reference_exendins.log
cat reference_exendins.log
cut -f1,2,8,9,10 reference_exendins.tsv
```
- **Expected genes.** Expect `LOC_00012324` (exendin-2) and `LOC_00012325` (exendin-4), as in
  the old step 4. `matched_text` is what names each gene in the plot.
- **Whether ANC is in.** Step 2 also draws genes NCBI calls glucagon-like (GCGL/ANC). Whether
  `LOC_00012323` (ANC) is in this set depends only on how EviAnn annotated it. If you want it
  in both figures and it is missing, add `--pattern 'exendin|glucagon-like'` and re-run 3a.
  Check that no real proglucagon got in.

## 2. One TOGA2 gene per exendin ortholog per query
`3b_collect_toga2_exendins.py` is `4b_collect_toga2_exendins.py` from the old repo, copied
unchanged (renamed only).
```
python3 ${SCRIPTDIR}/3b_collect_toga2_exendins.py \
    --run ${TOGADIR}/target_HelSus --gene-table reference_exendins.tsv \
    --status FI,I,PI --min-cds 150 --outdir toga2_exendins 2> toga2_exendins.log
cat toga2_exendins.log
cut -f1,3,4,7,9,10,11,16,18,19 toga2_exendins/tips.tsv | column -t

# optional: TOGA2's paralogous projections too (drawn dashed)
python3 ${SCRIPTDIR}/3b_collect_toga2_exendins.py \
    --run ${TOGADIR}/target_HelSus --gene-table reference_exendins.tsv \
    --include-paralogs --outdir toga2_exendins_paralogs 2> toga2_exendins_paralogs.log
```
See the old `4_toga2_exendin_tree.md` step 2 for how to read `tips.tsv`: `aa_len` far from
84–113 aa means a model problem, and many2many tips are copies TOGA2 could not split.

## 3. Plot
```
python3 ${SCRIPTDIR}/3c_plot_toga2_synteny.py \
    --run ${TOGADIR}/target_HelSus \
    --tips toga2_exendins/tips.tsv --reference-table reference_exendins.tsv \
    --boundaries SLC4A8,AQP5=LOC_00012326 --gff ${HELSUS_EVIANN_GFF} \
    --species-table ../species.tsv --species-tree ${SPTREE} \
    --svg toga2_synteny.svg --png toga2_synteny.png --layout toga2_synteny.tsv 2> toga2_synteny.log
cat toga2_synteny.log

python3 ${SCRIPTDIR}/3c_plot_toga2_synteny.py \
    --run ${TOGADIR}/target_HelSus \
    --tips toga2_exendins_paralogs/tips.tsv --reference-table reference_exendins.tsv \
    --boundaries SLC4A8,AQP5=LOC_00012326 --gff ${HELSUS_EVIANN_GFF} \
    --species-table ../species.tsv --species-tree ${SPTREE} \
    --svg toga2_synteny.paralogs.svg --png toga2_synteny.paralogs.png
```
- **The log.** It prints the name each reference gene got (`LOC_00012324: 'gene:Note=exendin-2'
  -> exendin-2`) and the boundary IDs. Check the SLC4A8 ID, as in the old 4d. If a name is
  wrong, override it with `--name-map` (TSV `gene_id<TAB>group`).
- **Drawing.** Colour = name of the reference exendin(s) a tip is orthologous to. A many2many
  tip gets the combined name (`exendin-2/4`), with its own colour. Filled = passed 3b's
  filters, hollow = failed them, dashed = paralogous projection. Labels are 3b tip names
  (`EX2`, `EX4b`, `EX2p`), so they read against `tips.tsv`.
- **Empty rows.** `--species-table` gives every species in `species.tsv` a row, in the same
  order as step 2. Species TOGA2 did not run on say "not a query in this TOGA2 run".
- **Contig names** are TOGA2's (e.g. HelSus `chr2`), not NCBI's (`NC_139459.1`). The two
  figures line up by species and gene order, not by coordinates.
