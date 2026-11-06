# Download files
Run the set-up block of `0_Environment_management.md` first.

## TOGA2 annotations
**Already downloaded on 2026-10-23.** Copy this section unchanged from
`1_download_files.md` ("## TOGA2" and its Readme) of
[Exendin_evolution](https://github.com/dannyjackson/Exendin_evolution). It fetches
Elise Lauterbur's TOGA2 run with rclone into `${TOGADIR}` (25 GB, UCSC browser and Nextflow
files excluded). Steps 3 and 4 use only `${TOGADIR}/target_HelSus`.

The HelSus EviAnn annotation that TOGA2 projected is also needed (step 3); it is the file
the old repo used:
```
HELSUS_EVIANN_GFF=${PROJDIR}/datasets/genomes/Heloderma_suspectum/Heloderma_suspectum.final.reference.gff3
ls -l ${HELSUS_EVIANN_GFF}
```

## Species table
One row per taxon: every tip of the Cactus tree, every TOGA2 query, and any outgroup you add
with `--extra`. Accessions come from `1_species_pinned.tsv` where listed, otherwise from NCBI
(RefSeq annotated > annotated > reference > best contig N50). Every NCBI candidate is
listed in `species_candidates.tsv`.
```
cd ${WORK}
python3 ${SCRIPTDIR}/1a_make_species_table.py \
    --tree ${SPTREE} \
    --toga2-run ${TOGADIR}/target_HelSus \
    --pinned ${SCRIPTDIR}/1_species_pinned.tsv \
    --out species.tsv --candidates species_candidates.tsv 2> species.log
cat species.log
```
Check:
- **Warnings about TOGA2 queries that are not pinned.** For those species the NCBI choice
  may not be the assembly TOGA2 ran on, so contig names will differ between the step-2 and
  step-3 figures. Look up the assembly each TOGA2 query used, and add a row for it to
  `1_species_pinned.tsv` (and commit it), then re-run.
- **Tree tips that are not `Genus_species`** are skipped unless pinned. A pinned row also
  supplies the species name.
- **Rows with `NA`** have no NCBI assembly. They stay in the table, so the plots show them
  as empty rows instead of leaving them out without saying so.
- To add outgroups such as zebra finch (whose GCGL gene sits between AQP5 and SLC4A8), use
  `--extra Taeniopygia_guttata,Gallus_gallus`. They are drawn below the tree.

`species.tsv` is the single list every later step reads. Edit it by hand if a choice is
wrong.

## Genomes, annotations and proteins
```
bash ${SCRIPTDIR}/1b_download_genomes.sh species.tsv 2>&1 | tee download.log
```
This downloads `genome,gff3,protein,seq-report` per accession into
`${GENOMEDIR}/<code>_<accession>/`, the layout the old `3a_scan_species.sh` used. So
genomes the old repo already fetched (`HelSus_GCF_057233095.1`, `HelCha_GCA_026122225.1`,
`AnoCar_...` etc.) are reused, and only a missing GFF or protein file is fetched. The
summary at the end lists which species have an NCBI annotation. Unannotated assemblies
(e.g. HelCha) are drawn as "no NCBI annotation" in step 2, and still searched in step 4.

On Biowulf, run it inside `tmux` or an `sinteractive` session: a dozen squamate genomes is
tens of GB.
