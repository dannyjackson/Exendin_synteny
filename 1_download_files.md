# Download files
Run the set-up block of `0_Environment_management.md` first.

## 1. TOGA2 run
**Already downloaded (2026-10-23).** Copy the "## TOGA2" section and its Readme unchanged
from `1_download_files.md` of
[Exendin_evolution](https://github.com/dannyjackson/Exendin_evolution): rclone of Elise
Lauterbur's run into `${TOGADIR}` (25 GB; UCSC browser and Nextflow files excluded). Steps 3
and 4 use `${TOGADIR}/target_HelSus`, and step 3 also uses the HelSus EviAnn annotation TOGA2
projected:
```
ls -d ${TOGADIR}/target_HelSus
ls -l ${HELSUS_EVIANN_GFF}
```

## 2. Species table
One row per taxon: every Cactus-tree tip, every TOGA2 query, and the outgroups given with
`--extra`. Chicken (`GalGal`) is the reference row of every figure and the root of the gene
trees. Accessions come from `1_species_pinned.tsv` where listed (chicken:
`GCF_016699485.2`, bGalGal1.mat.broiler.GRCg7b), otherwise from NCBI: RefSeq annotated >
annotated > reference > best contig N50.
```
cd ${WORK}
python3 ${SCRIPTDIR}/1a_make_species_table.py \
    --tree ${SPTREE} --toga2-run ${TOGADIR}/target_HelSus \
    --pinned ${SCRIPTDIR}/1_species_pinned.tsv --extra Gallus_gallus \
    --out species.tsv --candidates species_candidates.tsv 2> species.log
cat species.log
cut -f1 species.tsv | sort | uniq -d           # must print nothing: a duplicate row draws a species twice
```
- **TOGA2 queries that are not pinned** are warned about: NCBI's choice may not be the
  assembly TOGA2 ran on, and then contig names differ between steps 2 and 3. Pin the
  TOGA2 assembly in `1_species_pinned.tsv` and re-run.
- **Tree tips that are not `Genus_species`** are skipped unless pinned (the pinned row
  supplies the species name).
- **`NA` rows** have no NCBI assembly; they stay, so figures show them as empty rows.
- `species.tsv` is the one list every later step reads; edit it by hand if a choice is
  wrong. Every NCBI candidate is in `species_candidates.tsv`.

## 3. Genomes, annotations, proteins
```
bash ${SCRIPTDIR}/1b_download_genomes.sh species.tsv 2>&1 | tee download.log
```
Downloads `genome,gff3,protein,seq-report` per accession into
`${GENOMEDIR}/<code>_<accession>/` (the old `3a_scan_species.sh` layout, so its genomes are
reused and only a missing GFF or protein file is fetched), and indexes each FASTA. The
summary lists which species have an NCBI annotation; unannotated assemblies (e.g. HelCha)
are empty rows in step 2 and are still searched in step 4. Tens of GB: run in `tmux` or
`sinteractive`.
