# Download files
Run the set-up block of `0_Environment_management.md` first.

## 1. TOGA2 run
***downloaded on October 23rd, 2026***
```
ssh -L 53682:localhost:53682 helix.nih.gov

module load rclone

rclone config


rclone lsd gdrive: --drive-shared-with-me 

# Top-level contents (should show integrated_annotations, query_*, target_*)
rclone lsd gdrive: --drive-root-folder-id 1PqZeWAav5iu-Baj7kkXPoPnYT30n29ug

# Full recursive listing with sizes
rclone lsf -R --format "ps" --drive-root-folder-id 1PqZeWAav5iu-Baj7kkXPoPnYT30n29ug gdrive: > toga2_listing.txt
wc -l toga2_listing.txt
head toga2_listing.txt

# pull all files, excluding the genome browser files and the Nextflow working files
module load rclone

cat > ~/toga2_filters.txt << 'EOF'
- **/ucsc_browser_files/**
- **/nextflow/**
- .command.*
- .exitcode
+ **
EOF

# Check the size first; it should be about 26-27 GB
rclone size gdrive: --drive-root-folder-id 1PqZeWAav5iu-Baj7kkXPoPnYT30n29ug \
  --filter-from ~/toga2_filters.txt

# Then copy, inside tmux so it survives a disconnect
tmux new -s toga2
rclone copy gdrive: /data/Wilson_Lab/projects/squamates/Arizona_squamates/misc/Toxicofera_TOGA2_annotations/ \
  --drive-root-folder-id 1PqZeWAav5iu-Baj7kkXPoPnYT30n29ug \
  --filter-from ~/toga2_filters.txt \
  --progress --transfers 8 --checkers 16
```
Readme:
```
# TOGA2 annotations (Toxicofera / Anguimorpha)

- **Downloaded:** 2026-10-23, on helix.nih.gov
- **Source:** Google Drive folder shared by Elise Lauterbur, who ran the analysis
  https://drive.google.com/drive/folders/1PqZeWAav5iu-Baj7kkXPoPnYT30n29ug
- **Destination:** `/data/Wilson_Lab/projects/squamates/Arizona_squamates/misc/Toxicofera_TOGA2_annotations/`
- **Size:** 25G
- **Excluded:** UCSC browser files (`ucsc_browser_files/`, ~11.5 GB) and Nextflow
  working files (`nextflow/`, `.command.*`, `.exitcode`): not needed for analysis

Script for downloading:
https://github.com/dannyjackson/Exendin_evolution/blob/main/1_download_files.md
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

## 4. TOGA2 reference annotations
The non-HelSus TOGA2 runs (`target_HemCap`, `query_*/target_CroPor`) name genes by their
reference's annotation. `3c_plot_toga2_synteny.py` needs these GFFs to name genes. They must
be the same assemblies TOGA2 used (accessions from Elise). Genome FASTAs are not needed;
add `genome` to `--include` if you want them.

| code | species | accession |
|---|---|---|
| CroPor | *Crocodylus porosus* | GCF_001723895.1 |
| HemCap | *Hemicordylus capensis* | GCF_027244095.1 |

```
cd /data/Wilson_Lab/data/Exendin_Evolution/genomes

datasets download genome accession GCF_001723895.1 --include gff3,protein,seq-report --filename CroPor.zip
unzip -o -q CroPor.zip -d CroPor_GCF_001723895.1 && rm CroPor.zip

datasets download genome accession GCF_027244095.1 --include gff3,protein,seq-report --filename HemCap.zip
unzip -o -q HemCap.zip -d HemCap_GCF_027244095.1 && rm HemCap.zip

ls -lh CroPor_GCF_001723895.1/ncbi_dataset/data/GCF_001723895.1/genomic.gff HemCap_GCF_027244095.1/ncbi_dataset/data/GCF_027244095.1/genomic.gff
```