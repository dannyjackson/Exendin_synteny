# Exendin_synteny
Which exendin-family genes sit between **SLC4A8** and **AQP5** in each genome, and what
should they be called? Three lines of evidence are each drawn as the same synteny figure,
with the Cactus species tree as backbone, chicken as a reference row, and one colour per
gene name. A fourth step turns the sequence evidence into standard names and a named gene
tree.

| Step | File | What it does |
|---|---|---|
| 0 | `0_Environment_management.md` | paths, `exendin_env`, set-up block used by every step |
| 1 | `1_download_files.md` | TOGA2 run, species table, genomes / GFFs / proteins |
| 2 | `2_synteny_ncbi.md` | synteny from each assembly's NCBI annotation |
| 3 | `3_synteny_toga2.md` | synteny from TOGA2 projections of the HelSus reference |
| 4 | `4_reciprocal_blast.md` | reciprocal tblastn/blastx of every exendin vs every genome → curated orthologs/paralogs, models, synteny, gene tree |
| 5 | `5_nomenclature.md` | standard names (GCGL / positional letters / *Heloderma* published names), named tree, named synteny |
| 6 | `6_published_map.md` | The published exendins (helospectin I/II, helodermin, exendin-4) mapped codon by codon onto the HelSus genome; context figure with per-exendin pop-outs |
| 7 | `7_pseudogene_check.md` | HelSus_B / HelCha_C (no published match) against the helodermin precursor: exons, splice sites, frameshifts, stops, ORF, identity by region, dN/dS, TOGA2 status; intact EX2/EX1 as controls |

Built on [Exendin_evolution](https://github.com/dannyjackson/Exendin_evolution): its
drawing engine, interval scanner, tree drawer, rooting script, TOGA2 collectors and
nomenclature script are reused.

## Where each file comes from
**copy** = byte-identical to the old repo (renamed where the name differs).

| This repo | Origin | Used in | Role |
|---|---|---|---|
| `lib_synteny_engine.py` | **copy** of `4d_plot_exendin_synteny.py` | 2b 3c 4f 5d 6b 7b | drawing engine of every synteny figure; PNG writing for 6b, 7b |
| `lib_scan_interval.py` | **copy** of `2f_scan_interval.py` | 4c | locus calling from tblastn hits |
| `lib_rescue_locus.py` | **copy** of `2e_rescue_locus.py` | 4e `--method rescue` only | old cDNA-style rebuild, kept for comparison |
| `lib_draw_tree.py` | **copy** of `2c_draw_tree.py` | 4h | tree drawing |
| `lib_midpoint_root.py` | **copy** of `4c_midpoint_root.py` | 4, 5 | midpoint rooting |
| `3a_list_reference_exendins.py` | **copy** of `4a_list_reference_exendins.py` | 3 | |
| `3b_collect_toga2_exendins.py` | **copy** of `4b_collect_toga2_exendins.py` | 3 | |
| `4_reference_peptides.fa` | **copy** of `2_reference_peptides.fa` | 4a 4g 5a 7a | published peptides and precursors (7a: the helodermin precursor `EXENDIN2_precursor`) |
| `5a_assign_nomenclature.py` | copy of `2d_assign_nomenclature.py` **+ one fix** | 5 | the old file wrote `...` in place of the `differences` column |
| `lib_exendin_names.py` | new | all steps | name → colour group, fixed colours, `//` breaks, `--top` / `--inset`, species-table I/O |
| `lib_gff.py` | new | 6b 7a | GFF3 reader for NCBI and EviAnn GFFs (genes, transcripts, exons, CDS; genes by name or ID); TOGA2 query projections (BED12 + `loss_summary.tsv`) |
| `1_species_pinned.tsv`, `1a_make_species_table.py`, `1b_download_genomes.sh` | new | 1 | species table and downloads |
| `2a_ncbi_exendin_loci.py`, `2b_plot_ncbi_synteny.py` | new | 2 | NCBI loci and figure |
| `3c_plot_toga2_synteny.py` | new (old 4d's `main()` with name colours) | 3 | TOGA2 figure |
| `4a`–`4h` | new | 4 | query set, forward search, loci (with split-gene joining), reverse classification (with chromosome-aware placement, `--ignore-uncharacterized`), models (with missing-exon search), figure, tree inputs, tree drawing |
| `5b_name_genes.py`, `5c_rename_tips.py`, `5d_plot_named_synteny.py` | new | 5 | naming rules, tree relabelling, named synteny figure |
| `6_published_peptides.fa` | new | 6 (HelSus) | mature helospectin I, II, helodermin (Glu8,9), exendin-4 |
| `6_published_peptides_helcha.fa` | new | 6 (HelCha) | mature helospectin I, II, helodermin (Glu8,9), exendin-3 |
| `6a_map_published.py` | new | 6; functions reused by 7a | genome-wide tblastn, then a sensitive search in a window around the best loci; one codon per residue, ends extended over mismatches |
| `6b_plot_published_map.py` | new | 6; drawing helpers reused by 7b | panel A (SLC4A8–AQP5 from TOGA2, TOGA2 genes, curated reciprocal-BLAST genes), pop-outs B– per exendin |
| `7a_pseudogene_check.py` | new | 7 | the helodermin precursor mapped exon by exon onto chosen loci: splice sites, frameshifts, stops, ORF, identity by region, KR site, dN/dS, TOGA2 status; verdict intact / disrupted / partial |
| `7b_plot_pseudogene_check.py` | new | 7 | per-locus rows in precursor coordinates, summary table, dN/dS table |

## Conventions
- Every markdown file starts from the set-up block of `0_Environment_management.md`.
- Scripts are read-only on their inputs and write only to `--outdir` or named outputs.
- Scripts log `[info]` / `[warn]` to stderr: read the log before the figure.
- Species are 3+3 codes (`HelSus`); they match TOGA2 query names and, via
  `Genus_species` → `GenSpe`, the Cactus tree tips. `species.tsv` (step 1) is the one list
  every step reads.
- Python scripts use only the standard library; `blast` and `samtools` where they call them.
