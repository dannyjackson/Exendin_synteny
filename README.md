# Exendin_synteny
Where do the exendin-family genes sit between **SLC4A8** and **AQP5**, across taxa, by
three independent lines of evidence? Each line is drawn as the same synteny figure, with the
Cactus species tree as backbone and one colour per gene **name**:

1. what each assembly's own **NCBI annotation** calls an exendin (step 2);
2. what **TOGA2** projects from the HelSus reference (step 3);
3. what a **reciprocal BLAST** search of every exendin against every genome finds, curated
   into high-confidence orthologs and paralogs (step 4).

Built on [Exendin_evolution](https://github.com/dannyjackson/Exendin_evolution): its drawing
engine, interval scanner, locus rebuilder and TOGA2 collectors are reused unchanged.

| Step | File | What it does |
|---|---|---|
| 0 | `0_Environment_management.md`, `environment.yml` | Paths and the `exendin_env` environment |
| 1 | `1_download_files.md` | TOGA2 (as before), species table, genomes |
| 2 | `2_synteny_ncbi.md` | Synteny from NCBI GFFs |
| 3 | `3_synteny_toga2.md` | Synteny from TOGA2 |
| 4 | `4_reciprocal_blast.md` | Reciprocal BLAST, curated orthologs/paralogs, synteny of that set |

## Where each file comes from
Files marked **copy** are byte-identical to the old repo and can be copied from there
(renamed where the name differs). Everything else is new.

| This repo | Origin | Notes |
|---|---|---|
| `lib_synteny_engine.py` | **copy** of `4d_plot_exendin_synteny.py` | drawing engine for every figure (renamed: it is now a library, imported by 2b, 3c, 4f) |
| `lib_scan_interval.py` | **copy** of `2f_scan_interval.py` | locus calling from tblastn hits (4c runs it with `--blast-out`) |
| `lib_rescue_locus.py` | **copy** of `2e_rescue_locus.py` | protein/CDS models of curated loci (4e) |
| `3a_list_reference_exendins.py` | **copy** of `4a_list_reference_exendins.py` | renamed only |
| `3b_collect_toga2_exendins.py` | **copy** of `4b_collect_toga2_exendins.py` | renamed only |
| `4_reference_peptides.fa` | **copy** of `2_reference_peptides.fa` | published precursors: queries and reverse competitor |
| `.gitattributes` | **copy** | |
| `0_Environment_management.md` | adapted from `0_Environment_management.md` | new paths; fixes `hal_toga_env` / `soure` |
| `1_download_files.md` | adapted from `1_download_files.md` | TOGA2 section to be copied verbatim; new genome section |
| `1_species_pinned.tsv` | new | accessions fixed by hand (TOGA2 assemblies, old survey genomes) |
| `1a_make_species_table.py` | new | tree tips + TOGA2 queries → `species.tsv` with NCBI accessions |
| `1b_download_genomes.sh` | new (download loop modelled on `3a_scan_species.sh`) | genome, GFF3, proteins; reuses old downloads |
| `lib_exendin_names.py` | new | name → colour group (`exendin-2-like` = `exendin-2`), fixed colours, species table I/O |
| `2a_ncbi_exendin_loci.py` | new | family and boundary genes from each NCBI GFF |
| `2b_plot_ncbi_synteny.py` | new (driver in the style of `2g` / `3b`) | step-2 figure |
| `3c_plot_toga2_synteny.py` | new (`4d`'s `main()` with name colours) | step-3 figure |
| `4a_build_query_set.py` | new | forward queries, boundary queries, reverse database |
| `4b_forward_tblastn.sh` | new | tblastn against one whole genome |
| `4c_call_loci.py` | new | loci (via `lib_scan_interval.py`) and boundary clusters |
| `4d_reciprocal_classify.py` | new | blastx back; ortholog / paralog / rejected; confidence; placement |
| `4e_build_models.py` | new | models of curated loci (via `lib_rescue_locus.py`) |
| `4f_plot_curated_synteny.py` | new | step-4 figure |

Not carried over: the tree, nomenclature and copy-number scripts of the old repo
(`2a–2d`, `2g`, `3a`, `3b`, `4c`, `Working_Scripts/`). They answer other questions. The
step-4 `models/curated_aa.fa` and `curated_cds.fa` can go straight into the old
`2d_assign_nomenclature.py` and the alignment/tree steps.

## Conventions
As in the old repo:
- Every markdown file starts from the set-up block in `0_Environment_management.md`.
- Every script is read-only on its inputs and writes only to `--outdir` (or named outputs).
- Every script logs `[info]`/`[warn]` lines to stderr; read the log before the figure.
- Python scripts need only the standard library, plus `blast`/`samtools` where they run them.
- Species are referred to by a 3+3 code (`HelSus`). These codes match TOGA2 query names and,
  via `Genus_species` → `GenSpe`, the Cactus tree tips.
