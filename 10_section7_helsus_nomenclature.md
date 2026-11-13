## 7. One locus, three nomenclatures: the HelSus interval
One figure of the HelSus SLC4A8–AQP5 interval, drawn three times on the same axis:
| Panel | Annotation | Names shown |
|---|---|---|
| **A** | the genome's own annotations: EviAnn (the reference TOGA2 projected) and NCBI RefSeq | each annotation's gene ID and description (`LOC_00012325 · Exendin-4`, `LOC146345366 · exendin-4-like`) |
| **B** | TOGA2 with HelSus as both reference and query (`target_HelSus/query_HelSus`, the self-projection) | TOGA2's query-gene ID and the reference transcript it was projected from (`reg_6 · LOC_00012325-mRNA-1`) |
| **C** | step 5 (`names/names.tsv`) | your names: `GCGL`, `anc`, `A`, `B` … `Z`, `Y` … |

A faint vertical band, in the panel-C colour, runs from each named locus up through B and
A, so the three names of one locus line up. Where a band crosses empty track in A, that
annotation has no gene at the locus.

```
cd ${WORK}/toga2_exendin_tree
HELSUS_NCBI_GFF=${GENOMEDIR}/HelSus_GCF_057233095.1/ncbi_dataset/data/GCF_057233095.1/genomic.gff
python3 ${SCRIPTDIR}/10i_plot_helsus_nomenclature.py \
    --gff EviAnn=${HELSUS_EVIANN_GFF} --gff "NCBI RefSeq=${HELSUS_NCBI_GFF}" \
    --toga2 ${TOGADIR}/target_HelSus/query_HelSus --names names/names.tsv \
    --boundaries SLC4A8,AQP5=LOC_00012326 \
    --svg helsus_nomenclature.svg --pdf helsus_nomenclature.pdf --png helsus_nomenclature.png \
    --crosswalk names/HelSus_nomenclature_crosswalk.tsv 2> helsus_nomenclature.log
cat helsus_nomenclature.log
```
Read the log before the figure:
- **`view chr2:… (N kb), SLC4A8 … AQP5 …`**: the window, taken from the first `--gff`. This
  must be EviAnn, because TOGA2 and `names.tsv` use its coordinates. If SLC4A8 lies to the
  right in the assembly, the axis is flipped so that SLC4A8 is always on the left, as in
  every other figure.
- **`A/NCBI RefSeq: N gene(s) on NC_139459.1:…`**: NCBI's sequence was found through
  `sequence_report.jsonl` next to its GFF (chromosome `2` = `chr2`). A warning `no genes
  on …` means it was not. Add `--chrom-map chr2=NC_139459.1`.
- **`B: N TOGA2 query gene(s) in the view, M exendin-family`** and, if present, a
  warning about reference genes with no name. That warning means `--ref-gff` (default: the
  EviAnn GFF) does not name them.
- **`C: N named HelSus loci in the view; M TOGA2 exendin-family gene(s) without a name`**:
  the unnamed ones are drawn hollow grey in C. They are TOGA2 loci left out of the tree
  (not intact, or curated out in step 3), as in the step-5 figure. `not drawn: …` lists
  named loci outside the window or on another contig.
- **One line per locus** at the end of the log: its name, position and the ID in every
  annotation (`none` = no gene there). The same is in `--crosswalk` as a table.

Options:
- `--toga2-label query_gene | projection | both` (default `both`) sets what B shows.
- `--label-all` labels the grey genes too.
- `--pad-kb 150` widens the window past the boundaries. In the step-2 note,
  HelSus has copies 118 and 379 kb beyond GCGL: check whether they fall in the window.
- `--no-bands` removes the bands.
- `--width-mm 88` gives a single-column figure, and `--font-pt` sets every label size.

Reading it:
- **Colours differ by panel, and the legend says which:** A and B use the project's
  name colours (exendin-2, exendin-4, GCGL, as in steps 2 and 3); C uses step 5's class
  colours (GCGL green, Exendin_anc blue, lettered copies orange). The two schemes
  share hues, so the same hue means a different thing in A/B than in C. A blue
  `exendin-4` arrow in A can sit under an orange `Z` in C. That mismatch is the point of
  the figure: what each annotation calls a locus against where it sits in the tree.
- **Hollow** = NCBI pseudogene (A) or TOGA2 status other than FI/I (B); in C, a TOGA2
  exendin-family locus without a name. **Dashed** = a TOGA2 paralogous projection only.
- **Labels** are spread apart and joined to their gene by a leader line. Hover text in
  the SVG gives each gene's ID, description and coordinates.

`names/HelSus_nomenclature_crosswalk.tsv` has one row per named (or unnamed TOGA2) locus:
your name, class, position, and the overlapping IDs and labels in EviAnn, NCBI RefSeq and
TOGA2. It complements `names/HelSus_annotations.tsv` (step 6), which lists transcripts
for expression data.

| File | Content |
|---|---|
| `helsus_nomenclature.{svg,pdf,png}` | the three-panel figure |
| `names/HelSus_nomenclature_crosswalk.tsv` | every locus: your name ↔ EviAnn ID ↔ NCBI ID ↔ TOGA2 query gene |
