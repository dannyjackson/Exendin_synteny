# Synteny from TOGA2 annotations
The step-2 figure from TOGA2 instead of NCBI: same rows, layout and colours. TOGA2 projects
a reference annotation onto each query (five runs here: HelSus, HemCap, AnoCar, CroPor and
their integration), so it finds copies NCBI did not annotate, but a model is only as good
as the reference isoform behind it (HelCha exendin-3: TOGA2 used
isoform `mRNA-2`, whose first exon has no protein support, and produced a 198-aa model; see
step 4.1).

Steps 1–2 are the old repo's TOGA2 collection, unchanged; their `tips.tsv` feeds step 4. Step 3 draws the figures: every TOGA2 gene between SLC4A8 and AQP5, for each of the five TOGA2 runs, in the publication layout of step 2.

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

## 3. Plot: every gene in the interval, one figure per TOGA2 run (`3c`, new)
The same layout as the step-2 figure (`lib_pub_synteny.py`), drawn from TOGA2 instead of
NCBI. This time every gene TOGA2 annotated between SLC4A8 and AQP5 is drawn, not only the
exendins. There is one figure per run:

| Run | Reference | Query directories |
|---|---|---|
| HelSus | *Heloderma suspectum* (EviAnn) | `target_HelSus/query_*` |
| HemCap | *Hemicordylus capensis* | `target_HemCap/query_*` |
| AnoCar | *Anolis carolinensis* | `query_*/target_AnoCar/output` |
| CroPor | *Crocodylus porosus* | `query_*/target_CroPor/output` |
| integrated | AnoCar + CroPor (check: see below) | `integrated_annotations/query_*/integrated` |

```
cd ${WORK}/synteny_toga2
mkdir -p interval && cd interval
```

### 3.1 Reference gene names
TOGA2 names each projection `transcript#gene#chain`. The `gene` field is whatever the
reference annotation called it, so first check whether each run already uses gene
symbols:
```
for R in target_AnoCar target_CroPor; do echo "== $R"; cut -f1 ${TOGADIR}/query_HelCha/${R}/output/orthology_classification.tsv | sed -n '2,6p'; done
cut -f1 ${TOGADIR}/target_HemCap/query_HelCha/orthology_classification.tsv | sed -n '2,6p'
cut -f2 ${TOGADIR}/integrated_annotations/query_HelCha/integrated/query_genes.tsv | sed -n '2,6p'
ls ${TOGADIR}/integrated_annotations/query_HelCha/integrated/
```
- **Symbols (`AQP5`, `SLC4A8`):** no GFF is needed for that run.
- **IDs (`LOC_…`, `XM_…`, `rna-…`):** pass the reference's GFF with `--ref-gff`. It must be
  the annotation of the **same assembly** the run used. Any name the script cannot resolve
  is counted per query in the log (`... have no name`). Add the right GFF, or a
  `--name-map` (`ref_gene<TAB>symbol`), until those warnings are gone.
- **The integrated `ls`** shows which files TOGA2 integration wrote. The script needs
  `query_genes.bed` + `query_genes.tsv`, or else `query_annotation.bed`. The prefixes of the
  gene names (`cut -f2` above) show which references were merged.

```
HEMCAP_GFF=$(ls ${GENOMEDIR}/HemCap_*/ncbi_dataset/data/*/genomic.gff | head -1)
ANOCAR_GFF=$(ls ${GENOMEDIR}/AnoCar_*/ncbi_dataset/data/*/genomic.gff | head -1)
CROPOR_GFF=$(ls ${GENOMEDIR}/CroPor_*/ncbi_dataset/data/*/genomic.gff 2>/dev/null | head -1)
echo "HemCap ${HEMCAP_GFF:-MISSING}"; echo "AnoCar ${ANOCAR_GFF:-MISSING}"; echo "CroPor ${CROPOR_GFF:-MISSING}"
```
If a GFF is missing and its run uses IDs, download that reference's annotation into
`${GENOMEDIR}` the same way as step 1.3 (`datasets download genome accession <ACC> --include gff3`),
using the accession TOGA2 used. Ask Elise which assemblies those were.

### 3.2 One figure per run
`--near-kb 0` (the default here): in a split genome, each boundary line extends as far as
the longest complete interval among the species. `--split-genes all` (also the default)
draws every gene on that stretch, not only exendins.
```
python3 ${SCRIPTDIR}/3c_plot_toga2_synteny.py --dirs "${TOGADIR}/target_HelSus/query_{q}" \
    --ref-gff ${HELSUS_EVIANN_GFF} --boundaries SLC4A8,AQP5=LOC_00012326 \
    --species-table ../../species_toga2.tsv --species-tree ${SPTREE} \
    --svg interval_HelSus.svg --pdf interval_HelSus.pdf --png interval_HelSus.png \
    --layout interval_HelSus_layout.tsv --genes-out interval_HelSus_genes.tsv 2> interval_HelSus.log

python3 ${SCRIPTDIR}/3c_plot_toga2_synteny.py --dirs "${TOGADIR}/target_HemCap/query_{q}" \
    --ref-gff ${HEMCAP_GFF} \
    --species-table ../../species_toga2.tsv --species-tree ${SPTREE} \
    --svg interval_HemCap.svg --pdf interval_HemCap.pdf --png interval_HemCap.png \
    --layout interval_HemCap_layout.tsv --genes-out interval_HemCap_genes.tsv 2> interval_HemCap.log

python3 ${SCRIPTDIR}/3c_plot_toga2_synteny.py --dirs "${TOGADIR}/query_{q}/target_AnoCar/output" \
    --ref-gff ${ANOCAR_GFF} \
    --species-table ../../species_toga2.tsv --species-tree ${SPTREE} \
    --svg interval_AnoCar.svg --pdf interval_AnoCar.pdf --png interval_AnoCar.png \
    --layout interval_AnoCar_layout.tsv --genes-out interval_AnoCar_genes.tsv 2> interval_AnoCar.log

python3 ${SCRIPTDIR}/3c_plot_toga2_synteny.py --dirs "${TOGADIR}/query_{q}/target_CroPor/output" \
    ${CROPOR_GFF:+--ref-gff ${CROPOR_GFF}} \
    --species-table ../../species_toga2.tsv --species-tree ${SPTREE} \
    --svg interval_CroPor.svg --pdf interval_CroPor.pdf --png interval_CroPor.png \
    --layout interval_CroPor_layout.tsv --genes-out interval_CroPor_genes.tsv 2> interval_CroPor.log

python3 ${SCRIPTDIR}/3c_plot_toga2_synteny.py --dirs "${TOGADIR}/integrated_annotations/query_{q}/integrated" \
    --ref-gff ${ANOCAR_GFF} ${CROPOR_GFF:+--ref-gff ${CROPOR_GFF}} --ref-gff ${HEMCAP_GFF} --ref-gff ${HELSUS_EVIANN_GFF} \
    --boundaries SLC4A8,AQP5=LOC_00012326 \
    --species-table ../../species_toga2.tsv --species-tree ${SPTREE} \
    --svg interval_integrated.svg --pdf interval_integrated.pdf --png interval_integrated.png \
    --layout interval_integrated_layout.tsv --genes-out interval_integrated_genes.tsv 2> interval_integrated.log

grep -h "warn" interval_*.log | sort | uniq -c | sort -rn | head -40
```

### 3.3 Check each log before reading its figure
- **`N queries: ...`** should list all 14 TOGA2 genomes. A genome with no files in a run
  is drawn as "no TOGA2 annotation found".
- **`... have no name`**: reference IDs that the script could not resolve. Those genes are drawn as
  "uncharacterised". Fix this with `--ref-gff` / `--name-map` (3.1).
- **Boundaries:** each genome shows either a pair (`SLC4A8 … and AQP5 … on <contig>, N kb
  apart`) or a single gene (`(no pair on one contig)`). In the HelSus run, AQP5 is only
  found through `AQP5=LOC_00012326`. In the other runs, check that every genome found both.
  `no SLC4A8 gene` in most genomes means the symbol is spelled differently in that
  reference. Pass its ID as `SLC4A8=<ref gene ID>`.
- **`not drawn`** lists family genes outside the interval (as in step 2).
- **Boundary names:** the HelSus EviAnn annotation names genes `Similar to SLC4A8: ...`. The
  script reads `SLC4A8` from that, so every genome should now report an SLC4A8. If most
  genomes still say `no SLC4A8 gene`, pass the reference ID: `SLC4A8=<ref gene ID>`.
- `interval_<run>_genes.tsv` lists every gene in or near the interval, with its query gene
  ID, symbol, TOGA2 status and reference gene. Use it to check any gene on the figure.

### 3.4 Reading the figures
- **Colour:** exendins and GCGL keep the fixed colours used in every figure. **Every other
  gene is grey.** `--colour exendin` greys GCGL too. `--colour name` gives each gene symbol
  its own pale colour, for checking which genes are which. `interval_<run>_genes.tsv` names
  every gene, grey or not. SLC4A8 and AQP5 are both dark grey; the headers say which end is which.
- **Filled** = TOGA2 status FI/I. **Hollow** = PI, UL, M or L (not intact). **Dashed outline** =
  only paralogous or retrogene projections.
- **Reference genome rows:** a reference is only drawn if it is also a query. HelSus is (its
  self-projection); HemCap, AnoCar and CroPor are not.
- **Comparing runs** is the point of having five figures. A gene present in one run and
  absent in another is a reference-annotation difference, not biology. Only genes the
  reference itself annotates can be projected: the venom exendins can only come from the
  HelSus run.
