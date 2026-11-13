# WormSynteny-style riparian plot of the interval
**Question:** the step-9 question again (what each genome has between SLC4A8 and AQP5),
drawn with the method and look of **WormSynteny** (Zheng lab; https://wormsynteny.org,
https://github.com/Zhenglabhku/WormSynteny).

WormSynteny draws one region of a Cactus alignment across species in four stages:
1. **First wave:** `halLiftover` lifts the query region into every species.
2. **Stitch and filter:** fragments closer than a *gap value* are joined, and joined
   fragments shorter than a *filtering percentage* of the query are dropped. The
   **longest** fragment is the species' row; the others are kept as "(filtered)" extras.
3. **Second wave:** `halLiftover` lifts each species' fragments into the **next** species of
   the row order. These alignments are the **links**.
4. **Plot:** `gggenomes(seqs, genes, links)`, with genes filled by orthogroup.

**Why not the WormSynteny package itself.** The app is hard-wired to its bundled
*Caenorhabditis* HAL (`inst/ubuntu/evolver_Caenorhabditis.hal`). Its eleven species,
`caenorhabditis_` genome names, *C. elegans*-only queries and annotation tables are
fixed throughout the code. Pointed at this HAL, it does not run. Step 14 therefore runs the
same four stages on this project's HAL, with the same tools: `halLiftover` (Cactus module),
`bedtools merge` for stitching, and `gggenomes` for the figure. The only new code joins
them: `14a` builds gggenomes' three tables and `14b` draws them.

| WormSynteny | Here |
|---|---|
| query: a *C. elegans* gene or region | the HelSus SLC4A8–AQP5 interval (step 9, `region.tsv`) |
| species order: fixed phylogeny | the species-tree order of steps 9 and 13 (`genomes.tsv`) |
| genes: WormBase annotation of each species | each genome's own TOGA2 annotation (step 9, `target_genes.tsv`) |
| gene colour: orthogroup (Ma, Lau & Zheng 2024) | the HelSus gene TOGA2 projected from: the step-5 name for exendins (`A`, `B`, `anc` …, project colours), else its symbol (WormSynteny's palette) |
| gap value / filtering % (e.g. 200 bp / 5 % for one gene) | `GAP` / `FILTER` below, scaled up for a ~300 kb interval |

Run the set-up block of `0_Environment_management.md` first. Step 14 needs step 9's
`hal_lift` outputs and the patched 9a, which writes `ref_genes`. If step 13 has not been
run, only its `genomes` command (section 1) is needed, not halSynteny.

## 0. Software (once)
`halLiftover` and `bedtools` are already in use (steps 9 and 1). For gggenomes, add it to
the R environment of step 13:
```
source myconda
mamba activate genespace_env            # step 13, section 0
Rscript -e 'install.packages(c("gggenomes", "dplyr", "readr"), repos = "https://cloud.r-project.org")'
Rscript -e 'library(gggenomes); packageVersion("gggenomes")'     # >= 1.0
```

## 1. Region and row order
```
mkdir -p ${WORK}/wormsynteny && cd ${WORK}/wormsynteny
mamba activate exendin_env
module load cactus/3.3.0
HAL=$(ls ${CACTUSDIR}/*.hal | head -1); echo ${HAL}
LIFT=../synteny_cactus/hal_lift

# row order: reuse step 13's, or make it (no halSynteny needed)
GENOMES=../genespace_cactus/genomes.tsv
[ -s ${GENOMES} ] || python3 ${SCRIPTDIR}/13a_hal_to_blk.py genomes --indir ${LIFT} \
    --species-tree ${SPTREE} --species-table ../species_toga2.tsv --only-listed \
    --drop-species AnnSte,BarImb,ElgCoe --out genomes.tsv && GENOMES=../genespace_cactus/genomes.tsv
cat ${GENOMES}

# the query: the HelSus interval as BED4, as WormSynteny writes its query
REF=$(awk -F'\t' 'NR==2{print $1}' ${LIFT}/region.tsv)
awk -F'\t' 'NR==2{print $2"\t"$3-1"\t"$4"\tinterval"}' ${LIFT}/region.tsv > region.bed
LEN=$(awk '{print $3-$2}' region.bed); cat region.bed; echo "${LEN} bp"
```

## 2. First wave: lift, stitch, filter
```
GAP=5000        # WormSynteny "gap value", bp: fragments closer than this are joined
FILTER=5        # WormSynteny "filtering percentage": joined fragments shorter than this % of the interval are dropped
mkdir -p first frags
cut -f3 ${GENOMES} | tail -n +2 | while read G; do
    if [ "${G}" = "${REF}" ]; then cp region.bed frags/${G}.bed; continue; fi
    halLiftover ${HAL} ${REF} region.bed ${G} first/${G}.bed 2> first/${G}.log
    sort -k1,1 -k2,2n first/${G}.bed | bedtools merge -d ${GAP} -i - \
      | awk -v min=$(( LEN * FILTER / 100 )) '$3-$2 >= min {print $0"\tfrag"}' > frags/${G}.bed
    printf '%s\t%s raw\t%s kept\t%s kb kept\n' ${G} $(wc -l < first/${G}.bed) $(wc -l < frags/${G}.bed) \
      $(awk '{s+=$3-$2} END{printf "%.0f", s/1000}' frags/${G}.bed)
done | tee first_wave.tsv
```
- **`kept`** = fragments that survived stitching and filtering. A genome with **0 kept**
  has no row in the figure. Raise `GAP` or lower `FILTER` and re-run this block only.
- **Many kept fragments** = the interval is broken up in that genome. The longest becomes
  the row and the rest are `(filtered)` extras, which are drawn only with `--extras`.
  WormSynteny's advice applies: raising the gap value merges neighbouring fragments into
  the row, and a fragment far away is a genuine second site.
- **Choosing `GAP`:** WormSynteny's examples use 200–500 bp for single worm genes. Across
  squamates and a 300 kb interval, 5 kb joins most of a collinear region without bridging
  real rearrangements. Compare its rows with step 9, which uses 1 kb tiles chained within
  50 kb.

## 3. Second wave: links between neighbouring rows
```
mkdir -p pairs
cut -f3 ${GENOMES} | tail -n +2 > hal_order.txt
paste <(head -n -1 hal_order.txt) <(tail -n +2 hal_order.txt) | while read A B; do
    [ -s frags/${A}.bed ] || { echo "[warn] ${A}: no fragments, no links ${A} -> ${B}"; continue; }
    halLiftover --outPSL ${HAL} ${A} frags/${A}.bed ${B} pairs/${A}__${B}.psl 2> pairs/${A}__${B}.log
    printf '%s -> %s\t%s PSL lines\n' ${A} ${B} $(wc -l < pairs/${A}__${B}.psl)
done
```
As in WormSynteny, each row's fragments are lifted into the row below. Links therefore join
neighbours only, and something present in two non-adjacent rows but missing from the one
between them has no link.

## 4. Tables and plot
```
python3 ${SCRIPTDIR}/14a_wormsynteny_tables.py --genomes ${GENOMES} --lift ${LIFT} \
    --frag-dir frags --pair-dir pairs --species-table ../species_toga2.tsv \
    --gap ${GAP} --extras --out wormsynteny 2> wormsynteny_tables.log
cat wormsynteny_tables.log
 
mamba activate genespace_env
Rscript ${SCRIPTDIR}/14b_plot_wormsynteny.R --in wormsynteny --out wormsynteny_region \
    --width 10 --height 7
```
Read the log:
- **`<code>: row = <contig>:<start>-<end> (N kb); M other fragment(s)`**: the row drawn for
  each genome. Its contig should be step 9's `main_contig` (`hal_lift/targets.tsv`).
- **`no TOGA2 gene on its fragment(s)`**: the row is drawn bare. TOGA2 contigs are matched
  to HAL contigs with `.`, `_` and `__` treated alike (`NC_086153__1` = `NC_086153_1`), so
  this should only happen when a genome has no TOGA2 run or no annotation in the window.
- **`A -> B: N aligned block(s) … -> L link(s)`**: second-wave blocks between the two rows,
  chained into links with the same `GAP`. Zero links between two rows that both have
  fragments means their fragments do not align to each other directly. That can happen
  even when both align to HelSus.

Plot options: `--labels named` (default: step-5 names only) / `all` / `none`;
`--legend named` (default: the exendin array) / `all` (every orthogroup);
`--link-colour grey70`. Leave out `--extras` in 14a to hide the `(filtered)` rows, the
equivalent of WormSynteny's "Hide" button.

## Reading it
- **Rows** = one stitched fragment per genome, oriented by `gggenomes::sync()` to match
  its links. They are **not** to a common scale with the interval: each is as long as
  that genome's fragment.
- **Genes** = that genome's own TOGA2 predictions on the fragment, coloured by the HelSus
  gene they were projected from. The same colour in two rows means orthologous through
  HelSus. **Hollow** = TOGA2 status not intact.
- **Links** = sequence aligned between the two neighbouring rows, chained within `GAP`.
  A twisted link is an inversion; a gap in the links is sequence that does not align
  between those two genomes.
- **How it differs from step 9:** step 9 places every genome on HelSus coordinates via
  tiles, cuts long gaps with `//`, and can show several contigs per row. WormSynteny
  shows each genome's own fragment at its own length, with one fragment per row.

## Outputs
| File | Content |
|---|---|
| `region.bed` | the query (BED4) |
| `first/<hal>.bed`, `frags/<hal>.bed`, `first_wave.tsv` | first wave: raw lift; stitched and filtered fragments; counts |
| `pairs/<A>__<B>.psl` | second wave (`halLiftover --outPSL`) |
| `wormsynteny_seqs.tsv` | gggenomes `seqs`: one row per fragment (main or filtered), with its contig coordinates |
| `wormsynteny_genes.tsv` | gggenomes `genes`: TOGA2 genes on the fragments, orthogroup, TOGA2 status, FULL/PARTIAL coverage (as WormSynteny's "Aligned genes" table) |
| `wormsynteny_links.tsv` | gggenomes `links` |
| `wormsynteny_region.{pdf,png}` | the figure |
