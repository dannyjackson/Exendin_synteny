# GENESPACE riparian plots from the Cactus alignment
**Question:** what does the HelSus SLC4A8–AQP5 interval look like at chromosome scale, in the
GENESPACE style? Step 9 draws the interval itself, with TOGA2 genes on 1 kb tiles. This
step draws whole chromosomes and the few Mb around the interval, from synteny blocks
called directly on the HAL. No genes or orthogroups are involved.

**Approach.**
- **Blocks:** `halSynteny` (HAL tools, Krasheninnikova et al. 2020, *GigaScience*
  giaa047) calls synteny blocks for each pair of neighbouring rows of the figure,
  straight from the HAL, as PSL. It takes `--queryGenome` = upper row and `--targetGenome` = lower row.
- **Conversion:** `13a_hal_to_blk.py` reformats the PSL into GENESPACE's block table
  (`genome1, genome2, chr1, chr2, startBp1, endBp1, startBp2, endBp2, blkID`). It also does
  two things GENESPACE would otherwise do from genes:
  - it labels each block with the HelSus chromosome it traces back to, for the colours;
  - it projects the HelSus interval through the blocks, for the highlight.
- **Plot:** GENESPACE's own `riparian_engine()` (`13b_plot_riparian.R`). Its usual front end,
  `plot_riparian()`, needs a full GENESPACE run (OrthoFinder + MCScanX on proteins).
  `riparian_engine()` is the drawing core under it and takes a block table directly.
- **Row order** = the species-tree order of step 9, with HelSus as the reference for the colours.

Two figures:
| Figure | halSynteny run | Shows |
|---|---|---|
| `genespace_genome` | whole genomes, coarse blocks (`RES_GENOME`, 250 kb) | every chromosome, braids coloured by HelSus chromosome, the interval in purple |
| `genespace_region` | only the contigs carrying the interval (step 9), fine blocks (5 kb) | the interval ± `--window-kb` (2 Mb) in every genome |

Run the set-up block of `0_Environment_management.md` first. This step needs the outputs of
step 9 (`synteny_cactus/hal_lift`).

## 0. Software (once)
`halSynteny` is in the Cactus module, next to `halLiftover`. GENESPACE is an R package.
Only its plotting functions are used, so OrthoFinder and MCScanX are **not** needed.
```
module load cactus/3.3.0
which halSynteny && halSynteny --help 2>&1 | head -20

source myconda
mamba create -n genespace_env -c conda-forge -c bioconda --strict-channel-priority \
  r-base=4.3 r-data.table r-ggplot2 r-dbscan r-igraph r-r.utils r-devtools \
  bioconductor-biostrings bioconductor-rtracklayer
mamba activate genespace_env
Rscript -e 'devtools::install_github("jtlovell/GENESPACE", upgrade = "never")'
Rscript -e 'library(GENESPACE); packageVersion("GENESPACE"); exists("riparian_engine")'   # >= 1.3.1, TRUE
```
If `install_github` cannot reach GitHub from a compute node, run it on Helix
(`ssh helix.nih.gov`), where outbound access is allowed.

## 1. Genomes and row order
```
mkdir -p ${WORK}/genespace_cactus && cd ${WORK}/genespace_cactus
mamba activate exendin_env
HAL=$(ls ${CACTUSDIR}/*.hal | head -1); echo ${HAL}
LIFT=../synteny_cactus/hal_lift                      # step 9, 9a --outdir
DROP=AnnSte,BarImb,ElgCoe                            # as in step 9

python3 ${SCRIPTDIR}/13a_hal_to_blk.py genomes --indir ${LIFT} \
    --species-tree ${SPTREE} --species-table ../species_toga2.tsv --only-listed \
    --drop-species ${DROP} --out genomes.tsv 2> genomes.log
cat genomes.log; cat genomes.tsv
```
- `genomes.tsv`: one row per genome in figure order, with its HAL name and the HAL
  contig(s) the interval aligns to in step 9 (`interval_contigs`).
- Leave out `--only-listed` / `--drop-species` to draw every genome in the HAL.
- A genome flagged **`interval not aligned in step 9`** has no contigs to run the region
  halSynteny on, so the zoom figure breaks at that row (see 2.2).

## 2. halSynteny, one run per neighbouring pair
`halSynteny` traverses the whole alignment of the two genomes, so a whole-genome pair takes
anywhere from minutes to a few hours. Run the pairs in parallel with `swarm`.

### 2.1 Whole genomes (for `genespace_genome`)
```
RES_GENOME=250000        # --minBlockSize and --maxAnchorDistance, bp
mkdir -p psl_genome
awk -F'\t' -v hal=${HAL} -v r=${RES_GENOME} 'NR>1{h[NR]=$3} END{for(i=2;i<NR;i++)
    printf "halSynteny --queryGenome %s --targetGenome %s --minBlockSize %d --maxAnchorDistance %d %s psl_genome/%s__%s.psl\n",
    h[i], h[i+1], r, r, hal, h[i], h[i+1]}' genomes.tsv > halsynteny_genome.swarm
cat halsynteny_genome.swarm
swarm -f halsynteny_genome.swarm -g 32 -t 2 --time 24:00:00 --module cactus/3.3.0 \
    --logdir swarm_logs_genome --job-name halsyn_genome
```
- **Resolution:** `RES_GENOME` sets both the smallest block and the largest gap inside one.
  250 kb gives chromosome-scale braids that GENESPACE can draw in reasonable time. 1 Mb
  (`halSynteny`'s whole-genome example) is cleaner across snakes and lizards. 100 kb shows
  more rearrangements but produces tens of thousands of polygons.
- **Memory:** `-g 32` is a starting point. If a job fails with an out-of-memory error, re-run
  only that line with more.

### 2.2 Interval contigs only (for `genespace_region`)
`--queryChromosome` restricts each run to the contigs that carry the interval in the
upper genome (from step 9). These runs are quick, so fine blocks are affordable.
```
RES_REGION=5000           # halSynteny's default
mkdir -p psl_region
awk -F'\t' -v hal=${HAL} -v r=${RES_REGION} 'NR>1{h[NR]=$3; c[NR]=$4} END{for(i=2;i<NR;i++){
    n = split(c[i], cs, ",")
    if (!n) printf "echo [warn] no interval contigs for %s: region chain breaks here >&2\n", h[i]
    for (j = 1; j <= n; j++)
      printf "halSynteny --queryGenome %s --targetGenome %s --queryChromosome %s --minBlockSize %d --maxAnchorDistance %d %s psl_region/%s__%s.%s.psl\n",
      h[i], h[i+1], cs[j], r, r, hal, h[i], h[i+1], cs[j]}}' genomes.tsv > halsynteny_region.swarm
cat halsynteny_region.swarm
swarm -f halsynteny_region.swarm -g 16 -t 2 --time 4:00:00 --module cactus/3.3.0 \
    --logdir swarm_logs_region --job-name halsyn_region
```
Afterwards, check that every pair produced a non-empty PSL:
```
for f in psl_genome/*.psl psl_region/*.psl; do printf '%s\t%s\n' "$(wc -l < $f)" $f; done | sort -n | head
grep -il "error\|killed" swarm_logs_*/*.e 2>/dev/null
```

## 3. PSL → GENESPACE blocks
```
python3 ${SCRIPTDIR}/13a_hal_to_blk.py blk --genomes genomes.tsv --psl-dir psl_genome \
    --region ${LIFT}/region.tsv --min-chr-mb 1 --out blk_genome.csv 2> blk_genome.log
python3 ${SCRIPTDIR}/13a_hal_to_blk.py blk --genomes genomes.tsv --psl-dir psl_region \
    --region ${LIFT}/region.tsv --window-kb 2000 --out blk_region.csv 2> blk_region.log
cat blk_genome.log blk_region.log
```
Read the logs:
- **`<code>: interval projected to <contig>:<start>-<end>`**, one line per genome: where
  the HelSus interval lands after chaining through the blocks. **Compare these with step 9**
  (`hal_lift/targets.tsv`, `main_contig`): they should name the same contigs. If one
  doesn't, a block between that genome and its neighbour was mis-chained. Check it in
  `blk_genome.csv` before using the figure.
- **`A -> B: N block(s), X Mb of A, M traced to a HelSus chromosome`**: the per-pair
  coverage. A pair with few Mb is a poor alignment or an assembly-level mismatch. Blocks
  not traced back to HelSus are drawn grey.
- **`no PSL ... the river breaks here`**: a halSynteny job is missing (step 2).
- `--min-chr-mb 1` drops contigs whose blocks reach less than 1 Mb, which keeps unplaced
  scaffolds out of the whole-genome figure. Use 0 to keep everything.

## 4. Plot
```
mamba activate genespace_env
Rscript ${SCRIPTDIR}/13b_plot_riparian.R --blk blk_genome.csv --genomes genomes.tsv \
    --out genespace_genome --width 10 --height 7 --min-chr-len 1e6
Rscript ${SCRIPTDIR}/13b_plot_riparian.R --blk blk_region.csv --genomes genomes.tsv \
    --out genespace_region --width 7 --height 7 --min-chr-len 1e4

# only the interval coloured, everything else grey (plot_riparian's highlightBed look)
Rscript ${SCRIPTDIR}/13b_plot_riparian.R --blk blk_genome.csv --genomes genomes.tsv \
    --out genespace_genome_mono --mono --min-chr-len 1e6
```
Options: `--highlight "#7B3294"` (interval colour), `--inversion-colour "#A6611A"`
(inverted braids in their own colour), `--braid-alpha 0.6`, `--font 5` (chromosome
labels, pt), `--dpi 600`.

## Reading it
- **Rows** = genomes in species-tree order. **Bars** = chromosomes/contigs, as far as
  halSynteny blocks reach (physical scale, Mbp). The figure therefore shows the
  *aligned* extent of each chromosome, not its assembled length.
- **Braid colour** = the HelSus chromosome the block traces back to. Colours are chained
  through every row, so a HelSus chromosome can be followed down the tree, as in GENESPACE.
  Grey = no chain back to HelSus.
- **Purple** = the SLC4A8–AQP5 interval, projected block by block. Within a block the
  projection is linear. At 250 kb resolution it is accurate to the block, not to the gene;
  `genespace_region` and step 9 give the fine scale.
- **A twisted braid** = an inversion between the two neighbouring genomes.
- **Region figure:** each bar is a window of ±2 Mb around the projected interval, renamed
  `<contig>:<window start, Mb>` with its coordinates starting at 1. A genome with
  two bars has the interval split over two contigs or two distant places on one contig, as
  in step 9 (*A. sagrei*).
- **Ribbons only join neighbours**, as in step 9: something present in two non-adjacent
  genomes but missing between them has no braid.
- **Not GENESPACE synteny.** GENESPACE normally builds blocks from syntenic orthologs.
  Here the blocks come from the whole-genome alignment, so they include non-genic
  sequence, and no gene annotation is shown. For genes on the interval, use step 9.

## Outputs
| File | Content |
|---|---|
| `genomes.tsv` | row order, HAL names, interval contigs |
| `halsynteny_{genome,region}.swarm`, `psl_{genome,region}/` | halSynteny commands and their PSL |
| `blk_{genome,region}.csv` | GENESPACE block tables: one row per block, `refChr` = HelSus chromosome, `highlight` = 1 for the interval |
| `genespace_{genome,region}.{pdf,png}` | the riparian plots |
| `genespace_*_chromosomes.tsv` | riparian_engine's chromosome table (genome, chr, plotted length and position) |
