# Synteny from the Cactus alignment
**Question:** where does each part of the HelSus SLC4A8–AQP5 interval sit in every other
genome, by whole-genome alignment rather than by annotation? Steps 2 and 3 can only show
genes that an annotation (NCBI) or a projection (TOGA2) calls. The alignment also shows
what is in between: sequence that is conserved, lost, inverted, moved to another contig,
or duplicated, whether or not anyone annotated a gene there.

**Approach.**
- **Region:** the HelSus interval is SLC4A8 to AQP5 ± 25 kb.
- **Lift:** it is cut into 1 kb tiles. Each tile, and each HelSus gene in the interval, is
  lifted into every other genome of the HAL with `halLiftover`.
- **Blocks:** lifted tiles that continue in order on one contig and strand are chained into
  syntenic blocks.
- **Duplicates:** a tile that also lifts to a second locus is a secondary alignment, i.e.
  duplicated sequence.
- **Figure:** a river plot (as in GENESPACE / SVbyEye). One bar per genome in tree order,
  the lifted HelSus genes on each bar, and ribbons joining each segment's positions in
  neighbouring genomes.

Run the set-up block of `0_Environment_management.md` first.
```
module load cactus/3.3.0          # halLiftover, halStats
mkdir -p ${WORK}/synteny_cactus && cd ${WORK}/synteny_cactus
ls -lh ${CACTUSDIR}/*.hal
HAL=$(ls ${CACTUSDIR}/*.hal | head -1); echo ${HAL}
```

## 1. Check the HAL before lifting
```
halStats --tree ${HAL}
halStats --genomes ${HAL}
halStats --sequences HelSus ${HAL} | tr ',' '\n' | head
grep -v '^#' ${HELSUS_EVIANN_GFF} | cut -f1 | sort -u | head
```
- **Genome names:** the leaves of the tree are the genomes; names starting `Anc` are
  ancestral reconstructions and are skipped. Leaves may be codes (`HelSus`) or
  `Genus_species`; both are matched to the 3+3 codes of `species.tsv`. If the HelSus leaf
  has another name, pass it with `--ref`.
- **Chromosome names:** the HelSus chromosome names in the HAL must match the EviAnn GFF
  (`chr2`). If they differ (e.g. `NC_139459.1`), pass the HAL name with `--hal-chrom`.
  9a stops and lists the HAL names if it cannot find the chromosome.

## 2. Lift the interval into every genome
```
python3 ${SCRIPTDIR}/9a_hal_lift_region.py --hal ${HAL} --ref HelSus \
    --gff ${HELSUS_EVIANN_GFF} --boundaries SLC4A8,AQP5=LOC_00012326 \
    --pad-kb 25 --tile-kb 1 --threads 8 --outdir hal_lift 2> hal_lift.log
cat hal_lift.log
cat hal_lift/ref_genes.tsv | cut -f1-5,9
```
Read the log:
- **`region chr2:…`**: the HelSus interval and the SLC4A8 / AQP5 genes used. If several
  SLC4A8-like genes are listed, the pair closest together is used. To force one, pass
  `SLC4A8=<gene ID>`.
- **One line per genome:** `N kb of M kb aligned (X%), B primary block(s) on C
  contig(s), main <contig>; S secondary block(s)`.
  - Low % = little of the interval is alignable in that genome.
  - Many blocks or contigs = rearrangement or a fragmented assembly.
  - Secondary blocks = duplicated sequence.
- **`ref_genes.tsv`** lists every HelSus gene in the interval, with the name and family
  group used to colour it.

Each `halLiftover` call is quick for a region this size. The whole step takes minutes, but run it
in `sinteractive`, not on the login node.

## 3. River plot
```
python3 ${SCRIPTDIR}/9b_plot_hal_river.py --indir hal_lift \
    --species-table ../species.tsv --species-tree ${SPTREE} \
    --svg hal_river.svg --pdf hal_river.pdf --png hal_river.png \
    --layout hal_river_layout.tsv 2> hal_river.log
cat hal_river.log

# GENESPACE-style colours: every ribbon coloured by its position along HelSus
python3 ${SCRIPTDIR}/9b_plot_hal_river.py --indir hal_lift --ribbon-colour position \
    --species-table ../species.tsv --species-tree ${SPTREE} \
    --svg hal_river_position.svg --pdf hal_river_position.pdf --png hal_river_position.png 2> hal_river_position.log

# with duplicated (secondary) alignments as faded ribbons
python3 ${SCRIPTDIR}/9b_plot_hal_river.py --indir hal_lift --secondary \
    --species-table ../species.tsv --species-tree ${SPTREE} \
    --svg hal_river_dups.svg --pdf hal_river_dups.pdf --png hal_river_dups.png 2> hal_river_dups.log
```
Every genome in the alignment is drawn, not only the TOGA2 set. The plot uses the
publication layout of steps 2 and 3: 180 mm wide (`--width-mm 88` for one column), 7 pt
text, and a vector PDF. Genomes with no alignment to the interval are left out of the
river and named in the log (`not aligned to the interval`).

## Reading it
- **Bars** = one per genome, in species-tree order. A bar is the stretch of a contig that
  the HelSus interval aligns to, with HelSus at its place in the tree.
  - Several bars in one row = the interval is split over several contigs (assembly break
    or rearrangement).
  - Each genome is flipped to match HelSus and shifted so homologous sequence lines up
    roughly vertically.
  - `//` on a bar = an unaligned stretch longer than `--max-gap-kb` (100), cut short.
  - Contigs with less than `--min-frac` (5%) of a genome's aligned sequence are not drawn,
    and the log lists them.
- **Genes on the bars** = the HelSus genes lifted through the alignment.
  - Colours: exendins and GCGL keep their usual colours, every other gene is grey, and
    SLC4A8/AQP5 are dark grey.
  - Hollow = less than half of the gene aligned. Dashed = a second copy of the gene's
    sequence.
  - A missing arrow means the HelSus gene does not align there, not that the genome has
    no such gene. Steps 3 and 4 test for diverged copies.
- **Ribbons** connect each neighbouring pair of genomes: one ribbon per run of the HelSus
  interval that is collinear in both.
  - Grey = same orientation. Brown, twisting = inverted between the two genomes.
  - Exendin and GCGL ribbons are in the gene's colour, so the array can be followed down
    the tree.
  - Where ribbons stop, the sequence between them did not align in one of the two
    genomes.
  - `--ribbon-colour position` colours every ribbon by its position along HelSus (purple
    at SLC4A8 → yellow at AQP5), as in GENESPACE. That makes long-range moves easy to
    follow.
- **Ribbons only join neighbours.** Sequence present in two genomes that are not
  adjacent, but absent from the one between them, has no ribbon. Read the presence of a
  segment from the bars of each genome, and the ribbons for order and orientation.
- **`--secondary`** adds faded ribbons to the places where a HelSus segment aligns a second
  time in the next genome. Over the exendin array, this fan-out is how tandem copies show.

## Outputs
| File | Content |
|---|---|
| `hal_lift/region.tsv` | HelSus chromosome, region, boundary genes |
| `hal_lift/ref_genes.tsv` | HelSus genes in the region |
| `hal_lift/tile_hits.tsv` | where each 1 kb HelSus tile aligns in each genome (primary and secondary); the ribbons are drawn from this |
| `hal_lift/blocks.tsv` | every syntenic block: genome, contig, target and HelSus coordinates, strand, aligned bp, primary/secondary |
| `hal_lift/genes_lifted.tsv` | where each HelSus gene lands in each genome, and how much of it aligned |
| `hal_lift/targets.tsv` | per genome: aligned bp, fraction of the interval, main contig |
| `hal_lift/lifted/` | raw `halLiftover` output |
| `hal_river*.{svg,pdf,png}`, `hal_river_layout.tsv` | figures; every bar and ribbon with its HelSus coordinates |
