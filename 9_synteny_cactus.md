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
- **Figure:** each genome's blocks are coloured by their position along HelSus, and the
  lifted HelSus genes are drawn on top. HelSus is the colour key.

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

## 3. Plot
```
python3 ${SCRIPTDIR}/9b_plot_hal_synteny.py --indir hal_lift \
    --species-table ../species.tsv --species-tree ${SPTREE} \
    --svg hal_synteny.svg --pdf hal_synteny.pdf --png hal_synteny.png \
    --layout hal_synteny_layout.tsv 2> hal_synteny.log
cat hal_synteny.log
```
All genomes in the alignment are drawn, not only the TOGA2 set. The figure uses the
publication layout of steps 2 and 3: 180 mm wide (`--width-mm 88` for one column), 7 pt
text, and a vector PDF.

## Reading it
- **Colour band** = sequence aligned to the HelSus interval, coloured purple (SLC4A8 end)
  → yellow (AQP5 end).
  - Gradient running the same way as in HelSus: collinear.
  - Gradient reversed inside a block: an inversion relative to HelSus.
  - A colour missing from a genome: that part of the interval did not align there. It is
    lost, too diverged to align, or in an assembly gap.
- **Lines:** each genome is oriented to match HelSus and shifted so its blocks sit under
  the matching part of HelSus, so a vertical line through the figure crosses homologous
  sequence. Several lines for one genome = the interval spans several contigs (the
  assembly is broken, or the region is rearranged). Only contigs with at least
  `--min-frac` (5%) of the genome's aligned sequence are drawn, and the log lists the rest.
  `//` = an unaligned stretch longer than `--max-gap-kb` (100), cut short.
- **Faded band below** = secondary alignment: the same HelSus sequence aligns to a second
  place on that contig. Over the exendin array this is the signature of tandem copies.
  `--no-secondary` hides it.
- **Arrows** = HelSus genes lifted through the alignment, not the genome's own
  annotation.
  - Colours: exendins and GCGL keep their usual colours, every other gene is grey, and
    SLC4A8/AQP5 are dark grey.
  - Hollow = less than half of the gene aligned.
  - Dashed = the gene's secondary lift (a second copy).
  - Absence of an arrow means the HelSus gene does not align there. It does not mean the
    genome has no such gene: a diverged copy can be missed. Steps 3 and 4 test that.
- **Clipped lines:** the log warns if a line runs more than `--overhang` (30%) of the
  interval beyond HelSus and was clipped at the figure edge.

## Outputs
| File | Content |
|---|---|
| `hal_lift/region.tsv` | HelSus chromosome, region, boundary genes |
| `hal_lift/ref_genes.tsv` | HelSus genes in the region |
| `hal_lift/blocks.tsv` | every syntenic block: genome, contig, target and HelSus coordinates, strand, aligned bp, primary/secondary |
| `hal_lift/genes_lifted.tsv` | where each HelSus gene lands in each genome, and how much of it aligned |
| `hal_lift/targets.tsv` | per genome: aligned bp, fraction of the interval, main contig |
| `hal_lift/lifted/` | raw `halLiftover` output |
| `hal_synteny.{svg,pdf,png}`, `hal_synteny_layout.tsv` | figure and the drawn coordinates |
