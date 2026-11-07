# Synteny from NCBI annotations
**Question:** what does each assembly's own NCBI annotation put between SLC4A8 and AQP5?
This is the baseline; steps 3 and 4 redraw it from TOGA2 and from a sequence search in the
same layout and colours.

**Colour = gene name.** `name_group()` (`lib_exendin_names.py`) reduces each description to
a name: `exendin-2` and `exendin-2-like` share a colour; `glucagon-like` is GCGL. Colours
are fixed per name across all figures.

Run the set-up block of `0_Environment_management.md` first.
```
cd ${WORK}/synteny_ncbi
```

## 1. Family and boundary genes in every annotation
```
python3 ${SCRIPTDIR}/2a_ncbi_exendin_loci.py --species-table ../species.tsv --genomedir ${GENOMEDIR} \
    --outdir ncbi_loci 2> ncbi_loci.log
cat ncbi_loci.log
awk -F'\t' 'NR==1||$2!="boundary"' ncbi_loci/ncbi_loci.tsv | cut -f1,4-8,12,15 | column -t
```
A gene is a family gene if its symbol, description or product matches
`--pattern 'exendin|glucagon'` and not `--exclude receptor` (GLP1R/GCGR are left out). For
each species the log gives:
- **The boundary pair** used and how many family genes lie between. `[warn] ... chosen by
  description` = no gene is called SLC4A8 (or AQP5) exactly and a `-like` LOC was used:
  check it.
- **`[warn] ... on the same contig but N Mb apart`** = both boundaries on one chromosome
  but more than `--max-interval-kb` (5 Mb) apart with nothing between (AnoSag, ~57 Mb):
  a rearrangement or a mis-assigned gene. Each boundary is then used on its own.
- **No pair on one contig** (fragmented assembly): each boundary is chosen by the family
  genes near it; those genes are `near_boundary`.
- **Placement** of every family gene: `interval` (between the pair, ± `--pad-kb` 500, as
  HelSus has copies 118 and 379 kb beyond GCGL), `near_boundary`, `unplaced` (contig with
  no boundary; own track), `elsewhere` (not drawn: normally proglucagon). Check them:
  `awk -F'\t' '$15=="elsewhere"' ncbi_loci/ncbi_loci.tsv | cut -f1,4-8`
- **`[warn] ... gives no family group`**: matched the pattern but not a known name.

`ncbi_exendin_aa.fa` (family proteins) and `ncbi_boundary_aa.fa` feed step 4; chicken's
GCGL protein in it is the gene-tree outgroup.

## 2. Plot
```
python3 ${SCRIPTDIR}/2b_plot_ncbi_synteny.py --loci ncbi_loci/ncbi_loci.tsv \
    --species-tree ${SPTREE} --top GalGal \
    --svg ncbi_synteny.svg --png ncbi_synteny.png --layout ncbi_synteny.tsv 2> ncbi_synteny.log
cat ncbi_synteny.log
```
- **Drawing.** Filled = protein-coding; hollow = pseudogene or other biotype. Tracks are
  oriented with SLC4A8 on the left; a long boundary gene is clipped (`//` on the gene).
  Hover in the SVG for symbol, description, placement and protein.
- **Rows.** `--top GalGal` draws chicken first, above the tree. Other species not in the
  tree go below it (place them with `--tree-map`). A species with no NCBI annotation is an
  empty row saying so; that is not absence.
- **`//` on a track** = a stretch over `--max-gap-kb` (1000) with nothing drawn, shortened
  to `--break-kb` (20) and labelled with its true length (AnoSag). Separate tracks for one
  species = different contigs, i.e. an assembly break (AnoCar). Hover text and `--layout`
  keep true coordinates.
- `--placements interval,near_boundary,unplaced,elsewhere` draws every matched gene.
