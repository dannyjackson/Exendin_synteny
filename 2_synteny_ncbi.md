# Synteny of exendin-family genes from NCBI annotations
**Question:** what does each assembly's own NCBI annotation put between SLC4A8 and AQP5?
This figure is the baseline. Steps 3 and 4 redraw it from TOGA2 and from a sequence search,
in the same layout and colours, so the three can be compared row by row.

**Colour = gene name.** `name_group()` in `lib_exendin_names.py` reduces each NCBI
description to a name: `exendin-2` and `exendin-2-like` share a colour, as do `exendin-4`
and `exendin-4-like`. `glucagon-like` (GCGL, the ancestral gene of the array) gets its own
colour. The colours are fixed per name, so they are the same in all three figures.

Run the set-up block of `0_Environment_management.md` first.
```
cd ${WORK}/synteny_ncbi
```

## 1. Family and boundary genes in every annotation
```
python3 ${SCRIPTDIR}/2a_ncbi_exendin_loci.py \
    --species-table ../species.tsv --genomedir ${GENOMEDIR} \
    --outdir ncbi_loci 2> ncbi_loci.log
cat ncbi_loci.log
awk -F'\t' 'NR==1||$2!="boundary"' ncbi_loci/ncbi_loci.tsv | cut -f1,4-8,12,15 | column -t
```
Check the log for each species:
- **The boundary pair.** The line gives the two genes used and how many family genes lie
  between them. A `[warn] ... chosen by description` means no gene is called SLC4A8 (or AQP5)
  exactly, and a `-like` LOC was used. Check that it is the right gene.
- **No pair on one contig** (fragmented assembly): each boundary is chosen on its own, by how
  many family genes lie near it. Its family genes are `near_boundary`.
- **Placements.** `interval` genes are drawn on the main track. `unplaced` ones (exendin-named,
  on a contig with no boundary) get their own track. `elsewhere` genes are not drawn: this is
  normally proglucagon on its own chromosome. Check them anyway:
  ```
  awk -F'\t' '$15=="elsewhere"' ncbi_loci/ncbi_loci.tsv | cut -f1,4-8
  ```
- **`[warn] ... gives no family group`**: a gene matched `--pattern` but its name is not one
  the grouping knows. Look at it before deciding whether to widen `name_group()`.

Defaults worth knowing: `--pattern 'exendin|glucagon' --exclude receptor` (so GLP1R/GCGR are
left out), and `--pad-kb 500`. Genes up to 500 kb outside the boundary pair still count as
interval, because HelSus has copies 118 and 379 kb from ANC.

## 2. Plot
```
python3 ${SCRIPTDIR}/2b_plot_ncbi_synteny.py \
    --loci ncbi_loci/ncbi_loci.tsv --species-tree ${SPTREE} \
    --svg ncbi_synteny.svg --png ncbi_synteny.png --layout ncbi_synteny.tsv 2> ncbi_synteny.log
cat ncbi_synteny.log
```
- **Drawing.** Filled = protein-coding gene; hollow = pseudogene or other non-coding biotype.
  Hover in the SVG for the symbol, description, placement and protein. The drawing code is
  the old repo's `4d_plot_exendin_synteny.py`, copied unchanged as `lib_synteny_engine.py`.
  Orientation (SLC4A8 on the left), clipping of long boundary genes (`//`) and the tree work
  as in the old 4d / 2g / 3b figures. `--branch-lengths` draws the Cactus tree with lengths.
- **Tree matching.** Tips are matched to species codes as `Genus_species` → `GenSpe`.
  Species not in the tree (outgroups) are drawn below it; place them with `--tree-map`.
- **Empty rows.** A species with no NCBI annotation is drawn as an empty row saying so.
  An empty row does not mean the gene is absent: steps 3 and 4 are how to find out.
- To see every matched gene, including `elsewhere` ones:
  `--placements interval,near_boundary,unplaced,elsewhere`.
