# Synteny from NCBI annotations (TOGA2 genomes only)
**Question:** what does each assembly's own NCBI annotation put between SLC4A8 and AQP5?
This is the baseline; steps 3 and 4 redraw it from TOGA2 and from a sequence search in the
same layout and colours.

**Taxa.** Only the genomes in the TOGA2 run (`query_*` of `target_HelSus`), plus chicken as
the reference row. `species.tsv` (step 1) also holds every Cactus-tree tip; those without a
TOGA2 projection are left out here, so this figure and the step-3 TOGA2 figure have the same
rows. The Cactus tree is still the backbone; the engine prunes it to the species drawn.

**Colour = gene name.** `name_group()` (`lib_exendin_names.py`) reduces each description to
a name: `exendin-2` and `exendin-2-like` share a colour; `glucagon-like` is GCGL. Colours
are fixed per name across all figures.

Run the set-up block of `0_Environment_management.md` first.
```
cd ${SCRIPTDIR} && git pull && cd ${WORK}
```

## 1. Restrict the species table to the TOGA2 genomes
```
KEEP=GalGal     # non-TOGA2 rows to keep (chicken: reference row, step-4 outgroup). KEEP= for none

awk -F'\t' -v keep="${KEEP}" '
  BEGIN{n=split(keep,k,","); for(i=1;i<=n;i++) K[k[i]]=1}
  /^#/{next}
  !h{for(i=1;i<=NF;i++) c[$i]=i; h=1; print; next}
  $c["toga2_query"]!="" || ($c["code"] in K)' species.tsv > species_toga2.tsv

# every TOGA2 query has a row, and no row claims a query TOGA2 does not have
ls -d ${TOGADIR}/target_HelSus/query_* | sed 's#.*/query_##' | sort > toga2_queries.txt
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) c[$i]=i; next} $c["toga2_query"]!=""{print $c["toga2_query"]}' \
    species_toga2.tsv | sort > table_queries.txt
comm -3 toga2_queries.txt table_queries.txt      # must print nothing

# every TOGA2 genome is pinned to the assembly TOGA2 ran on
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) c[$i]=i; next}
  $c["toga2_query"]!="" && $c["source"]!="pinned"{print $c["code"]"\t"$c["accession"]"\t"$c["source"]}' \
    species_toga2.tsv                            # must print nothing

wc -l < toga2_queries.txt                        # number of TOGA2 genomes
cut -f1-5 species_toga2.tsv | column -t
```
- **`comm` prints a code** in the left column = a TOGA2 query with no row (`species.tsv` was
  built without `--toga2-run`, or the query's code does not follow `Genus_species` → `GenSpe`):
  add a row to `1_species_pinned.tsv` and re-run step 1.2. In the right column = a row whose
  `toga2_query` is not in this run.
- **The pinning check prints a code** = its accession was picked by NCBI ranking, not taken
  from the TOGA2 run. If it is not the assembly TOGA2 used, contig names differ from step 3.
  Pin the TOGA2 assembly in `1_species_pinned.tsv`, then re-run steps 1.2–1.3 and this step.
- `species.tsv` is unchanged; only this step reads `species_toga2.tsv`.

## 2. Family and boundary genes in every annotation
```
cd ${WORK}/synteny_ncbi
python3 ${SCRIPTDIR}/2a_ncbi_exendin_loci.py --species-table ../species_toga2.tsv --genomedir ${GENOMEDIR} \
    --outdir ncbi_loci 2> ncbi_loci.log
cat ncbi_loci.log
cut -f1 ncbi_loci/ncbi_loci.tsv | sort -u | grep -vx species | wc -l    # = TOGA2 genomes + KEEP rows
awk -F'\t' 'NR==1||$2!="boundary"' ncbi_loci/ncbi_loci.tsv | cut -f1,4-8,12,15 | column -t
```
A gene is a family gene if its symbol, description or product matches
`--pattern 'exendin|glucagon'` and not `--exclude receptor` (GLP1R/GCGR are left out). For
each species the log gives:
- **The boundary pair** used and how many family genes lie between. `[warn] ... chosen by
  description` = no gene is called SLC4A8 (or AQP5) exactly and a `-like` LOC was used:
  check it.
- **`[warn] ... on the same contig but N Mb apart`** = both boundaries on one chromosome
  but more than `--max-interval-kb` (5 Mb) apart with nothing between: a rearrangement or a
  mis-assigned gene. Each boundary is then used on its own.
- **No pair on one contig** (fragmented assembly): each boundary is chosen by the family
  genes near it; those genes are `near_boundary`.
- **Placement** of every family gene: `interval` (between the pair, ± `--pad-kb` 500, as
  HelSus has copies 118 and 379 kb beyond GCGL), `near_boundary`, `unplaced` (contig with
  no boundary; own track), `elsewhere` (not drawn: normally proglucagon). Check them:
  `awk -F'\t' '$15=="elsewhere"' ncbi_loci/ncbi_loci.tsv | cut -f1,4-8`
- **`no NCBI annotation`** (e.g. HelCha, the JANEZZ contigs): expected for TOGA2 genomes that
  NCBI never annotated. These are empty rows here and are covered by steps 3 and 4.
- **`[warn] ... gives no family group`**: matched the pattern but not a known name.

`ncbi_exendin_aa.fa` (family proteins) and `ncbi_boundary_aa.fa` now hold the TOGA2 genomes
and chicken only; chicken's GCGL protein in it is the gene-tree outgroup.

## 3. Plot (publication figure)
Only the SLC4A8–AQP5 interval is drawn. Text is 7 pt at the final size; there is no title,
no contig names and no notes, since those go in the caption.
```
python3 ${SCRIPTDIR}/2b_plot_ncbi_synteny.py --loci ncbi_loci/ncbi_loci.tsv \
    --species-table ../species_toga2.tsv --species-tree ${SPTREE} --top GalGal \
    --svg ncbi_synteny.svg --pdf ncbi_synteny.pdf --png ncbi_synteny.png \
    --layout ncbi_synteny.tsv 2> ncbi_synteny.log
cat ncbi_synteny.log
grep "not drawn" ncbi_synteny.log          # every family gene left out, and why
```
- **One contig carries both SLC4A8 and AQP5:** the species gets one line, from SLC4A8 to AQP5.
  Genes that 2a counted as `interval` only because of its `--pad-kb` 500 margin are **not
  drawn**. They are listed as `outside the SLC4A8-AQP5 interval`. Check HelSus: its copies
  118 and 379 kb from GCGL may fall here.
- **SLC4A8 and AQP5 are on different contigs:** the species gets two lines. The SLC4A8 contig
  is left-aligned and the AQP5 contig right-aligned, so each sits at its end of the
  interval. Each line shows only the family genes on its interval side, within `--near-kb`
  (1000). Family genes on the far side are logged as `on the far side of ...`. A species
  with only one boundary annotated gets one such line.
- **Orientation:** SLC4A8 is always on the left. A contig with only one boundary is oriented
  by that gene's strand in HelSus (`--reference-query`). The first log line gives the
  strands used.
- **Family genes on a contig with no boundary** are left out (`--unplaced` draws them).
- **Tree:** the `pruned` line lists the Cactus taxa without a TOGA2 projection, which is
  expected. `not in the tree` should name only GalGal, which `--top` draws above the tree.
  Any TOGA2 genome named there needs a `--tree-map` row (`tree_tip<TAB>code`).
- **Reading the figure:** filled = protein-coding, hollow = pseudogene. Arrows give the
  strand after orientation. `//` = a stretch longer than `--max-gap-kb` (1000) with nothing
  drawn, cut to `--break-kb` (20). Boundary genes are clipped to `--clip-kb` (15) on their
  outer side. "no NCBI annotation" means the assembly was never annotated, not that the
  genes are absent.
- **Size:** `--width-mm 180` (double column, default) or `--width-mm 88` (single column).
  `--font-pt` sets every label. `--full-names` writes *Heloderma suspectum* instead of
  *H. suspectum*, and `--labels` puts a label on each gene. The PDF is vector, for
  submission. The PNG is `--dpi` 600. Hover text in the SVG and the `--layout` table keep
  the true coordinates.
