# Standard names, the named gene tree and the named synteny figure
**Why.** After step 4 each gene carries the name of whichever annotation its best reverse
hit came from (NCBI "exendin-3-like", TOGA2 "exendin-2/3", published). Those names are
neither monophyletic nor consistent between sources. The published names exendin-1 to -4
are **mature-peptide** names, given to venom peptides in 1984–1992 before any genome
existed. This step replaces the inherited names with names that state only what the data
support.

**Working hypothesis.** One ancestral, single-copy exendin-family gene, **GCGL** (Wang et
al. 2012; Irwin's "exendin"), sits between SLC4A8 and AQP5 and has duplicated several times
in squamates.

**Rules** (`5b_name_genes.py`)
| Case | Name | What it claims |
|---|---|---|
| one family gene in the genome | `<species>_GCGL` | the ancestral single-copy state |
| several genes in the genome | `<species>_A`, `_B`, `_C` … from SLC4A8 to AQP5 | position only: gene A of one species need not be orthologous to gene A of another |
| *Heloderma* (HelSus, HelCha) | `_EX1` … `_EX4` or `_GCGL`: the published mature peptide the gene encodes at ≥ 80% (`5a`) | which published peptide the gene encodes; no match ≥ 80% → its positional letter |
| outgroup (chicken NCBI protein) | `GalGal_GCGL` if GCGL-type | not counted |
| published precursor | `Published_EX2/3/4`, unchanged | the published sequence itself |

- **Copies counted**: every curated locus of the species, remnants and pseudogenes
  included, so "single copy" means no second family gene at all.
- **Two *Heloderma* copies matching one peptide**: the better match keeps the plain name,
  the other gets its letter appended (`HelSus_EX4`, `HelSus_EX4_D`); a tie gives both letters.
- **"SLC4A8 to AQP5"**: on the contig with both boundaries, position measured from SLC4A8
  towards AQP5, whichever way the assembly runs. Fragmented assemblies: SLC4A8 contig, then
  unplaced contigs, then AQP5 contig; `order_basis` in `names_names.tsv` records which.

**Scripts.** `5a_assign_nomenclature.py` is the old `2d_assign_nomenclature.py` with one
fix: the old file wrote a literal `...` where `differences` belongs, so its rows were one
column short of the header (old `nomenclature_assignments.tsv` files are shifted from
`differences` on). `5b` applies the rules, `5c` relabels trees and FASTA files, `5d` draws
the named synteny figure; the tree is drawn with `4h`.

Inputs are step 4's curated loci, models and tree (`4_reciprocal_blast.md` steps 4–7).
Run the set-up block of `0_Environment_management.md` first.
```
cd ${WORK}/nomenclature
R4=${WORK}/reciprocal_blast
GG=$(awk -F'\t' '$3=="GalGal"{print $1}' ${R4}/reciprocal/curated_loci.tsv | paste -sd,)
echo "chicken's curated locus: ${GG:-none}"        # the same gene as the outgroup protein
```

## 1. Mature peptide of every gene vs the published peptides
Every curated model (not only tree tips), the outgroup and the published precursors.
```
cat ${R4}/models/curated_aa.fa ${R4}/tree_outgroup.faa \
    <(seqkit grep -r -p '^Published_' ${R4}/tree/tree_aa.fa) > genes_aa.fa
python3 ${SCRIPTDIR}/5a_assign_nomenclature.py --query genes_aa.fa \
    --reference ${SCRIPTDIR}/4_reference_peptides.fa --out-prefix nomenclature 2> nomenclature.log
awk -F'\t' 'NR==1||$1~/^(Hel|Published)/' nomenclature_assignments.tsv | cut -f1-3,5-10 | column -t -s$'\t'
```
- **Verdicts.** "encodes published EXENDIN*i* exactly" = the gene behind that peptide;
  "EXENDIN*i*-type, N differences" = same class (≥ 80%); "no published match" = never
  described as a peptide. Most squamate copies outside *Heloderma* are ~50% to everything
  published (HydCya GCGL 52% to zebra finch GCGL; ShiCro 52% to GLP-1), which is why rule 2
  uses position, not peptide.
- **`differences`.** Only deletions at the C-terminal end (`S29-,R30-,…`) = a truncated
  model, not a different gene: check its `models.tsv` row (`filled`, `q_cover`).
- **Expected *Heloderma* names**: HelSus EX2, EX4 and GCGL (89.7% to zebra finch GCGL);
  HelCha EX1, EX3 and GCGL. HelSus/HelCha `Ex2_3` (`HADGTFTSDLSTYLENQ…`, 58% to GLP-1) match
  nothing published and keep letters, as do the pseudogene copies.
- **`how_found`** ≠ `motif` (cleavage_site / best_window) = the precursor breaks the
  strict motif; check by eye. `n_hormone_domains` counts strict-motif matches only (1 =
  exendin/GCGL, 3 = real proglucagon); diverged squamate peptides can show 0.

## 2. Apply the naming rules
```
python3 ${SCRIPTDIR}/5b_name_genes.py \
    --curated ${R4}/reciprocal/curated_loci.tsv --boundaries ${R4}/reciprocal/boundaries_used.tsv \
    --assignments nomenclature_assignments.tsv ${GG:+--drop $GG} \
    --out-prefix names 2> names.log
cat names.log
column -t -s$'\t' names_names.tsv | less -S
```
Outputs: `names_names.tsv` (the record: old name, new name, rule, position, evidence; keep
it with the tree), `names_map.tsv` (old → new; everything downstream reads it),
`names_groups.tsv` (colour class), `names_rules.tsv` + `names_legend.txt` (tip codes).
- **`--drop`** leaves chicken's partial curated locus out of the counts (it is the outgroup
  protein's gene).
- **Log**: the *Heloderma* lines list each gene's new name (compare with step 1);
  `[warn] ... unplaced or out-of-interval` = that species' letters rest on an assumed
  position (see `order_basis`); `[warn] new name X given to ...` = resolve by hand.
- **Overriding** a name: edit column 2 of `names_map.tsv` and note why in `names_names.tsv`.

## 3. Relabel the trees
Branch lengths and supports are untouched; nothing is rebuilt.
```
python3 ${SCRIPTDIR}/5c_rename_tips.py --map names_map.tsv --outdir renamed \
    ${R4}/tree/exendin_aa.rooted.nwk ${R4}/tree/exendin_aa.aln.fa ${R4}/tree/tree_aa.fa

# codon tree, if built:  ... --outdir renamed ${R4}/tree/exendin_codon.rooted.nwk
python3 ${SCRIPTDIR}/5c_rename_tips.py --map names_map.tsv --outdir renamed \
    ${R4}/tree/exendin_codon.rooted.nwk ${R4}/tree/exendin_aa.aln.fa ${R4}/tree/tree_aa.fa
```
The log lists any tip not in the map; there should be none.

## 4. Named gene tree
```
OG=$(awk -F'\t' 'NR==FNR{m[$1]=$2; next} {n=split($0,a,","); s=""; for(i=1;i<=n;i++) s=s (i>1?",":"") (a[i] in m ? m[a[i]] : a[i]); print s}' \
       names_map.tsv ${R4}/tree/outgroup.txt)          # the outgroup's new name
python3 ${SCRIPTDIR}/4h_draw_tree.py --tree renamed/exendin_aa.rooted.nwk \
    --groups names_groups.tsv --sources names_rules.tsv --legend names_legend.txt \
    ${OG:+--outgroup $OG} --row 18 \
    --title "Exendin-family genes, SLC4A8–AQP5; colour = published peptide matched (≥80%)" \
    --png exendin_named.png > exendin_named.svg
```
- **Colour** = published peptide encoded at ≥ 80% (5a); grey = no close published match.
  It is a sequence statement, independent of the name.
- **Code after each tip** = the rule that named it: `[S]` only gene in its genome (GCGL),
  `[L]` positional letter, `[H nn%]` *Heloderma* published match at that identity, `[O]`
  outgroup, `[P]` published precursor. Hover for the old name and the evidence.
- Supports SH-aLRT/UFBoot (trust ≥ 80 / 95); bold = *Heloderma*.

## 5. Named synteny figure
```
python3 ${SCRIPTDIR}/5d_plot_named_synteny.py --species-table ../species.tsv \
    --loci ${R4}/reciprocal/reciprocal_loci.tsv --boundaries-used ${R4}/reciprocal/boundaries_used.tsv \
    --models ${R4}/models/models.tsv --names names_names.tsv --groups names_groups.tsv \
    ${GG:+--rename $GG=GalGal_GCGL} --top GalGal --species-tree ${SPTREE} \
    --svg named_synteny.svg --png named_synteny.png --layout named_synteny.tsv 2> named_synteny.log
cat named_synteny.log
```
- **Labels** = new names without the species prefix; **colours** as in the named tree.
- **Outlines from the 4e models**: filled = intact; hollow = pseudogene (stop codon);
  dashed = partial or flagged model.
- **`--top GalGal`** = chicken as the first row above the tree; `--rename` names chicken's
  locus, which `5b --drop` left out. A "drawn with their old names" warning lists loci that
  need a `--rename`.
- Letters run SLC4A8 → AQP5 along each track, so the figure checks 5b's ordering directly.

## Reading the result
Put `names_names.tsv` beside the tree.
1. **`_GCGL` tips** (single-copy genomes, *Heloderma* GCGL, chicken) should form a clade that
   follows the species tree if they are one ancestral gene. A single-copy gene outside it
   may be a retained duplicate: revisit rule 1 for that species rather than renaming by the
   tree.
2. **Letters** are positions; clades gathering the same letter across species would show
   conserved order, but copy numbers differ, so letters shift.
3. **Published names** tie venom peptides to genes. Genes that pair in the tree but carry
   different published names (HelSus_EX2 with HelCha_EX1, as in the old repo) are
   candidates for one gene named twice; `Published_EX*` tips show which gene each published
   precursor sits with.
4. **Root**: placed on chicken by construction; it orients the tree, it does not test where
   the duplications happened.
