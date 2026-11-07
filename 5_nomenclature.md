# Standard names for the exendin-family genes, and the named gene tree
**Why.** The names the curated genes carry after step 4 are inherited from whichever
annotation their best reverse hit came from (NCBI "exendin-3-like", TOGA2 "exendin-2/3",
the published peptides). In the step-4 tree they are neither monophyletic nor consistent
between sources. The published names exendin-1 to -4 are **mature-peptide** names, given to
peptides purified from *Heloderma* venom in 1984–1992, before any genome existed (see the
table at the top of the old repo's `2_exendin_origin.md`). This step replaces the inherited
names with names that say only what the data support.

**Working hypothesis behind the scheme.** There is one ancestral exendin-family gene, the
single-copy **GCGL** of Wang et al. 2012 (Irwin's "exendin"), which sits between SLC4A8 and
AQP5 and has duplicated several times in squamates.

**The rules** (`5b_name_genes.py`):

| Case | Name | What it claims |
|---|---|---|
| one exendin-family gene in the genome | `<species>_GCGL` | the ancestral single-copy state |
| several genes in the genome | `<species>_A`, `_B`, `_C` … in order from SLC4A8 to AQP5 | **position only**: gene A of one species need not be orthologous to gene A of another |
| *Heloderma* (HelSus, HelCha) | `<species>_EX1` … `_EX4`, or `_GCGL`, by the published mature peptide the gene encodes (≥ 80% identity) | which published peptide the gene encodes; a copy matching nothing published keeps its positional letter |
| outgroup tip (chicken) | `<species>_GCGL` if GCGL-type | not counted in any genome |

Details of the rules:
- **What counts as a copy.** Every curated (high-confidence) locus of a species counts,
  remnants and disabled copies included. So "single copy" means a genome with no second
  family gene at all.
- **Two *Heloderma* copies with the same peptide.** The better match (e.g. the exact one)
  keeps the plain name, and the other gets its positional letter appended (`HelSus_EX4`,
  `HelSus_EX4_D`).
- **How "SLC4A8 to AQP5" is read.** On the contig carrying both boundary genes, position is
  measured from SLC4A8 towards AQP5, whichever way round the assembly has them. In a
  fragmented assembly the order is: SLC4A8 contig, then unplaced contigs, then AQP5 contig.
  The `order_basis` column of `names_names.tsv` says which applied to each gene.

**Scripts**

| File | Origin |
|---|---|
| `5a_assign_nomenclature.py` | **copy** of the old repo's `2d_assign_nomenclature.py`, with one fix (below) |
| `5b_name_genes.py` | new: applies the rules |
| `5c_rename_tips.py` | new: relabels trees and FASTA files with the new names |
| `4h_draw_tree.py` | from step 4 (`--legend` added): draws the named tree |
| `4_reference_peptides.fa` | **copy** of the old `2_reference_peptides.fa`: the published mature peptides |

The fix in `5a`: the old `2d` wrote a literal `...` where the `differences` column belongs,
so every row had one column fewer than the header, with `differences` and `second_match`
run together. Old `nomenclature_assignments.tsv` files are shifted by one column from
`differences` on. Nothing else in the script changed.

Run the set-up block of `0_Environment_management.md` first.
```
mkdir -p ${WORK}/nomenclature
cd ${WORK}/nomenclature
R4=${WORK}/reciprocal_blast
```

## 1. Rebuild the step-4 tree with every gene in the region
The tree holds every curated exendin-family locus between SLC4A8 and AQP5. Remnants,
flagged models and pseudogenes are all included; `4g` no longer filters on model quality.
Stop codons are masked (`*` → `X` in the protein, stop codons → `NNN` in the CDS) so
MAFFT and IQ-TREE accept them. `tree/tips.tsv` says which tips were masked or flagged.
The only loci left out are those with no 4e model, and those 4d placed `elsewhere`,
i.e. outside the region.

The earlier step-4 tree had no squamate GCGL, because the old 4g dropped every locus 4d
called "disabled". 4d's flag includes `domain_starts_X` from the old interval scanner
(`lib_scan_interval.py`), which finds the hormone domain with the strict motif
`[HY][ASG][DEQ][AGP][ATVIGSM]F[TSI]`. Squamate GCGL breaks that motif (`YSEGTLAS`), so the
scanner assumed a default domain start and flagged every GCGL locus. With no filtering in
4g, that flag no longer matters. Re-run 4g and the tree, as in `4_reciprocal_blast.md`
step 7:
```
cd ${R4}
python3 ${SCRIPTDIR}/4g_tree_inputs.py --curated reciprocal/curated_loci.tsv --models-dir models \
    --outgroup-fasta tree_outgroup.faa --outgroup-group "GCGL (glucagon-like)" \
    --reverse-table queries/reverse_table.tsv --outdir tree 2> tree/tree_inputs.log
cat tree/tree_inputs.log                         # "left out" lines: only no-model and outside-region loci
grep -c GCGL tree/groups.tsv                     # squamate GCGL tips are in
awk -F'\t' '$6=="yes" && $7!="-"' tree/tips.tsv  # tips kept with masked stops or flagged models

OG=$(cat tree/outgroup.txt)
IQTREE=$(command -v iqtree3 || command -v iqtree2 || command -v iqtree)
cd tree
mafft --localpair --maxiterate 1000 tree_aa.fa > exendin_aa.aln.fa
$IQTREE -s exendin_aa.aln.fa -m MFP -B 1000 -alrt 1000 -T AUTO ${OG:+-o $OG} --prefix exendin_aa -redo
if [[ -n "$OG" ]]; then cp exendin_aa.treefile exendin_aa.rooted.nwk
else python3 ${SCRIPTDIR}/lib_midpoint_root.py exendin_aa.treefile > exendin_aa.rooted.nwk; fi
cd ${WORK}/nomenclature
```
- **Short remnants.** Single-exon remnants of 20–60 aa are kept. They are mostly gaps in
  the alignment and tend to sit on long, poorly supported branches, so read their placement
  with that in mind.
- **Leaving tips out.** Drop single tips with `--exclude NAME,...`. For a selection
  analysis, where pseudogenes and partial models do not belong, use `--exclude-disabled`
  and `--exclude-flagged`.

## 2. Compare every gene's mature peptide with the published peptides
Run on every curated model, not only the tree tips, so loci left out of the tree are named
too. Include the outgroup.
```
cat ${R4}/models/curated_aa.fa ${R4}/tree_outgroup.faa > genes_aa.fa
python3 ${SCRIPTDIR}/5a_assign_nomenclature.py --query genes_aa.fa \
    --reference ${SCRIPTDIR}/4_reference_peptides.fa --out-prefix nomenclature 2> nomenclature.log

awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i}
  {print $h["gene"]"\t"$h["mature_peptide"]"\t"$h["n_hormone_domains"]"\t"$h["how_found"]"\t"\
         $h["best_match"]"\t"$h["identity"]"\t"$h["verdict"]}' nomenclature_assignments.tsv | column -t -s$'\t' | less -S
```
Read it as in the old repo:
- **Verdicts.**
  - **"encodes published EXENDIN<i> exactly"**: this locus is the gene behind that peptide.
  - **"EXENDIN<i>-type, N difference(s)"**: same peptide class (≥ 80%), another copy,
    individual or species.
  - **"no published match"**: a copy never described as a peptide. Expect this for most
    squamate copies. HydCya GCGL is ~52% identical to the zebra finch GCGL peptide and
    ShiCro's copy ~52% to GLP-1, so the published panel cannot classify squamate genes
    outside *Heloderma*. That is why rule 2 uses position, not peptide.
- **`how_found`** says how the mature peptide was located. Check any row that is not
  `motif` by eye: `cleavage_site` and `best_window` mean the precursor's hormone domain
  does not match the strict motif.
- **`n_hormone_domains`**: 1 for an exendin or GCGL, 3 for real proglucagon. It counts
  strict-motif matches only, so diverged squamate peptides can show 0. The 155-aa GCGL
  models of ColCon, HydCya and ThaEle carry a second and third glucagon-like segment;
  look at their full precursors before treating that as anything.
- **Expected *Heloderma* matches.** From the old repo's curation: HelSus `EXENDIN2` and
  `EXENDIN4` exactly, HelCha `EXENDIN1` and `EXENDIN3` exactly. If one of those comes back
  lower, its 4e model is wrong (check `models/models.tsv`), not the name.
- ***Heloderma* GCGL may not reach 80%** to the zebra finch GCGL peptide. If so it keeps a
  positional letter. The old repo also noted that the *Heloderma* GCGL peptides lack the
  dibasic site after the hormone domain (`ERYS..` where other squamates have `KR`), so their
  mature peptide runs on and scores lower. Decide by eye whether to call it GCGL (step 3).

## 3. Apply the naming rules
```
python3 ${SCRIPTDIR}/5b_name_genes.py \
    --curated ${R4}/reciprocal/curated_loci.tsv \
    --boundaries ${R4}/reciprocal/boundaries_used.tsv \
    --assignments nomenclature_assignments.tsv \
    --drop GalGal_Ex3 \
    --out-prefix names 2> names.log
cat names.log
column -t -s$'\t' names_names.tsv | less -S
```
- **`--drop GalGal_Ex3`.** Chicken's curated locus is the same gene as its NCBI protein
  used as the outgroup. Dropping the partial locus keeps chicken's GCGL from being counted
  twice. Use the name your `curated_loci.tsv` has, or leave `--drop` out if chicken has no
  curated locus.
- **The *Heloderma* lines** in the log list each gene's new name. Check them against the
  expected matches in step 2.
- **`[warn] ... rest on an unplaced or out-of-interval locus`**: that species' letters
  depend on an assumed position. Look at its `order_basis` column before using the
  letters in text.
- **`[warn] new name X given to ...`**: two genes would get the same name. Fix it by hand
  in `names_map.tsv`; no automatic rule should decide this.
- **Overriding a name.** Edit `names_map.tsv` (column 2); everything downstream reads it.
  The *Heloderma* GCGL of step 2 is the expected case: if you judge it GCGL by eye, set its
  new name to `<species>_GCGL` there, and note why in `names_names.tsv`.

`names_names.tsv` is the record of the renaming: old name, new name, rule, position, and
the published-peptide evidence. Keep it with the tree.

## 4. Relabel the trees
Branch lengths and support values are untouched, so the step-4 trees are relabelled, not
rebuilt.
```
python3 ${SCRIPTDIR}/5c_rename_tips.py --map names_map.tsv --outdir renamed \
    ${R4}/tree/exendin_aa.rooted.nwk ${R4}/tree/exendin_aa.aln.fa ${R4}/tree/tree_aa.fa
# and the codon tree, if built:
# python3 ${SCRIPTDIR}/5c_rename_tips.py --map names_map.tsv --outdir renamed ${R4}/tree/exendin_codon.rooted.nwk
```
Every name not in the map is listed in the log line. There should be none.

## 5. Draw the named tree
```
OG=$(awk -F'\t' 'NR==FNR{m[$1]=$2; next} {n=split($0,a,","); s=""; for(i=1;i<=n;i++) s=s (i>1?",":"") (a[i] in m ? m[a[i]] : a[i]); print s}' \
       names_map.tsv ${R4}/tree/outgroup.txt)
echo "outgroup: ${OG:-none (midpoint-rooted)}"
python3 ${SCRIPTDIR}/4h_draw_tree.py --tree renamed/exendin_aa.rooted.nwk \
    --groups names_groups.tsv --sources names_rules.tsv --legend names_legend.txt \
    ${OG:+--outgroup $OG} --row 18 \
    --title "Exendin-family genes, SLC4A8–AQP5; colour = published peptide matched (≥80%)" \
    --png exendin_named.png > exendin_named.svg
```
The `OG=` line translates the outgroup's old name into its new one. Then:
- **Colour** says which published peptide a gene encodes, at ≥ 80% identity (5a). Grey is
  "no close published match", which is most copies outside *Heloderma*. Colour is a
  sequence statement, independent of the name.
- **The code after each tip** says which rule named it:
  - `[S]` the only gene in its genome, so GCGL;
  - `[L]` a positional letter;
  - `[H 100%]` *Heloderma*, named by the published peptide it matches, at that identity;
  - `[O]` outgroup.

  Hover over a code in the SVG for the old name and the evidence.
- **Supports and bold.** Node labels are SH-aLRT/UFBoot; trust a node at ≥ 80 / ≥ 95.
  Bold = *Heloderma*.

## Reading the result
Put `names_names.tsv` beside the tree.
1. **`_GCGL` tips** (single-copy genomes, plus *Heloderma* GCGL-type and the outgroup) should
   form a clade following the species tree if they are one ancestral gene. Where they do
   not, the single-copy gene of that genome may be a retained duplicate rather than GCGL.
   That is a reason to revisit rule 1 for that species, not to rename by the tree.
2. **Letters** are positions. To ask whether position is conserved, look for clades that
   gather the same letter across species. Don't expect it: different species have different
   numbers of copies, so letters shift.
3. ***Heloderma* published names** tie the venom peptides to genes. Two genes that pair in
   the tree and carry different published names (e.g. HelSus_EX2 and HelCha_EX1) are
   candidates for one gene named twice. The old repo found exactly this pair (HelCha
   exendin-1 pairing with HelSus exendin-2).
4. **The root.** It is placed on the chicken branch by construction: it orients the tree,
   it does not test where the duplications happened. See the old `2_exendin_origin.md`
   step 5 for the outgroup set that would test it.


cd ${WORK}/nomenclature
R4=${WORK}/reciprocal_blast
python3 ${SCRIPTDIR}/5d_plot_named_synteny.py --species-table ../species.tsv \
    --loci ${R4}/reciprocal/reciprocal_loci.tsv --boundaries-used ${R4}/reciprocal/boundaries_used.tsv \
    --models ${R4}/models/models.tsv --names names_names.tsv --groups names_groups.tsv \
    --rename GalGal_Ex3=GalGal_GCGL --top GalGal \
    --species-tree ${SPTREE} \
    --svg named_synteny.svg --png named_synteny.png --layout named_synteny.tsv 2> named_synteny.log
cat named_synteny.log


cd ${WORK}/reciprocal_blast
grep -c '^GalGal' ../species.tsv          # must be 1: delete a duplicate row first (the doubled chicken track)
rm -rf loci reciprocal models tree
python3 ${SCRIPTDIR}/4c_call_loci.py --species-table ../species.tsv --forward-dir forward \
    --queries queries/forward_queries.faa --query-table queries/query_table.tsv \
    --boundary-queries queries/boundary_queries.faa --outdir loci 2> loci.log
grep -E "joined|not joined" loci.log
python3 ${SCRIPTDIR}/4d_reciprocal_classify.py --species-table ../species.tsv --genomedir ${GENOMEDIR} \
    --loci-dir loci --queries-dir queries --threads 8 --outdir reciprocal 2> reciprocal.log
python3 ${SCRIPTDIR}/4e_build_models.py --species-table ../species.tsv --genomedir ${GENOMEDIR} \
    --curated reciprocal/curated_loci.tsv --loci-dir loci \
    --queries queries/forward_queries.faa --query-table queries/query_table.tsv \
    --outdir models 2> models.log
seqkit grep -r -p '^Hel' models/curated_aa.fa      # HelSus/HelCha GCGL should now start with M...


mkdir -p tree
python3 ${SCRIPTDIR}/4g_tree_inputs.py --curated reciprocal/curated_loci.tsv --models-dir models \
    --outgroup-fasta tree_outgroup.faa --outgroup-group "GCGL (glucagon-like)" \
    --reverse-table queries/reverse_table.tsv --exclude GalGal_Ex3 --outdir tree 2> tree/tree_inputs.log