# Exendin gene tree from TOGA2: collection, curated alignment
**Goal:** one alignment of every exendin gene TOGA2 annotated, across all four TOGA2 runs.
Every curation decision is recorded in a table, so the alignment can be rebuilt from the
TOGA2 files at any time. The tree is built from the curated alignment (step 4), rooted on GIP, with a topology test of where the exendins come from.

**What goes in:**
- **Gene set:** every query gene whose *reference* gene is named exendin in that run's
  reference annotation (exendin-N, exendin-N-like, unnumbered exendin).
- **Runs:** HelSus, HemCap, AnoCar and CroPor. Only references that annotate an exendin
  contribute; the log says which genes each run used.
- **Ortholog projections:** kept only if TOGA2 calls them intact (status FI or I).
- **The HelSus exendin-4 isoform `mRNA-2` is never used** (`--exclude-projection`). Its 419-bp
  first exon is an annotation artifact (no protein support). It was the only model TOGA2
  produced for HelCha exendin-3, which is replaced by the rebuilt model in step 3. In every
  other genome, `mRNA-2` scored L/M/UL and failed TOGA2's orthology step, so excluding it
  changes nothing there.
- **Unclassified projections are left out:** these are in TOGA2's annotation, but their
  transcript did not pass TOGA2's orthology step (typically a second isoform of an
  ortholog). Where they overlap an ortholog, the ortholog's projection represents the
  locus anyway. `--keep-unclassified` keeps them.
- **Paralogous projections:** TOGA2's own `#paralog` projections, i.e. extra copies found
  through a paralogous chain. Held to the same rule: kept only if intact (FI or I), so
  lost (L), uncertain (UL), missing (M) and partial (PI) paralogs are dropped.
  `--paralog-status` sets a different list (e.g. `FI,I,PI`), or `any`.
- **Left out:** retrogene projections, and CDS that are missing, shorter than 150 bp, or
  not a multiple of 3.

**One locus, one sequence.**
- **Within a run:** projections that overlap on one strand (isoforms, or one locus reached
  from two reference exendins) are one locus. The best status represents it, then the
  lowest isoform number, then the longest CDS.
- **Across runs:** the same locus found by several runs is kept once. The order is best
  status, then the order of `--run` (HelSus first), then the longest CDS. The other runs
  are recorded in `also_in`, which is independent support for the locus.

`--include-gcgl` also collects the genes named glucagon-like (GCGL, the ancestral gene).
Use it when you want an outgroup to root the tree.

Run the set-up block of `0_Environment_management.md` first. `10a` reads its helpers from
`3b_collect_toga2_exendins.py` and `3c_plot_toga2_synteny.py`, so keep all three in `${SCRIPTDIR}`.
```
mkdir -p ${WORK}/toga2_exendin_tree && cd ${WORK}/toga2_exendin_tree
```

## 1. Collect
```
python3 ${SCRIPTDIR}/10a_collect_toga2_exendins.py \
    --run HelSus="${TOGADIR}/target_HelSus/query_{q}" \
    --run HemCap="${TOGADIR}/target_HemCap/query_{q}" \
    --run AnoCar="${TOGADIR}/query_{q}/target_AnoCar/output" \
    --run CroPor="${TOGADIR}/query_{q}/target_CroPor/output" \
    --ref-gff HelSus=${HELSUS_EVIANN_GFF} \
    --ref-gff HemCap=${GENOMEDIR}/HemCap_GCF_027244095.1/ncbi_dataset/data/GCF_027244095.1/genomic.gff \
    --ref-gff CroPor=${GENOMEDIR}/CroPor_GCF_001723895.1/ncbi_dataset/data/GCF_001723895.1/genomic.gff \
    --ref-gff AnoCar=$(ls ${GENOMEDIR}/AnoCar_*/ncbi_dataset/data/*/genomic.gff | head -1) \
    --exclude-projection 'LOC_00012325-mRNA-2#' \
    --outdir collect 2> collect.log
cat collect.log
cut -f1-8,15-19 collect/candidates.tsv | column -t
```
Read the log:
- **`run X: N reference gene(s) named exendin: ...`** = the reference exendins each run
  projected. For HelSus, expect `LOC_00012324` (exendin-2) and `LOC_00012325` (exendin-4).
  `NONE` with a warning = that reference annotates no exendin, or its names were not
  resolved (wrong or missing GFF; check `3c`'s `have no name` warnings).
- **`kept: orthologs with status FI,I; paralogs with status FI,I`** = the filters in force.
- **`left out: ...`** = how many loci each filter removed (`paralog_status_L` and so on). `candidates.tsv` keeps every one,
  with its reason, so nothing disappears silently.
- **One line per genome** = the kept tips with status and length. The precursors are 84–113 aa.
  A tip well outside that range is a model problem. The known case is HelCha exendin-3:
  TOGA2 used reference isoform `mRNA-2`, so expect ~198 aa.

## 1b. Rebuild HelCha exendin-3 (no TOGA2 model)
With `mRNA-2` excluded, TOGA2 has no model for HelCha's exendin-3 locus: its `mRNA-1`
never reached TOGA2's orthology step (no usable chain over the canonical exon 1). The gene is
rebuilt from the published *H. horridum* exendin-3 cDNA (AJ580309) aligned to the HelCha
contig, as in Exendin_evolution, and added in step 3 with `--add-cds`.
```
mkdir -p rescued
cat << 'EOF' > rescued/AJ580309.fa
>AJ580309.1 Heloderma horridum mRNA for exendin 3
GGTGCCTGTGCAAAGGAGAAGATGAAAATCATCCTGTGGCTGTGTGTTTTCGGGCTGTTCCTTGCAACTT
TATTCCCTGTCAGCTGGCAAATGCCTGTTGAATCTGGGTTGTCTTCTGAGGATTCTGCAAGCTCAGAAAG
CTTTGCTTCGAAGATTAAGCGACATAGTGATGGAACATTTACCAGTGACTTGTCAAAACAGATGGAAGAG
GAGGCAGTGCGGTTATTTATTGAGTGGCTTAAGAACGGAGGACCAAGTAGCGGGGCACCTCCGCCATCGG
GTTAAGGTCTTTCAATTGTGGAACAAGACACACACCTGATGTTTGACGACCATTTTAAAGAAATGTTTCC
AGCAATAAGTCACATGTCTTTGTTTCCAATTAGTGAGCGACACAGCC
EOF
HELCHA_FASTA=$(ls ${GENOMEDIR}/HelCha_GCA_026122225.1/ncbi_dataset/data/GCA_026122225.1/*_genomic.fna)
grep -c JANEZZ010000386.1 ${HELCHA_FASTA}.fai          # must be 1
python3 ${SCRIPTDIR}/lib_rescue_locus.py --genome ${HELCHA_FASTA} --contig JANEZZ010000386.1 \
    --query rescued/AJ580309.fa --name HelCha_Ex3 --outdir rescued
cat rescued/HelCha_Ex3.exons.tsv rescued/HelCha_Ex3.aa.fa
```
Expect 3 exons on the minus strand covering the cDNA end to end, exon 1 near
4,399,218–4,399,326, and an 87-aa protein identical to the published exendin-3 precursor
(P20394): `MKIILWLCVFGLFLATLFPVSWQMPVESGLSSEDSASSESFASKIKRHSDGTFTSDLSKQMEEEAVRLFIEWLKNGGPSSGAPPPSG`.

## 1c. Outgroups: GCG, GIP and chicken GCGL
The family's other ligands root the tree and show where the exendin copies attach:
proglucagon (GCG) and GIP from chicken, human and HelSus, and chicken GCGL, the ortholog
used to define GCGL (Wang et al. 2012).
- **Chicken and HelSus GCG/GIP** come from their NCBI annotations, already downloaded in step 1
  (GalGal `GCF_016699485.2`, HelSus `GCF_057233095.1`). The gene is found by symbol, then by
  description. Its MANE-tagged transcript is used if there is one, else the longest CDS.
- **Human GCG/GIP** come from their RefSeq mRNAs (`NM_002054`, `NM_004123`), and **chicken
  GCGL** from EU718628, fetched as CDS records.
```
mkdir -p outgroups
command -v efetch || module load edirect
efetch -db nuccore -id NM_002054,NM_004123 -format fasta_cds_na > outgroups/human_cds.fa
efetch -db nuccore -id EU718628 -format fasta_cds_na > outgroups/EU718628_cds.fa
grep '>' outgroups/*.fa

python3 ${SCRIPTDIR}/10d_outgroup_cds.py --species-table ../species.tsv --genomedir ${GENOMEDIR} \
    --gene GalGal:GCG,GIP --gene HelSus:GCG,GIP \
    --cds-fasta HomSap=outgroups/human_cds.fa --cds-fasta GalGal=outgroups/EU718628_cds.fa:GCGL \
    --out outgroups/outgroups_cds.fa 2> outgroups/outgroups.log
cat outgroups/outgroups.log
cut -f1-7 outgroups/outgroups_cds.fa.tsv | column -t
```
Expect 7 tips: `GalGal_GCG`, `GalGal_GIP`, `HelSus_GCG`, `HelSus_GIP`, `HomSap_GCG`, `HomSap_GIP`,
`GalGal_GCGL`, each `ok` in the checks column.
- **`no gene X`**: that annotation names the gene differently. Find it with
  `grep -i "glucagon\|inhibitory polypeptide" <genomic.gff> | awk -F'\t' '$3=="gene"'`, then pass
  its symbol (`--gene HelSus:LOC...`), or a description regex (`--describe GIP='...'`).
- **`N genes match`**: check that the one used is right (e.g. `GCG` and not a `GCG`-like LOC).
- **`no ATG start` / `no final stop`**: a partial record. EU718628 should be the complete
  chicken GCGL CDS: check it says `ok`.
- **Proteins:** `HomSap_GCG` should be the 180-aa human proglucagon and `HelSus_GCG` should
  match the HelSus proglucagon `GCG_precursor` in `4_reference_peptides.fa` (O12956).

**Trim the outgroups to the part that matches an exendin precursor.** Proglucagon carries
three hormone domains (glucagon, GLP-1, GLP-2), and GIP has a long C-terminal peptide.
Aligned whole, they open long gap blocks that hold no information about the exendins.
Keep signal peptide, propeptide, the first hormone domain and its cleavage site.
- **Human proglucagon:** signal 1–20, GRPP 21–50, glucagon 53–81, IP-1 to ~89. Keep `1-89`.
- **Human GIP precursor:** signal 1–21, propeptide 22–51, GIP 52–93, then the C-terminal
  peptide. Keep `1-95`.
- **Chicken and HelSus:** check where the same domains end in the viewer of step 2.
Add the `trim` rows in step 3, e.g. `HomSap_GCG	trim	1-89	signal, GRPP, glucagon, IP-1; GLP-1 and GLP-2 removed`.

## 2. First alignment
Codon-aware alignment with MACSE, as in step 4. Frameshifts and internal stops become `NNN`/`X`
instead of breaking the frame.
```
mkdir -p aln0
cat collect/exendin_cds.fa rescued/HelCha_Ex3.cds.fa outgroups/outgroups_cds.fa > aln0/input_cds.fa
macse -prog alignSequences -seq aln0/input_cds.fa -out_NT aln0/exendin.macse_NT.fa -out_AA aln0/exendin.macse_AA.fa
macse -prog exportAlignment -align aln0/exendin.macse_NT.fa -codonForFinalStop --- -codonForInternalStop NNN -codonForExternalFS --- -codonForInternalFS NNN -out_NT aln0/exendin.codon.fa -out_AA aln0/exendin.aa.fa
```
For viewing only, add the published precursors (exendin-2, -3, -4) as guide rows. `--keeplength`
leaves the tips' alignment unchanged:
```
seqkit grep -r -p "^EXENDIN[0-9]_precursor" ${SCRIPTDIR}/4_reference_peptides.fa > published_precursors.faa
mafft --add published_precursors.faa --keeplength aln0/exendin.aa.fa > aln0/exendin.aa.with_published.fa
python3 ${SCRIPTDIR}/10c_plot_alignment.py --aln aln0/exendin.aa.with_published.fa \
    --tips collect/candidates.tsv --species-tree ${SPTREE} \
    --svg aln0/exendin.aa.svg --pdf aln0/exendin.aa.pdf --png aln0/exendin.aa.png
```
Open `aln0/exendin.aa.pdf`.
- **Rows:** grouped by reference exendin (Ex2, Ex4, …) in species-tree order, with the
  published precursors in grey at the top.
- **Annotation columns:** TOGA2 status (red if not FI/I), o/p (ortholog / paralogous
  projection) and length (red outside 80–120 aa).
- **Residues:** coloured by class, with `X` (a stop or unknown codon) in red.
- **Shading:** to shade the mature peptide, find its columns from the published rows, then
  re-run with e.g. `--region 48-86`.

## 3. Curate
Copy the template that 10a wrote, then edit the copy. Re-running 10a overwrites
`collect/curation.tsv`, so never edit the template itself:
```
cp collect/curation.tsv curation.tsv
```
It has one row per tip, all `keep`, with automatic notes. Change `action` and give the
reason in `note`:

| action | value | use for |
|---|---|---|
| `keep` | | as is |
| `drop` | | a broken model, a pseudogene, a duplicate of another tip, a copy you cannot place |
| `trim` | `first-last` amino acids to keep (1-based, unaligned) | a model with a spurious extension, e.g. `42-128` |
| `rename` | new tip name | a tip whose name should reflect what it is, e.g. `HelCha_Ex4` → `HelCha_Ex3` |

Several rows for one tip are applied in order (e.g. a `trim` row, then a `rename` row).

Things to look for:
- **N- or C-terminal extensions** that the published precursors and the other tips lack:
  a mis-annotated reference isoform, as for HelCha exendin-3. Trim them.
- **Rows that only align over the peptide or only over the signal peptide:** a partial
  projection.
- **Paralogous projections (p)** that are identical to an ortholog tip from the same genome:
  the same locus reached twice. Drop one. Check `contig`/`start` in `candidates.tsv`.
- **`X` inside the hormone domain** or a domain that does not start His/Tyr: a pseudogene.
  Keep it if you want it in the tree as such, or drop it.

Then apply the table and realign:
```
mkdir -p curated aln1
python3 ${SCRIPTDIR}/10b_apply_curation.py --cds collect/exendin_cds.fa --curation curation.tsv \
    --add-cds rescued/HelCha_Ex3.cds.fa --add-cds outgroups/outgroups_cds.fa \
    --out-prefix curated/exendin 2> curated.log
cat curated.log                        # must list your drops / renames / trims, and HelCha_Ex3 + 7 outgroups 'added'
grep -c '>' curated/exendin_cds.fa     # stop here if this fails: nothing below can run
macse -prog alignSequences -seq curated/exendin_cds.fa -out_NT aln1/exendin.macse_NT.fa -out_AA aln1/exendin.macse_AA.fa
macse -prog exportAlignment -align aln1/exendin.macse_NT.fa -codonForFinalStop --- -codonForInternalStop NNN -codonForExternalFS --- -codonForInternalFS NNN -out_NT aln1/exendin.codon.fa -out_AA aln1/exendin.aa.fa
mafft --add published_precursors.faa --keeplength aln1/exendin.aa.fa > aln1/exendin.aa.with_published.fa
python3 ${SCRIPTDIR}/10c_plot_alignment.py --aln aln1/exendin.aa.with_published.fa \
    --tips collect/candidates.tsv --curation-log curated/exendin_log.tsv --species-tree ${SPTREE} \
    --svg aln1/exendin.aa.svg --pdf aln1/exendin.aa.pdf --png aln1/exendin.aa.png
```
Repeat (edit `curation.tsv`, run this block, which writes `aln1/` again) until the alignment
is clean. `curated/exendin_log.tsv` is the record of every decision and its reason.

Sequences from outside TOGA2 (e.g. a model rebuilt with `lib_rescue_locus.py`) can be added
with `--add-cds rebuilt.fa`. They are logged as `added`.

## 4. Tree
Built on the curated **protein** alignment: GIP, GCG and GCGL are too divergent from the
exendins at the codon level for a codon model to be reliable on the deep nodes. Set `ALN` to
the alignment you settled on in step 3 (MACSE, the L-INS-i version, or the C-terminally
trimmed one).
```
cd ${WORK}/toga2_exendin_tree
ALN=aln1/exendin.aa.fa
IQTREE=$(command -v iqtree3 || command -v iqtree2 || command -v iqtree)
mkdir -p tree
python3 ${SCRIPTDIR}/10e_tree_inputs.py --aln ${ALN} --outdir tree 2> tree/tree_inputs.log
cat tree/tree_inputs.log
OG=$(cat tree/outgroup.txt); echo ${OG}

$IQTREE -s ${ALN} -m MFP -B 1000 -alrt 1000 -o ${OG} -T AUTO --prefix tree/exendin_aa -redo
grep "Best-fit model according to" tree/exendin_aa.iqtree
python3 ${SCRIPTDIR}/4h_draw_tree.py --tree tree/exendin_aa.treefile --groups tree/groups.tsv --outgroup ${OG} \
    --title "Exendin-family genes from TOGA2, curated; rooted on GIP" --png tree/exendin_aa.png > tree/exendin_aa.svg
```
Check `tree_inputs.log` first. It counts the tips in each class: 3 GIP, 3 GCG, 1 chicken
GCGL, the squamate GCGL, and the exendins. A tip listed under the wrong class has a name
the script does not recognise. Rename it in `curation.tsv` (`<Code>_GIP`, `<Code>_GCG`,
`<Code>_GCGL`, anything else = exendin).

**Rooting on GIP.** GIP is the family's single-domain ligand outside the glucagon /
GCGL / exendin group, so the root falls on its branch. If the three GIP tips are not a
clade, IQ-TREE says so in `tree/exendin_aa.log` (`grep -i outgroup tree/exendin_aa.log`).
Then the root is unreliable, and the GIP alignment needs checking first.

### Where do the exendins come from? Topology test
The ML tree shows one answer. Its support at the deep nodes is often weak with ~90 columns,
so test the three origins directly. For each hypothesis, find the best tree under that
constraint, then compare all of them with the unconstrained tree in an AU test:

| Hypothesis | Constraint | Meaning |
|---|---|---|
| H1 | exendins + squamate GCGL form a clade, excluding chicken GCGL | the exendins arose by duplicating GCGL within squamates |
| H2 | exendins sister to all GCGL (chicken included), together excluding GCG | the duplication predates the bird–squamate split |
| H3 | exendins + GCG form a clade, excluding GCGL | the exendins derive from proglucagon, not GCGL |
```
MODEL=$(grep "Best-fit model according to" tree/exendin_aa.iqtree | awk '{print $NF}'); echo ${MODEL}
for H in H1_squamate_GCGL H2_sister_GCGL H3_GCG; do $IQTREE -s ${ALN} -m ${MODEL} -g tree/${H}.constraint.nwk -o ${OG} -T AUTO --prefix tree/${H} -redo; done
cat tree/exendin_aa.treefile tree/H1_squamate_GCGL.treefile tree/H2_sister_GCGL.treefile tree/H3_GCG.treefile > tree/candidates.nwk
$IQTREE -s ${ALN} -m ${MODEL} -z tree/candidates.nwk -n 0 -zb 10000 -au -T AUTO --prefix tree/AUtest -redo
grep -A8 "USER TREES" tree/AUtest.iqtree
```
The trees are listed in that order: 1 = ML, 2 = H1, 3 = H2, 4 = H3. A hypothesis with
`p-AU < 0.05` is rejected. More than one hypothesis not rejected means this alignment cannot
choose between them, and that should be reported as such. The answer is not settled by the
ML tree alone.

### Reading the tree
- **Support:** node labels are SH-aLRT/UFBoot. Trust a node at ≥ 80 / ≥ 95. Bold tips are
  *Heloderma*.
- **Expected backbone:** GIP outside; then proglucagon (chicken, human, HelSus) and GCGL
  (chicken + squamates) as separate clades; then the exendins. Where the exendins attach is
  the question tested above.
- **Within the exendins:**
  - Orthologs across species (e.g. `HelSus_Ex2_Ex4c` = exendin-4 with `HelCha_Ex3` =
    exendin-3) should pair.
  - Copies within one genome that group together are recent tandem duplicates.
  - The anguid 76–78 aa copies form their own lineage.
- **Colours** follow the tip names, which are still TOGA2's provisional labels. Name the
  genes from this tree with step 5. Until then, read groups from the topology, not the
  colours.

### Codon tree of the exendins (optional)
For relationships among recent copies (within *Heloderma*, within anguids), the codon
alignment has more signal. Drop GCG and GIP, then midpoint-root. A GCGL root would assume
the answer to H1/H2.
```
seqkit grep -v -r -p '_GCG$' -p '_GIP$' aln1/exendin.codon.fa > tree/exendin_codon.ingroup.fa
$IQTREE -s tree/exendin_codon.ingroup.fa -st CODON -m MFP -B 1000 -alrt 1000 -T AUTO --prefix tree/exendin_codon -redo
python3 ${SCRIPTDIR}/lib_midpoint_root.py tree/exendin_codon.treefile > tree/exendin_codon.rooted.nwk
python3 ${SCRIPTDIR}/4h_draw_tree.py --tree tree/exendin_codon.rooted.nwk --groups tree/groups.tsv \
    --title "Exendins and GCGL, codon model, midpoint-rooted" --png tree/exendin_codon.png > tree/exendin_codon.svg
```
If you trimmed the protein alignment in step 3, cut the codon alignment to match
(`seqkit subseq -r 1:<3 x protein columns>`).

## 5. Name the genes and draw the figure
Names come from the tree and from each gene's position in its genome:

| Name | Rule |
|---|---|
| `<Sp>_GCGL` | every tip in the GCGL clade |
| `<Sp>_Exendin_anc` | the ancestral anguimorph exendin: one gene per genome in the clade that holds the Shinisaurus, Heloderma and anguid orthologues |
| `<Sp>_A`, `_B`, `_C` … | other exendin genes, in order from that genome's Exendin_anc **toward AQP5** |
| `<Sp>_Z`, `_Y`, `_X` … | other exendin genes, in order from Exendin_anc **toward SLC4A8** |
| `<Sp>_U1`, `_U2` … | exendin genes on a different contig from Exendin_anc, which cannot be ordered |
| `*_GCG`, `*_GIP` | outgroups, unchanged |

Only exendin genes are counted for the letters; GCGL and other genes are skipped.

**Choosing the two clades.** Pass tips whose smallest common clade *is* that gene:
- **Exendin_anc:** one Heloderma, one Shinisaurus and one anguid gene of the
  one-per-genome clade. In the current tree that is `HelSus_Ex2_Ex4`, `ShiCro_Ex2_Ex4` and
  `ElgMul_Ex2_Ex4`. The log lists every tip in the clade: check it contains no GCGL and
  nothing you would call a copy.
- **GCGL:** chicken GCGL and one squamate GCGL, e.g. `GalGal_GCGL,HelSus_GCGL`. Use the old
  name, e.g. `HelSus_Ex3_Ex4`, if the GCGL rows were not renamed in step 3.

```
cd ${WORK}/toga2_exendin_tree
GENES=../synteny_toga2/interval/interval_HelSus_genes.tsv      # from step 3 (3c --genes-out)
ls ${GENES}
python3 ${SCRIPTDIR}/10f_name_genes_figure.py --tree tree/exendin_aa.treefile --outgroup "$(cat tree/outgroup.txt)" \
    --anc-tips HelSus_Ex2_Ex4,ShiCro_Ex2_Ex4,ElgMul_Ex2_Ex4 --gcgl-tips GalGal_GCGL,HelSus_Ex3_Ex4 \
    --candidates collect/candidates.tsv --curation-log curated/exendin_log.tsv \
    --rescued HelCha_Ex3=rescued/HelCha_Ex3.exons.tsv --genes ${GENES} \
    --species-table ../species_toga2.tsv --species-tree ${SPTREE} \
    --aln ${ALN} --aln-region "48-86:hormone domain" --drop-species BarImb \
    --outdir names --svg exendin_figure.svg --pdf exendin_figure.pdf --png exendin_figure.png 2> names.log
cat names.log
column -t -s$'\t' names/names.tsv
```
Read the log before the figure:
- **`N genes in the Exendin_anc clade`** (e.g. ElgCoe, ShiCro): a genome has more than one
  gene there. The one closest to the clade's base is called Exendin_anc, and the others
  become lettered copies. If synteny or the tree says another is the orthologue, override it.
- **`no Exendin_anc; copies counted from …GCGL`**: e.g. HelCha, where TOGA2 projected no
  Exendin_anc orthologue. Its copies are lettered from GCGL instead, so its letters are
  not directly comparable with HelSus's. Note this in the legend, or override the names.
- **`no boundary on the anchor's contig`**: the AQP5 direction was assumed, so check those
  letters by hand.
- **`no coordinates`**: a tip that isn't in `candidates.tsv`. Give it with `--coords`
  (`tip contig start end strand`) or `--rescued`.

**Overrides.** Write a table of `old_tip<TAB>new_name` and pass it with `--names overrides.tsv`.
Each override is marked in `names.tsv`. `names.tsv` records every name with the rule and
distance behind it, and `names/renamed.treefile` is the tree under the new names.

**The figure** is 180 mm wide (`--width-mm`; panel A is `--left-mm` 105), with 7 pt text, and
written as a vector PDF.
- **Panel A:** each genome's SLC4A8–AQP5 interval from the TOGA2 run, with every exendin-family
  gene labelled by its new name.
  - Colours: GCGL green, Exendin_anc blue, lettered copies orange, other genes grey.
  - Hollow genes are exendin loci in the TOGA2 annotation that are not in the tree (not
    intact, or curated out).
  - `--no-context` hides the grey genes. `--label-tiers 2` makes it shorter when labels
    are sparse.
- **Panel B:** the GIP-rooted tree with the same names and colours.
  - Filled dot = UFBoot ≥ 95 and SH-aLRT ≥ 80. Open dot = UFBoot ≥ 70.
- **Panel C** (`--aln`): the protein alignment the tree was built from, filtered to the
  exendin genes only (Exendin_anc and the lettered copies). GCGL, GCG and GIP rows are
  left out, and so are the columns that are gaps in every exendin.
  - Rows follow panel B's order, with the same names and colours. Residues are coloured
    by class, with `X` in red.
  - It runs the full width under A and B, always on one line: if the alignment does not fit
    at 7 pt, its text shrinks until it does. The log gives the size used (`text 6.4 pt`) and
    warns below 5 pt. `--aln-font-pt` sets the size directly.
  - `--aln-region START-END:LABEL` draws a labelled box over those columns. The numbering
    is **panel C's own**, after the gap columns are removed, so read the column numbers
    off a first draw and then set it (e.g. the mature peptide, `--aln-region "48-86:hormone
    domain"`). It can be given more than once, e.g. for the signal peptide.

**Leaving a genome out of the figure.** `--drop-species BarImb` removes *Barisia imbricata*
from all three panels: its row in A, its tips in B and its rows in C. Its tips are pruned
from the tree drawing (the tree is not re-inferred), and the other genes keep their names
and supports. To remove it from the analysis itself, drop its tips in `curation.tsv` and
re-run steps 3–5.

To carry the names into the alignment as well:
```
cut -f1,2 names/names.tsv | tail -n +2 > names/map.tsv
seqkit replace -p '^(\S+)' -r '{kv}' -k names/map.tsv aln1/exendin.aa.fa > aln1/exendin.aa.named.fa
```

## 6. Genome annotations for each locus (to pair with expression data)
For every exendin-family locus of one genome (GCGL, Exendin_anc, lettered copies): which
gene and transcript IDs of that genome's own annotations cover it. These are the IDs
expression data are quantified against. The name from step 5 is only the join key.
```
HELSUS_NCBI_GFF=${GENOMEDIR}/HelSus_GCF_057233095.1/ncbi_dataset/data/GCF_057233095.1/genomic.gff
python3 ${SCRIPTDIR}/10g_gene_table.py --species HelSus --names names/names.tsv \
    --gff EviAnn=${HELSUS_EVIANN_GFF} --gff NCBI=${HELSUS_NCBI_GFF} \
    --out names/HelSus_annotations.tsv
```
The screen gives one line per locus: its position, then each annotation's gene(s) there,
with NCBI GeneID and number of transcripts. `names/HelSus_annotations.tsv` has one row per
locus × overlapping gene, with these columns:
- **Gene:** ID, NCBI GeneID, biotype, strand and span, and the description.
- **Overlap:** in bp and as % of the locus.
- **Transcripts:** every transcript ID (EviAnn `LOC_…-mRNA-n`, NCBI `XM_…`/`XR_…`) and its
  product.

Read it like this:
- **Which IDs to use:** use the column of whichever annotation the expression data were
  quantified against. For a gene-level count matrix, that's `gene_id` / `ncbi_geneid`. For a
  transcript-level one (salmon, kallisto), use `transcripts`.
- **`none`:** that annotation has no gene at the locus. It is unannotated there, and an
  expression pipeline based on that annotation has no feature for it. Earlier work found
  most HelSus copies unannotated. To quantify them, add the locus to the annotation (e.g. a
  GTF built from the TOGA2 projection, `query_annotation.gtf.gz` of `target_HelSus/query_HelSus`).
- **`(spans locus: host gene?)`:** a gene more than 5× longer than the locus overlaps it,
  typically on the other strand (check `gene_strand`). The locus sits inside that gene. Its
  expression is not the exendin's.
- **Several genes at one locus:** check `overlap_pct_of_locus` and the strand. The one that
  covers the locus on the same strand is the exendin's model.
- **Coordinates:** loci use TOGA2's sequence names (`chr2`). NCBI's `NC_139459.1` is
  matched through the NCBI sequence report next to its GFF (chromosome `2`). If that file is
  missing, give the name with `--chrom-map chr2=NC_139459.1`.

Use `--species HelCha` (or any other code with an annotation) for another genome.

## Outputs
| File | Content |
|---|---|
| `collect/candidates.tsv` | every exendin locus found in any run: tip, run, other runs that found it, ortholog/paralog, status, coordinates, length, internal stops, kept or why not |
| `collect/exendin_cds.fa`, `collect/exendin_aa.fa` | the uncurated sequences |
| `collect/curation.tsv` | the curation template |
| `curation.tsv` | your curation decisions |
| `curated/exendin_{cds,aa}.fa`, `curated/exendin_log.tsv` | curated sequences and the log of what was done to each tip |
| `aln1/exendin.codon.fa`, `aln1/exendin.aa.fa` | the curated codon and protein alignments, for the tree |
| `aln*/exendin.aa.{svg,pdf,png}` | the alignment view |
| `tree/groups.tsv`, `tree/outgroup.txt`, `tree/*.constraint.nwk` | tip classes and colours, the GIP root, constraint trees for H1–H3 |
| `tree/exendin_aa.{treefile,iqtree,svg,png}` | protein ML tree, rooted on GIP, and its drawing |
| `tree/AUtest.iqtree` | AU test of the ML tree against H1–H3 |
| `tree/exendin_codon.*` | optional codon tree of the exendins and GCGL |
| `names/names.tsv`, `names/renamed.treefile` | every new name with its rule; the tree under the new names |
| `exendin_figure.{svg,pdf,png}` | the figure: named interval synteny (A), gene tree (B), exendin alignment (C) |
| `names/<Sp>_annotations.tsv` | per locus: the genome annotations' gene and transcript IDs overlapping it (EviAnn, NCBI), for expression data |
