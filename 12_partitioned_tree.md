# Partitioned tree: topology from the slow block, rates from the fast block
**Problem.** The step-4 tree puts the exendins sister to *all* GCGL, chicken included, with
weak support — which would mean the duplication predates the bird–squamate split and every
non-anguimorph lineage lost the copy independently. The synteny says otherwise: GCGL is
universal across squamates, exendin_anc is anguimorph-restricted. One of the two is wrong.

**Likely cause.** Irwin & Prentice (2011) found that on the *Heloderma* lineage the exendin
toxin peptide sequences, but not their pro or signal peptides, evolved very rapidly. A
fast-evolving block attached to a slow backbone is the standard way to drag a clade toward
the root. Panel C of the step-5 figure shows it directly: near-identical signal peptides,
wildly divergent hormone domains.

**Fix.** Split the alignment at the mature-peptide boundary. Build the topology from the
**slow** block (signal + pro), where the phylogenetic signal is. Then map the **fast** block
(mature peptide) onto that fixed topology and read the rates off it. The same split is the
right design for the selection tests, and it turns Irwin & Prentice's rate claim into a
hypothesis you can test with 14 genomes instead of the handful they had.

Everything below is in the **step-5 nomenclature** (`<Sp>_GCGL`, `<Sp>_Exendin_anc`,
`<Sp>_A` / `_B` / `_Z`, and the unchanged `*_GCG` / `*_GIP` outgroups), so run
`10h_apply_names.py` (end of step 5) first — the `aln1/` files still carry the provisional
TOGA2 labels and `--highlight HelSus_A` would not match anything in them.
```
cd ${WORK}/toga2_exendin_tree
ls named/exendin.aa.fa || echo "run 10h_apply_names.py first (step 5)"
ALN=named/exendin.aa.fa
CODON=named/exendin.codon.fa
IQTREE=$(command -v iqtree3 || command -v iqtree2 || command -v iqtree)
OG=$(cat named/outgroup.txt); echo ${OG}
```

## 1. Split the alignment
```
python3 ${SCRIPTDIR}/12a_partition_alignment.py --aln ${ALN} --codon-aln ${CODON} \
    --reference ${SCRIPTDIR}/4_reference_peptides.fa --outdir partitioned 2> partitioned.log
cat partitioned.log
cut -f1,2,4,6 partitioned/boundaries.tsv | column -t
```
The boundary is the alignment column where the mature peptide starts, found per sequence by
the same rule as `5a_assign_nomenclature.py` and taken as the consensus.
- **`N/N sequences agree`** — if the range spans more than a few columns the alignment is
  ragged over the cleavage site. Fix the alignment, or set `--boundary` by eye from panel C,
  before trusting either partition.
- **`found by 'cleavage_site', not the motif`** — that sequence's hormone domain does not
  match the canonical motif (expect this for the divergent anguid copies and for anything
  lacking His1). Check its row in `boundaries.tsv`.
- **GCG and GIP** carry more than one hormone domain, so their boundary may land elsewhere.
  If they disagree with the rest, drop them from the consensus with
  `--exclude HomSap_GCG,GalGal_GCG,HelSus_GCG,HomSap_GIP,GalGal_GIP,HelSus_GIP` — they stay
  in both partitions, they just stop voting on where the boundary is.

## 2. Topology from the slow block
```
$IQTREE -s partitioned/slow.aa.fa -m MFP -B 1000 -alrt 1000 -o ${OG} -T AUTO \
    --prefix partitioned/slow -redo
grep "Best-fit model" partitioned/slow.iqtree
python3 ${SCRIPTDIR}/4h_draw_tree.py --tree partitioned/slow.treefile --groups named/groups.tsv \
    --outgroup ${OG} --title "Signal peptide + pro-region only" --png partitioned/slow.png \
    > partitioned/slow.svg
```
Compare this with `named/exendin_aa.treefile` from step 4. The question is whether the
exendins move from sister-to-all-GCGL to nested within squamate GCGL once the fast block is
removed. If they do, step 4's topology was an artifact and this is the tree to report.

**If the slow block cannot resolve the deep nodes either** — support still below 70 — that is
a real result, not a failure: ~50 columns of signal peptide carry little information. Say so,
and fix it by adding taxa rather than sites: GCGL orthologs from more vertebrates (turtle,
crocodile, *Xenopus*, teleosts) break the long branches that cause the problem in the first
place.

## 3. Re-run the origin test on the slow block
The H1/H2/H3 hypotheses from step 4 are exactly what this partition is meant to decide.
```
python3 ${SCRIPTDIR}/10e_tree_inputs.py --aln partitioned/slow.aa.fa --outdir partitioned/hyp 2>&1 | tail -5
M=$(grep "Best-fit model" partitioned/slow.iqtree | awk '{print $NF}')
for H in H1_squamate_GCGL H2_sister_GCGL H3_GCG; do
  $IQTREE -s partitioned/slow.aa.fa -m ${M} -g partitioned/hyp/${H}.constraint.nwk -o ${OG} \
      -T AUTO --prefix partitioned/${H} -redo
done
cat partitioned/slow.treefile partitioned/H1_squamate_GCGL.treefile \
    partitioned/H2_sister_GCGL.treefile partitioned/H3_GCG.treefile > partitioned/candidates.nwk
$IQTREE -s partitioned/slow.aa.fa -m ${M} -z partitioned/candidates.nwk -n 0 -zb 10000 -au \
    -T AUTO --prefix partitioned/AUtest -redo
grep -A8 "USER TREES" partitioned/AUtest.iqtree
```
Report this AU test alongside the step-4 one. If the two partitions reject different
hypotheses, the fast block was driving step 4's answer and the slow block's is the one to
believe.

## 4. Map the fast block onto the slow topology
`-te` fixes the topology and re-estimates only branch lengths, so the two trees differ in
lengths alone.
```
$IQTREE -s partitioned/fast.aa.fa -m MFP -te partitioned/slow.treefile -T AUTO \
    --prefix partitioned/fast_on_slow -redo

python3 ${SCRIPTDIR}/12b_branch_rates.py \
    --slow-tree partitioned/slow.treefile --fast-tree partitioned/fast_on_slow.treefile \
    --highlight HelSus_A,HelSus_B,HelCha_A,HelCha_B \
    --highlight-name "Heloderma venom copies" \
    --out-prefix partitioned/rates 2> partitioned/rates.log
cat partitioned/rates.log
```
`--highlight` takes tip names in the step-5 nomenclature; the test set is their terminal
branches plus any internal branch lying wholly inside them. Set it to whichever copies your
naming calls the venom genes — `cut -f2,4 names/names.tsv | grep copy` lists the candidates,
and `names/names.tsv` records which letter each published exendin ended up as.

Output: a per-branch table (slow length, fast length, ratio, log2 shift against the tree
median), a permutation test of the highlighted branches against the rest, and a tree with
branches coloured by rate shift.

**Reading it.** A positive median shift on the *Heloderma* branches with a small p-value
confirms Irwin & Prentice with more taxa. Three other outcomes are worth as much:
- **The anguid copies (`ElgCoe_Y`, `DopGra_Z`) also show a large shift.** Then rapid mature-peptide
  evolution is not a venom phenomenon, and the venom story needs rethinking.
- **The shift is on the stem branch, not the tips.** The burst happened once, at recruitment,
  and the copies have been conserved since — consistent with helodermatid venom composition
  being conserved over ~30 My.
- **No shift anywhere.** Then the deep-node problem is not rate heterogeneity and you should
  look at composition bias or alignment error instead.

## 5. One rate per partition, and selection
A partitioned run gives the overall rate ratio in a single number:
```
$IQTREE -s ${CODON} -p partitioned/partitions.nex -m MFP -B 1000 \
    -T AUTO --prefix partitioned/both -redo
grep -A6 "Partition information" partitioned/both.iqtree       # relative rate per partition
```
Then run the selection tests on the **fast codon partition**, where the signal is, with the
venom branches as foreground. Testing the whole precursor dilutes a mature-peptide signal
across ~50 invariant codons:
```
hyphy relax --alignment partitioned/fast.codon.fa --tree partitioned/slow.treefile \
    --test "Heloderma" --output partitioned/relax_fast.json
hyphy absrel --alignment partitioned/fast.codon.fa --tree partitioned/slow.treefile \
    --output partitioned/absrel_fast.json
```
RELAX needs the foreground branches labelled in the tree (`{Heloderma}` tags); add them by
hand or with `hyphy label-tree`. Run the same two on `slow.codon.fa` as the control — the
prediction is that the signal and pro-regions show nothing.

## Outputs
| File | Content |
|---|---|
| `partitioned/boundaries.tsv` | per sequence: how the mature peptide was found and which column it starts at |
| `partitioned/{slow,fast}.{aa,codon}.fa` | the two partitions |
| `partitioned/partitions.nex`, `partitions_aa.nex` | partition files for IQ-TREE / RAxML |
| `partitioned/slow.{treefile,iqtree,png}` | topology from the signal + pro-region |
| `partitioned/AUtest.iqtree` | origin test (H1–H3) on the slow block |
| `partitioned/fast_on_slow.treefile` | mature-peptide branch lengths on that topology |
| `partitioned/rates_branches.tsv`, `rates.{svg,pdf,png}` | per-branch rate shift and the permutation test |
| `partitioned/both.iqtree` | relative rate of each partition in one model |
