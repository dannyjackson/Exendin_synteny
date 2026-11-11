# Exendin gene tree from TOGA2: collection, curated alignment
**Goal:** one alignment of every exendin gene TOGA2 annotated, across all four TOGA2 runs.
Every curation decision is recorded in a table, so the alignment can be rebuilt from the
TOGA2 files at any time. The tree is built from the curated alignment (next section, to come).

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

## 2. First alignment
Codon-aware alignment with MACSE, as in step 4. Frameshifts and internal stops become `NNN`/`X`
instead of breaking the frame.
```
mkdir -p aln0
macse -prog alignSequences -seq collect/exendin_cds.fa -out_NT aln0/exendin.macse_NT.fa -out_AA aln0/exendin.macse_AA.fa
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
Edit `collect/curation.tsv` (copy it first: `cp collect/curation.tsv curation.tsv`). It has one
row per tip, all `keep`, with automatic notes. Change `action` and give the reason in `note`:

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
    --add-cds rescued/HelCha_Ex3.cds.fa --out-prefix curated/exendin 2> curated.log
cat curated.log
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
