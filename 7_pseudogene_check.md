# Is the unmatched Heloderma copy (HelSus_B / HelCha_C) a gene or a pseudogene?
**Why.** Both *Heloderma* genomes carry one exendin-family copy that matches no published
peptide (step 5): **HelSus_B** and **HelCha_C**.
- **Where it sits.** In HelSus it lies between GCGL and the exendin-2 gene (step 6,
  panel A).
- **Tree.** The two copies pair in the step-5 tree (93.5/99) as orthologs. They group
  with DopGra_A among the anguimorph array copies, not with the HelSus_EX2 / HelCha_EX1
  pair.

Matching no venom peptide could mean a gene that is expressed but was never purified, or
a gene that has decayed. This step tells them apart by comparing each copy, exon by exon,
with the helodermin (exendin-2) precursor. The intact exendin-2 genes (HelSus_EX2 and its
HelCha ortholog HelCha_EX1) are run alongside as controls.

**Method** (`7a_pseudogene_check.py`)
1. **Map the precursor.** The 84-aa helodermin precursor (UniProt P04204: signal peptide
   1–22, propeptide 23–46, mature helodermin 47–84) is searched with tblastn against each
   locus ± 10 kb. The hits are chained into exons as in step 6: one codon per precursor
   residue, ends extended over mismatches. Any stretch of 5 or more residues still
   uncovered is searched again on its own, between its neighbouring blocks (step 4e's
   gap fill), so a short diverged exon is not missed.
2. **Read the junctions** on the gene's strand.
   - **Gap < 60 bp (inside one exon).** The genomic length minus 3 × the residues it
     should hold is an indel. Not a multiple of 3 = **frameshift**.
   - **Longer gap (intron).** A donor GT (or GC) and acceptor AG that keep the reading
     frame are looked for near the block ends. None found = **no splice site**.
3. **Assemble the CDS.** The exons are joined, the first extended back to the start
   codon and the last on to its stop (only if ≤ 6 residues are missing at that end; more
   means the exon was not found). The CDS is translated.
4. **Check the coding sequence.** It reports:
   - the start codon (ATG or not);
   - in-frame **stop codons**, and where the ORF ends;
   - identity to the precursor, overall and per region, counted codon by codon so a
     frameshift does not hide how similar the rest is;
   - the **KR** processing site before the mature peptide.
5. **Verdict.**
   - **intact**: ATG; no stop, frameshift or splice loss; ≥ 90% of the precursor found;
     a stop where the precursor ends.
   - **disrupted**: any stop, frameshift, lost splice site, or an ORF that ends early;
     `reasons` lists them.
   - **partial**: anything else, e.g. an exon not found but no disruption.
6. **Selection (dN/dS).** Nei-Gojobori dN/dS for each `--pair`, over codons both loci
   place on the same precursor residue.
   - A pseudogene evolves neutrally (dN/dS near 1); a functional gene is constrained
     (well below 1).
   - The key comparisons are:
     - **HelSus_B vs HelCha_C**: the same copy in two species, dated by their split;
     - the control pair **HelSus_EX2 vs HelCha_EX1**;
     - **HelSus_B vs HelSus_EX2**: the paralogs.
7. **TOGA2 status (`--toga2-run`).** TOGA2's own status of the projection at each locus
   is added: FI/I = intact, PI = partly intact, L/UL = lost or uncertain loss, M =
   missing, PG = pseudogene. It is an independent opinion on the same question.

**Figure** (`7b_plot_pseudogene_check.py`)
- **A. One row per locus, in precursor coordinates, under the precursor.**
  - **Cells.** Each cell is what the locus's codon for that residue encodes: tinted =
    identical, red letter = different, red box = stop, dot = not found.
  - **Marks.** Exon boundaries (intron length, donor..acceptor; red ✗ if no splice
    site), frameshifts (purple ▲ with the indel), and where the ORF ends.
  - **Verdict** is on the left.
- **B.** The summary table.
- **C.** dN/dS.

**Files** (all new):
- `7a_pseudogene_check.py`
- `7b_plot_pseudogene_check.py`

The reference precursor comes from `4_reference_peptides.fa` (`EXENDIN2_precursor`,
`EXENDIN2_mature`).

Run the set-up block of `0_Environment_management.md` first.
```
mkdir -p ${WORK}/pseudogene_check
cd ${WORK}/pseudogene_check

python3 ${SCRIPTDIR}/7a_pseudogene_check.py --species-table ../species.tsv --genomedir ${GENOMEDIR} \
    --curated ${WORK}/reciprocal_blast/reciprocal/curated_loci.tsv \
    --names ${WORK}/nomenclature/names_names.tsv \
    --locus HelSus_B --locus HelCha_C --control HelSus_EX2 --control HelCha_EX1 \
    --pair HelSus_B,HelCha_C --pair HelSus_EX2,HelCha_EX1 --pair HelSus_B,HelSus_EX2 \
    --toga2-run ${TOGADIR}/target_HelSus --outdir check 2> check.log
cat check.log
python3 ${SCRIPTDIR}/7b_plot_pseudogene_check.py --check-dir check --out pseudogene_check
```
- **Names.** `--locus` / `--control` take the step-5 names (or step-4 names without
  `--names`). Check `nomenclature/names_names.tsv` that HelSus_B and HelCha_C are the
  copies you mean: positional letters change if step 5 is re-run after the step-4 fixes.
- **Speed.** Each locus is searched only in its own ±10 kb window, so this takes seconds;
  no genome-wide search.

**Outputs** (`check/`)
- `summary.tsv`: per locus, the verdict and its reasons, exons, intron lengths, splice
  dinucleotides, start codon, ORF length and first stop, identities (overall / signal /
  propeptide / mature), KR site, TOGA2 status, and the aligned proteins.
- `exons.tsv`: per exon, the precursor residues, genomic coordinates, and the
  donor/acceptor and intron after it.
- `events.tsv`: every disruption: stop, frameshift, in-frame indel, splice, start, ORF end.
- `residues.tsv`: codon by codon, as in step 6.
- `cds.fna`, `protein.faa`: the assembled CDSs and translations, with the precursor.
- `pairs.tsv`: codons compared, synonymous and nonsynonymous differences, dS, dN, dN/dS.
- `pseudogene_check.png` / `.svg`: the figure.

**Reading it**
- **Controls first.** HelSus_EX2 and HelCha_EX1 must come out **intact**. If they don't,
  the mapping is wrong for this gene structure (e.g. the ~10-residue last exon was
  missed), not the gene. Fix that before reading the test loci.
- **Disruptions shared by both species** (same stop or frameshift in HelSus_B and
  HelCha_C) arose before *H. suspectum* and *H. horridum* split, and the gene has been
  dead since. Disruptions in only one species are recent.
- **No disruption but a diverged mature peptide** (low `mature` identity, KR intact, ATG,
  stop in place) is an intact gene that encodes a different, unpurified peptide. Then
  dN/dS between HelSus_B and HelCha_C should be well below 1.
- **dN/dS needs differences.** HelSus and HelCha are close relatives, so a pair may
  differ at only a handful of sites (`Sd`, `Nd` in `pairs.tsv`). With fewer than ~5 of
  each, treat dN/dS as a direction, not a number.
- **A missing first exon** ("precursor 1–20 not found") is common for diverged
  signal-peptide exons and is not evidence of decay by itself. It counts against the
  gene only if the remaining exons also carry disruptions.
- **Compare with TOGA2's status.** An intact verdict where TOGA2 says L or PG (or the
  reverse) is worth opening in the alignment (panel A).
