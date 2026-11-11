# Which species was each published exendin purified from?
**Why.** The four published peptides were purified from Gila monster or beaded lizard
venom, but the species labels are uncertain:
- **Exendin-1 (helospectin).** Named from *H. suspectum*, yet Vandermeers et al. 1987
  found helospectins typical of *H. horridum* venom.
- **The other three.** Helodermin (exendin-2) and exendin-4 are reported from
  *H. suspectum* and exendin-3 from *H. horridum*. The venom came from animals sold under
  those names before *H. horridum* was split into four species.

A genome resolves this: a peptide encoded **exactly** by a gene of one species, and only
approximately by every gene of the other, was most likely purified from the first.
Step 6 asked where each peptide maps best within one genome. This step maps every
peptide onto **every** curated exendin-family gene of both genomes, so the best match in
each genome and the gap between the two can be read side by side.

**Method** (`8a_published_vs_loci.py`)
- **Peptides.** `8_published_peptides.fa` holds the mature exendin-1 (helospectin I),
  exendin-2 (helodermin, 35 aa, Glu8,9), exendin-3 and exendin-4. Helospectin II is
  helospectin I without Ser38, so it is left out.
- **Loci.** Every curated locus (step 4, with its step-5 name) of each genome in
  `--species` (default HelSus, HelCha).
- **Mapping, per peptide and locus:**
  - a sensitive tblastn of the peptide against the locus ± 5 kb;
  - chained through the best hit inside the locus, with 7a's splice-aware chaining and
    gap fill (the mature exendin-3/-4 span two exons; the last is ~10 residues);
  - one codon per residue, with ends extended over mismatches;
  - the genome residues just before residue 1 (the processing site, KR) and just after
    the last residue (Gly = amidation signal) are read from the genome.
- **Ranking.** Within each genome, loci are ranked per peptide. The best genome for a
  peptide is the one whose best locus matches most residues; the margin is how many more.

**Figure** (`8b_plot_published_vs_loci.py`), one plot:
- **Summary (top).** One row per peptide: who reported it from which venom; the best gene
  of each genome with its identity (shaded by identity, ✓ for an exact match); the best
  genome and its margin.
- **One block per peptide.**
  - **Rows.** The published sequence, then the loci of each genome (those the peptide
    hits, best first; the rest listed on one line).
  - **Cells.** What each locus's codon for each residue encodes: tinted = identical, red
    letter = different, red box = stop, dot = not found, | = exon boundary.
  - **Margins.** The genome flanks in grey and the identity on the right. ★ marks the
    best gene of each genome.

**Changes to earlier scripts in this round** (re-run steps 6 and 7):
- `7a_pseudogene_check.py`. Two fixes:
  - In the chaining, an aligned block now earns score only for the precursor residues it
    **adds**. Before, a second Pro-rich block re-covering residues 76–84 earned its full
    score, so the controls' last exon was taken from an off-target block (HelSus_EX2:
    NSPHPPAP, 2,905 bp, in place of the RTSPPP exon 2,137 bp away), and they came out
    "partial".
  - A splice junction must now keep an **open** reading frame. The coding bases it adds
    may contain no stop, and a long diverged stretch may differ by up to half its length
    in codons. Before, the B/C copies' diverged propeptide forced a rare GC donor that
    ran into a stop ("ORF ends" after residue 28).
- `6a_map_published.py`. The same "residues it adds" rule in its chaining; it was behind
  helospectin's three-block mapping in step 6.

Run the set-up block of `0_Environment_management.md` first.
```
mkdir -p ${WORK}/published_vs_loci
cd ${WORK}/published_vs_loci

python3 ${SCRIPTDIR}/8a_published_vs_loci.py --species-table ../species.tsv --genomedir ${GENOMEDIR} \
    --curated ${WORK}/reciprocal_blast/reciprocal/curated_loci.tsv \
    --names ${WORK}/nomenclature/names_names.tsv --species HelSus,HelCha --outdir match 2> match.log
cat match.log
python3 ${SCRIPTDIR}/8b_plot_published_vs_loci.py --match-dir match --out published_vs_loci
```
Each locus is searched only in its own window, so this takes about a minute.

**Outputs** (`match/`)
- `best.tsv`: per peptide, the best locus and identity in each genome, the genome
  flanks there, the best genome and the margin.
- `matches.tsv`: every peptide × locus: identity, exons, flanks, coordinates, and the
  locus's rank within its genome.
- `residues.tsv`: codon by codon.
- `published_vs_loci.png` / `.svg`: the figure.

**Reading it**
- **One genome exact, the other not.** That is the evidence for the species of origin.
  - The margin matters most for exendin-3 vs exendin-4. They differ only at residues
    2–3, so each will match the other's gene at 37/39 and its own at 39/39.
  - A margin of 2 there is the whole signal, so check those two residues in the blocks.
- **Exact in neither genome.** The peptide may come from a population or species not
  sequenced here (*H. horridum* sensu lato is four species; HelCha is one), from a gene
  missing from the assembly, or from a sequencing error in the 1980s protein data.
  - Compare the best gene's differences with known sequencing ambiguities: Glu/Gln at
    helodermin 8–9 was such a case.
- **Exact in both.** The peptide is conserved between the species and gives no
  provenance. Report it as such.
- **Flanks.** KR before residue 1 and G after the last residue are the processing and
  amidation signals. Exendin-3/-4 are amidated; helodermin and helospectin are not
  reported as amidated. A best gene with the wrong flanks is worth a second look.
