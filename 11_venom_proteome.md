# Is each exendin gene's protein in the venom?
**Question:** which of the named exendin proteins (step 5) have been seen in a *Heloderma*
venom proteome, and — where one has not — would it have been seen if it were there?

The paralogs share their signal peptide and pro-region almost exactly, so most tryptic
peptides cannot tell `<Sp>_Exendin_anc` from `<Sp>_B`. Only peptides over the divergent regions
are **diagnostic**. A negative result is worth reporting only if diagnostic peptides existed
and were in the size range an LC-MS/MS run reports. `11a_search_venom_proteome.py` works that
out and then asks whether any were observed.

## 1. Get a peptide list
The published gila monster venom proteome is ProteomeXchange **PXD001343** (Sanggaard et al.
2015, *J Proteomics* 117:1-11, with the companion *Data in Brief* 3:137-42 listing the
identified proteins).
```
mkdir -p ${WORK}/venom_proteome && cd ${WORK}/venom_proteome
# identified proteins and peptides: supplementary tables of the Data in Brief paper (PMC4510068)
# raw/processed MS: ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2015/01/PXD001343
```
The script takes a plain list of peptide sequences, or any TSV/CSV in which it finds the
peptide column itself (PEAKS, MaxQuant `peptides.txt`, an mzIdentML export converted to
text). `K.PEPTIDER.S` flanks and `(+15.99)` modifications are stripped. Save it as
`PXD001343_peptides.tsv`.

If you only have the list of *identified proteins* rather than peptides, that is much weaker:
identification was by homology, so a divergent paralog could have been reported under the
name of its closest known relative. Prefer the peptide-level data.

## 2. Build a background
Any protein that was in the sample competes for a peptide. A peptide shared with a kallikrein
is not evidence for an exendin.
```
# Heloderma proteins from UniProt, as a background set
curl -s "https://rest.uniprot.org/uniprotkb/stream?query=taxonomy_id:34736&format=fasta" \
    > heloderma_uniprot.fa
grep -c '>' heloderma_uniprot.fa
```

## 3. Run
```
python3 ${SCRIPTDIR}/11a_search_venom_proteome.py \
    --proteins ${WORK}/toga2_exendin_tree/named/exendin.aa.fa \
    --peptides PXD001343_peptides.tsv \
    --background heloderma_uniprot.fa \
    --outdir search 2> search.log
cat search.log
```
`named/exendin.aa.fa` is the step-5-renamed alignment (10h); gaps are stripped on reading,
so the aligned file is fine as a source of protein sequences and every result comes out in
the `<Sp>_Exendin_anc` / `<Sp>_A` nomenclature. Pass **every** paralog in `--proteins`, not a subset: a peptide shared with a paralog you left
out would look diagnostic. `--only` limits what is printed without changing that.

## 4. Verdicts
| Verdict | Meaning |
|---|---|
| `DETECTED` | a peptide matching this protein and no other was observed |
| `AMBIGUOUS` | peptides were observed, but all are shared with another paralog — no evidence either way for this copy |
| `ABSENT - informative` | not observed, and diagnostic peptides were in the detectable range: the protein was looked for and not found |
| `ABSENT - uninformative` | not observed, and no diagnostic peptide was detectable: this dataset could never have distinguished it. Report it as untested, not absent |

`search/observed_matches.tsv` gives each observed peptide and the proteins it matches;
`search/peptides.tsv` gives every in-silico peptide with its length, mass, whether it is
diagnostic, and which paralogs share it.

## 5. Caveats specific to this dataset
These belong in the text, not just your notes:
- **2D-PAGE, major spots only.** The venom was run on a 2D gel and only major
  Coomassie-stained spots were excised for digestion. Anything of low abundance was never
  analysed.
- **Small peptides.** A 2D gel retains and stains ~4-10 kDa peptides poorly. The exendin
  precursors are ~10 kDa and the mature peptides ~4 kDa, so the method is biased against
  exactly this size class. The `note` column flags proteins under 12 kDa for this reason.
- **De novo sequencing.** Identification used PEAKS de novo sequencing with cross-species
  homology search. De novo cannot separate Ile from Leu, so the script treats them as
  identical by default (`--exact` turns that off). Homology-based identification also means
  a divergent paralog may have been reported under a known relative's name.
- **Species and individual.** The dataset is *H. suspectum suspectum* from one source.
  Venom composition varies; an absence in one animal is not an absence in the species.

A negative result here therefore supports "not a major venom component", not "not in venom".
The decisive experiment is venom-gland RNA-seq (step 6) or targeted MS (PRM/MRM) against the
diagnostic peptides this script lists — `peptides.tsv` filtered to `diagnostic=yes` and
`detectable=yes` is the transition list to order.
