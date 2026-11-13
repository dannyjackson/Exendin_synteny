## 8. HelCha and the bad reference transcript
TOGA2 copies the reference transcript's exon structure onto each query. The HelSus EviAnn
reference has two isoforms of exendin-4 (`LOC_00012325`):
- **`mRNA-1`** starts with the canonical first exon.
- **`mRNA-2`** starts with a 419-bp alternative first exon whose CDS has no protein support.

At HelCha's exendin-3 locus TOGA2 projected only `mRNA-2`. The model therefore took the
wrong first exon, gaining ~110 residues (198 aa), and missed the real exon 1 ~5.9 kb away.
Section 1 excludes `mRNA-2` and 1b rebuilds the locus from the published *H. horridum* cDNA
(AJ580309). This figure documents the problem and the fix:

| Panel | Shows |
|---|---|
| **A** (the error) | The HelCha exendin-3 locus, exon by exon, 5′ → 3′: every TOGA2 projection over the curated exons, with status and protein length, above the curated model with its name from `names.tsv`. TOGA2 exons the curated model lacks are outlined red, with their size. Dashed lines carry the curated exon ends up through the TOGA2 rows, so every boundary difference shows. A line says so when `mRNA-1` has no projection there. |
| **B** (cause) | Every isoform of `LOC_00012325` in the HelSus reference, exon by exon, 5′ → 3′, on the **same bp scale as A**. The exon only `mRNA-2` has is outlined red, with its size. Its position matches the red exon of A directly above. |
| **C** (evidence) | The HelCha **genome itself**, read at the curated exons, spliced, translated, and aligned residue by residue to the published exendin-3 precursor (P20394). Differences are red and the mature exendin-3 is bracketed. Each residue sits on the colour of its exon, and codons split by an intron are underlined. Below the alignment, each exon's residues and genomic coordinates are listed. |

The curated model's CDS in A comes from the same genome translation as panel C.
Panel C translates the genome, not the 1b output file. It reads the exon coordinates of
`rescued/HelCha_Ex3.exons.tsv` from the genome FASTA (through its `.fai` index; samtools is
not needed), so every residue in the figure traces to a base on `JANEZZ010000386.1`.

```
cd ${WORK}/toga2_exendin_tree
HELCHA_FASTA=$(ls ${GENOMEDIR}/HelCha_GCA_026122225.1/ncbi_dataset/data/GCA_026122225.1/*_genomic.fna)
python3 ${SCRIPTDIR}/10j_plot_helcha_transcript_issue.py --species HelCha \
    --toga2 ${TOGADIR}/target_HelSus/query_HelCha \
    --ref-gff ${HELSUS_EVIANN_GFF} --names names/names.tsv \
    --ref-gene LOC_00012325 --bad-tx LOC_00012325-mRNA-2 \
    --rescued rescued/HelCha_Ex3.exons.tsv --rescued-cds rescued/HelCha_Ex3.cds.fa \
    --genome ${HELCHA_FASTA} --genome-contig JANEZZ010000386.1 \
    --published ${SCRIPTDIR}/4_reference_peptides.fa \
    --svg helcha_transcript_issue.svg --pdf helcha_transcript_issue.pdf \
    --png helcha_transcript_issue.png 2> helcha_transcript_issue.log
cat helcha_transcript_issue.log
```
The log reads like this, and every number in it should match sections 1 and 1b:
- **`A: TOGA2 LOC_00012325-mRNA-2#… I …, N aa; exon(s) not in the curated model: … (419 bp)`**:
  the projection behind the old HelCha model, with ~198 aa as in section 1's log, and
  its wrong first exon.
- **`A: other LOC_00012325 isoforms projected onto HelCha: …; at this locus: none`**: where
  `mRNA-1` went, if anywhere.
- **`B: LOC_00012325-mRNA-2: …; exon(s) no other isoform has: … (419 bp)`**: the cause,
  read straight from the GFF. If no exon is listed, `--bad-tx` is spelled differently
  in the GFF.
- **`C: JANEZZ010000386.1 spliced at 3 exon(s) = N bp; the CDS of HelCha_Ex3.cds.fa at base …:
  261 bp = 87 aa`**: the genome read and translated. If the script stops with `the CDS … is
  not in the genome spliced at …`, then `--genome` or `--genome-contig` is wrong, or the exon
  table and CDS come from different 1b runs. Without `--rescued-cds`, the longest ORF is
  used, as 1b does.
- **`C: P20394 exendin-3 precursor: 87 aa; identical at 87/87 positions (100.0%)`**: the
  result. Section 1b expects an identical 87-aa protein. Anything less is drawn in red in
  panel C, and is worth a look before the figure is used.
- **`C: exon N: residues …, coding bases JANEZZ010000386.1:…`**: the same table as under
  panel C. Exon 1 should include 4,399,218–4,399,326 (section 1b).

Options:
- `--genome-contig` is the contig's name in the FASTA. It also picks TOGA2's contig
  (`JANEZZ010000386__1` = `JANEZZ010000386.1`). Without it, the TOGA2 contig with
  projections over the curated exons is used, and the script stops if there are several.
- `--published-id` / `--mature-id` (`EXENDIN3_precursor` / `EXENDIN3_mature`) choose the
  records in `4_reference_peptides.fa`. Use `--mature-id ''` for no bracket.
- `--flank-kb` (1) sets how far A and B extend past the transcripts.
- `--width-mm` and `--font-pt` work as in the other figures. The alignment wraps onto
  further lines if it does not fit.

For another genome or gene, change `--species`, `--toga2`, `--ref-gene`, `--bad-tx`,
`--rescued*`, `--genome*` and `--published-id`.

| File | Content |
|---|---|
| `helcha_transcript_issue.{svg,pdf,png}` | the three-panel figure |

The full interval for HelCha, with every TOGA2 gene and the curated names, is in the
step-5 figure (`exendin_figure.pdf`, panel A).
