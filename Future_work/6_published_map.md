# Where the published exendins are encoded in the Gila monster genome
**Why.** Exendin-1 (helospectin I and II), exendin-2 (helodermin) and exendin-4 were named
from peptides purified from *Heloderma* venom between 1984 and 1992. Step 5 names genes by
comparing each gene's predicted peptide with these. This step goes the other way: it
starts from each purified peptide and finds the codons that encode each of its residues in
the HelSus genome, with no gene model in between. That shows which gene, and which exon,
carries each published name.

**Helospectin I and II.** They are the same peptide except for one C-terminal residue
(helospectin II = helospectin I without Ser38). Both are mapped independently, and panel B
colours each residue by whether the two peptides use the same codon for it:
- If both come from one gene, all 37 shared residues match and only Ser38 is
  highlighted. Helospectin II's genome flank then shows the codon after its end: Ser
  means the genome encodes the residue it lacks (C-terminal trimming); anything else
  means a second gene.
- If they map to different loci, the panel shows both loci, one above the other.

**Method** (`6a_map_published.py`):
1. **Genome-wide search.** tblastn against the genome, with masking and composition
   adjustment off (`-seg no -comp_based_stats 0`) and word size 3, as short queries need.
2. **Candidate loci.** Hits are grouped into loci. Each locus is scored by its best
   collinear chain of HSPs around its strongest HSP, with at most 4 exons, introns of at
   most 10 kb, and every other HSP at least 14 bits.
3. **Window re-search.** The best 5 loci are searched again with a sensitive tblastn,
   each in a ±5 kb window. This finds the ~10-residue last exon (GPSSGAPPPS), which a
   genome-wide search misses (as in 4e's gap fill).
4. **Choice and codons.** The locus with the most identical residues wins. Each residue
   gets one codon; where two exons overlap in the query, a single cut point is chosen.
5. **Flanks.** The genome residues just before residue 1 (the processing site, e.g. KR)
   and just after the last residue (e.g. G, the amidation signal) are translated from
   the same exon.
6. **Genes.** Overlapping NCBI genes are reported with the exon number. A published
   peptide carried by an "uncharacterized protein LOC…" is named here.

**Figure** (`6b_plot_published_map.py`):
- **A.** The whole chromosome with the region boxed, and below it the region with every
  NCBI gene in it, out to SLC4A8 and AQP5 if they are within 3 Mb. Exendin-family genes
  are in their name colour; genes carrying a mapped peptide are outlined. A marker under
  each mapped locus is joined to its pop-out.
- **B–D.** One pop-out per published exendin, each showing:
  - the NCBI gene models of the locus, with each peptide's exons beneath them;
    minus-strand loci are drawn reversed, so every peptide reads left to right;
  - the peptide above the genome's translation of its codons: identical residues tinted,
    mismatches in red, stops filled red, the intron length between exons, and the
    genome flanks in grey;
  - identity, the overlapping gene, and the runner-up locus.

**Files**

| File | Origin |
|---|---|
| `6_published_peptides.fa` | new: the mature peptides (helodermin as corrected by Vandermeers et al. 1987: Glu8,9, 35 aa) |
| `6a_map_published.py` | new: mapping |
| `6b_plot_published_map.py` | new: figure |
| `lib_gff.py` | new: GFF3 reader (genes, transcripts, exons, CDS of a region) |

Exendin-3 is from *H. horridum*, so it is left out. To map it anyway (e.g. onto HelCha),
add a line to `6_published_peptides.fa`; it gets its own pop-out automatically.

Run the set-up block of `0_Environment_management.md` first.
```
mkdir -p ${WORK}/published_map
cd ${WORK}/published_map

python3 ${SCRIPTDIR}/6a_map_published.py --species-table ../species.tsv --genomedir ${GENOMEDIR} \
    --code HelSus --peptides ${SCRIPTDIR}/6_published_peptides.fa --threads 8 --outdir map

python3 ${SCRIPTDIR}/6b_plot_published_map.py --map-dir map --out published_map \
    --toga2-gff ${HELSUS_EVIANN_GFF} --boundaries SLC4A8,AQP5=LOC_00012326 \
    --curated ${WORK}/reciprocal_blast/reciprocal/curated_loci.tsv \
    --names ${WORK}/nomenclature/names_names.tsv 2> published_map.log

cat published_map.log


####### repeat with HelCha

HELSUS_EVIANN_GFF=${PROJDIR}/datasets/genomes/Heloderma_suspectum/Heloderma_suspectum.final.reference.gff3
mkdir -p ${WORK}/published_map_helcha && cd ${WORK}/published_map_helcha

python3 ${SCRIPTDIR}/6a_map_published.py --species-table ../species.tsv --genomedir ${GENOMEDIR} \
    --code HelCha --peptides ${SCRIPTDIR}/6_published_peptides_helcha.fa --threads 8 --outdir map

python3 ${SCRIPTDIR}/6b_plot_published_map.py --map-dir map --out published_map_helcha \
    --toga2-run ${TOGADIR}/target_HelSus --toga2-query HelCha \
    --toga2-gff ${HELSUS_EVIANN_GFF} --boundaries SLC4A8,AQP5=LOC_00012326 \
    --curated ${WORK}/reciprocal_blast/reciprocal/curated_loci.tsv \
    --names ${WORK}/nomenclature/names_names.tsv 2> published_map_helcha.log

cat published_map_helcha.log

```
The genome-wide search reuses step 4's BLAST database (`${GENOMEDIR}/blastdb/HelSus_<accession>`)
and takes a few minutes. Its hits are kept in `map/genome_hits.tsv`, so changing a
chaining option and re-running 6a does not search again; delete that file to force a new
search.

**Outputs**
- `map/map.tsv`: one row per peptide, with
  - the locus, strand, exon count and intron lengths;
  - identical residues out of the peptide length;
  - the genome residues before and after the peptide;
  - the overlapping NCBI gene(s) and exon(s);
  - the runner-up locus with its identity.
- `map/residues.tsv`: one row per residue, with its codon coordinates, the genome
  residue, and a status (identical / mismatch / stop / gap / unmapped).
- `map/hsps.tsv`, `map/candidates.tsv`: the exons of each mapping, and every candidate
  locus with its scores.
- `published_map.png` / `.svg`: the figure.

**Reading it**
- **Strength of the match.** Check `identical` against the runner-up. *Heloderma* has
  several similar exendin copies, so a clear winner (e.g. 39/39 against 33/39) means the
  name belongs to that gene. Near-ties mean the peptide could come from more than one.
- **Gene name.** The `genes` column gives the NCBI name of the gene that encodes each
  published peptide. If it is "uncharacterized protein LOC146345366", that is the HelSus
  protein that rejected EryReg's interval gene in step 4d (see `--ignore-uncharacterized`).
