# Reciprocal BLAST: curated exendin-family genes, their models, synteny and gene tree
**Why.** Steps 2–3 show what the annotations say, but the array is under-annotated (HelSus:
four of six copies had no gene model; HemCap's GFF has one exendin where a scan found four
loci) and NCBI and TOGA2 name the same gene differently. Here every exendin from steps 2–3
searches every genome, and a locus is kept only if it searches back to an exendin-family
protein.

**Method**
| Step | Script | What happens |
|---|---|---|
| 1 | `4a` | forward queries (NCBI + TOGA2 + published exendin precursors, deduplicated), the HelSus SLC4A8/AQP5 proteins, and a reverse database (HelSus RefSeq proteome + all queries) |
| 2 | `4b` | tblastn of all queries against each whole genome |
| 3 | `4c` | hits → loci with the old interval scanner (`lib_scan_interval.py`); genes split at an intron are joined; boundary genes located |
| 4 | `4d` | blastx of each locus back to the reverse database → ortholog / paralog / rejected, confidence, placement |
| 5 | `4e` | protein + CDS model of each curated locus from its tblastn alignments, with a search for exons the genome-wide search missed |
| 6 | `4f` | synteny figure of the curated set |
| 7 | `4g`, `4h` | gene tree of every curated gene in the region (+ chicken, + published precursors) |

**Classes** (4d). *Ortholog of group G*: best reverse hit is a G protein and this locus is
the genome's best forward hit for the G queries (reciprocal best hit). *Paralog*: best
reverse hit is a G protein, but another locus is the best forward hit for G. *Rejected*: best
reverse hit is a non-family protein (proglucagon, GIP, VIP, PACAP …). **High confidence**
(the curated set): reverse margin ≥ 1.2 over the best non-family protein, forward e-value ≤
1e-5, and both exon 1 and hormone domain found. The names here (`Ex4`, `Ex2p1`, `Ex3_4`) are
inherited from the annotation of each locus's best reverse hit; step 5 replaces them.

Run the set-up block of `0_Environment_management.md` first.
```
cd ${WORK}/reciprocal_blast
HELSUS_PROT=${GENOMEDIR}/HelSus_GCF_057233095.1/ncbi_dataset/data/GCF_057233095.1/protein.faa
```

## 1. Query set and reverse database
```
python3 ${SCRIPTDIR}/4a_build_query_set.py \
    --ncbi-dir ../synteny_ncbi/ncbi_loci \
    --toga2-tips ../synteny_toga2/toga2_exendins/tips.tsv \
    --toga2-reference ../synteny_toga2/reference_exendins.tsv \
    --published ${SCRIPTDIR}/4_reference_peptides.fa \
    --reverse-proteome HelSus=${HELSUS_PROT} \
    --outdir queries_raw 2> queries.log
cat queries.log
awk -F'\t' 'NR>1 && ($6<70 || $6>150)' queries_raw/query_table.tsv | cut -f1-6   # outside 70-150 aa
```
- **Log.** How many HelSus proteome proteins were named exendin / GCGL / proglucagon (if
  GCGL gets a name `name_group()` does not know, GCGL loci are rejected); and "identical
  sequence named X and Y" warnings, which are the NCBI/TOGA2 naming conflicts (the NCBI or
  published name is kept).
- **Remove broken models before the search.** Precursors are 84–113 aa. The 198-aa
  `HelCha__toga2__Ex3_4__1` is the TOGA2 HelCha exendin-3 model built on the faulty
  `mRNA-2` isoform; its unsupported first exon is a *Heloderma*-specific repeat that made
  ~47,000 loci in HelSus and ~35,000 in HelCha. Drop it (and any other query step 3's check
  shows behaving like it); the published exendin-3 precursor and other species' copies still
  find the gene. `CanAsp__ncbi__GCGL__1` (161 aa) is long but behaves normally: keep it.
```
BAD=HelCha__toga2__Ex3_4__1                       # comma-separated if several
mkdir -p queries
for F in forward_queries.faa reverse_db.faa; do seqkit grep -v -p ${BAD//,/ -p } queries_raw/$F > queries/$F; done
for T in query_table reverse_table; do
  awk -F'\t' -v B="$BAD" 'BEGIN{split(B,b,","); for(i in b) x[b[i]]=1} !($1 in x)' queries_raw/$T.tsv > queries/$T.tsv
done
cp queries_raw/boundary_queries.faa queries/
```
  Everything downstream reads `queries/`. Re-running `4a` writes `queries_raw/` again, so
  repeat this filter after it.
- **More competitors** for distant taxa: add `--reverse-proteome AnoCar=...`.

## 2. Forward tblastn against every genome
Submit as a Biowulf swarm (makeblastdb + tblastn ≈ 10 min per squamate genome):
```
awk -F'\t' -v S=${SCRIPTDIR} -v W=$PWD -v G=${GENOMEDIR} 'NR>1 && $3!="NA"{
  printf "source myconda; mamba activate exendin_env; cd %s; SPECIES_TABLE=../species.tsv GENOMEDIR=%s QUERYDIR=queries OUT=forward THREADS=8 bash %s/4b_forward_tblastn.sh %s\n", W, G, S, $1}' \
  ../species.tsv > forward.swarm
swarm -f forward.swarm -g 16 -t 8 --time 8:00:00 --logdir swarm_logs
```
Writes `forward/<code>.tblastn.tsv`. BLAST databases are built once into
`${GENOMEDIR}/blastdb` (`BLASTDB_DIR`) and reused. Progress of a running job: the last query
with hits in `forward/<code>.tblastn.tsv.tmp` against the query order in
`forward/<code>.queries.faa`. A single new species can be run alone and appended in step 3.

## 3. Call loci
```
python3 ${SCRIPTDIR}/4c_call_loci.py --species-table ../species.tsv --forward-dir forward \
    --queries queries/forward_queries.faa --query-table queries/query_table.tsv \
    --boundary-queries queries/boundary_queries.faa --outdir loci 2> loci.log
grep -v "joined" loci.log
grep -E "joined|not joined" loci.log
awk -F'\t' 'NR>1{n[$2]++} END{for(s in n) print s, n[s]}' loci/forward_loci.tsv | sort -k2,2nr | head
awk -F'\t' 'NR>1{n[$2"\t"$10]++} END{for(k in n) if(n[k]>200) print n[k]"\t"k}' loci/forward_loci.tsv
```
- **Loci per species** should be tens (family copies plus their glucagon-superfamily
  relatives, which step 4 rejects). Thousands in one species = a broken query; the last
  command names it (best query making >200 loci). Remove it as in step 1 and re-run.
- **Split genes are joined.** The scanner can call one gene as two loci: the signal-peptide
  exon (often modelled running on into the intron, with a false stop) and the rest from the
  hormone domain on (HelSus and HelCha GCGL: two pieces 712 bp apart). Neighbours on one
  strand within `--merge-kb` (5) are joined when the queries hitting both cover consecutive
  parts of themselves (the downstream part starts and ends ≥ `--merge-min-offset` (10)
  residues later). Tandem copies cover the same part twice and stay separate
  (ElgMul_Ex2p4/p5: "not joined"). Joined loci carry `;merged` in their state.
- **Positive control.** HelSus must recover the loci of the old copy-number survey (GCGL,
  the Frag1 copy, exendin-2, exendin-4, the two distal copies). If any are missing, lower
  `--min-identity` (35) before trusting an absence elsewhere.
- **Boundaries** (`boundary_hits.tsv`): the best SLC4A8 / AQP5 cluster should have
  `query_cov` near 1; paralogs (SLC4A10, other aquaporins) follow at lower `frac_of_top`.
- **Adding one species later** (after its step-2 job): run `4c` with a one-row species
  table into `loci_<code>`, then append its rows to the three `loci/*.tsv` tables and move
  `loci_<code>/scan/<code>` into `loci/scan/`.

## 4. Reverse blastx and classification
```
python3 ${SCRIPTDIR}/4d_reciprocal_classify.py --species-table ../species.tsv --genomedir ${GENOMEDIR} \
    --loci-dir loci --queries-dir queries --threads 8 --outdir reciprocal 2> reciprocal.log
cat reciprocal.log
column -t reciprocal/summary.tsv
cut -f1,3,8-12,23,25,27,28 reciprocal/reciprocal_loci.tsv | column -t | less -S
```
Minutes for a few hundred loci (one blastx per species). Use a new `--outdir` (or delete
`reciprocal/`) after changing the queries: an existing `reverse_db` is reused.
- **`rejected`**: `rev_best_desc` says what it is (proglucagon, GIP, VIP …).
- **`paralog`**: `reason` names the ortholog locus of the same group.
- **`medium`**: `reason` gives the failed test (margin, e-value, remnant).
- **Cut-offs** without re-running blastx: `--reuse-blastx --min-margin 1.5` (or `--max-evalue`).
- **`curated_loci.tsv`** = the high-confidence set used from here on.
- **`disabled`** = a premature stop, *or* a hormone domain not starting His/Tyr according to
  the old scanner's strict motif. The second also flags every squamate GCGL (its peptide
  `YSEGTLAS…` breaks the motif), so `disabled` is not used to filter anything downstream:
  4e/4g judge stops from the models.

## 5. Protein and CDS models of the curated loci
```
python3 ${SCRIPTDIR}/4e_build_models.py --species-table ../species.tsv --genomedir ${GENOMEDIR} \
    --curated reciprocal/curated_loci.tsv --loci-dir loci \
    --queries queries/forward_queries.faa --query-table queries/query_table.tsv \
    --outdir models 2> models.log
tail -1 models.log
awk -F'\t' 'NR==1 || $10!="ok"' models/models.tsv | column -t -s$'\t'
seqkit grep -r -p '^Hel' models/curated_aa.fa        # Heloderma exendins should end ...PPPS
```
How a model is built:
- **From the tblastn alignments**, for every query of the locus's own name group; the best
  model is kept (hormone domain > coverage ≥ 0.7 > 60–150 aa > collinear > coverage). A
  long but faulty query can win on bitscore and still give a worse model.
- **Exon junctions.** Where two hits overlap in the query, the cut is placed where the two
  together match the query best (a stop counts against its hit), so intron read-through is
  dropped. The CDS is cut codon by codon under the same alignment and always translates to
  the protein.
- **Missing exons.** A short exon cannot pass the genome-wide e-value: the *Heloderma*
  exendins' ~10-residue last exon (`…GPSSGAPPPS`) was missing from every model, cutting each
  mature peptide ~30 residues in. Every stretch of ≥ 5 query residues no hit covers is
  searched again with tblastn in its own small window (≤ `--fill-extend-bp` 10 kb beyond the
  first/last hit, or between two hits); `filled` counts the exons added. `--no-fill` turns
  it off.
- **Hormone domain**: the old motif `[HY][ASG][DEQ][AGP][ATVIGSM]F[TSI]` or the relaxed
  `[HY][ASGE][DEQ][AGP].[FYLI]` squamate copies need (GCGL `YSEGTLAS`, PhrGui `YSEADLGS`,
  ShiCro `HSDGTYTS`); `domain_start` says `(relaxed)`. Residues before the first Met (within
  `--trim-to-met` 20) are cut.

Status other than `ok`: `internal_stop` = pseudogene; `low_cover` / short `length` = remnant
or diverged end; `no_hormone_motif` = check by eye. `--method rescue` runs the old 2e
approach (built for cDNA queries; fragments protein-query models) for comparison only.

## 6. Synteny figure of the curated set
```
python3 ${SCRIPTDIR}/4f_plot_curated_synteny.py --species-table ../species.tsv \
    --loci reciprocal/reciprocal_loci.tsv --boundaries-used reciprocal/boundaries_used.tsv \
    --species-tree ${SPTREE} --top GalGal \
    --svg curated_synteny.svg --png curated_synteny.png --layout curated_synteny.tsv
```
`--confidence high,medium` adds the medium loci (labels end in `?`). Colour = name group of
the best reverse hit; filled = ortholog, hollow = paralog, dashed = 4d `disabled` (which
includes the false GCGL flag of step 4; step 5's named figure draws outlines from the
models instead). Boundaries are tblastn's, so unannotated genomes have them too.
Reading the three figures together: a locus here but not in step 2 = not annotated by NCBI;
not in step 3 = not projected by TOGA2; a colour that differs = a naming conflict (4a's
duplicate warnings).

## 7. Gene tree
Every curated gene in the region, plus chicken's NCBI GCGL protein as outgroup and the
published exendin-2/-3/-4 precursors. Chicken's own genome model is partial (~half the
precursor) and is the same gene as the outgroup, so it is left out.
```
cd ${WORK}/reciprocal_blast
seqkit grep -r -p '^GalGal_' ../synteny_ncbi/ncbi_loci/ncbi_exendin_aa.fa > tree_outgroup.faa
grep '>' tree_outgroup.faa                          # keep one record: the interval GCGL
GG=$(awk -F'\t' '$3=="GalGal"{print $1}' reciprocal/curated_loci.tsv | paste -sd,)   # chicken's curated locus
echo "chicken locus left out: ${GG:-none}"

mkdir -p tree
python3 ${SCRIPTDIR}/4g_tree_inputs.py --curated reciprocal/curated_loci.tsv --models-dir models \
    --outgroup-fasta tree_outgroup.faa --outgroup-group "GCGL (glucagon-like)" \
    --published ${SCRIPTDIR}/4_reference_peptides.fa \
    --reverse-table queries/reverse_table.tsv ${GG:+--exclude $GG} --outdir tree 2> tree/tree_inputs.log
cat tree/tree_inputs.log
awk -F'\t' '$6=="yes" && $7!="-"' tree/tips.tsv     # tips kept with masked stops or flagged models

OG=$(cat tree/outgroup.txt)
IQTREE=$(command -v iqtree3 || command -v iqtree2 || command -v iqtree)
cd tree
mafft --localpair --maxiterate 1000 tree_aa.fa > exendin_aa.aln.fa
$IQTREE -s exendin_aa.aln.fa -m MFP -B 1000 -alrt 1000 -T AUTO ${OG:+-o $OG} --prefix exendin_aa -redo
if [[ -n "$OG" ]]; then cp exendin_aa.treefile exendin_aa.rooted.nwk
else python3 ${SCRIPTDIR}/lib_midpoint_root.py exendin_aa.treefile > exendin_aa.rooted.nwk; fi
python3 ${SCRIPTDIR}/4h_draw_tree.py --tree exendin_aa.rooted.nwk --groups groups.tsv \
    --sources name_sources.tsv ${OG:+--outgroup $OG} \
    --title "Curated exendin-family genes (reciprocal BLAST)" --png exendin_aa.png > exendin_aa.svg
cd ..
```
- **Tips.** Every curated locus in the SLC4A8–AQP5 region (`--placements`): remnants,
  flagged models and pseudogenes included. Stops are masked (`*` → `X`; stop codons → `NNN`
  in the CDS). Left out only: loci with no model, `elsewhere` loci, and `--exclude`.
  `tree/tips.tsv` lists every locus and why. `--exclude-disabled --exclude-flagged` give the
  strict set (e.g. for selection analyses). Single-exon remnants (20–60 aa) sit on long,
  weakly supported branches.
- **Published precursors** `Published_EX2/3/4` (P04204, P20394, P26349); exendin-1 is
  published only as a mature peptide and is not added. A *Heloderma* model identical to one
  sits with it as an identical sequence (IQ-TREE keeps both).
- **Root** on chicken GCGL, by construction: it orients the tree, it does not test where
  the duplications happened (that needs GCGL from several vertebrates; old
  `2_exendin_origin.md` step 5).
- **Drawing.** Colour = inherited name group; `[R]/[N]/[T]/[P]/[O]` after each tip = where
  that name came from (reference RefSeq proteome / another species' NCBI annotation / TOGA2 /
  published / outgroup); hover for the protein. Supports SH-aLRT/UFBoot (trust ≥ 80 / 95);
  bold = *Heloderma*.

### Codon tree (optional)
The protein alignment back-translated, so gaps never split a codon. Outgroup and published
tips have no genomic CDS: drop them and midpoint-root.
```
cd tree
seqkit grep -v -r -p '^Published_' -f <(tr ',' '\n' < outgroup.txt) exendin_aa.aln.fa > exendin_aa.aln.noextra.fa
pal2nal.pl exendin_aa.aln.noextra.fa tree_cds.fa -output fasta -nogap > exendin_codon.aln.fa 2> pal2nal.log
cat pal2nal.log                                     # a mismatch names the tip to fix
$IQTREE -s exendin_codon.aln.fa -st CODON -m MFP -B 1000 -alrt 1000 -T AUTO --prefix exendin_codon -redo
python3 ${SCRIPTDIR}/lib_midpoint_root.py exendin_codon.treefile > exendin_codon.rooted.nwk
python3 ${SCRIPTDIR}/4h_draw_tree.py --tree exendin_codon.rooted.nwk --groups groups.tsv \
    --sources name_sources.tsv --title "Curated exendin-family genes, codon model, midpoint-rooted" \
    --png exendin_codon.png > exendin_codon.svg
cd ..
```
`-nogap` drops gapped columns; remove it to keep them as missing data. Trust deep nodes from
the protein tree and relationships among recent copies (e.g. within *Heloderma*) from the
codon tree.
