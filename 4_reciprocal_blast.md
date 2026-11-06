# Reciprocal BLAST: a curated set of exendin orthologs and paralogs
**Why.** Steps 2 and 3 show what the annotations say. The array is systematically
under-annotated, though: in HelSus four of six copies had no gene model, and HemCap's GFF
shows one exendin where a scan found four loci. NCBI and TOGA2 also name the same gene
differently. This step asks the genomes directly, with every exendin from steps 2–3 as a
query, and keeps a locus only if it searches back to an exendin-family protein.

**Method**
1. **Forward.** tblastn of all exendin-family proteins (NCBI, TOGA2 and the published
   precursors) plus the HelSus SLC4A8 and AQP5 proteins against every whole genome.
2. **Loci.** Hits are clustered into loci by the old repo's `2f_scan_interval.py`, copied
   unchanged as `lib_scan_interval.py`. Each locus is classified `candidate_gene` (exon 1 and
   hormone domain both found), `N_term_only` or `domain_only`. Premature stops and the first
   residue of the hormone domain are checked.
3. **Reverse.** blastx of each locus against a whole reference proteome (HelSus RefSeq) plus
   all forward queries. The reference proteome supplies the competitors (proglucagon, GIP,
   VIP, PACAP, secretin, GHRH): a locus whose best reverse hit is one of them is rejected.
   Query proteins from the locus's own species are ignored, so a locus is never named just by
   its own annotation.
4. **Orthologs and paralogs.**
   - **Ortholog of group G:** the best reverse hit is a G protein, and this locus is the
     genome's best forward hit for the G queries (a reciprocal best hit).
   - **Paralog named after G:** the best reverse hit is a G protein, but another locus is the
     best forward hit for G.
   - **High confidence:** reverse margin ≥ 1.2 over the best non-family protein, forward
     e-value ≤ 1e-5, and a `candidate_gene`. Everything else that passes the reverse test is
     medium.
5. **Plot** the high-confidence set in the step-2/3 layout and colours, with SLC4A8/AQP5
   located by tblastn, so unannotated genomes get boundary genes too.

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
    --outdir queries 2> queries.log
cat queries.log
cut -f1-6 queries/query_table.tsv | column -t
```
- `4_reference_peptides.fa` is `2_reference_peptides.fa` from the old repo, copied unchanged.
  Only its `*_precursor` records are used: exendin-2/-3/-4 as queries, and proglucagon as a
  reverse competitor.
- **Reverse proteome line.** Check how many proteins in the reference proteome were named
  exendin / GCGL / proglucagon. If HelSus calls ANC something the grouping does not know,
  ANC loci will be rejected. Fix `name_group()`, or rename it in a copy of the proteome.
- **Duplicate warnings.** "identical sequence named X and Y" means NCBI and TOGA2 call the
  same protein different things. The NCBI (or published) name is kept. Read these warnings:
  they are the nomenclature conflicts.
- **More competitors.** Add another proteome (e.g. `--reverse-proteome AnoCar=...`) to
  strengthen the reverse test in distant taxa.

## 2. Forward tblastn against every genome
One job per genome. Interactively:
```
export SPECIES_TABLE=../species.tsv GENOMEDIR QUERYDIR=queries OUT=forward THREADS=8
for C in $(awk -F'\t' 'NR>1 && $3!="NA"{print $1}' ../species.tsv); do bash ${SCRIPTDIR}/4b_forward_tblastn.sh $C; done
```
Or as a Biowulf swarm (makeblastdb plus tblastn is about an hour per squamate genome):
```
awk -F'\t' -v S=${SCRIPTDIR} -v W=$PWD -v G=${GENOMEDIR} 'NR>1 && $3!="NA"{
  printf "source myconda; mamba activate exendin_env; cd %s; SPECIES_TABLE=../species.tsv GENOMEDIR=%s QUERYDIR=queries OUT=forward THREADS=8 bash %s/4b_forward_tblastn.sh %s\n", W, G, S, $1}' \
  ../species.tsv > forward.swarm
swarm -f forward.swarm -g 16 -t 8 --time 8:00:00 --logdir swarm_logs
```
Databases go to `${GENOMEDIR}/blastdb` and are reused. Change this with `BLASTDB_DIR`.

## 3. Call loci
```
python3 ${SCRIPTDIR}/4c_call_loci.py --species-table ../species.tsv --forward-dir forward \
    --queries queries/forward_queries.faa --query-table queries/query_table.tsv \
    --boundary-queries queries/boundary_queries.faa --outdir loci 2> loci.log
cat loci.log
cut -f1-9 loci/forward_loci.tsv | column -t | less -S
cut -f1-8 loci/boundary_hits.tsv | column -t
```
- **Positive control.** HelSus must recover the loci the old copy-number survey found: ANC,
  the Frag1 copy, exendin-2, exendin-4 and the two distal copies. If any are missing, lower
  `--min-identity` before trusting an absence elsewhere.
- **Expected extras.** The forward set also finds proglucagon and its relatives. That is
  expected; the reverse step removes them.
- **Boundary hits.** In `boundary_hits.tsv`, the best SLC4A8 / AQP5 cluster should have
  `query_cov` near 1. Paralogs (SLC4A10, other aquaporins) follow at lower `frac_of_top`.

## 4. Reverse blastx and classification
```
python3 ${SCRIPTDIR}/4d_reciprocal_classify.py --species-table ../species.tsv --genomedir ${GENOMEDIR} \
    --loci-dir loci --queries-dir queries --threads 8 --outdir reciprocal 2> reciprocal.log
cat reciprocal.log
column -t reciprocal/summary.tsv
cut -f1,3,8-12,23,25,27,28 reciprocal/reciprocal_loci.tsv | column -t | less -S
```
- **Classes.**
  - **`rejected`:** the reverse hit names what it is (`rev_best_desc`); expect proglucagon,
    GIP, VIP etc.
  - **`paralog`:** `reason` names the locus that is the ortholog for the same group.
  - **`medium`:** `reason` says which test failed (margin, e-value, or a remnant locus).
- **Changing the cut-offs.** To try others without re-running blastx, re-run with
  `--reuse-blastx --min-margin 1.5` (or `--max-evalue`).
- **`curated_loci.tsv`** is the high-confidence set; the step-4 dataset is this file.
- **Disabled copies.** `disabled = yes` marks a premature stop, or a hormone domain not
  starting His/Tyr. These are kept as paralogs, as presence evidence, but should not go into
  a selection analysis.

## 5. Protein and CDS models of the curated loci
```
python3 ${SCRIPTDIR}/4e_build_models.py --species-table ../species.tsv --genomedir ${GENOMEDIR} \
    --curated reciprocal/curated_loci.tsv --queries queries/forward_queries.faa \
    --outdir models 2> models.log
column -t models/models.tsv
```
- **Method.** Each locus is rebuilt by the old repo's `2e_rescue_locus.py`, copied unchanged
  as `lib_rescue_locus.py`. Its best forward query of the same group is spliced over the
  locus, and the longest ORF is reported.
- **What to check.** Anything with status other than `ok` (no hormone motif, or a length far
  from 84–113 aa) needs a look at its `<name>.exons.tsv`, exactly as in the old step-2 QC.
- **Output.** `models/curated_aa.fa` and `models/curated_cds.fa` are ready for the old repo's
  alignment and tree steps.

## 6. Plot
```
python3 ${SCRIPTDIR}/4f_plot_curated_synteny.py --species-table ../species.tsv \
    --loci reciprocal/reciprocal_loci.tsv --boundaries-used reciprocal/boundaries_used.tsv \
    --species-tree ${SPTREE} \
    --svg curated_synteny.svg --png curated_synteny.png --layout curated_synteny.tsv

# with medium-confidence loci too (labels end in '?')
python3 ${SCRIPTDIR}/4f_plot_curated_synteny.py --species-table ../species.tsv \
    --loci reciprocal/reciprocal_loci.tsv --boundaries-used reciprocal/boundaries_used.tsv \
    --species-tree ${SPTREE} --confidence high,medium \
    --svg curated_synteny.all.svg --png curated_synteny.all.png
```
- **Drawing.** Colour = name group of the best reverse hit, with the same colours as steps 2
  and 3. Filled = ortholog, hollow = paralog, dashed = disabled copy. Labels are the names
  in `curated_loci.tsv` (`Ex4`, `Ex2p1` …).
- **Reading the three figures together.**
  - A locus filled here but missing from step 2 is a gene NCBI did not annotate.
  - One here but missing from step 3 is a copy TOGA2 did not project as an ortholog.
  - A colour that differs between figures is a naming conflict between the annotation and
    the reciprocal search; `4a`'s duplicate warnings list them.
