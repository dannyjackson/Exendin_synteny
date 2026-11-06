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

# check the output
cd ${WORK}/reciprocal_blast
# which queries make the Heloderma loci
awk -F'\t' '$2=="HelSus"||$2=="HelCha"{n[$2"\t"$10]++} END{for(k in n) print n[k]"\t"k}' loci/forward_loci.tsv | sort -k1,1nr | head
# what kind of loci they are (expect almost all N_term_only / domain_only)
awk -F'\t' '$2=="HelSus"{split($9,s,";"); n[s[1]]++} END{for(k in n) print n[k], k}' loci/forward_loci.tsv
# length of every forward query: the precursors are 84-113 aa
awk -F'\t' 'NR>1 && ($6<70 || $6>150)' queries/query_table.tsv | cut -f1-6


# filter out faulty m-RNA annotation
BAD=HelCha__toga2__Ex3_4__1
mkdir -p queries_v2 forward_v2
seqkit grep -v -p $BAD queries/forward_queries.faa > queries_v2/forward_queries.faa
seqkit grep -v -p $BAD queries/reverse_db.faa      > queries_v2/reverse_db.faa
cp queries/boundary_queries.faa queries_v2/
for T in query_table reverse_table; do awk -F'\t' -v B=$BAD '$1!=B' queries/$T.tsv > queries_v2/$T.tsv; done
for F in forward/*.tblastn.tsv; do awk -F'\t' -v B=$BAD '$2!=B' $F > forward_v2/$(basename $F); done

python3 ${SCRIPTDIR}/4c_call_loci.py --species-table ../species.tsv --forward-dir forward_v2 \
    --queries queries_v2/forward_queries.faa --query-table queries_v2/query_table.tsv \
    --boundary-queries queries_v2/boundary_queries.faa --outdir loci_v2 2> loci_v2.log
awk -F'\t' 'NR>1{n[$2]++} END{for(s in n) print s, n[s]}' loci_v2/forward_loci.tsv | sort -k2,2nr | head

cd ${WORK}/reciprocal_blast
rm -rf forward loci queries reciprocal reciprocal.log
mv forward_v2 forward
mv loci_v2 loci
mv queries_v2 queries
mv loci_v2.log loci.log
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
# classify
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

### Add chicken outgroup

```
cd ${WORK}
printf 'GalGal\tGallus_gallus\tGCF_016699485.2\t\t\tpinned\toutgroup, not on the tree\n' >> species.tsv
bash ${SCRIPTDIR}/1b_download_genomes.sh species.tsv 2>&1 | tee download.log     # fetches only what is missing

cd synteny_ncbi
python3 ${SCRIPTDIR}/2a_ncbi_exendin_loci.py --species-table ../species.tsv --genomedir ${GENOMEDIR} --outdir ncbi_loci 2> ncbi_loci.log
grep GalGal ncbi_loci.log
python3 ${SCRIPTDIR}/2b_plot_ncbi_synteny.py --loci ncbi_loci/ncbi_loci.tsv --species-tree ${SPTREE} \
    --inset GalGal --inset-title "Chicken (not on the tree)" \
    --svg ncbi_synteny.svg --png ncbi_synteny.png --layout ncbi_synteny.tsv 2> ncbi_synteny.log

cd ${WORK}/reciprocal_blast
SPECIES_TABLE=../species.tsv GENOMEDIR=${GENOMEDIR} QUERYDIR=queries OUT=forward THREADS=8 \
    bash ${SCRIPTDIR}/4b_forward_tblastn.sh GalGal
```
#### Rerun 4b on chicken
```
cd ${WORK}/reciprocal_blast
ls ${GENOMEDIR}/GalGal_GCF_016699485.2/ncbi_dataset/data/GCF_016699485.2/   # needs *_genomic.fna (1b)

SPECIES_TABLE=../species.tsv GENOMEDIR=${GENOMEDIR} QUERYDIR=queries OUT=forward THREADS=8 \
    bash ${SCRIPTDIR}/4b_forward_tblastn.sh GalGal

rm -rf loci_GalGal
```
#### Rerun 4c on chicken
```
cd ${WORK}/reciprocal_blast
awk -F'\t' 'NR==1 || $1=="GalGal"' ../species.tsv > species_GalGal.tsv

python3 ${SCRIPTDIR}/4c_call_loci.py --species-table species_GalGal.tsv --forward-dir forward \
    --queries queries/forward_queries.faa --query-table queries/query_table.tsv \
    --boundary-queries queries/boundary_queries.faa --outdir loci_GalGal 2> loci_GalGal.log
cat loci_GalGal.log

for T in forward_loci forward_scores boundary_hits; do
  grep -q '^GalGal' loci/$T.tsv 2>/dev/null && echo "!! GalGal already in loci/$T.tsv" \
    || tail -n +2 loci_GalGal/$T.tsv >> loci/$T.tsv
done
mv loci_GalGal/scan/GalGal loci/scan/ && rm -rf loci_GalGal species_GalGal.tsv
```
#### Rerun 4d, 4e
```
# classify
python3 ${SCRIPTDIR}/4d_reciprocal_classify.py --species-table ../species.tsv --genomedir ${GENOMEDIR} \
    --loci-dir loci --queries-dir queries --threads 8 --outdir reciprocal 2> reciprocal.log
cat reciprocal.log
column -t reciprocal/summary.tsv
```
python3 ${SCRIPTDIR}/4f_plot_curated_synteny.py --species-table ../species.tsv \
    --loci reciprocal/reciprocal_loci.tsv --boundaries-used reciprocal/boundaries_used.tsv \
    --species-tree ${SPTREE} \
    --svg curated_synteny.svg --png curated_synteny.png --layout curated_synteny.tsv --inset GalGal

##### with medium-confidence loci too (labels end in '?')
python3 ${SCRIPTDIR}/4f_plot_curated_synteny.py --species-table ../species.tsv \
    --loci reciprocal/reciprocal_loci.tsv --boundaries-used reciprocal/boundaries_used.tsv \
    --species-tree ${SPTREE} --confidence high,medium \
    --svg curated_synteny.all.svg --png curated_synteny.all.png --inset GalGal
  
## 7. Gene tree of the curated set
```
mkdir -p tree
python3 ${SCRIPTDIR}/4g_tree_inputs.py --curated reciprocal/curated_loci.tsv --models-dir models \
    --outgroup-species GalGal --outdir tree 2> tree/tree_inputs.log
cat tree/tree_inputs.log
OG=$(cat tree/outgroup.txt)
IQTREE=$(command -v iqtree3 || command -v iqtree2 || command -v iqtree)

cd tree
seqkit fx2tab -n -l tree_aa.fa | sort -k2,2n | sed -n '1p;$p'      # shortest and longest tip
mafft --localpair --maxiterate 1000 tree_aa.fa > exendin_aa.aln.fa
$IQTREE -s exendin_aa.aln.fa -m MFP -B 1000 -alrt 1000 -T AUTO ${OG:+-o $OG} --prefix exendin_aa
grep "Best-fit model according to" exendin_aa.iqtree

if [[ -n "$OG" ]]; then cp exendin_aa.treefile exendin_aa.rooted.nwk
else python3 ${SCRIPTDIR}/lib_midpoint_root.py exendin_aa.treefile > exendin_aa.rooted.nwk; fi

python3 ${SCRIPTDIR}/4h_draw_tree.py --tree exendin_aa.rooted.nwk --groups groups.tsv ${OG:+--outgroup $OG} \
    --title "Curated exendin-family genes (reciprocal BLAST), ${OG:+rooted on $OG}${OG:-midpoint-rooted}" \
    --png exendin_aa.png > exendin_aa.svg
cd ..
```
- **Tips.** `tree/tips.tsv` lists every curated locus, whether it is in the tree, and why not.
  By default the tree takes only 4e models with status `ok` (hormone-domain motif, 60–150
  aa) and no disabled copies. Override with `--keep-flagged`, `--include-disabled`, or
  `--exclude NAME,...`, e.g. after looking at the alignment.
- **Root.** Chicken GCGL (`--outgroup-species GalGal`) when it passed 4d/4e; otherwise
  midpoint. With one outgroup sequence the root is placed on its branch by construction: it
  orients the tree, it does not test where the duplications happened. To test that, add GCGL
  orthologs from several vertebrates, as described in the old `2_exendin_origin.md` step 5.
- **Reading it.** Tip colours are the name groups of the synteny figures, so a tip whose
  colour disagrees with its clade is a naming conflict to trace (see `rev_best` in
  `reciprocal_loci.tsv`). Node labels are SH-aLRT/UFBoot; trust a node at ≥80 / ≥95. Bold =
  *Heloderma*.
- **Optional codon tree.** Same steps as the old `4_toga2_exendin_tree.md` step 3: MACSE on
  `tree/tree_cds.fa`, then `$IQTREE -st CODON` on the exported codon alignment.