#!/usr/bin/env bash
# Forward search of step 4: tblastn of every exendin-family query and the two boundary
# proteins (4a) against ONE whole genome.
#
#   bash 4b_forward_tblastn.sh CODE
#
# Run once per species (see 4_reciprocal_blast.md for a loop and a swarm file).
# Environment (defaults in brackets):
#   SPECIES_TABLE  species.tsv from 1a                       [species.tsv]
#   GENOMEDIR      genomes from 1b                           [/data/Wilson_Lab/data/Exendin_Evolution/genomes]
#   BLASTDB_DIR    where nucleotide BLAST databases are kept [${GENOMEDIR}/blastdb]
#   QUERYDIR       4a output                                 [queries]
#   OUT            output directory                          [forward]
#   THREADS        [8];  EVALUE [1e-3]
# Writes ${OUT}/${CODE}.tblastn.tsv with columns
#   saccver qseqid qstart qend sstart send pident length evalue bitscore sframe qseq sseq
# (2f's -outfmt 6 fields with the contig added in front; 4c splits it per contig for 2f).
# -seg no, as in 2f/2e: low-complexity masking hides the short single-exon hits of these
# small precursors. The database is built once per genome and reused.

set -euo pipefail
CODE=${1:?usage: bash 4b_forward_tblastn.sh CODE}
SPECIES_TABLE=${SPECIES_TABLE:-species.tsv}
GENOMEDIR=${GENOMEDIR:-/data/Wilson_Lab/data/Exendin_Evolution/genomes}
BLASTDB_DIR=${BLASTDB_DIR:-${GENOMEDIR}/blastdb}
QUERYDIR=${QUERYDIR:-queries}
OUT=${OUT:-forward}
THREADS=${THREADS:-8}
EVALUE=${EVALUE:-1e-3}

ACC=$(awk -F'\t' -v C="$CODE" '$1==C{print $3}' "$SPECIES_TABLE")
[[ -n "$ACC" && "$ACC" != "NA" ]] || { echo "!! $CODE: no accession in $SPECIES_TABLE"; exit 1; }
FA=$(ls ${GENOMEDIR}/${CODE}_${ACC}/ncbi_dataset/data/${ACC}/*_genomic.fna 2>/dev/null | head -1)
[[ -n "$FA" ]] || { echo "!! $CODE: no FASTA under ${GENOMEDIR}/${CODE}_${ACC} (run 1b)"; exit 1; }
for F in forward_queries.faa boundary_queries.faa; do
  [[ -s ${QUERYDIR}/$F ]] || { echo "!! missing ${QUERYDIR}/$F (run 4a)"; exit 1; }
done
mkdir -p "$BLASTDB_DIR" "$OUT"

DB=${BLASTDB_DIR}/${CODE}_${ACC}
if [[ ! -f ${DB}.ndb && ! -f ${DB}.nal && ! -f ${DB}.nin && ! -f ${DB}.00.nin ]]; then
  echo "=== $CODE: makeblastdb $(basename "$FA")"
  makeblastdb -in "$FA" -dbtype nucl -parse_seqids -out "$DB" -title "${CODE}_${ACC}" > "${DB}.makeblastdb.log"
fi

cat ${QUERYDIR}/forward_queries.faa ${QUERYDIR}/boundary_queries.faa > ${OUT}/${CODE}.queries.faa
echo "=== $CODE: tblastn $(grep -c '>' ${OUT}/${CODE}.queries.faa) queries, ${THREADS} threads"
tblastn -query ${OUT}/${CODE}.queries.faa -db "$DB" -evalue "$EVALUE" -seg no \
  -max_target_seqs 1000 -num_threads "$THREADS" \
  -outfmt "6 saccver qseqid qstart qend sstart send pident length evalue bitscore sframe qseq sseq" \
  > ${OUT}/${CODE}.tblastn.tsv.tmp
mv ${OUT}/${CODE}.tblastn.tsv.tmp ${OUT}/${CODE}.tblastn.tsv
rm -f ${OUT}/${CODE}.queries.faa
echo "=== $CODE: $(wc -l < ${OUT}/${CODE}.tblastn.tsv) HSPs on $(cut -f1 ${OUT}/${CODE}.tblastn.tsv | sort -u | wc -l) contigs"
