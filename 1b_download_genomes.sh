#!/usr/bin/env bash
# Download genome, NCBI annotation (GFF3), proteins and sequence report for every row of
# species.tsv (from 1a_make_species_table.py), and index each FASTA.
#
#   bash 1b_download_genomes.sh species.tsv
#
# Layout (same as the old repo's 3a_scan_species.sh, so genomes it already fetched are reused):
#   ${GENOMEDIR}/<code>_<accession>/ncbi_dataset/data/<accession>/{*_genomic.fna,genomic.gff,protein.faa}
# A directory that already has the FASTA is not downloaded again; a missing GFF or protein
# file is fetched on its own and unzipped into the same directory. Assemblies NCBI has not
# annotated have no GFF/protein: step 2 draws them as "no NCBI annotation", step 4 still
# searches them.

set -u
GENOMEDIR=${GENOMEDIR:-/data/Wilson_Lab/data/Exendin_Evolution/genomes}
TABLE=${1:-species.tsv}
[[ -f $TABLE ]] || { echo "!! not found: $TABLE"; exit 1; }
command -v datasets >/dev/null || { echo "!! datasets not on PATH (mamba activate exendin_env)"; exit 1; }
mkdir -p "$GENOMEDIR"

fetch() {   # fetch ACC INCLUDE DEST
  local ACC=$1 INC=$2 DEST=$3 ZIP
  ZIP=$(mktemp -u "${GENOMEDIR}/.dl_XXXXXX").zip
  datasets download genome accession "$ACC" --include "$INC" --filename "$ZIP" >/dev/null 2>&1 \
    && unzip -o -q "$ZIP" -d "$DEST" && rm -f "$ZIP" && return 0
  rm -f "$ZIP"; return 1
}

# only the first three columns; empty fields become placeholders (read collapses empty tab fields)
grep -v '^#' "$TABLE" | awk -F'\t' 'NR>1{print $1"\t"($2==""?"-":$2)"\t"($3==""?"NA":$3)}' \
| while IFS=$'\t' read -r CODE SPECIES ACC; do
  [[ -z "${CODE:-}" ]] && continue
  if [[ -z "${ACC:-}" || "$ACC" == "NA" ]]; then echo "=== $CODE: no accession, skipped"; continue; fi
  DEST=${GENOMEDIR}/${CODE}_${ACC}
  D=${DEST}/ncbi_dataset/data/${ACC}
  echo "=== $CODE $SPECIES $ACC"
  if ! ls ${D}/*_genomic.fna >/dev/null 2>&1; then
    echo "    downloading genome, gff3, protein, seq-report ..."
    fetch "$ACC" genome,gff3,protein,seq-report "$DEST" || { echo "    !! download failed"; continue; }
  fi
  [[ -f ${D}/genomic.gff ]]  || fetch "$ACC" gff3 "$DEST"    || true
  [[ -f ${D}/protein.faa ]]  || fetch "$ACC" protein "$DEST" || true
  FA=$(ls ${D}/*_genomic.fna 2>/dev/null | head -1)
  [[ -n "$FA" ]] || { echo "    !! no FASTA under $D"; continue; }
  [[ -f ${FA}.fai ]] || samtools faidx "$FA"
  printf '    fasta %s\n    gff   %s\n    prot  %s\n' "$(basename "$FA")" \
    "$([[ -f ${D}/genomic.gff ]] && echo yes || echo 'none (not annotated by NCBI)')" \
    "$([[ -f ${D}/protein.faa ]] && echo yes || echo none)"
done

echo
echo "=== summary"
{ echo -e "code\taccession\tfasta\tgff\tprotein"
  grep -v '^#' "$TABLE" | awk -F'\t' 'NR>1 && $3!="NA" && $3!=""{print $1"\t"$3}' | while read -r C A; do
    D=${GENOMEDIR}/${C}_${A}/ncbi_dataset/data/${A}
    printf '%s\t%s\t%s\t%s\t%s\n' "$C" "$A" \
      "$(ls ${D}/*_genomic.fna >/dev/null 2>&1 && echo yes || echo NO)" \
      "$([[ -f ${D}/genomic.gff ]] && echo yes || echo no)" \
      "$([[ -f ${D}/protein.faa ]] && echo yes || echo no)"
  done
} | { command -v column >/dev/null && column -t || cat; }
