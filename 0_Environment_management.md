# Environment management
Adapted from `0_Environment_management.md` of
[Exendin_evolution](https://github.com/dannyjackson/Exendin_evolution), with two fixes: its
activation block named `hal_toga_env` (the environment is `exendin_env`), and `soure` was a
typo for `source`.

## Working directory (once)
```
PROJDIR=/data/Wilson_Lab/projects/GLP1_evolution/
mkdir -p ${PROJDIR}/scripts && cd ${PROJDIR}/scripts
git clone https://github.com/dannyjackson/Exendin_synteny

mkdir -p ${PROJDIR}/analyses/exendin_synteny
cd ${PROJDIR}/analyses/exendin_synteny
mkdir -p synteny_ncbi synteny_toga2 reciprocal_blast nomenclature
```

## Software (once)
The old `exendin_env` has everything except an SVG→PNG converter. Add it, or create the
environment from scratch:
```
source myconda
mamba install -n exendin_env -c conda-forge librsvg          # existing environment

mamba create -n exendin_env -c conda-forge -c bioconda --strict-channel-priority \
  python=3.11 pandas biopython ncbi-datasets-cli \
  samtools bedtools seqkit gffread \
  ucsc-twobittofa ucsc-fatotwobit ucsc-bigbedtobed ucsc-liftover \
  blast miniprot exonerate mafft macse trimal pal2nal iqtree paml hyphy librsvg
```
A minimal environment for this repo alone is `environment.yml`
(`mamba env create -f ${SCRIPTDIR}/environment.yml`).

| Purpose | Tools | Steps |
|---|---|---|
| choose and download NCBI assemblies | `datasets`, `unzip` | 1 |
| FASTA index, region extraction | samtools | 1, 4 |
| forward / reverse / missing-exon searches | blast (makeblastdb, tblastn, blastx) | 4 |
| alignment and trees | mafft, iqtree, pal2nal | 4, 5 |
| PNG copies of SVG figures | rsvg-convert (else inkscape / ImageMagick / cairosvg) | 2–5 |
| everything else | Python 3 standard library | all |

Only the species tree in `${CACTUSDIR}` is used; no HAL tools are needed.

## Set-up block
Run at the top of every step.
```
source myconda
mamba activate exendin_env

PROJDIR=/data/Wilson_Lab/projects/GLP1_evolution/
TOGADIR=/data/Wilson_Lab/projects/squamates/Arizona_squamates/misc/Toxicofera_TOGA2_annotations/
CACTUSDIR=/data/Wilson_Lab/projects/squamates/Arizona_squamates/misc/Toxicofera_cactus_alignment_2026
GENOMEDIR=/data/Wilson_Lab/data/Exendin_Evolution/genomes
SCRIPTDIR=${PROJDIR}/scripts/Exendin_synteny
SPTREE=${CACTUSDIR}/Anguimorpha_species_tree_iteration-01_5-23-26.nwk
WORK=${PROJDIR}/analyses/exendin_synteny
HELSUS_EVIANN_GFF=${PROJDIR}/datasets/genomes/Heloderma_suspectum/Heloderma_suspectum.final.reference.gff3

cd ${SCRIPTDIR} && git pull && cd ${WORK}
```
