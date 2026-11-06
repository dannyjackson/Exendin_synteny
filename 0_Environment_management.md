# Environment management
Adapted from `0_Environment_management.md` of
[Exendin_evolution](https://github.com/dannyjackson/Exendin_evolution). Two fixes against
that file: the activation block named `hal_toga_env` (the environment is `exendin_env`), and
`soure` was a typo for `source`.

## Make working directory
```
PROJDIR=/data/Wilson_Lab/projects/GLP1_evolution/
TOGADIR=/data/Wilson_Lab/projects/squamates/Arizona_squamates/misc/Toxicofera_TOGA2_annotations/
CACTUSDIR=/data/Wilson_Lab/projects/squamates/Arizona_squamates/misc/Toxicofera_cactus_alignment_2026
GENOMEDIR=/data/Wilson_Lab/data/Exendin_Evolution/genomes

mkdir -p ${PROJDIR}/scripts
cd ${PROJDIR}/scripts
git clone https://github.com/dannyjackson/Exendin_synteny

SCRIPTDIR=${PROJDIR}/scripts/Exendin_synteny

mkdir -p ${PROJDIR}/analyses/exendin_synteny
cd ${PROJDIR}/analyses/exendin_synteny
mkdir -p synteny_ncbi synteny_toga2 reciprocal_blast
```

## Software
The old `exendin_env` already has everything this pipeline needs except an SVG-to-PNG
converter. If it exists, add one and reuse it:
```
source myconda
mamba install -n exendin_env -c conda-forge librsvg      # provides rsvg-convert
```
Otherwise create it from scratch (same as the old repo, plus `librsvg`):
```
source myconda
mamba create -n exendin_env -c conda-forge -c bioconda --strict-channel-priority \
  python=3.11 pandas biopython \
  ncbi-datasets-cli \
  samtools bedtools seqkit gffread \
  ucsc-twobittofa ucsc-fatotwobit ucsc-bigbedtobed ucsc-liftover \
  blast miniprot exonerate \
  mafft macse trimal pal2nal \
  iqtree paml hyphy \
  librsvg

mamba activate exendin_env
mamba env export --no-builds > exendin_env.yml
```
A minimal environment for this repo alone is in `environment.yml`
(`mamba env create -f ${SCRIPTDIR}/environment.yml`).

| Purpose | Tools | Steps |
|---|---|---|
| Pick and download NCBI assemblies | `datasets` (ncbi-datasets-cli), `unzip` | 1 |
| FASTA indexing and region extraction | samtools | 1, 4 |
| Forward / reverse searches | blast (makeblastdb, tblastn, blastx) | 4 |
| PNG copies of the SVG figures | rsvg-convert (librsvg), else inkscape / ImageMagick / cairosvg | 2, 3, 4 |
| Everything else | Python 3 standard library only | all |

Cactus/HAL tools are not needed: only the species tree in `${CACTUSDIR}` is used.

## Set up environment
Run this at the top of every markdown file
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

cd ${SCRIPTDIR}
git pull
cd ${WORK}
```
