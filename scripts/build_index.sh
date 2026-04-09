#!bin/bash
set -euo pipefail

# Use args
# FASTA_GENOME=$1
# ANNOTATION_GTF=$2
# OUTPUT_DIR=$3
# FASTQ_DIR=$4
# GENOME_NAME=${5:-hagan2021.fa}
NTHREADS=${1:-8}
MEMORY=${2:-32}

# GENOME_NAME="Hu_PGMD_10x_ref"
GENOME_NAME="PitayaGenomic"
GENOME_GTF="/work/TALC/mdsc519_2026w/students/jamie/Dragon/data/reference/PitayaGenomic.gtf"
GENOME_FASTA="/work/TALC/mdsc519_2026w/students/jamie/Dragon/data/reference/PitayaGenomic.fa"
FASTQ_DIR="/work/TALC/mdsc519_2026w/students/jamie/Dragon/data/fastq"
REF_OUTPUT_DIR="/work/TALC/mdsc519_2026w/students/jamie/Dragon/data/reference/PitayaGenomic_cellranger_ref"

# Sample IDs to process
declare -A group_ids
group_ids["SRR24952453"]="Post"
group_ids["SRR24952454"]="CK"

# mkdir -p "$REF_OUTPUT_DIR"

# CELLRANGER to generate the reference genome index
# cellranger version is 8.0.1, which is compatible with reference genome version 2.0.0

if [[ -d "$REF_OUTPUT_DIR" ]]; then
    echo "Reference genome already exists at $REF_OUTPUT_DIR, skipping mkref..."
else
    echo "Building Cell Ranger reference genome index for $GENOME_NAME..."
    mkdir "$REF_OUTPUT_DIR"  # Clean up any stale files
    cellranger mkref \
        --genome=$GENOME_NAME \
        --fasta=$GENOME_FASTA \
        --genes=$GENOME_GTF \
        --nthreads=$NTHREADS \
        --mempercore=$MEMORY \
        --ref-version=2.0.0 \
        --output-dir=$REF_OUTPUT_DIR
fi
    
    
##########################################################
# Rename and gzip FASTQ files to cellranger format
##########################################################
echo "Preparing FASTQ files for cellranger..."
for sample in "${!group_ids[@]}"; do
    if [[ -f "${FASTQ_DIR}/${sample}_1.fastq" && -f "${FASTQ_DIR}/${sample}_2.fastq" && ! -f "${FASTQ_DIR}/${sample}_S1_L001_R1_001.fastq.gz" && ! -f "${FASTQ_DIR}/${sample}_S1_L001_R2_001.fastq.gz" ]]; then
        echo "  Renaming and compressing ${sample}..."
        pigz -p "${NTHREADS}" -c "${FASTQ_DIR}/${sample}_1.fastq" > "${FASTQ_DIR}/${sample}_S1_L001_R1_001.fastq.gz"
        pigz -p "${NTHREADS}" -c "${FASTQ_DIR}/${sample}_2.fastq" > "${FASTQ_DIR}/${sample}_S1_L001_R2_001.fastq.gz"
    else
        echo "  FASTQ files for ${sample} are already prepared, skipping..."
    fi
done
echo "  FASTQ files are ready."

