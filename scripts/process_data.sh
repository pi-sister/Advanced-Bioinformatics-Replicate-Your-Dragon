#!/bin/bash
set -euo pipefail

REF = "/work/TALC/mdsc519_2026w/students/jamie/Dragon/data/reference/PitayaGenomic_cellranger_ref"
FASTQ_DIR="/work/TALC/mdsc519_2026w/students/jamie/Dragon/data/fastq"
OUTDIR="/work/TALC/mdsc519_2026w/students/jamie/Dragon/data/output/cellranger"
# IMAGE=

mkdir -p "${OUTDIR}"
cd "${OUTDIR}"

# Reference genome parameters
declare -A group_ids
# group_ids["SRR24952452"]="Trypsin"
group_ids["SRR24952453"]="Post"
group_ids["SRR24952454"]="CK"

# Visium sample
VIS_SAMPLE_ID="SRR25533465"

# Run parameters
THREADS=${1:-8}
MEM=${2:-32}



##########################################################
# Run Cell Ranger count for each sample
##########################################################
for sample in "${!group_ids[@]}"; do
    sample_name="${group_ids[$sample]}"
    echo "Processing sample: $sample ($sample_name)"
    cellranger count \
        --id="${sample_name}" \
        --transcriptome="${REF}" \
        --fastqs="${FASTQ_DIR}" \
        --sample="${sample}" \
        --localcores="${THREADS}" \
        --localmem="${MEM}" \
        --create-bam=true
done

###############################################
# Visium
################################################
# echo "Processing Visium sample: $VIS_SAMPLE_ID"
# cellranger count \
#     --id="Visium_${VIS_SAMPLE_ID}" \
#     --transcriptome="${REF}" \