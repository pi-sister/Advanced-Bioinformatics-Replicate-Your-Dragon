#!/bin/bash
set -euo pipefail

REF="/work/TALC/mdsc519_2026w/students/jamie/Dragon/data/reference/PitayaGenomic_cellranger_ref"
FASTQ_DIR="/work/TALC/mdsc519_2026w/students/jamie/Dragon/data/fastq"
OUTDIR="/work/TALC/mdsc519_2026w/students/jamie/Dragon/data/output/cellranger"
# IMAGE=

mkdir -p "${OUTDIR}"
cd "${OUTDIR}"

# Reference genome parameters
declare -A group_ids
# group_ids["SRR24952452"]="Trypsin"
group_ids["SRR24952453"]="post"
group_ids["SRR24952454"]="ck"
# group_ids["SRR25533465"]="visium"

# Visium sample
VIS_SAMPLE_ID="SRR25533465"
VIS_AREA=""
VIS_IMAGE="/work/TALC/mdsc519_2026w/students/jamie/Dragon/data/slides/tissue_lowres_image.png"
VIS_SLIDE_ID=""

# Run parameters
THREADS=${1:-8}
MEM=${2:-32}



##########################################################
# Run Cell Ranger count for each sample (in parallel)
##########################################################
for sample in "${!group_ids[@]}"; do
    sample_name="${group_ids[$sample]}"
    echo "Processing sample: $sample ($sample_name)"
    
    if [ -f "${OUTDIR}/${sample_name}/outs/web_summary.html" ]; then
        echo "  Already completed, skipping..."
    else
        cellranger count \
            --id="${sample_name}" \
            --transcriptome="${REF}" \
            --fastqs="${FASTQ_DIR}" \
            --sample="${sample}" \
            --localcores="${THREADS}" \
            --localmem="${MEM}" \
            --create-bam=true &
    fi
done
wait
echo "All samples completed!"

###############################################
# Visium
################################################
echo "Processing Visium sample: $VIS_SAMPLE_ID"
if [ -f "${OUTDIR}/visium/outs/web_summary.html" ]; then
    echo "  Already completed, skipping..."
else
    if [ -z "$VIS_AREA" ]; then
        echo "  Warning: VIS_AREA or VIS_SLIDE_ID not set, running with --unknown-slide..."
        spaceranger count \
            --id="visium" \
            --transcriptome="${REF}" \
            --fastqs="${FASTQ_DIR}" \
            --sample="${VIS_SAMPLE_ID}" \
            --localcores="${THREADS}" \
            --localmem="${MEM}" \
            --image="${VIS_IMAGE}" \
            --unknown-slide="visium-1" \
            --create-bam=true
    else
        spaceranger count \
            --id="visium" \
            --transcriptome="${REF}" \
            --fastqs="${FASTQ_DIR}" \
            --sample="${VIS_SAMPLE_ID}" \
            --localcores="${THREADS}" \
            --localmem="${MEM}" \
            --chemistry=auto \
            --create-bam=true \
            --image="${VIS_IMAGE}" \
            --slide="${VIS_SLIDE_ID}" \
            --area="${VIS_AREA}"
    fi
fi