#!/bin/bash
set -euo pipefail

# Script to convert sra to fastq
# Usage: bash convert_sra_to_fastq.sh  [output_dir]

# if no arguments provided, print usage and exit
# if [ $# -lt 1 ]; then
#     echo "Usage: bash convert_sra_to_fastq.sh <sra_file> [output_dir]"
#     echo "**Example:"
#     echo "  bash convert_sra_to_fastq.sh output_fastq/"
#     exit 1
# fi

# SAMPLES=("SRR24952452" "SRR24952453" "SRR24952454" "SRR25533465")
# SAMPLES=("SRR24952453_Post" "SRR24952454")
SAMPLES=("SRR25533465") # Visium sample, only one run
# Output dir optional, defaults to data/
OUTPUT_DIR="${1:-data/}"
FASTQ_DIR="$OUTPUT_DIR/fastq"
# FASTQ_DIR="$OUTPUT_DIR/fastq2"

# Create output dir if it doesn't exist
[ -d "$OUTPUT_DIR" ] || mkdir -p "$OUTPUT_DIR"
[ -d "$FASTQ_DIR" ] || mkdir -p "$FASTQ_DIR"

echo ""
echo "Converting .sra files to FASTQ"

CONVERTED_COUNT=0
SKIPPED_COUNT=0

# Loop through each sample and convert .sra to FASTQ
# "${SAMPLES[@]}" is the array of sample IDs defined at the top, "$SAMPLES" points to the first element of the array
for SAMPLE in "${SAMPLES[@]}"; do
    echo "========== Processing sample: $SAMPLE =========="
    SRA_FILE="$OUTPUT_DIR/${SAMPLE}/${SAMPLE}.sra"
    FASTQ1_FILE="$FASTQ_DIR/${SAMPLE}_1.fastq"
    FASTQ2_FILE="$FASTQ_DIR/${SAMPLE}_2.fastq"
    
    # echo "DEBUG: SRA_FILE=$SRA_FILE"
    # echo "DEBUG: FASTQ1_FILE=$FASTQ1_FILE"
    # echo "DEBUG: FASTQ2_FILE=$FASTQ2_FILE"

    # Check if sra file exists
    if [ ! -f "$SRA_FILE" ]; then
        echo "Error: SRA file not found: $SRA_FILE"
        continue
    fi

    # Check if FASTQ already exists
    if [ -f "$FASTQ1_FILE" ] && [ -f "$FASTQ2_FILE" ]; then
        echo "Skipping $SAMPLE (already converted)"
        ((SKIPPED_COUNT++)) || true
    else
        echo "Converting $SRA_FILE"
        if fasterq-dump "$SRA_FILE" \
            --split-files \
            --threads 4 \
            --outdir "$FASTQ_DIR"; then
            echo "Successfully converted $SAMPLE"
            ((CONVERTED_COUNT++)) || true
        else
            echo "ERROR: Failed to convert $SAMPLE"
        fi

    if 
    fi
    echo "========== Done with $SAMPLE =========="
done
    
echo ""
echo "Conversion complete: $CONVERTED_COUNT new files converted, $SKIPPED_COUNT files skipped"

# echo "Done."

# Get run ID from filename
# FILENAME=$(basename "$FASTQ_FILE" .fastq)
# FASTA_FILE="$OUTPUT_DIR/${FILENAME}.fasta"

# echo "Converting $FASTQ_FILE to FASTA..."
# echo "Output: $FASTA_FILE"

# # Skip if FASTA already exists
# if [ -f "$FASTA_FILE" ]; then
#     echo "FASTA file already exists. Skipping..."
#     exit 0
# fi

# Convert FASTQ to FASTA by extracting header and sequence lines
# FASTQ format has 4 lines per record: header, sequence, optional header, quality
# We want to convert it to FASTA format, which has 2 lines per record: header (starting with >) and sequence
# 
# awk 'BEGIN {FS = "\n" ; RS = "" } { print ">" substr($1, 2) "\n" $2 }' \
#     "$FASTQ_FILE" > "$FASTA_FILE"

# echo "Conversion complete!"
# ls -lh "$FASTA_FILE"

