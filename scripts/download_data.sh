#!/bin/bash
set -euo pipefail

# This script downloads sequencing data from the SRA for specified project IDs, saves metadata, and converts .sra files to FASTQ format.
OUTDIR="${1:-data}"
PROJECT_IDS=("PRJNA974579" "PRJNA1002459")

# Create directories only if they don't exist
[ -d "$OUTDIR/metadata" ] || mkdir -p "$OUTDIR/metadata"
[ -d "$OUTDIR/fastq" ] || mkdir -p "$OUTDIR/fastq"
cd "$OUTDIR"

# Loop through each project ID and download metadata and sequencing data
for PROJECT_ID in "${PROJECT_IDS[@]}"; do
    echo "=========================================="
    echo "Processing project: $PROJECT_ID"
    echo "=========================================="
    
    echo "Saving RunInfo metadata for project: $PROJECT_ID"
    
    # Save sequencing metadata (always refresh to capture any updates)
    esearch -db sra -query "$PROJECT_ID" \
    | efetch -format runinfo > metadata/"${PROJECT_ID}_runinfo.csv"
    
    # Check if any runs are missing before downloading
    SRA_CACHE="${HOME}/.ncbi/public/sra"
    
    # Get list of run accessions from metadata
    # Skip the header line and extract the first column (Run), then get unique values
    RUN_IDS=$(tail -n +2 metadata/"${PROJECT_ID}_runinfo.csv" | cut -d',' -f1 | sort -u)
    
    # Check for missing .sra files. not sure if they are in the .ncbi cache or the current data directory
    MISSING_RUNS=()
    for RUN_ID in $RUN_IDS; do
        if [ ! -f "$SRA_CACHE/$RUN_ID/$RUN_ID.sra" ] && [ ! -f "/work/TALC/mdsc519_2026w/students/jamie/Dragon/data/$RUN_ID/$RUN_ID.sra" ]; then
            MISSING_RUNS+=("$RUN_ID")
        fi
    done
    
    # If there are no missing runs, skip download
    if [ ${#MISSING_RUNS[@]} -eq 0 ]; then
        echo "All runs for $PROJECT_ID already downloaded. Skipping..."
    else
        echo "Downloading ${#MISSING_RUNS[@]} missing sequencing runs for $PROJECT_ID"
        prefetch "$PROJECT_ID" --max-size 200GB
    fi
done