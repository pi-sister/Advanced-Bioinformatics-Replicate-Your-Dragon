#!/bin/bash
#SBATCH --job-name=run_dragon
#SBATCH --partition=bigmem
#SBATCH --cpus-per-task=32
#SBATCH --output=slurm_outputs/i-p_%j.out
#SBATCH --error=slurm_outputs/i-p_%j.err
#SBATCH --mem=256G
#SBATCH --time=12:00:00

source /home/jamie.macdonald2/miniforge3/bin/activate /work/TALC/mdsc519_2026w/students/jamie/Dragon/Env_files

# Build the Cell Ranger reference genome index
module load cellranger/v8.0.1
bash /work/TALC/mdsc519_2026w/students/jamie/Dragon/scripts/build_index.sh 32 128

# Process the sequencing data with Cell Ranger count (parallel)
bash /work/TALC/mdsc519_2026w/students/jamie/Dragon/scripts/process_data.sh 16 128

# scRNA-seq data pipeline
