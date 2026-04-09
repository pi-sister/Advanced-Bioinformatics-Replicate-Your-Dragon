#!/bin/bash
#SBATCH --job-name=data_grab
#SBATCH --cpus-per-task=8
#SBATCH --output=slurm_outputs/download_%j.out
#SBATCH --error=slurm_outputs/download_%j.err
#SBATCH --mem=16G
#SBATCH --time=12:00:00

source /home/jamie.macdonald2/miniforge3/bin/activate /work/TALC/mdsc519_2026w/students/jamie/Dragon/Env_download
bash /work/TALC/mdsc519_2026w/students/jamie/Dragon/scripts/download_data.sh data
bash /work/TALC/mdsc519_2026w/students/jamie/Dragon/scripts/convert_sra_to_fastq.sh
