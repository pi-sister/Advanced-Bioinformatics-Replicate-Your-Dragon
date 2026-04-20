#!/bin/bash

#SBATCH --job-name=sc_drag
#SBATCH --cpus-per-task=4
#SBATCH --output=slurm_outputs/sc_%j.out
#SBATCH --error=slurm_outputs/sc_%j.err
#SBATCH --mem=32G
#SBATCH --time=3:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jamie.macdonald2@ucalgary.ca


source /home/jamie.macdonald2/miniforge3/bin/activate /work/TALC/mdsc519_2026w/students/jamie/Dragon/Env_files

# python /work/TALC/mdsc519_2026w/students/jamie/Dragon/scripts/scrna-seq.py

python /work/TALC/mdsc519_2026w/students/jamie/Dragon/scripts/scrna-seq_by_condition.py