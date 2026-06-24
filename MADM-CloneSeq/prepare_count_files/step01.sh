#!/bin/bash
#SBATCH --ntasks=8
#SBATCH --job-name=step01.ID
#SBATCH --output=step01.out
#SBATCH --time=12:00:00
#SBATCH --mem=50G
#SBATCH --no-requeue

module load R/4.4.0

Rscript ./step_01_prepare_count_files.R
