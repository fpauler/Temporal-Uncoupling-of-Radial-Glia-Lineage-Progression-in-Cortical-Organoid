#!/bin/bash
#SBATCH --ntasks=10
#SBATCH --job-name=prep_pb
#SBATCH --output=prep_pb.out
#SBATCH --time=24:00:00
#SBATCH --mem=100G
#SBATCH --no-requeue

module load R/4.4.1
Rscript prepare_pb.R
