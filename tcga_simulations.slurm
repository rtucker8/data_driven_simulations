#!/bin/bash

#SBATCH --job-name=data_based_simulations
#SBATCH --mail-user=rbtucker@umich.edu
#SBATCH --mail-type=BEGIN,END
#SBATCH --cpus-per-task=4
#SBATCH --mem=10GB
#SBATCH --time=02:00:00
#SBATCH --account=wdem2
#SBATCH --partition=standard
#SBATCH --output=logs/%x-%A_%a.log
#SBATCH --array=1-192

cd /home/rbtucker/tcga_simulations
module load R

R run_cluster_sum.R $SLURM_ARRAY_TASK_ID