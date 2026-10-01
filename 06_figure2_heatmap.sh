#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=100G
#SBATCH --job-name=06_figure2
#SBATCH --time=12:00:00
#SBATCH --partition=general
#SBATCH --constraint=epyc3
#SBATCH --account=a_imb_cpdg
#SBATCH -o logs/%x.out
#SBATCH -e logs/%x.err

# 100G: the backfill loads the full plink genotype matrix
cd "${SLURM_SUBMIT_DIR}"
module load r/4.4.0-combo-EPYC3-only
Rscript 06_figure2_heatmap.R
