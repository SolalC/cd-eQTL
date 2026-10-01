#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=100G
#SBATCH --job-name=21_chenComparison
#SBATCH --time=04:00:00
#SBATCH --partition=general
#SBATCH --constraint=epyc3
#SBATCH --account=a_imb_cpdg
#SBATCH -o logs/%x.out
#SBATCH -e logs/%x.err

# Submit from this folder after 01 and 02: sbatch 21_chenComparison.sh
# (memory is for the GTEx variant lookup table; genotypes are not loaded)
cd "${SLURM_SUBMIT_DIR}"
module load r/4.4.0-combo-EPYC3-only
Rscript 21_chenComparison.R
