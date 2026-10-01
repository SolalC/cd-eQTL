#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=100G
#SBATCH --job-name=14_multiTissue
#SBATCH --time=24:00:00
#SBATCH --partition=general
#SBATCH --constraint=epyc3
#SBATCH --account=a_imb_cpdg
#SBATCH -o logs/%x.out
#SBATCH -e logs/%x.err

cd "${SLURM_SUBMIT_DIR}"
module load r/4.4.0-combo-EPYC3-only
Rscript 14_multiTissueRegression.R
