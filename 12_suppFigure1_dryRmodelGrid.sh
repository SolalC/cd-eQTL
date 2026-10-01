#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --job-name=12_suppFig1
#SBATCH --time=00:20:00
#SBATCH --partition=general
#SBATCH --constraint=epyc3
#SBATCH --account=a_imb_cpdg
#SBATCH -o logs/%x.out
#SBATCH -e logs/%x.err

# r/4.4.2: patchwork >= 1.2 is needed with ggplot2 3.5 (the 4.4.0-combo module
# has patchwork 1.1.3, which fails in ggsave)
cd "${SLURM_SUBMIT_DIR}"
module load r/4.4.2
Rscript 12_suppFigure1_dryRmodelGrid.R
