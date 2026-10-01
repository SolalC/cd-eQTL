#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=400G
#SBATCH --job-name=01_rhythmicGenes
#SBATCH --time=2:00:00
#SBATCH --partition=general
#SBATCH --array 1-50
#SBATCH --constraint=epyc3
#SBATCH --account=a_imb_cpdg
#SBATCH -o logs/%x_%a.out
#SBATCH -e logs/%x_%a.err

# Submit from this folder: sbatch 01_rhythmicGenes.sh
cd "${SLURM_SUBMIT_DIR}"
module load r/4.4.0-combo-EPYC3-only
Rscript 01_rhythmicGenes.R $SLURM_ARRAY_TASK_ID $SLURM_CPUS_PER_TASK
