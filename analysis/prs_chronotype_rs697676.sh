#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --job-name=prs336
#SBATCH --time=00:30:00
#SBATCH --partition=general
#SBATCH --account=a_imb_cpdg
#SBATCH -o /scratch/project_mnt/S0007/solal/circadian/Results/PRS/%x_%j.log

# Download the scoring file on the login node first (see the R script header);
# the script reuses it from the output directory.
cd "${SLURM_SUBMIT_DIR}"
module load r/4.4.2-heavy
Rscript prs_chronotype_rs697676.R /scratch/project_mnt/S0007/solal/circadian/Results/PRS
