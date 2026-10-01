#!/bin/bash
#SBATCH --job-name=mortality_qiime2
#SBATCH --partition=cu-1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=2-00:00:00
#SBATCH --output=slurm_logs/qiime2_%j.out
#SBATCH --error=slurm_logs/qiime2_%j.err

bash "${PIPELINE_SCRIPT_DIR:-${SLURM_SUBMIT_DIR:?Submit this script with sbatch from the repository root}/sequencing/Mor_scripts}/run_mortality_qiime2.sh"
