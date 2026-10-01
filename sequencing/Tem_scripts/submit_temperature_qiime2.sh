#!/bin/bash
#SBATCH --job-name=temperature_qiime2
#SBATCH --partition=cu-1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=7-00:00:00
#SBATCH --output=slurm_logs/qiime2_%j.out
#SBATCH --error=slurm_logs/qiime2_%j.err

bash "${PIPELINE_SCRIPT_DIR:-${SLURM_SUBMIT_DIR:?Submit this script with sbatch from the repository root}/sequencing/Tem_scripts}/run_temperature_qiime2.sh"
