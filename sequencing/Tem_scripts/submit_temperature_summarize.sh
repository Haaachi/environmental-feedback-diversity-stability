#!/bin/bash
#SBATCH --job-name=temp_summarize
#SBATCH --partition=cu-1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --time=01:00:00
#SBATCH --output=slurm_logs/summarize_%j.out
#SBATCH --error=slurm_logs/summarize_%j.err

set -euo pipefail

source "${CONDA_SH:-$HOME/anaconda3/etc/profile.d/conda.sh}"
conda activate "${QIIME2_ENV:-qiime2-amplicon-2023.9}"

WORK_DIR="${WORK_DIR:-/lustre/home/zfhu/16s_project/analysis}"
cd "${WORK_DIR}"

echo "Start: $(date)"

qiime demux summarize \
    --i-data "${WORK_DIR}/demux_old.qza" \
    --o-visualization "${WORK_DIR}/demux_old.qzv"

qiime demux summarize \
    --i-data "${WORK_DIR}/demux_new.qza" \
    --o-visualization "${WORK_DIR}/demux_new.qzv"

echo "End: $(date)"
echo "Done."
