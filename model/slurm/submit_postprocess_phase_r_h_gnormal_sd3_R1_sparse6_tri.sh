#!/bin/bash
set -e
mkdir -p slurm_logs

CONFIG="configs/phase_r_h_gnormal_sd3_R1_sparse6_tri.yaml"

echo "=== Submitting postprocess for sparse-resource gamma-normal sd3 R1 phase_r_h ==="
echo "  config: ${CONFIG}"
echo ""

CONDA_INIT="source ${CONDA_SH:-$HOME/anaconda3/etc/profile.d/conda.sh} && conda activate ${CRMODEL_ENV:-crph-model}"

JOB=$(sbatch --parsable \
  --cpus-per-task=4 \
  --mem=32G \
  --time=2:00:00 \
  --job-name=crph_rh_gn3_R1_pp \
  --output=slurm_logs/rh_gn3_R1_postprocess.out \
  --error=slurm_logs/rh_gn3_R1_postprocess.err \
  --wrap="$CONDA_INIT && python run.py postprocess --config ${CONFIG} --strict")

echo "Postprocess job submitted: $JOB"
echo "Monitor: squeue -u \$USER"
echo "Log:"
echo "  tail -n 80 slurm_logs/rh_gn3_R1_postprocess.out"
echo "  tail -n 80 slurm_logs/rh_gn3_R1_postprocess.err"
