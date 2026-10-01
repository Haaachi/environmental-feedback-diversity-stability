#!/bin/bash
set -e
mkdir -p slurm_logs

CONDA_INIT="source ${CONDA_SH:-$HOME/anaconda3/etc/profile.d/conda.sh} && conda activate ${CRMODEL_ENV:-crph-model}"

CONFIGS=(
  "configs/phase_r_h_gnormal_sd3_R1_sparse6_tri_S1_h015_e500.yaml"
  "configs/phase_r_h_gnormal_sd3_R1_sparse6_tri_S12_h015_e5000.yaml"
  "configs/phase_r_h_gnormal_sd3_R1_sparse6_tri_S24_h015_e500.yaml"
)

NAMES=(
  "S1_h015_e500"
  "S12_h015_e5000"
  "S24_h015_e500"
)

echo "=== Submitting postprocess jobs for sd3/R1 h015 S-sensitivity ==="
echo ""

for idx in "${!CONFIGS[@]}"; do
  CONFIG="${CONFIGS[$idx]}"
  NAME="${NAMES[$idx]}"
  JOB=$(sbatch --parsable \
    --cpus-per-task=4 \
    --mem=64G \
    --time=4:00:00 \
    --job-name="crph_pp_${NAME}" \
    --output="slurm_logs/rh_gn3_R1_${NAME}_postprocess.out" \
    --error="slurm_logs/rh_gn3_R1_${NAME}_postprocess.err" \
    --wrap="$CONDA_INIT && python run.py postprocess --config ${CONFIG} --strict")
  echo "Submitted postprocess ${NAME}: ${JOB}"
done
