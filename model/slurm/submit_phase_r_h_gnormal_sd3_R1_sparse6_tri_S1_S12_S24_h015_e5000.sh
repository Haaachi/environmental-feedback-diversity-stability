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

echo "=== Submitting sd3/R1 sparse-6 h015 S-sensitivity phase_r_h ==="
echo "  gamma : Normal(0, 3^2), quenched independently per community"
echo "  R0    : 1.0"
echo "  h     : [0, 0.15]"
echo "  r     : [0, 1]"
echo "  grid  : 100x100"
echo "  ens   : S=1 and S=24 use 500; S=12 uses 5000"
echo ""

for idx in "${!CONFIGS[@]}"; do
  CONFIG="${CONFIGS[$idx]}"
  NAME="${NAMES[$idx]}"
  JOB=$(sbatch --parsable \
    --array=1-100 \
    --cpus-per-task=32 \
    --mem=32G \
    --time=8:00:00 \
    --job-name="crph_${NAME}" \
    --output="slurm_logs/rh_gn3_R1_${NAME}_%a.out" \
    --error="slurm_logs/rh_gn3_R1_${NAME}_%a.err" \
    --wrap="$CONDA_INIT && python run.py simulate --config ${CONFIG} --row-index-env SLURM_ARRAY_TASK_ID --workers \$SLURM_CPUS_PER_TASK")
  echo "Submitted ${NAME}: ${JOB}"
  echo "  Status: python run.py status --config ${CONFIG}"
done
