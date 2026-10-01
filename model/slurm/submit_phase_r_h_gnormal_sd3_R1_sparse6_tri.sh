#!/bin/bash
set -e
mkdir -p slurm_logs

CONFIG="configs/phase_r_h_gnormal_sd3_R1_sparse6_tri.yaml"

echo "=== Submitting sparse-resource gamma-normal sd3 R1 phase_r_h ==="
echo "  config : ${CONFIG}"
echo "  R0     : 1.0 for all resources"
echo "  gamma  : Normal(0, 3^2), quenched independently per community"
echo "  traits : sparse_random_links, resource_degree=6"
echo "  scan   : r in [0,1], h in [0,0.2], 100x100"
echo ""

CONDA_INIT="source ${CONDA_SH:-$HOME/anaconda3/etc/profile.d/conda.sh} && conda activate ${CRMODEL_ENV:-crph-model}"

JOB=$(sbatch --parsable \
  --array=1-100 \
  --cpus-per-task=32 \
  --mem=16G \
  --time=3:00:00 \
  --job-name=crph_rh_gn3_R1 \
  --output=slurm_logs/rh_gn3_R1_%a.out \
  --error=slurm_logs/rh_gn3_R1_%a.err \
  --wrap="$CONDA_INIT && python run.py simulate --config ${CONFIG} --row-index-env SLURM_ARRAY_TASK_ID --workers \$SLURM_CPUS_PER_TASK")

echo "Job submitted: $JOB"
echo "Monitor:  squeue -u \$USER"
echo "Status:   python run.py status --config ${CONFIG}"
