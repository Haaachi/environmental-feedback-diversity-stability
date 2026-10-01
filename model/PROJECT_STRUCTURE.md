# Project structure

This repository contains simulation and analysis code for the CR-pH feedback
model used in the manuscript.

## Core files

```text
run.py
MODEL_SPEC.md
requirements.txt
environment.yml
pyproject.toml
configs/
src/crmodel/model/
src/crmodel/io/
src/crmodel/analysis/
scripts/plot/plot_main_stability_diversity_phase.py
scripts/plot/plot_all_stability_diversity_phase.py
scripts/check_data_catalog.py
data/runs/datasets.yaml
```

## Main manuscript configuration

```text
configs/phase_r_h_gnormal_sd3_R1_sparse6_tri.yaml
configs/phase_r_h_gnormal_sd3_R1_sparse6_tri_S12_h015_e5000.yaml
```

Dataset id:

```text
rh_gnormal_sd3_R1_sparse6_tri_S12_h015_g100_e5000
```

This is the main Fig. 3 model result: `gamma ~ Normal(0, 3^2)`,
`R0 = 1.0`, sparse 6-resource preference vectors, triangular pH response,
`S = 12`, `M = 24`, `100 x 100` r-h grid, `h_stress <= 0.15`, and 5000
random communities per grid point.

## Optional analyses

```text
configs/prelim/
configs/sensitivity/
configs/phase_r_h_gnormal_sd3_R1_sparse6_tri.yaml
configs/phase_r_h_gnormal_sd5_R05_sparse6_tri.yaml
configs/phase_r_h_tophat*.yaml
scripts/run_s_scaling.py
scripts/plot_community_cv_hist.py
scripts/plot_oscillating_communities.py
experiments/emergent_oscillation_v3/
slurm/
```

These support comparison runs, sensitivity checks, S-scaling diagnostics,
representative trajectory generation, or server execution.

## Configuration scope

The public repository should prioritize configurations that reproduce the
reported figures and the robustness checks discussed in the manuscript. Older
exploratory settings can be documented in `data/runs/datasets.yaml`, but should
not be presented as equal-weight manuscript analyses unless they are used in the
text or supplement.

## Excluded from this repository

Large generated outputs are not committed:

```text
raw HDF5 outputs
processed NPZ outputs
generated figures
job logs
Python bytecode and cache files
server-specific packaging snapshots
older exploratory scripts
```

Restore large simulation outputs under `data/runs/` using the paths recorded in
`data/runs/datasets.yaml`.
