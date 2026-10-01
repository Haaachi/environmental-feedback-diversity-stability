# CR-pH feedback model

This repository contains simulation and analysis code for a daily-passage
consumer-resource model with pH feedback. The model was used to generate and
analyze the manuscript phase diagrams connecting growth, self-generated
environmental stress, community fluctuation, and diversity-stability
relationships.

The manuscript main simulation is:

```text
rh_gnormal_sd3_R1_sparse6_tri_S12_h015_g100_e5000
```

It corresponds to:

```text
S = 12 species
M = 24 resources
R0_j = 1.0
gamma_j ~ Normal(0, 3^2), drawn independently per community and quenched
6 sparse resources per species
triangular Phi(p)
r_mean x h_stress scan, 100 x 100 grid
h_stress in [0, 0.15]
5000 random communities per grid point
30 passage cycles, final 10 saved
f = 1/30
```

## Repository contents

```text
run.py                         main simulate/postprocess/status entry point
MODEL_SPEC.md                  model equations and parameter convention
REPRODUCE.md                   commands for reproducing model analyses
PROJECT_STRUCTURE.md           file inventory and included analyses
FIGURE_OUTPUTS.md              phase-map output convention
configs/                       manuscript, comparison, and sensitivity configs
src/crmodel/model/             ODE, serial dilution, resource and trait draws
src/crmodel/analysis/          last-3 and last-10 metric computation
src/crmodel/io/                config, manifest, and dataset registry helpers
scripts/plot/                  main stability-diversity phase-map plotting
scripts/run_s_scaling.py       optional S-scaling simulation helper
experiments/emergent_oscillation_v3/
                                representative trajectory experiment code
data/runs/datasets.yaml        lightweight dataset inventory; no large data
```

Large raw HDF5/NPZ outputs are not committed. Restore them under `data/runs/`
following `data/runs/datasets.yaml`.

## Environment

Conda users can create the recommended environment with:

```bash
conda env create -f environment.yml
conda activate crph-model
```

Alternatively, install the Python dependencies with:

```bash
python -m pip install -r requirements.txt
python -m pip install "matplotlib>=3.7"
```

## Quick check

Run a small smoke test before launching large scans:

```bash
python scripts/smoke_test.py
```

This writes a small output under `data/runs/smoke_test/`.

## Main commands

Simulate the manuscript main r-h phase plane:

```bash
python run.py simulate \
  --config configs/phase_r_h_gnormal_sd3_R1_sparse6_tri_S12_h015_e5000.yaml \
  --workers 32
```

Postprocess raw row files:

```bash
python run.py postprocess \
  --config configs/phase_r_h_gnormal_sd3_R1_sparse6_tri_S12_h015_e5000.yaml \
  --strict
```

Check completion:

```bash
python run.py status \
  --config configs/phase_r_h_gnormal_sd3_R1_sparse6_tri_S12_h015_e5000.yaml
```

Generate the main stability-diversity maps after restoring raw rows:

```bash
python scripts/plot/plot_main_stability_diversity_phase.py \
  --run-dir data/runs/rh_gnormal_sd3_R1_sparse6_tri_S12_h015_g100_e5000 \
  --out-dir figures/main_stability_diversity/rh_gnormal_sd3_R1_sparse6_tri_S12_h015_g100_e5000_cv010_oscmin010 \
  --osc-frac-min 0.10
```

## Main readout and classification

The manuscript phase-map analysis uses the B-based readout.

The raw simulation stores the final 10 passage endpoints:

```text
N_last10
B_last10
p_end_last10
seed
gamma, for gamma-normal runs
```

Stability and diversity labels are computed from the final 3 endpoints. A
community is treated as collapsed when mean total B is below `1e-3`, and as
fluctuating when it is non-collapsed and `community_CV_B > 0.1`.

The main stability-diversity effect-size map uses the pairwise probability
difference between fluctuating and stable communities:

```text
Delta_P = P(H_fluctuating > H_stable) - P(H_fluctuating < H_stable)
```

## Notes

The sd5/R0.5 configuration is retained as an alternative/sensitivity run, not
as the manuscript main result. The e1000, h <= 0.2 run is retained as a
supporting broad-scan reference. Historical comparison and sensitivity datasets
are recorded in `data/runs/datasets.yaml`.
