# Reproducing the model analyses

This document gives the command sequence used for the manuscript model
analyses. Large simulation outputs are not committed to the repository; restore
them under `data/runs/` or regenerate them with `run.py`.

## 1. Create the environment

```bash
conda env create -f environment.yml
conda activate crph-model
```

## 2. Run a smoke test

```bash
python scripts/smoke_test.py
```

The smoke test uses:

```text
configs/smoke_test.yaml
```

It runs a small `2 x 2` r-h grid with two random communities per grid point,
then postprocesses the output. This validates the solver, config loader, HDF5
output contract, and metric pipeline.

## 3. Main manuscript r-h scan

The main manuscript configuration is:

```text
configs/phase_r_h_gnormal_sd3_R1_sparse6_tri_S12_h015_e5000.yaml
```

Run the full scan:

```bash
python run.py simulate \
  --config configs/phase_r_h_gnormal_sd3_R1_sparse6_tri_S12_h015_e5000.yaml \
  --workers 32
```

For a SLURM array, submit one row per task using the scripts in `slurm/`.

Check completion:

```bash
python run.py status \
  --config configs/phase_r_h_gnormal_sd3_R1_sparse6_tri_S12_h015_e5000.yaml
```

Postprocess:

```bash
python run.py postprocess \
  --config configs/phase_r_h_gnormal_sd3_R1_sparse6_tri_S12_h015_e5000.yaml \
  --strict
```

Expected output layout:

```text
data/runs/rh_gnormal_sd3_R1_sparse6_tri_S12_h015_g100_e5000/
  phase_r_h_meta.h5
  raw/phase_r_h_rows/phase_r_h_row_%04d.h5
  processed/phase_r_h_metrics_last3.npz
  processed/phase_r_h_metrics_last10.npz
  snapshot/run_manifest.json
```

## 4. Main stability-diversity phase maps

After raw rows are present:

```bash
python scripts/plot/plot_main_stability_diversity_phase.py \
  --run-dir data/runs/rh_gnormal_sd3_R1_sparse6_tri_S12_h015_g100_e5000 \
  --out-dir figures/main_stability_diversity/rh_gnormal_sd3_R1_sparse6_tri_S12_h015_g100_e5000_cv010_oscmin010 \
  --osc-frac-min 0.10
```

The plotting script recomputes B-based labels and writes individual phase-map
panels plus a compressed metric file. The manuscript PR uses unthresholded B
over the final three passages, with sample standard deviations and division by
the total species-pool size S. Both the reorganized and legacy phase2 input
formats use this convention. Last-ten PR diagnostics are optional analyses.

After updating from an earlier repository version, rerun this plotting command
from the raw model outputs to replace the previous PR maps and metric file:

```text
figures/main_stability_diversity/rh_gnormal_sd3_R1_sparse6_tri_S12_h015_g100_e5000_cv010_oscmin010/
  stability_diversity_metrics_cv010_absB.npz
  *.png
  *.pdf
```

## 5. Additional analyses

The repository also includes configurations for:

```text
gamma controls:       gamma = 0, gamma [-10,+10], gamma [-5,+10]
historical gm5p5:     gamma [-5,+5] comparison runs
top-hat Phi(p):       pH-response sensitivity
S sensitivity:        S = 1, 24, 36, 48
sd5/R0.5:             alternative gamma/resource-supply sensitivity
S-scaling:            fixed-parameter species-pool-size scans
h-f scan:             dilution-related scan, if restored externally
```

Registered dataset names and expected file locations are listed in:

```text
data/runs/datasets.yaml
```
