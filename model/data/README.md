# Data layout

Large raw simulation outputs are not included in this repository.

Restore or download data under:

```text
data/runs/
```

The lightweight registry is:

```text
data/runs/datasets.yaml
```

New `run.py`-style outputs should follow:

```text
data/runs/<dataset_id>/
  <scan>_meta.h5
  raw/<scan>_rows/<scan>_row_%04d.h5
  processed/<scan>_metrics_last3.npz
  processed/<scan>_metrics_last10.npz
  snapshot/run_manifest.json
```

Legacy old-format data, if restored, should keep its original contract:

```text
data/runs/legacy_data_Final/
  CR_phase1_meta.h5
  CR_phase1_processed.npz
  CR_phase2_meta.h5
  CR_phase2_processed.npz
  phase1_rows/
  phase2_rows/
```
