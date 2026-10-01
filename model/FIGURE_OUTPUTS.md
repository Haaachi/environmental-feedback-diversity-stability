# Stability-diversity phase-map outputs

Unified phase-map outputs use:

```text
fluctuation definition: community_CV_B > 0.1
diversity activity threshold: absolute B_i > 1e-3
readout: B-based metrics
minimum oscillating fraction for effect-size display: 0.10
directory suffix: cv010_oscmin010
```

The `oscmin010` convention is a visualization/reporting gate for
diversity-effect and participation maps. The underlying oscillation label is
still defined by `community_CV_B > 0.1`.

Each output directory contains:

```text
01_fluctuation_fraction_cv010.*
01b_collapse_fraction_B.*
02_delta_mean/
02_rank_biserial/
02_cohens_d/
02_spearman_cv_diversity/
03_std_participation_ratio_B_osc_mean.*
03b_std_participation_ratio_B_osc_norm_mean.*
03c_one_minus_max_std_share_B_osc_mean.*
stability_diversity_metrics_cv010_absB.npz
```

The diversity variants are:

```text
richness_last_absB
richness_mean_last3_absB
richness_cum3_absB
shannon_last_absB
shannon_mean_last3_absB
shannon_cum3_absB
```

## Manuscript main output

```text
figures/main_stability_diversity/
  rh_gnormal_sd3_R1_sparse6_tri_S12_h015_g100_e5000_cv010_oscmin010/
```

Data source:

```text
data/runs/rh_gnormal_sd3_R1_sparse6_tri_S12_h015_g100_e5000
```

Parameters:

```text
R0 = 1.0
gamma_j ~ Normal(0, 3^2), independently drawn per community and quenched
S = 12
M = 24
sparse 6-resource preferences
triangular Phi(p)
r_mean in [0, 1]
h_stress in [0, 0.15]
100 x 100 grid
5000 communities per grid point
```

## Reproduction

```bash
python scripts/plot/plot_main_stability_diversity_phase.py \
  --run-dir data/runs/rh_gnormal_sd3_R1_sparse6_tri_S12_h015_g100_e5000 \
  --out-dir figures/main_stability_diversity/rh_gnormal_sd3_R1_sparse6_tri_S12_h015_g100_e5000_cv010_oscmin010 \
  --cv-thresh 0.1 \
  --osc-frac-min 0.10
```

Batch plotting for available registered datasets:

```bash
python scripts/plot/plot_all_stability_diversity_phase.py
```
