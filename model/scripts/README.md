# Scripts

Main manuscript plotting:

```text
scripts/plot/plot_main_stability_diversity_phase.py
```

Default output convention: `community_CV_B > 0.1`, `osc_frac_min = 0.10`,
directory suffix `cv010_oscmin010`.

Batch plotting across registered/available datasets:

```text
scripts/plot/plot_all_stability_diversity_phase.py
```

Dataset availability check:

```text
scripts/check_data_catalog.py
```

End-to-end smoke test:

```text
scripts/smoke_test.py
```

Optional S-scaling diagnostics:

```text
scripts/run_s_scaling.py
scripts/plot_community_cv_hist.py
scripts/plot_oscillating_communities.py
```

Generated figures should go under `figures/`, which is ignored by default.
