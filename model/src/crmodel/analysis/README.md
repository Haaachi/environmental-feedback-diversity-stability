# Analysis modules

This package contains metric functions used by `run.py postprocess`.

```text
metrics.py                  current last-3 and last-10 metric implementation
metrics_last3.py            retained last-3 helper
metrics_turnover_last10.py  retained last-10 turnover helper
```

The manuscript phase-map script uses B-based endpoints. Stability, collapse,
community CV, and primary fluctuation participation are computed from the final
three passage endpoints; endpoint diversity uses the final passage.
The manuscript PR is calculated from unthresholded B, with sample standard
deviations, and normalized by the total species-pool size S.

The separate last-ten postprocessing output contains optional cumulative
diversity, synchrony, participation and turnover diagnostics. Its participation
ratio is not the final-three-passage PR reported in the manuscript.
