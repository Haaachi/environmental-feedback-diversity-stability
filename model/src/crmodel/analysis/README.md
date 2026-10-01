# Analysis modules

This package contains metric functions used by `run.py postprocess`.

```text
metrics.py                  current last-3 and last-10 metric implementation
metrics_last3.py            retained last-3 helper
metrics_turnover_last10.py  retained last-10 turnover helper
```

The main phase-map analysis uses B-based endpoints. Stability, collapse,
community CV, and last-day diversity are computed from the final three passage
endpoints. Cumulative diversity, synchrony, participation ratio, and turnover
diagnostics use the final ten endpoints.
