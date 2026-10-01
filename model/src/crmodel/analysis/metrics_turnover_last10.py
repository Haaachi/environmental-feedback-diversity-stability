"""Last-10-day cumulative diversity and turnover metrics.

Only cumulative diversity / turnover diagnostics should use the final
ten days. The core stability and diversity metrics remain last-3-day
metrics in :mod:`crmodel.analysis.metrics_last3`.
"""

from __future__ import annotations

import numpy as np


def shannon_from_values(x: np.ndarray, axis: int = -1) -> np.ndarray:
    """Return Shannon diversity along ``axis`` for non-negative values."""
    total = x.sum(axis=axis, keepdims=True)
    p = np.divide(x, total, out=np.zeros_like(x, dtype=float), where=total > 0)
    log_p = np.where(p > 0, np.log(p), 0.0)
    return -(p * log_p).sum(axis=axis)


def cumulative_diversity_gap_last10(x: np.ndarray) -> dict[str, np.ndarray]:
    """Compute cumulative-last10 diversity diagnostics.

    Parameters
    ----------
    x
        Endpoint tensor with shape ``(..., S, 10)``. Use B for the main
        paper-facing readout unless explicitly diagnosing live abundance.

    Returns
    -------
    dict
        ``shannon_last``, ``shannon_cum10``, ``delta_shannon_cum10_last``,
        ``hill_last``, ``hill_cum10``, and ``delta_hill_cum10_last``.
    """
    x_last = x[..., -1]
    x_cum = x.sum(axis=-1)

    h_last = shannon_from_values(x_last, axis=-1)
    h_cum = shannon_from_values(x_cum, axis=-1)
    hill_last = np.exp(h_last)
    hill_cum = np.exp(h_cum)

    return {
        "shannon_last": h_last,
        "shannon_cum10": h_cum,
        "delta_shannon_cum10_last": h_cum - h_last,
        "hill_last": hill_last,
        "hill_cum10": hill_cum,
        "delta_hill_cum10_last": hill_cum - hill_last,
    }

