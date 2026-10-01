"""Last-3-day metrics used by the original phase-diagram analysis.

These functions intentionally preserve the existing convention:
all stability, collapse, survival, and last-day diversity metrics are
computed from the final three daily endpoints.
"""

from __future__ import annotations

import numpy as np


COLLAPSE_THRESH = 1e-3
REL_AB_THRESH = 1e-3
CV_THRESH = 1e-1


def compute_metrics_last3(
    x: np.ndarray,
    rel_thresh: float = REL_AB_THRESH,
    collapse_thresh: float = COLLAPSE_THRESH,
    cv_thresh: float = CV_THRESH,
) -> dict[str, np.ndarray]:
    """Compute per-ensemble metrics from an endpoint tensor.

    Parameters
    ----------
    x
        Array with shape ``(..., S, 3)``. The leading dimensions are
        preserved, typically ``(n_axis2, n_ensemble)`` for one row.

    Returns
    -------
    dict
        Arrays over the leading dimensions: total, is_collapsed,
        surv_frac, shannon, community_CV, total_CV, is_oscillating.
    """
    per_cycle_total = x.sum(axis=-2)
    total = per_cycle_total.mean(axis=-1)
    is_collapsed = total < collapse_thresh

    total_CV = np.where(
        is_collapsed,
        0.0,
        per_cycle_total.std(axis=-1, ddof=1) / np.where(is_collapsed, 1.0, total),
    )

    x_last = x[..., -1]
    total_last = x_last.sum(axis=-1)
    safe_total_last = np.where(total_last > 0, total_last, 1.0)[..., np.newaxis]
    rel_ab = x_last / safe_total_last
    alive = rel_ab > rel_thresh
    surv_frac = np.where(is_collapsed, 0.0, alive.sum(axis=-1) / x.shape[-2])

    pos_mask = x_last > 0
    pos_sum = x_last.sum(axis=-1, where=pos_mask, keepdims=True)
    safe_pos_sum = np.where(pos_sum > 0, pos_sum, 1.0)
    p = np.where(pos_mask, x_last / safe_pos_sum, 0.0)
    log_p = np.where(p > 0, np.log(p), 0.0)
    shannon = np.where(is_collapsed, 0.0, -(p * log_p).sum(axis=-1))

    sp_std = x.std(axis=-1, ddof=1)
    sum_all_std = sp_std.sum(axis=-1)
    safe_total = np.where(total > 0, total, 1.0)
    community_CV = np.where(is_collapsed, 0.0, sum_all_std / safe_total)

    is_oscillating = (~is_collapsed) & (community_CV >= cv_thresh)

    return {
        "total": total,
        "is_collapsed": is_collapsed,
        "surv_frac": surv_frac,
        "shannon": shannon,
        "community_CV": community_CV,
        "total_CV": total_CV,
        "is_oscillating": is_oscillating,
    }

