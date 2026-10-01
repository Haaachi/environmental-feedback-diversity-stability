"""Comprehensive metrics for CR-pH model phase diagram analysis.

All functions are vectorized over arbitrary leading dimensions.
Typical input shape: (n_axis2, n_ensemble, S, T).
"""

from __future__ import annotations

import numpy as np

COLLAPSE_THRESH = 1e-3
REL_AB_THRESH = 0.01
CV_THRESH = 0.1


def _shannon(x: np.ndarray, axis: int = -1) -> np.ndarray:
    total = x.sum(axis=axis, keepdims=True)
    p = np.divide(x, total, out=np.zeros_like(x, dtype=float), where=total > 0)
    with np.errstate(divide="ignore"):
        log_p = np.where(p > 0, np.log(p), 0.0)
    return -(p * log_p).sum(axis=axis)


def _richness(x: np.ndarray, rel_thresh: float = REL_AB_THRESH,
              axis: int = -1) -> np.ndarray:
    total = x.sum(axis=axis, keepdims=True)
    safe_total = np.where(total > 0, total, 1.0)
    rel = x / safe_total
    return (rel > rel_thresh).sum(axis=axis)


def _bray_curtis_consecutive(x: np.ndarray) -> np.ndarray:
    """Mean Bray-Curtis dissimilarity between consecutive time points.

    x: (..., S, T). Returns (...,).
    """
    T = x.shape[-1]
    if T < 2:
        return np.zeros(x.shape[:-2])
    bc_sum = np.zeros(x.shape[:-2])
    for t in range(T - 1):
        a = x[..., t]
        b = x[..., t + 1]
        num = np.abs(a - b).sum(axis=-1)
        den = (a + b).sum(axis=-1)
        with np.errstate(invalid="ignore"):
            bc_sum += np.where(den > 0, num / den, 0.0)
    return bc_sum / (T - 1)


def _synchrony_phi(x: np.ndarray) -> np.ndarray:
    """Loreau & de Mazancourt synchrony: Var(total) / [sum SD_i]^2.

    x: (..., S, T). Returns (...,). NaN when denominator is zero.
    """
    total_ts = x.sum(axis=-2)
    var_total = np.var(total_ts, axis=-1, ddof=1)
    std_per_sp = np.std(x, axis=-1, ddof=1)
    sum_std = std_per_sp.sum(axis=-1)
    denom = sum_std ** 2
    with np.errstate(invalid="ignore"):
        return np.where(denom > 0, var_total / denom, np.nan)


def _std_participation_ratio(x: np.ndarray) -> np.ndarray:
    """Standard-deviation participation ratio.

    x: (..., S, T). Returns (...,).

    PR_sigma = (sum sigma_i)^2 / sum(sigma_i^2)
      Range [1, S]. Higher = more species contribute to total fluctuation.
      Equals S when all species fluctuate equally; equals 1 when only one does.

    Returns NaN where all sigma_i == 0 (no variance).
    """
    sigma = np.std(x, axis=-1, ddof=1)  # (..., S)
    sum_sigma = sigma.sum(axis=-1)       # (...,)
    sum_sigma_sq = (sigma ** 2).sum(axis=-1)  # (...,)
    with np.errstate(invalid="ignore"):
        return np.where(sum_sigma_sq > 0,
                        sum_sigma ** 2 / sum_sigma_sq, np.nan)


def _eigenvalue_dominance(x: np.ndarray) -> np.ndarray:
    """Ratio of largest eigenvalue to trace of covariance matrix.

    x: (..., S, T). Returns (...,).

    eig_dom = lambda_max / tr(C)
      Range [0, 1]. Higher = single coherent mode dominates total variance.

    Returns NaN where trace(C) == 0.
    For S=1, always returns 1.0 (trivially, the only eigenvalue IS the trace).
    """
    leading = x.shape[:-2]
    S, T = x.shape[-2], x.shape[-1]

    # S=1 short-circuit: covariance matrix is 1x1, eig_dom = 1 trivially
    if S == 1:
        var = np.var(x[..., 0, :], axis=-1, ddof=1)
        result = np.where(var > 0, 1.0, np.nan)
        return result

    # Vectorized with chunking to avoid OOM for large S
    flat = x.reshape(-1, S, T)
    n = flat.shape[0]
    result = np.full(n, np.nan, dtype=np.float64)

    # Process in chunks of ~10000 to keep memory bounded
    CHUNK = 10000
    for start in range(0, n, CHUNK):
        end = min(start + CHUNK, n)
        chunk = flat[start:end]  # (chunk_size, S, T)

        mean = chunk.mean(axis=-1, keepdims=True)
        centered = chunk - mean
        C = np.einsum("nsi,nri->nsr", centered, centered) / (T - 1)

        trace = np.trace(C, axis1=-2, axis2=-1)
        eigvals = np.linalg.eigvalsh(C)
        lam_max = eigvals[:, -1]

        with np.errstate(invalid="ignore"):
            result[start:end] = np.where(trace > 0, lam_max / trace, np.nan)

    return result.reshape(leading)


def compute_last3_metrics(
    x: np.ndarray,
    collapse_thresh: float = COLLAPSE_THRESH,
    rel_thresh: float = REL_AB_THRESH,
    cv_thresh: float = CV_THRESH,
) -> dict[str, np.ndarray]:
    """Stability and diversity metrics from the last 3 daily endpoints.

    Parameters
    ----------
    x : ndarray, shape (..., S, 3)

    Returns
    -------
    dict — arrays shaped (...,)
    """
    per_cycle_total = x.sum(axis=-2)
    total = per_cycle_total.mean(axis=-1)
    is_collapsed = total < collapse_thresh

    safe_total = np.where(is_collapsed, 1.0, total)
    total_CV = np.where(
        is_collapsed, 0.0,
        per_cycle_total.std(axis=-1, ddof=1) / safe_total)

    sp_std = x.std(axis=-1, ddof=1)
    community_CV = np.where(
        is_collapsed, 0.0, sp_std.sum(axis=-1) / safe_total)

    is_oscillating = (~is_collapsed) & (community_CV >= cv_thresh)

    x_last = x[..., -1]
    shannon = np.where(is_collapsed, 0.0, _shannon(x_last, axis=-1))
    richness = np.where(is_collapsed, 0, _richness(x_last, rel_thresh, axis=-1))
    S_count = x.shape[-2]
    surv_frac = np.where(is_collapsed, 0.0, richness / S_count)
    bray_curtis = np.where(is_collapsed, 0.0, _bray_curtis_consecutive(x))

    return {
        "total": total,
        "is_collapsed": is_collapsed,
        "surv_frac": surv_frac,
        "shannon": shannon,
        "richness": richness,
        "community_CV": community_CV,
        "total_CV": total_CV,
        "bray_curtis": bray_curtis,
        "is_oscillating": is_oscillating,
    }


def compute_last10_metrics(
    x: np.ndarray,
    collapse_thresh: float = COLLAPSE_THRESH,
    rel_thresh: float = REL_AB_THRESH,
    cv_thresh: float = CV_THRESH,
) -> dict[str, np.ndarray]:
    """Turnover and temporal diversity metrics from last 10 endpoints.

    Parameters
    ----------
    x : ndarray, shape (..., S, T)  where T >= 2

    Returns
    -------
    dict — arrays shaped (...,)
    """
    per_cycle_total = x.sum(axis=-2)
    total = per_cycle_total.mean(axis=-1)
    is_collapsed = total < collapse_thresh

    safe_total = np.where(is_collapsed, 1.0, total)
    sp_std = x.std(axis=-1, ddof=1)
    community_CV = sp_std.sum(axis=-1) / safe_total
    is_oscillating = (~is_collapsed) & (community_CV >= cv_thresh)

    x_last = x[..., -1]
    x_cum = x.sum(axis=-1)

    shannon_last = _shannon(x_last, axis=-1)
    shannon_cum = _shannon(x_cum, axis=-1)
    delta_shannon = shannon_cum - shannon_last

    richness_last = _richness(x_last, rel_thresh, axis=-1)
    richness_cum = _richness(x_cum, rel_thresh, axis=-1)
    delta_richness = richness_cum - richness_last

    phi = _synchrony_phi(x)
    phi = np.where(is_collapsed | (~is_oscillating), np.nan, phi)

    std_pr = _std_participation_ratio(x)
    std_pr = np.where(is_collapsed | (~is_oscillating), np.nan, std_pr)

    eig_dom = _eigenvalue_dominance(x)
    eig_dom = np.where(is_collapsed | (~is_oscillating), np.nan, eig_dom)

    return {
        "shannon_last": shannon_last,
        "shannon_cum": shannon_cum,
        "delta_shannon": delta_shannon,
        "richness_last": richness_last,
        "richness_cum": richness_cum,
        "delta_richness": delta_richness,
        "synchrony_phi": phi,
        "std_participation_ratio": std_pr,
        "eigenvalue_dominance": eig_dom,
    }
