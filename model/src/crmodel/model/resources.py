"""Resource vector construction from YAML config."""

from __future__ import annotations

import numpy as np


def build_resource_vectors(
    res_cfg: dict, M: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Build R0 and gamma vectors from resource config.

    Supports three modes controlled by ``res_cfg``:

    * **linear_gradient** (``gamma_mode: linear_gradient``): gamma is a
      linear gradient from ``gamma_min`` to ``gamma_max`` across M
      resources.  R0 is uniform (scalar ``R0``).
    * **normal** (``gamma_mode: normal``): gamma is a fixed random draw from
      ``Normal(gamma_mean, gamma_sd)`` using ``gamma_seed``.  R0 is uniform.
    * **normal_per_community**: each community draws its own quenched gamma
      vector from ``Normal(gamma_mean, gamma_sd)``.  The returned gamma vector
      is a NaN sentinel; draw concrete vectors with
      ``draw_community_gamma``.
    * **grouped** (has ``groups`` key): each group specifies count, R0,
      and gamma.  ``R0`` can be per-group or a single scalar applied to
      every resource (uniform mode).
    * **flat** (has ``R0`` and ``gamma`` as lists): used directly.
    """
    gamma_mode = res_cfg.get("gamma_mode", "")

    if gamma_mode == "linear_gradient":
        gamma_min = float(res_cfg["gamma_min"])
        gamma_max = float(res_cfg["gamma_max"])
        gamma = np.linspace(gamma_min, gamma_max, M)
        r0_val = float(res_cfg.get("R0", 1.0))
        R0 = np.full(M, r0_val, dtype=np.float64)
        return R0, gamma

    if gamma_mode == "normal":
        gamma_mean = float(res_cfg.get("gamma_mean", 0.0))
        gamma_sd = float(res_cfg["gamma_sd"])
        gamma_seed = int(res_cfg.get("gamma_seed", 0))
        rng = np.random.default_rng(gamma_seed)
        gamma = rng.normal(gamma_mean, gamma_sd, size=M).astype(np.float64)
        r0_val = float(res_cfg.get("R0", 1.0))
        R0 = np.full(M, r0_val, dtype=np.float64)
        return R0, gamma

    if is_gamma_per_community(res_cfg):
        r0_val = float(res_cfg.get("R0", 1.0))
        R0 = np.full(M, r0_val, dtype=np.float64)
        gamma = np.full(M, np.nan, dtype=np.float64)
        return R0, gamma

    if "groups" in res_cfg:
        R0_list: list[float] = []
        gamma_list: list[float] = []
        uniform_R0 = res_cfg.get("R0")
        for grp in res_cfg["groups"]:
            n = int(grp["count"])
            r0 = float(uniform_R0) if uniform_R0 is not None else float(grp["R0"])
            g = float(grp["gamma"])
            R0_list.extend([r0] * n)
            gamma_list.extend([g] * n)
        if len(R0_list) != M:
            raise ValueError(
                f"Resource groups sum to {len(R0_list)}, expected M={M}"
            )
        return (
            np.array(R0_list, dtype=np.float64),
            np.array(gamma_list, dtype=np.float64),
        )

    R0 = np.asarray(res_cfg["R0"], dtype=np.float64)
    gamma = np.asarray(res_cfg["gamma"], dtype=np.float64)
    if len(R0) != M or len(gamma) != M:
        raise ValueError(f"R0/gamma length mismatch with M={M}")
    return R0, gamma


def is_gamma_per_community(res_cfg: dict) -> bool:
    """Return whether gamma should be drawn once per community."""
    return str(res_cfg.get("gamma_mode", "")).lower() in {
        "normal_per_community",
        "normal_community",
        "quenched_normal",
    }


def draw_community_gamma(res_cfg: dict, M: int, community_seed: int) -> np.ndarray:
    """Draw a community-specific quenched gamma vector.

    A separate seed stream is derived from the community seed so gamma disorder
    follows the community identity without reusing the trait RNG stream.
    """
    gamma_mean = float(res_cfg.get("gamma_mean", 0.0))
    gamma_sd = float(res_cfg["gamma_sd"])
    seed_offset = int(res_cfg.get("gamma_seed_offset", 730_000_000))
    rng = np.random.default_rng(int(community_seed) + seed_offset)
    return rng.normal(gamma_mean, gamma_sd, size=M).astype(np.float64)
