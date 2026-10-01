"""Trait generation for CR-pH model species."""

from __future__ import annotations

import numpy as np


def draw_sparse_preferences(
    rng: np.random.Generator,
    n_species: int,
    n_resources: int,
    degree: int,
) -> np.ndarray:
    """Draw sparse resource preferences with rows summing to one.

    Each species randomly selects ``degree`` out of ``n_resources``
    resources, with Dirichlet-distributed weights.
    """
    degree = int(min(max(degree, 1), n_resources))
    pref = np.zeros((n_species, n_resources), dtype=np.float64)
    for i in range(n_species):
        cols = rng.choice(n_resources, size=degree, replace=False)
        weights = rng.dirichlet(np.ones(degree))
        pref[i, cols] = weights
    return pref


def draw_dense_preferences(
    rng: np.random.Generator,
    n_species: int,
    n_resources: int,
) -> np.ndarray:
    """Draw dense resource preferences with rows summing to one.

    Every species can consume all resources. Weights are drawn from a
    symmetric Dirichlet distribution, producing generalist consumers
    with random relative preferences.
    """
    pref = rng.dirichlet(np.ones(n_resources), size=n_species)
    return pref.astype(np.float64)


def draw_quenched_traits(
    seed: int,
    n_species: int,
    n_resources: int,
    r_mean: float,
    growth_cv: float,
    resource_degree: int,
    p_opt_min: float,
    p_opt_max: float,
    preference_model: str = "sparse_random_links",
) -> tuple[np.ndarray, np.ndarray]:
    """Draw resource preferences, pH optima, and max uptake rates.

    Parameters
    ----------
    preference_model : str
        "sparse_random_links" (default) — each species uses
        ``resource_degree`` resources.
        "dense_dirichlet" — each species can consume all resources.
    """
    rng = np.random.default_rng(seed)

    if preference_model == "dense_dirichlet":
        v_pref = draw_dense_preferences(
            rng=rng,
            n_species=n_species,
            n_resources=n_resources,
        )
    else:
        v_pref = draw_sparse_preferences(
            rng=rng,
            n_species=n_species,
            n_resources=n_resources,
            degree=resource_degree,
        )

    p_opt = rng.uniform(float(p_opt_min), float(p_opt_max), size=n_species)
    if growth_cv <= 0:
        r_i = np.full(n_species, float(r_mean), dtype=np.float64)
    else:
        z = rng.standard_normal(n_species)
        r_i = float(r_mean) * np.maximum(1.0 + float(growth_cv) * z, 0.0)
    v_max = v_pref * r_i[:, None]
    return v_max.astype(np.float64), p_opt.astype(np.float64)
