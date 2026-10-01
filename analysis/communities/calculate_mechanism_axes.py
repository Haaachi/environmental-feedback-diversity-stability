#!/usr/bin/env python3
# ===========================================================
# calculate_mechanism_axes.py
#
# Mechanism-oriented metrics for separating dominant-species-driven
# fluctuation from distributed/complexity-driven fluctuation.
#
# Output:
#   processed/mechanism_axes/mechanism_metrics.csv
# ===========================================================

import glob
import os
from itertools import combinations

import numpy as np
import pandas as pd


BASE_DIR = os.environ.get("COMMUNITY_WORKSPACE", os.path.dirname(os.path.abspath(__file__)))
PROC_DIR = os.path.join(BASE_DIR, "processed")
OUT_DIR = os.path.join(PROC_DIR, "mechanism_axes")

COLLAPSE_THRESHOLD = 0.05
MIN_TOTAL_BIOMASS_FOR_DIVERSITY = 0.05

COMMUNITY_CV_THRESHOLD = 0.25

WINDOW_CONFIG = {
    "last3": {
        "temperature": {"days": [8, 9, 10], "replicas": [1, 2, 3]},
        "mortality": {"days": [4, 5, 6], "replicas": [1]},
    },
    "last4": {
        "temperature": {"days": [7, 8, 9, 10], "replicas": [1]},
        "mortality": {"days": [3, 4, 5, 6], "replicas": [1]},
    },
    "long": {
        "temperature": {"days": list(range(3, 11)), "replicas": [1]},
        "mortality": {"days": [3, 4, 5, 6], "replicas": [1]},
    },
}


def load_abs_abundance():
    frames = []
    for experiment in ("mortality", "temperature"):
        pattern = os.path.join(PROC_DIR, experiment, "W*", "abs_abundance.csv")
        for path in sorted(glob.glob(pattern)):
            frames.append(pd.read_csv(path))

    if not frames:
        raise FileNotFoundError(f"No abs_abundance.csv files found under {PROC_DIR}")

    df = pd.concat(frames, ignore_index=True)
    df["experiment"] = df["experiment"].astype(str)
    df["condition"] = df["condition"].astype(str)
    df["community"] = df["community"].astype(int)
    df["replica"] = df["replica"].astype(int)
    df["day"] = df["day"].astype(int)
    df["taxon"] = df["taxon"].astype(str)
    return df


def bray_curtis(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    denom = np.sum(x) + np.sum(y)
    if denom <= 0:
        return np.nan
    return float(np.sum(np.abs(x - y)) / denom)


def shannon(vals):
    vals = np.asarray(vals, dtype=float)
    total = np.sum(vals)
    if total <= 0:
        return np.nan
    p = vals[vals > 0] / total
    if len(p) == 0:
        return np.nan
    return float(-np.sum(p * np.log(p)))


def participation_metrics(weights):
    weights = np.asarray(weights, dtype=float)
    total = np.sum(weights)
    if total <= 0:
        return {
            "ipr": np.nan,
            "pr": np.nan,
            "max_share": np.nan,
        }
    shares = weights / total
    ipr = float(np.sum(shares ** 2))
    return {
        "ipr": ipr,
        "pr": float(1.0 / ipr) if ipr > 0 else np.nan,
        "max_share": float(np.max(shares)),
    }


def eig_metrics(matrix, standardize):
    x = np.asarray(matrix, dtype=float)
    if x.ndim != 2 or x.shape[0] < 3 or x.shape[1] < 2:
        return {
            "rho1": np.nan,
            "mode_dispersion": np.nan,
            "eigen_ipr": np.nan,
            "eigen_pr": np.nan,
            "n_modes_input": 0,
        }

    sd = np.std(x, axis=0, ddof=1)
    keep = sd > 0
    x = x[:, keep]
    sd = sd[keep]
    if x.shape[1] < 2:
        return {
            "rho1": np.nan,
            "mode_dispersion": np.nan,
            "eigen_ipr": np.nan,
            "eigen_pr": np.nan,
            "n_modes_input": int(x.shape[1]),
        }

    x = x - np.mean(x, axis=0)
    if standardize:
        x = x / sd

    cov = np.cov(x, rowvar=False, ddof=1)
    eigvals = np.linalg.eigvalsh(cov)
    eigvals = np.clip(eigvals, 0, None)
    total = float(np.sum(eigvals))
    if total <= 0:
        return {
            "rho1": np.nan,
            "mode_dispersion": np.nan,
            "eigen_ipr": np.nan,
            "eigen_pr": np.nan,
            "n_modes_input": int(x.shape[1]),
        }

    shares = eigvals / total
    rho1 = float(np.max(shares))
    eigen_ipr = float(np.sum(shares ** 2))
    return {
        "rho1": rho1,
        "mode_dispersion": float(1.0 - rho1),
        "eigen_ipr": eigen_ipr,
        "eigen_pr": float(1.0 / eigen_ipr) if eigen_ipr > 0 else np.nan,
        "n_modes_input": int(x.shape[1]),
    }


def clean_absolute_matrix(df_rep, days):
    df = df_rep.copy()
    day_totals = df.groupby("day")["rel_abund"].sum().rename("day_rel_total")
    df = df.merge(day_totals, on="day", how="left")
    df["rel_abund_clean"] = np.where(
        df["day_rel_total"] > 0,
        df["rel_abund"] / df["day_rel_total"] * 100.0,
        np.nan,
    )
    df["abs_abund_clean"] = np.where(
        df["rel_abund_clean"].isna() | df["OD"].isna(),
        np.nan,
        df["rel_abund_clean"] / 100.0 * df["OD"],
    )

    taxa = sorted(df["taxon"].dropna().unique())
    abs_mat = (
        df.pivot_table(index="day", columns="taxon",
                       values="abs_abund_clean", aggfunc="sum")
        .reindex(index=days, columns=taxa)
        .fillna(0.0)
    )
    rel_mat = (
        df.pivot_table(index="day", columns="taxon",
                       values="rel_abund_clean", aggfunc="sum")
        .reindex(index=days, columns=taxa)
        .fillna(0.0)
    )
    return abs_mat, rel_mat


def compute_temporal_bc(rel_mat):
    vals = []
    arr = rel_mat.to_numpy(dtype=float)
    for i, j in combinations(range(arr.shape[0]), 2):
        bc = bray_curtis(arr[i, :], arr[j, :])
        if not np.isnan(bc):
            vals.append(bc)
    return float(np.mean(vals)) if vals else np.nan


def compute_diversity_gap(abs_mat):
    arr = abs_mat.to_numpy(dtype=float)
    totals = arr.sum(axis=1)

    daily_h = [
        shannon(arr[i, :])
        for i in range(arr.shape[0])
        if totals[i] >= MIN_TOTAL_BIOMASS_FOR_DIVERSITY
    ]
    mean_daily = float(np.mean(daily_h)) if daily_h else np.nan
    cumulative = shannon(arr.sum(axis=0))
    gap = cumulative - mean_daily if not np.isnan(cumulative) and not np.isnan(mean_daily) else np.nan
    return cumulative, mean_daily, gap


def compute_leave_one_out(abs_mat):
    arr = abs_mat.to_numpy(dtype=float)
    if arr.shape[0] < 2:
        return {
            "leave_one_out_total_std_max_effect": np.nan,
            "leave_one_out_total_std_min_effect": np.nan,
            "leave_one_out_total_std_max_abs_effect": np.nan,
            "leave_one_out_total_std_driver_taxon": "",
            "leave_one_out_sum_abs_std_max_effect": np.nan,
            "leave_one_out_sum_abs_std_min_effect": np.nan,
            "leave_one_out_sum_abs_std_max_abs_effect": np.nan,
            "leave_one_out_sum_abs_std_driver_taxon": "",
        }

    sd = np.std(arr, axis=0, ddof=1)
    full_sum_abs_std = float(np.nansum(sd))
    total = arr.sum(axis=1)
    full_std = float(np.std(total, ddof=1))

    if full_std <= 0:
        total_effects = []
    else:
        total_effects = []
        taxa = list(abs_mat.columns)
        for idx, taxon in enumerate(taxa):
            minus_std = float(np.std(total - arr[:, idx], ddof=1))
            total_effects.append((taxon, (full_std - minus_std) / full_std))

    if full_sum_abs_std <= 0:
        sum_effects = []
    else:
        sum_effects = []
        taxa = list(abs_mat.columns)
        for idx, taxon in enumerate(taxa):
            remaining = np.delete(arr, idx, axis=1)
            if remaining.shape[1] == 0:
                minus_sum = 0.0
            else:
                minus_sum = float(np.nansum(np.std(remaining, axis=0, ddof=1)))
            sum_effects.append(
                (taxon, (full_sum_abs_std - minus_sum) / full_sum_abs_std)
            )

    def summarize_effects(effects, prefix):
        if not effects:
            return {
                f"{prefix}_max_effect": np.nan,
                f"{prefix}_min_effect": np.nan,
                f"{prefix}_max_abs_effect": np.nan,
                f"{prefix}_driver_taxon": "",
            }
        values = np.array([x[1] for x in effects], dtype=float)
        max_i = int(np.nanargmax(values))
        return {
            f"{prefix}_max_effect": float(np.nanmax(values)),
            f"{prefix}_min_effect": float(np.nanmin(values)),
            f"{prefix}_max_abs_effect": float(np.nanmax(np.abs(values))),
            f"{prefix}_driver_taxon": effects[max_i][0],
        }

    return {
        **summarize_effects(total_effects, "leave_one_out_total_std"),
        **summarize_effects(sum_effects, "leave_one_out_sum_abs_std"),
    }


def collapsed_record(window, experiment, condition, community, replica, days, mean_total):
    return {
        "window": window,
        "experiment": experiment,
        "condition": condition,
        "community": community,
        "replica": replica,
        "days": ",".join(map(str, days)),
        "n_days": len(days),
        "n_species": 0,
        "collapsed": True,
        "mean_total_abs": mean_total,
        "community_cv": 0.0,
        "sum_abs_std": 0.0,
        "temporal_bc_mean": 0.0,
        "total_biomass_std": 0.0,
        "synchrony_phi": np.nan,
        "asynchrony": np.nan,
        "species_ipr_var": np.nan,
        "species_pr_var": np.nan,
        "max_var_share": np.nan,
        "species_ipr_std": np.nan,
        "species_pr_std": np.nan,
        "max_std_share": np.nan,
        "rho1_corr": np.nan,
        "mode_dispersion_corr": np.nan,
        "eigen_ipr_corr": np.nan,
        "eigen_pr_corr": np.nan,
        "rho1_cov": np.nan,
        "mode_dispersion_cov": np.nan,
        "eigen_ipr_cov": np.nan,
        "eigen_pr_cov": np.nan,
        "leave_one_out_total_std_max_effect": np.nan,
        "leave_one_out_total_std_min_effect": np.nan,
        "leave_one_out_total_std_max_abs_effect": np.nan,
        "leave_one_out_total_std_driver_taxon": "",
        "leave_one_out_sum_abs_std_max_effect": np.nan,
        "leave_one_out_sum_abs_std_min_effect": np.nan,
        "leave_one_out_sum_abs_std_max_abs_effect": np.nan,
        "leave_one_out_sum_abs_std_driver_taxon": "",
        "cumulative_shannon": np.nan,
        "mean_daily_shannon": np.nan,
        "cumulative_excess_shannon": np.nan,
        "community_cv_threshold": COMMUNITY_CV_THRESHOLD,
        "fluctuating_species_instability": False,
    }


def compute_group(window, experiment, condition, community, replica, days, df_rep):
    mean_total_abs = float(df_rep.groupby("day")["abs_abund"].sum().mean())
    if np.isnan(mean_total_abs) or mean_total_abs < COLLAPSE_THRESHOLD:
        return collapsed_record(
            window, experiment, condition, community, replica, days, mean_total_abs
        )

    abs_mat, rel_mat = clean_absolute_matrix(df_rep, days)
    arr = abs_mat.to_numpy(dtype=float)
    n_days, n_species = arr.shape

    sd = np.std(arr, axis=0, ddof=1) if n_days >= 2 else np.full(n_species, np.nan)
    var = sd ** 2
    sum_abs_std = float(np.nansum(sd))
    community_cv = (
        sum_abs_std / mean_total_abs
        if mean_total_abs > 0 and not np.isnan(sum_abs_std)
        else np.nan
    )
    total_series = arr.sum(axis=1)
    total_biomass_std = (
        float(np.std(total_series, ddof=1)) if n_days >= 2 else np.nan
    )

    std_part = participation_metrics(sd)
    var_part = participation_metrics(var)
    corr_eig = eig_metrics(arr, standardize=True)
    cov_eig = eig_metrics(arr, standardize=False)
    loo = compute_leave_one_out(abs_mat)

    temporal_bc = compute_temporal_bc(rel_mat)
    cumulative_h, mean_daily_h, diversity_gap = compute_diversity_gap(abs_mat)

    phi = (
        (total_biomass_std ** 2) / (sum_abs_std ** 2)
        if sum_abs_std > 0 and not np.isnan(total_biomass_std)
        else np.nan
    )
    asynchrony = 1.0 - phi if not np.isnan(phi) else np.nan

    fluctuating = (
        not np.isnan(community_cv)
        and community_cv >= COMMUNITY_CV_THRESHOLD
    )

    return {
        "window": window,
        "experiment": experiment,
        "condition": condition,
        "community": community,
        "replica": replica,
        "days": ",".join(map(str, days)),
        "n_days": n_days,
        "n_species": int(n_species),
        "collapsed": False,
        "mean_total_abs": mean_total_abs,
        "community_cv": community_cv,
        "sum_abs_std": sum_abs_std,
        "temporal_bc_mean": temporal_bc,
        "total_biomass_std": total_biomass_std,
        "synchrony_phi": phi,
        "asynchrony": asynchrony,
        "species_ipr_var": var_part["ipr"],
        "species_pr_var": var_part["pr"],
        "max_var_share": var_part["max_share"],
        "species_ipr_std": std_part["ipr"],
        "species_pr_std": std_part["pr"],
        "max_std_share": std_part["max_share"],
        "rho1_corr": corr_eig["rho1"],
        "mode_dispersion_corr": corr_eig["mode_dispersion"],
        "eigen_ipr_corr": corr_eig["eigen_ipr"],
        "eigen_pr_corr": corr_eig["eigen_pr"],
        "rho1_cov": cov_eig["rho1"],
        "mode_dispersion_cov": cov_eig["mode_dispersion"],
        "eigen_ipr_cov": cov_eig["eigen_ipr"],
        "eigen_pr_cov": cov_eig["eigen_pr"],
        **loo,
        "cumulative_shannon": cumulative_h,
        "mean_daily_shannon": mean_daily_h,
        "cumulative_excess_shannon": diversity_gap,
        "community_cv_threshold": COMMUNITY_CV_THRESHOLD,
        "fluctuating_species_instability": bool(fluctuating),
    }


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    df = load_abs_abundance()
    records = []

    for window, exp_cfg in WINDOW_CONFIG.items():
        for experiment, cfg in exp_cfg.items():
            days = cfg["days"]
            replicas = set(cfg["replicas"])

            df_exp = df[df["experiment"] == experiment]
            for condition in sorted(df_exp["condition"].unique()):
                df_cond = df_exp[df_exp["condition"] == condition]
                for community in sorted(df_cond["community"].unique()):
                    df_comm = df_cond[df_cond["community"] == community]
                    for replica in sorted(df_comm["replica"].unique()):
                        if int(replica) not in replicas:
                            continue
                        df_rep = df_comm[
                            (df_comm["replica"] == replica)
                            & (df_comm["day"].isin(days))
                        ].copy()
                        if df_rep.empty:
                            continue
                        records.append(
                            compute_group(
                                window,
                                experiment,
                                condition,
                                int(community),
                                int(replica),
                                days,
                                df_rep,
                            )
                        )

    out = pd.DataFrame.from_records(records)
    out = out.sort_values(
        ["window", "experiment", "condition", "community", "replica"]
    )
    out_path = os.path.join(OUT_DIR, "mechanism_metrics.csv")
    out.to_csv(out_path, index=False)
    print(f"Wrote {len(out)} rows to {out_path}")


if __name__ == "__main__":
    main()
