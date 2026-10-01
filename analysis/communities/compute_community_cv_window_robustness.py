"""
compute_community_cv_window_robustness.py
=========================================

Robustness analysis for the community_CV time-window choice.

Definitions follow compute_fluctuations.py:
  sum_abs_std = sum_i std(abs_abund_clean_i)

For the window-start robustness curve, the window always ends at the final
experimental day and the first day is swept while keeping at least 3 days:
  mortality:   first day 1..4, end day 6,  R1 only
  temperature: first day 1..8, end day 10, R1 only

Classification is fixed by the last-3-day R1 community_CV criterion
(community_cv >= 0.25):
  mortality:   days 4..6,  R1
  temperature: days 8..10, R1

Temperature replicate-distance analysis uses the last 3 days and all three
replicates. W5 collapsed communities other than C3/C7/C12 are excluded from
that scatter, giving 12*5 - 9 = 51 communities. If a sample-day is collapsed,
its relative-abundance vector is set to all zeros before Bray-Curtis distances
are computed between replicates.
"""

import os
from itertools import combinations

import numpy as np
import pandas as pd


BASE_DIR = os.environ.get("COMMUNITY_WORKSPACE", os.path.dirname(os.path.abspath(__file__)))
PROC_DIR = os.path.join(BASE_DIR, "processed")
OUT_DIR = os.path.join(PROC_DIR, "community_cv_window_robustness")

COLLAPSE_THRESHOLD = 0.05
COMMUNITY_CV_THRESHOLD = 0.25

EXPERIMENT_CONFIG = {
    "mortality": {
        "final_day": 6,
        "start_days": list(range(1, 5)),
        "last3_days": [4, 5, 6],
        "window_replicas": [1],
        "conditions": [f"W{i}" for i in range(1, 6)],
    },
    "temperature": {
        "final_day": 10,
        "start_days": list(range(1, 9)),
        "last3_days": [8, 9, 10],
        "window_replicas": [1],
        "conditions": [f"W{i}" for i in range(1, 6)],
    },
}

TEMPERATURE_W5_ACTIVE_COMMUNITIES = {3, 7, 12}


def read_condition_abundance(experiment, condition):
    path = os.path.join(PROC_DIR, experiment, condition, "abs_abundance.csv")
    if not os.path.exists(path):
        return pd.DataFrame()
    return pd.read_csv(path)


def community_universe(df, experiment, condition):
    communities = set(int(x) for x in df["community"].dropna().unique())
    if experiment == "temperature" and condition == "W5":
        communities.update(range(1, 13))
    return sorted(communities)


def calc_bray_curtis(vec_a, vec_b):
    a = np.array(vec_a, dtype=float)
    b = np.array(vec_b, dtype=float)
    denom = a.sum() + b.sum()
    if denom == 0:
        return 0.0
    return float(np.abs(a - b).sum() / denom)


def calc_temporal_bc_mean(df_win, days):
    if df_win.empty or len(days) < 2:
        return 0.0

    taxa = sorted(df_win["taxon"].dropna().unique())
    if len(taxa) == 0:
        return 0.0

    day_vecs = {}
    for day in days:
        df_day = df_win[df_win["day"] == day]
        if df_day.empty or df_day["abs_abund"].sum() < COLLAPSE_THRESHOLD:
            day_vecs[day] = np.zeros(len(taxa), dtype=float)
            continue

        vals = []
        for taxon in taxa:
            row = df_day[df_day["taxon"] == taxon]
            if len(row) == 0:
                vals.append(0.0)
            else:
                val = row["rel_abund"].values[0]
                vals.append(float(val) if not pd.isna(val) else 0.0)
        day_vecs[day] = np.array(vals, dtype=float)

    bc_vals = [
        calc_bray_curtis(day_vecs[day_a], day_vecs[day_b])
        for day_a, day_b in combinations(days, 2)
    ]
    return round(float(np.mean(bc_vals)), 6) if len(bc_vals) > 0 else 0.0


def calc_community_cv(df_comm_rep, days):
    if df_comm_rep.empty:
        return {
            "community_cv": 0.0,
            "temporal_bc_mean": 0.0,
            "sum_abs_std": 0.0,
            "mean_total_abs": 0.0,
            "n_days_valid": 0,
            "n_taxa_clean": 0,
            "collapsed": True,
        }

    df_win = df_comm_rep[df_comm_rep["day"].isin(days)].copy()
    if df_win.empty:
        return {
            "community_cv": 0.0,
            "temporal_bc_mean": 0.0,
            "sum_abs_std": 0.0,
            "mean_total_abs": 0.0,
            "n_days_valid": 0,
            "n_taxa_clean": 0,
            "collapsed": True,
        }

    total_abs_by_day = df_win.groupby("day")["abs_abund"].sum()
    mean_total_abs = float(total_abs_by_day.mean())
    temporal_bc_mean = calc_temporal_bc_mean(df_win, days)
    if pd.isna(mean_total_abs) or mean_total_abs < COLLAPSE_THRESHOLD:
        return {
            "community_cv": 0.0,
            "temporal_bc_mean": 0.0,
            "sum_abs_std": 0.0,
            "mean_total_abs": mean_total_abs,
            "n_days_valid": int(total_abs_by_day.gt(0).sum()),
            "n_taxa_clean": 0,
            "collapsed": True,
        }

    df_core = df_win.copy()
    day_totals = df_core.groupby("day")["rel_abund"].sum().rename("day_total")
    df_core = df_core.merge(day_totals, on="day", how="left")
    df_core["rel_abund_clean"] = np.where(
        df_core["day_total"] > 0,
        df_core["rel_abund"] / df_core["day_total"] * 100,
        np.nan,
    )
    df_core["abs_abund_clean"] = np.where(
        df_core["rel_abund_clean"].isna() | df_core["OD"].isna(),
        np.nan,
        df_core["rel_abund_clean"] / 100.0 * df_core["OD"],
    )

    taxon_stds = []
    for taxon in sorted(df_core["taxon"].unique()):
        df_taxon = df_core[df_core["taxon"] == taxon]
        vals = []
        for day in days:
            row = df_taxon[df_taxon["day"] == day]
            vals.append(row["abs_abund_clean"].values[0]
                        if len(row) > 0 else np.nan)
        vals = np.array(vals, dtype=float)
        valid = vals[~np.isnan(vals)]
        if len(valid) >= 2:
            taxon_stds.append(float(np.std(valid, ddof=1)))

    if len(taxon_stds) == 0 or mean_total_abs <= 0:
        community_cv = np.nan
        sum_abs_std = np.nan
    else:
        sum_abs_std = float(np.sum(taxon_stds))
        community_cv = float(sum_abs_std / mean_total_abs)

    return {
        "community_cv": round(community_cv, 6),
        "temporal_bc_mean": temporal_bc_mean,
        "sum_abs_std": round(sum_abs_std, 6),
        "mean_total_abs": round(mean_total_abs, 6),
        "n_days_valid": int(total_abs_by_day.gt(0).sum()),
        "n_taxa_clean": int(df_core["taxon"].nunique()),
        "collapsed": False,
    }


def compute_window_cv_table():
    rows = []
    for experiment, cfg in EXPERIMENT_CONFIG.items():
        final_day = cfg["final_day"]
        for condition in cfg["conditions"]:
            df = read_condition_abundance(experiment, condition)
            if df.empty:
                continue

            for community in community_universe(df, experiment, condition):
                df_comm = df[df["community"] == community]
                for replica in cfg["window_replicas"]:
                    df_rep = df_comm[df_comm["replica"] == replica]
                    for first_day in cfg["start_days"]:
                        days = list(range(first_day, final_day + 1))
                        metrics = calc_community_cv(df_rep, days)
                        rows.append({
                            "experiment": experiment,
                            "condition": condition,
                            "community": community,
                            "replica": replica,
                            "first_day": first_day,
                            "final_day": final_day,
                            "n_days": len(days),
                            "days": ",".join(str(x) for x in days),
                            **metrics,
                        })

    cv_df = pd.DataFrame(rows)
    last3_rows = []
    for experiment, cfg in EXPERIMENT_CONFIG.items():
        first_day = cfg["last3_days"][0]
        sub = cv_df[
            (cv_df["experiment"] == experiment)
            & (cv_df["replica"] == 1)
            & (cv_df["first_day"] == first_day)
        ].copy()
        sub = sub[[
            "experiment", "condition", "community", "replica",
            "community_cv", "sum_abs_std", "temporal_bc_mean", "collapsed",
        ]].rename(columns={
            "community_cv": "classification_community_cv",
            "sum_abs_std": "classification_sum_abs_std",
            "temporal_bc_mean": "classification_temporal_bc_mean",
            "collapsed": "classification_collapsed",
        })
        last3_rows.append(sub)

    class_df = pd.concat(last3_rows, ignore_index=True)
    class_df["classification_community_cv_threshold"] = COMMUNITY_CV_THRESHOLD
    class_df["fluctuation_class"] = np.where(
        class_df["classification_community_cv"] >= COMMUNITY_CV_THRESHOLD,
        "Fluctuation",
        "Stable",
    )

    cv_df = cv_df.merge(
        class_df,
        on=["experiment", "condition", "community", "replica"],
        how="left",
    )
    cv_df["include_in_window_average"] = ~(
        (cv_df["experiment"] == "temperature")
        & (cv_df["condition"] == "W5")
    )
    return cv_df


def sample_day_is_collapsed(df_sample_day):
    if df_sample_day.empty:
        return True
    total_abs = df_sample_day["abs_abund"].sum()
    return pd.isna(total_abs) or float(total_abs) < COLLAPSE_THRESHOLD


def rel_vector_for_sample_day(df_comm, replica, day, taxa):
    df_sample = df_comm[
        (df_comm["replica"] == replica)
        & (df_comm["day"] == day)
    ]
    if sample_day_is_collapsed(df_sample):
        return np.zeros(len(taxa), dtype=float)

    vals = []
    for taxon in taxa:
        row = df_sample[df_sample["taxon"] == taxon]
        if len(row) == 0:
            vals.append(0.0)
        else:
            val = row["rel_abund"].values[0]
            vals.append(float(val) if not pd.isna(val) else 0.0)
    return np.array(vals, dtype=float)


def compute_temperature_replicate_distance(cv_df):
    cv_rows = []
    for condition in EXPERIMENT_CONFIG["temperature"]["conditions"]:
        df = read_condition_abundance("temperature", condition)
        if df.empty:
            continue
        for community in community_universe(df, "temperature", condition):
            df_comm = df[df["community"] == community]
            for replica in [1, 2, 3]:
                metrics = calc_community_cv(
                    df_comm[df_comm["replica"] == replica],
                    EXPERIMENT_CONFIG["temperature"]["last3_days"],
                )
                cv_rows.append({
                    "experiment": "temperature",
                    "condition": condition,
                    "community": community,
                    "replica": replica,
                    "community_cv_last3": metrics["community_cv"],
                    "sum_abs_std_last3": metrics["sum_abs_std"],
                    "temporal_bc_mean_last3": metrics["temporal_bc_mean"],
                    "collapsed_last3": metrics["collapsed"],
                })

    temp_cv = pd.DataFrame(cv_rows)
    temp_cv.to_csv(
        os.path.join(OUT_DIR, "temperature_last3_replicate_cv.csv"),
        index=False,
    )

    max_cv = (
        temp_cv.groupby(["experiment", "condition", "community"], as_index=False)
        .agg(
            community_cv_max_last3=("community_cv_last3", "max"),
            sum_abs_std_max_last3=("sum_abs_std_last3", "max"),
            temporal_bc_mean_max_last3=("temporal_bc_mean_last3", "max"),
            n_replicates_cv=("replica", "nunique"),
        )
    )
    max_cv["fluctuation_class"] = np.where(
        max_cv["community_cv_max_last3"] >= COMMUNITY_CV_THRESHOLD,
        "Fluctuation",
        "Stable",
    )

    rows = []
    for condition in EXPERIMENT_CONFIG["temperature"]["conditions"]:
        df = read_condition_abundance("temperature", condition)
        if df.empty:
            continue
        for community in sorted(df["community"].dropna().unique()):
            community = int(community)
            if (condition == "W5"
                    and community not in TEMPERATURE_W5_ACTIVE_COMMUNITIES):
                continue

            df_comm = df[df["community"] == community]
            taxa = sorted(df_comm["taxon"].unique())
            bc_vals = []
            for day in EXPERIMENT_CONFIG["temperature"]["last3_days"]:
                vectors = {
                    replica: rel_vector_for_sample_day(
                        df_comm, replica, day, taxa)
                    for replica in [1, 2, 3]
                }
                for rep_a, rep_b in combinations([1, 2, 3], 2):
                    bc_vals.append(calc_bray_curtis(vectors[rep_a],
                                                    vectors[rep_b]))

            rows.append({
                "experiment": "temperature",
                "condition": condition,
                "community": community,
                "mean_replicate_bc_last3": round(float(np.mean(bc_vals)), 6),
                "n_distances": len(bc_vals),
            })

    dist_df = pd.DataFrame(rows)
    dist_df = dist_df.merge(
        max_cv,
        on=["experiment", "condition", "community"],
        how="left",
    )
    return dist_df


def save_summary(cv_df):
    cv_for_pooled_average = cv_df[cv_df["include_in_window_average"]].copy()
    cv_df[~cv_df["include_in_window_average"]].to_csv(
        os.path.join(OUT_DIR, "window_cv_excluded_from_average.csv"),
        index=False,
    )

    summary = (
        cv_for_pooled_average.groupby(["experiment", "first_day",
                                       "fluctuation_class"], as_index=False)
        .agg(
            n=("community_cv", "count"),
            mean_community_cv=("community_cv", "mean"),
            sd_community_cv=("community_cv", "std"),
        )
    )
    summary["sem_community_cv"] = (
        summary["sd_community_cv"] / np.sqrt(summary["n"])
    )

    by_condition = (
        cv_df.groupby(["experiment", "condition", "first_day",
                       "fluctuation_class"], as_index=False)
        .agg(
            n=("community_cv", "count"),
            mean_community_cv=("community_cv", "mean"),
            sd_community_cv=("community_cv", "std"),
        )
    )
    by_condition["sem_community_cv"] = (
        by_condition["sd_community_cv"] / np.sqrt(by_condition["n"])
    )

    summary.to_csv(os.path.join(OUT_DIR, "window_cv_summary.csv"),
                   index=False)
    by_condition.to_csv(
        os.path.join(OUT_DIR, "window_cv_summary_by_condition.csv"),
        index=False,
    )


def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    cv_df = compute_window_cv_table()
    cv_df.to_csv(os.path.join(OUT_DIR, "community_cv_by_window_start.csv"),
                 index=False)
    save_summary(cv_df)

    dist_df = compute_temperature_replicate_distance(cv_df)
    dist_df.to_csv(
        os.path.join(OUT_DIR, "temperature_replicate_bc_vs_max_cv.csv"),
        index=False,
    )

    print(f"Output directory: {OUT_DIR}")
    print(f"community_cv_by_window_start rows: {len(cv_df)}")
    print(f"temperature replicate-distance rows: {len(dist_df)}")
    print("Temperature window starts: 1..8, ending at day 10")
    print("Mortality window starts: 1..4, ending at day 6")
    print(f"Community CV threshold: {COMMUNITY_CV_THRESHOLD}")


if __name__ == "__main__":
    main()
