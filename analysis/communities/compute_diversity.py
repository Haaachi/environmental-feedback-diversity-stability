"""
Compute diversity for each experiment x condition x community x replicate.

Taxon presence: rel_abund > REL_THRESHOLD (1%).
Community persistence: mean daily total abs_abund >= COLLAPSE_THRESHOLD (0.05).

Metrics for endpoint and pooled windows:
1. richness: number of taxa exceeding 1% relative abundance.
2. survival_fraction: richness / theoretical_n, read from the Excel header.
3. shannon: Shannon H from the full abundance distribution.
4. simpson: 1 - sum(p_i ** 2).
5. evenness: H / ln(richness).

Endpoint (alpha_diversity.csv):
- Richness uses the last-day presence rule.
- Shannon/Simpson normalize all taxa on that day, including taxa below 1%.
  They do not renormalize within the present subset.
- Evenness is shannon / ln(richness).

Pooled window (gamma_diversity.csv):
- Sum each taxon's absolute abundance across the analysis window.
- Richness counts taxa above 1% on any included day.
- Shannon/Simpson normalize pooled abundance across all taxa, including taxa
  that never exceed 1%. Evenness is gamma_shannon / ln(gamma_richness).
- Pooled metric column names have a gamma_ prefix.

Presence filtering is separate from probability-distribution construction.
Earlier subset-only normalization underestimated diversity by dropping
low-abundance taxa. The threshold now affects richness and the evenness
normalizer; Shannon/Simpson use the full distribution.

Outputs: processed/diversity/{full,early}/, including alpha_diversity.csv,
gamma_diversity.csv and mean_daily_diversity.csv.
"""

import os
import re
import numpy as np
import pandas as pd
import openpyxl

BASE_DIR = os.environ.get("COMMUNITY_WORKSPACE", os.path.dirname(os.path.abspath(__file__)))
PROC_DIR = os.path.join(BASE_DIR, "processed")
OUT_DIR  = os.path.join(PROC_DIR, "diversity")

# Parameters
REL_THRESHOLD      = 1      # Taxon presence: rel_abund > 1%.
COLLAPSE_THRESHOLD = 0.05   # Community persistence: mean total abs_abund >= this value.

WINDOWS = {
    "full": {
        "experiments": {
            "temperature": {
                "last_day":       10,
                "gamma_days":     [8, 9, 10],
                "use_replicates": [1, 2, 3],
            },
        },
    },
    "early": {
        "experiments": {
            "mortality": {
                "last_day":       6,
                "gamma_days":     [4, 5, 6],
                "use_replicates": [1],
            },
            "temperature": {
                "last_day":       6,
                "gamma_days":     [4, 5, 6],
                "use_replicates": [1],
            },
        },
    },
}

def include_community(experiment, condition, community):
    # All communities in the unified Excel workbooks are part of the
    # statistical sample space. Temperature W5 communities that collapsed are
    # retained here and recorded as richness/diversity = 0 by the collapse rule.
    return True


# Initial species count: read n=XX from the Excel Community header.
def get_theoretical_n(path):
    wb   = openpyxl.load_workbook(path)
    ws   = wb[wb.sheetnames[0]]
    theo = {}
    for row in ws.iter_rows(values_only=True):
        c0 = row[0]
        if c0 and str(c0).startswith("Community"):
            m_comm = re.search(r"Community\s+(\d+)", str(c0))
            m_n    = re.search(r"n\s*=\s*(\d+)", str(c0))
            if m_comm and m_n:
                theo[int(m_comm.group(1))] = int(m_n.group(1))
    return theo


# Community persistence check
def is_community_alive(df_rep):
    mean_total_abs = df_rep.groupby("day")["abs_abund"].sum().mean()
    if pd.isna(mean_total_abs):
        return False, None
    return (mean_total_abs >= COLLAPSE_THRESHOLD,
            round(float(mean_total_abs), 6))


def day_total_abs(df_day):
    if df_day is None or len(df_day) == 0:
        return np.nan
    return float(pd.to_numeric(df_day["abs_abund"], errors="coerce").sum())


def is_day_alive(df_day):
    total = day_total_abs(df_day)
    if pd.isna(total):
        return False, np.nan
    return total >= COLLAPSE_THRESHOLD, round(total, 6)


def valid_days_by_od(df_window):
    totals = (
        df_window.groupby("day")["abs_abund"]
        .sum()
        .rename("day_total_abs")
        .reset_index()
    )
    totals["day_od_valid"] = totals["day_total_abs"] >= COLLAPSE_THRESHOLD
    return totals


# Diversity calculations with separate presence and distribution rules
def _shannon_simpson(p):
    """
    Return Shannon and Simpson diversity from a normalized probability vector p.
    Zero probabilities are valid; sum(p) = 1.
    """
    p = np.asarray(p, dtype=float)
    p_nz = p[p > 0]                   # Exclude zeros from the logarithm; the conventional value of 0 * ln(0) is zero.
    shannon = float(-np.sum(p_nz * np.log(p_nz))) if len(p_nz) > 0 else 0.0
    simpson = float(1.0 - np.sum(p ** 2))
    return shannon, simpson


def calc_diversity(weights_all, alive_mask, theoretical_n):
    """
    Shared alpha/gamma diversity calculation.

    Parameters
    ----------
    weights_all : 1D array
        All taxon abundance weights: daily rel_abund for alpha, pooled abs_abund
        for gamma. Zero weights are allowed.
    alive_mask : 1D bool array of the same length
        Taxa counted as present for richness.
    theoretical_n : int or float
        Initial species count used for survival_fraction.

    Returns
    -------
    dict
        Richness, survival fraction, Shannon, Simpson and evenness.

    Richness is sum(alive_mask). Shannon/Simpson normalize all weights and do
    not apply the presence threshold to the probability distribution.
    """
    weights_all = np.asarray(weights_all, dtype=float)
    alive_mask  = np.asarray(alive_mask,  dtype=bool)

    s = int(alive_mask.sum())
    total = float(weights_all.sum())

    if s == 0 or total <= 0:
        return dict(richness=0, survival_fraction=0.0,
                    shannon=0.0, effective_shannon=0.0,
                    simpson=0.0, evenness=0.0)

    # Use the full distribution without dropping taxa.
    p = weights_all / total
    shannon, simpson = _shannon_simpson(p)

    # Use richness, the number of present species, in the evenness denominator.
    evenness = float(shannon / np.log(s)) if s > 1 else 0.0

    surv = (round(s / theoretical_n, 4)
            if (not np.isnan(theoretical_n) and theoretical_n > 0)
            else np.nan)

    return dict(
        richness          = s,
        survival_fraction = surv,
        shannon           = round(shannon,  4),
        effective_shannon = round(float(np.exp(shannon)), 4),
        simpson           = round(simpson,  4),
        evenness          = round(evenness, 4),
    )


# Empty record template for collapsed communities
DIVERSITY_COLS = ["richness", "survival_fraction",
                  "shannon", "effective_shannon", "simpson", "evenness"]

def collapsed_record(base, theoretical_n, mean_total_abs, prefix=""):
    return {**base,
            **{f"{prefix}{k}": 0 if k == "richness" else 0.0
               for k in DIVERSITY_COLS},
            "theoretical_n":  theoretical_n,
            "mean_total_abs": mean_total_abs,
            "collapsed":      True}


# Main workflow
def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    theo_n = {
        "mortality":   get_theoretical_n(
            os.path.join(BASE_DIR, "data", "Mortality_Unified.xlsx")),
        "temperature": get_theoretical_n(
            os.path.join(BASE_DIR, "data", "Temperature_Unified.xlsx")),
    }

    for window_name, cfg in WINDOWS.items():
        print(f"\n=== {window_name} ===")

        alpha_records = []
        gamma_records = []
        mean_daily_records = []

        for experiment, exp_cfg in cfg["experiments"].items():
            last_day   = exp_cfg["last_day"]
            gamma_days = exp_cfg["gamma_days"]
            use_reps   = exp_cfg["use_replicates"]

            for cond in [f"W{i}" for i in range(1, 6)]:
                csv_path = os.path.join(PROC_DIR, experiment,
                                        cond, "abs_abundance.csv")
                if not os.path.exists(csv_path):
                    continue
                df = pd.read_csv(csv_path)

                for comm in sorted(df["community"].unique()):
                    if not include_community(experiment, cond, comm):
                        continue
                    df_comm = df[df["community"] == comm]
                    t_n     = theo_n[experiment].get(comm, np.nan)

                    for rep in use_reps:
                        base = dict(experiment=experiment, condition=cond,
                                    community=comm, replica=rep)

                        # Evaluate community persistence over the analysis window.
                        df_window = df_comm[
                            (df_comm["replica"] == rep) &
                            (df_comm["day"].isin(gamma_days))
                        ]
                        if len(df_window) == 0:
                            continue

                        alive_comm, mean_total_abs = is_community_alive(df_window)
                        day_validity = valid_days_by_od(df_window)
                        valid_days = set(
                            day_validity.loc[
                                day_validity["day_od_valid"], "day"
                            ].astype(int)
                        )
                        n_valid_days = len(valid_days)

                        # Alpha diversity: endpoint day.
                        df_rd = df_comm[
                            (df_comm["replica"] == rep) &
                            (df_comm["day"]     == last_day)
                        ]
                        alpha_day_alive, alpha_day_total_abs = is_day_alive(df_rd)

                        if not alpha_day_alive or len(df_rd) == 0:
                            alpha_records.append(
                                collapsed_record(
                                    {**base, "day": last_day},
                                    t_n, alpha_day_total_abs, prefix=""))
                        else:
                            # Use all relative abundances; the presence mask uses >1%.
                            rel_vals   = df_rd["rel_abund"].values
                            alive_mask = rel_vals > REL_THRESHOLD

                            div = calc_diversity(rel_vals, alive_mask, t_n)
                            alpha_records.append({
                                **base,
                                "day":            last_day,
                                **div,
                                "theoretical_n":  t_n,
                                "mean_total_abs": alpha_day_total_abs,
                                "collapsed":      False,
                            })

                        # Mean daily diversity across the final/gamma window.
                        daily_divs = []
                        for day_i in gamma_days:
                            df_day = df_window[df_window["day"] == day_i]
                            if len(df_day) == 0:
                                continue
                            day_alive, _day_total_abs = is_day_alive(df_day)
                            if not day_alive:
                                daily_divs.append(
                                    calc_diversity([], [], t_n)
                                )
                                continue
                            rel_vals_day = df_day["rel_abund"].values
                            alive_mask_day = rel_vals_day > REL_THRESHOLD
                            daily_divs.append(
                                calc_diversity(rel_vals_day, alive_mask_day, t_n)
                            )

                        if len(daily_divs) == 0:
                            mean_daily_records.append(
                                collapsed_record(
                                    {**base, "days": str(gamma_days)},
                                    t_n, mean_total_abs,
                                    prefix="mean_daily_"))
                        else:
                            mean_vals = {}
                            for metric in DIVERSITY_COLS:
                                vals = [
                                    d[metric] for d in daily_divs
                                    if not pd.isna(d[metric])
                                ]
                                mean_vals[f"mean_daily_{metric}"] = (
                                    round(float(np.mean(vals)), 4)
                                    if len(vals) > 0 else np.nan
                                )
                            mean_daily_records.append({
                                **base,
                                "days":           str(gamma_days),
                                **mean_vals,
                                "theoretical_n":  t_n,
                                "mean_total_abs": mean_total_abs,
                                "n_valid_days":   n_valid_days,
                                "collapsed":      n_valid_days == 0,
                            })

                        # Gamma: pool only OD-valid days so collapsed-day
                        # sequencing noise cannot create artificial diversity.
                        df_window_valid = df_window[
                            df_window["day"].isin(valid_days)
                        ]
                        if len(df_window_valid) == 0:
                            gamma_records.append(
                                collapsed_record(
                                    {**base, "days": str(gamma_days)},
                                    t_n, mean_total_abs, prefix="gamma_"))
                        else:
                            # Sum each taxon's abs_abund over the window, including all taxa.
                            pooled_abs = (
                                df_window_valid.groupby("taxon")["abs_abund"]
                                .sum()
                                .sort_index()
                            )

                            # Presence: rel_abund > 1% on any included day.
                            alive_taxa = set(
                                df_window_valid.loc[
                                    df_window_valid["rel_abund"] > REL_THRESHOLD,
                                    "taxon"
                                ]
                            )
                            alive_mask = pooled_abs.index.isin(alive_taxa)

                            div = calc_diversity(
                                pooled_abs.values,
                                alive_mask,
                                t_n,
                            )
                            gamma_records.append({
                                **base,
                                "days":           str(gamma_days),
                                **{f"gamma_{k}": v for k, v in div.items()},
                                "theoretical_n":  t_n,
                                "mean_total_abs": mean_total_abs,
                                "n_valid_days":   n_valid_days,
                                "collapsed":      False,
                            })

        # Save outputs.
        out_subdir = os.path.join(OUT_DIR, window_name)
        os.makedirs(out_subdir, exist_ok=True)

        pd.DataFrame(alpha_records).to_csv(
            os.path.join(out_subdir, "alpha_diversity.csv"), index=False)
        pd.DataFrame(gamma_records).to_csv(
            os.path.join(out_subdir, "gamma_diversity.csv"), index=False)
        pd.DataFrame(mean_daily_records).to_csv(
            os.path.join(out_subdir, "mean_daily_diversity.csv"), index=False)

        print(f"  alpha (last_day={last_day}):    {len(alpha_records)} rows")
        print(f"  gamma (pooled {gamma_days}): {len(gamma_records)} rows")
        print(f"  mean daily ({gamma_days}):   {len(mean_daily_records)} rows")

    print(f"\nDone! Output directory: {OUT_DIR}")
    print(f"  Taxon presence:     rel_abund > {REL_THRESHOLD}%")
    print(f"  Community persistence:       mean_total_abs >= {COLLAPSE_THRESHOLD}")
    print(f"  Richness:       number of present taxa (threshold {REL_THRESHOLD}%)")
    print(f"  Shannon/Simpson: full taxon distribution (no present-subset renormalization)")
    print(f"  Evenness denominator:  ln(richness)")
    print(f"  Gamma pooling:  summed abs_abund (biomass-weighted)")


if __name__ == "__main__":
    main()
