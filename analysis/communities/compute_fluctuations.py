"""
compute_fluctuations.py
=======================
璁＄畻姣忎釜 experiment 脳 condition 脳 community 脳 replica 鐨勯渿鑽℃寚鏍?
淇濈暀涓変釜鎸囨爣锛?  1. community_cv
       = 危岬?蟽岬abs / mean_total_abs
       鍒嗗瓙锛氭牳蹇僼axon鐨刟bs_std涔嬪拰锛堝叏閮ㄦ湁鏁堢偣鍚?璁＄畻std锛?       鍒嗘瘝锛氭椂闂寸獥鍐呯兢钀芥€荤敓鐗╅噺鍧囧€硷紙= mean OD锛?       鐗╃悊鎰忎箟锛氱兢钀芥暣浣撴尝鍔ㄩ噺 / 缇よ惤骞冲潎鐢熺墿閲?
  2. total_biomass_cv
       = std(OD_t) / mean(OD_t)锛屼粎鐢ㄦ椂闂寸獥鍐呯殑OD
       绾疧D灞傞潰鐨勭敓鐗╅噺娉㈠姩锛屼笉渚濊禆娴嬪簭鏁版嵁

  3. temporal_bc_mean
       鏃堕棿绐楀唴鎵€鏈夊ぉ涓や袱涔嬮棿 Bray-Curtis 璺濈鐨勫潎鍊?       锛坋.g. Day 8/9/10 鈫?3瀵癸細8vs9, 8vs10, 9vs10锛?       鍙叧娉ㄧ兢钀界粨鏋勬湰韬殑闇囪崱锛屼笌缁濆OD鏃犲叧
       涓板害鍚戦噺锛氭瘡澶╁瓨娲籺axon鐨?rel_abund锛堟浜＄疆0锛夛紝
       涓嶄娇鐢?abs_abund锛岀‘淇濅笌OD閲忕骇瑙ｈ€?
鏍稿績璁捐锛氶噸褰掍竴鍖?鈫?娑堥櫎娴嬪簭鎵规鏁堝簲
  Step 1  鍣煶杩囨护锛圢OISE_THRESHOLD=0锛岄粯璁や笉杩囨护锛?  Step 2  淇濈暀taxon姣忓ぉ rel_abund 閲嶅綊涓€鍖栧埌100%
          abs_abund_clean = rel_abund_clean / 100 脳 OD
  Step 3  璁＄畻浠ヤ笂涓変釜鎸囨爣

宕╂簝鍒ゅ畾锛氭椂闂寸獥鍐呮瘡澶?total abs_abund 鍧囧€?< COLLAPSE_THRESHOLD (0.05)
  鈫?鎵€鏈夋寚鏍囩疆0锛宑ollapsed = True

杈撳嚭锛堜繚鎸佷笌鍘熺増鐩稿悓璺緞缁撴瀯锛夛細
  processed/fluctuations/
  鈹溾攢鈹€ full/
  鈹?  鈹溾攢鈹€ taxon_level.csv
  鈹?  鈹斺攢鈹€ community_level.csv
  鈹斺攢鈹€ early/
      鈹溾攢鈹€ taxon_level.csv
      鈹斺攢鈹€ community_level.csv
"""

import os
import numpy as np
import pandas as pd
from itertools import combinations

BASE_DIR = os.environ.get("COMMUNITY_WORKSPACE", os.path.dirname(os.path.abspath(__file__)))
PROC_DIR = os.path.join(BASE_DIR, "processed")
OUT_DIR  = os.path.join(PROC_DIR, "fluctuations")

# 鈹€鈹€ 鍙傛暟 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
COLLAPSE_THRESHOLD = 0.05  # mean total abs_abund 浣庝簬姝ゅ€艰涓篶ollapsed
COMMUNITY_CV_THRESHOLD = 0.265
LOG_PSEUDOCOUNT_FRAC = 1e-4
LOG_PSEUDOCOUNT_FLOOR = 1e-12

TIME_WINDOWS = {
    "full": {
        "temperature": {"days": [8, 9, 10], "replicas": [1, 2, 3]},
    },
    "early": {
        "mortality": {"days": [4, 5, 6], "replicas": [1]},
        "temperature": {"days": [4, 5, 6], "replicas": [1]},
    },
    "last4": {
        "mortality": {"days": [3, 4, 5, 6], "replicas": [1]},
        "temperature": {"days": [7, 8, 9, 10], "replicas": [1]},
    },
}

def include_community(experiment, condition, community):
    # All communities in the unified Excel workbooks are part of the
    # statistical sample space. Temperature W5 communities that collapsed are
    # retained here and recorded as stable collapsed communities (CV = 0).
    return True

def force_missing_as_collapsed(experiment, condition):
    # Temperature W5 has collapsed communities whose sequencing profiles are
    # not biologically meaningful. They still belong in denominators as stable
    # collapsed communities.
    return experiment == "temperature" and condition == "W5"

def condition_community_universe(df, experiment, condition):
    communities = set(int(x) for x in df["community"].dropna().unique())
    if force_missing_as_collapsed(experiment, condition):
        communities.update(range(1, 13))
    return sorted(communities)

# 鈹€鈹€ collapsed缇よ惤鎸囨爣妯℃澘 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
def make_collapsed_record(base, mean_total_abs):
    return {
        **base,
        "community_cv":      0.0,
        "community_cv_rel":  0.0,
        "sum_abs_std":       0.0,
        "sum_rel_std":       0.0,
        "weighted_log_sd":   0.0,
        "log_pseudocount":   np.nan,
        "total_biomass_cv":  0.0,
        "total_biomass_std": 0.0,
        "temporal_bc_mean":  0.0,
        "temporal_bc_max":   0.0,
        "n_taxa_clean":      0,
        "mean_total_abs":    (round(float(mean_total_abs), 6)
                              if mean_total_abs is not None
                              and not np.isnan(mean_total_abs)
                              else None),
        "n_days_valid":      0,
        "collapsed":         True,
    }

# 鈹€鈹€ Bray-Curtis璺濈 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
def calc_bray_curtis(p1, p2):
    p1    = np.array(p1, dtype=float)
    p2    = np.array(p2, dtype=float)
    denom = p1.sum() + p2.sum()
    if denom == 0:
        return np.nan
    return float(np.sum(np.abs(p1 - p2)) / denom)

def log_pseudocount_for(mean_total_abs):
    if mean_total_abs is None or pd.isna(mean_total_abs) or mean_total_abs <= 0:
        return LOG_PSEUDOCOUNT_FLOOR
    return max(float(mean_total_abs) * LOG_PSEUDOCOUNT_FRAC,
               LOG_PSEUDOCOUNT_FLOOR)

def weighted_log_sd_from_taxa(df_taxon):
    if df_taxon.empty or "log_abs_sd" not in df_taxon.columns:
        return np.nan

    mean_abs = (
        df_taxon["abs_mean"]
        .fillna(0.0)
        .clip(lower=0.0)
        .to_numpy(dtype=float)
    )
    total_mean_abs = float(np.sum(mean_abs))
    if total_mean_abs <= 0:
        df_taxon["mean_abs_weight"] = 0.0
        return np.nan

    weights = mean_abs / total_mean_abs
    df_taxon["mean_abs_weight"] = weights

    log_sd = df_taxon["log_abs_sd"].to_numpy(dtype=float)
    valid = np.isfinite(weights) & np.isfinite(log_sd)
    valid_weight = float(np.sum(weights[valid]))
    if valid_weight <= 0:
        return np.nan

    return float(np.sum(weights[valid] * log_sd[valid]) / valid_weight)

# 鈹€鈹€ taxon灞傞潰鎸囨爣锛堜粎淇濈暀 community_cv 鎵€闇€鐨?abs_std锛夆攢鈹€鈹€鈹€鈹€鈹€鈹€鈹€
def compute_taxon_metrics(vals_abs_clean, log_pseudocount):
    """
    vals_abs_clean锛氶噸褰掍竴鍖栧悗鐨刟bs涓板害鏃堕棿搴忓垪锛堝惈0锛岀己澶卞ぉ涓簄an锛?    杩斿洖锛歛bs_std锛堝惈0璁＄畻锛夊拰 abs_mean锛堝惈0璁＄畻锛?    鐢ㄤ簬 community_cv 鐨勫垎瀛愮疮鍔?    """
    vals = np.array(vals_abs_clean, dtype=float)
    valid = vals[~np.isnan(vals)]   # 淇濈暀0
    if len(valid) >= 2:
        abs_std  = float(np.std(valid, ddof=1))
        abs_mean = float(np.mean(valid))
        log_abs_sd = float(
            np.std(np.log(np.maximum(valid, 0.0) + log_pseudocount),
                   ddof=1)
        )
    else:
        abs_std  = np.nan
        abs_mean = float(np.nanmean(vals)) if not np.all(np.isnan(vals)) else np.nan
        log_abs_sd = np.nan
    return {"abs_std": abs_std, "abs_mean": abs_mean,
            "log_abs_sd": log_abs_sd}

# 鈹€鈹€ 涓绘祦绋?鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    for window_name, exp_cfg in TIME_WINDOWS.items():
        days = "configured per experiment"
        print(f"\n=== 鏃堕棿绐楀彛: {window_name} (Day {days}) ===")

        taxon_records     = []
        community_records = []

        for experiment, cfg in exp_cfg.items():
            days = cfg["days"]
            use_reps = set(cfg["replicas"])
            print(f"  {experiment}: days={days}, replicas={sorted(use_reps)}")

            for cond in [f"W{i}" for i in range(1, 6)]:
                csv_path = os.path.join(PROC_DIR, experiment,
                                        cond, "abs_abundance.csv")
                if not os.path.exists(csv_path):
                    continue

                df = pd.read_csv(csv_path)

                for comm in condition_community_universe(df, experiment, cond):
                    if not include_community(experiment, cond, comm):
                        continue
                    df_comm = df[df["community"] == comm]

                    if force_missing_as_collapsed(experiment, cond):
                        reps_to_check = sorted(use_reps)
                    else:
                        reps_to_check = sorted(
                            int(rep) for rep in df_comm["replica"].unique()
                            if int(rep) in use_reps
                        )

                    for rep in reps_to_check:
                        df_rep = df_comm[
                            (df_comm["replica"] == rep) &
                            (df_comm["day"].isin(days))
                        ]

                        base = dict(
                            experiment=experiment,
                            condition=cond,
                            community=comm,
                            replica=rep,
                            n_days=len(days),
                        )

                        if len(df_rep) == 0:
                            if force_missing_as_collapsed(experiment, cond):
                                community_records.append(
                                    make_collapsed_record(base, 0.0))
                                print(f"  [collapsed-missing] {experiment} "
                                      f"{cond} C{comm} R{rep}")
                            continue

                        # 鈹€鈹€ 宕╂簝妫€鏌?鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
                        mean_total_abs = (
                            df_rep.groupby("day")["abs_abund"]
                            .sum().mean()
                        )
                        if pd.isna(mean_total_abs) or \
                                mean_total_abs < COLLAPSE_THRESHOLD:
                            community_records.append(
                                make_collapsed_record(base, mean_total_abs))
                            print(f"  [collapsed] {experiment} {cond} "
                                  f"C{comm} R{rep} "
                                  f"mean_total_abs="
                                  f"{mean_total_abs:.4f}" if not pd.isna(
                                      mean_total_abs) else
                                  f"  [collapsed] {experiment} {cond} "
                                  f"C{comm} R{rep} mean_total_abs=NaN")
                            continue

                        n_days_valid = int(
                            df_rep.groupby("day")["abs_abund"]
                            .sum().gt(0).sum()
                        )

                        # 鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲
                        # 鎸囨爣2锛歵otal_biomass_cv
                        # 鐩存帴鐢ㄦ椂闂寸獥鍐呮瘡澶╃殑 total OD
                        # 鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲
                        od_by_day = (
                            df_rep.groupby("day")["OD"].mean()
                        )
                        od_vals = od_by_day.dropna().values
                        if len(od_vals) >= 2 and np.mean(od_vals) > 0:
                            total_biomass_std = round(
                                float(np.std(od_vals, ddof=1)), 6)
                            total_biomass_cv = round(
                                float(total_biomass_std / np.mean(od_vals)), 6)
                        else:
                            total_biomass_std = np.nan
                            total_biomass_cv = np.nan

                        # 鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲
                        # 鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲
                        core_taxa = sorted(df_rep["taxon"].unique())

                        if len(core_taxa) == 0:
                            log_pseudocount = log_pseudocount_for(mean_total_abs)
                            community_records.append({
                                **base,
                                "community_cv":     np.nan,
                                "community_cv_rel": np.nan,
                                "sum_abs_std":      np.nan,
                                "sum_rel_std":      np.nan,
                                "weighted_log_sd":  np.nan,
                                "log_pseudocount":  log_pseudocount,
                                "total_biomass_cv": total_biomass_cv,
                                "total_biomass_std": total_biomass_std,
                                "temporal_bc_mean": np.nan,
                                "temporal_bc_max":  np.nan,
                                "n_taxa_clean":     0,
                                "mean_total_abs":   round(float(mean_total_abs), 6),
                                "n_days_valid":     n_days_valid,
                                "collapsed":        False,
                            })
                            continue

                        df_core = df_rep[
                            df_rep["taxon"].isin(core_taxa)
                        ].copy()

                        # 鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲
                        # Step 2锛氶噸鏂板綊涓€鍖?                        # 鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲
                        day_totals = (
                            df_core.groupby("day")["rel_abund"]
                            .sum().rename("day_total")
                        )
                        df_core = df_core.merge(
                            day_totals, on="day", how="left")

                        df_core["rel_abund_clean"] = np.where(
                            df_core["day_total"] > 0,
                            df_core["rel_abund"]
                            / df_core["day_total"] * 100,
                            np.nan
                        )
                        df_core["abs_abund_clean"] = np.where(
                            df_core["rel_abund_clean"].isna() |
                            df_core["OD"].isna(),
                            np.nan,
                            df_core["rel_abund_clean"]
                            / 100.0 * df_core["OD"]
                        )

                        # 鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲
                        # Step 3a锛歵axon灞傞潰 鈫?community_cv 鍒嗗瓙
                        # 鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲
                        taxon_rows = []
                        log_pseudocount = log_pseudocount_for(mean_total_abs)
                        for taxon in core_taxa:
                            df_t = (df_core[df_core["taxon"] == taxon]
                                    .sort_values("day"))

                            vals_abs_clean = []
                            for d in days:
                                row = df_t[df_t["day"] == d]
                                vals_abs_clean.append(
                                    row["abs_abund_clean"].values[0]
                                    if len(row) > 0 else np.nan
                                )

                            # 闈為浂鏈夋晥鐐?< 2 鈫?璺宠繃
                            m = compute_taxon_metrics(vals_abs_clean,
                                                      log_pseudocount)
                            row_data = {
                                **base,
                                "taxon":    taxon,
                                "abs_std":  m["abs_std"],
                                "abs_mean": m["abs_mean"],
                                "log_abs_sd": m["log_abs_sd"],
                            }
                            taxon_rows.append(row_data)

                        # 鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲
                        # community_cv = 危蟽岬abs / mean_total_abs
                        # 鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲
                        if len(taxon_rows) > 0:
                            df_taxon   = pd.DataFrame(taxon_rows)
                            weighted_log_sd = round(
                                weighted_log_sd_from_taxa(df_taxon), 6)
                            taxon_records.extend(
                                df_taxon.to_dict(orient="records"))
                            valid_stds = df_taxon["abs_std"].dropna()
                            sum_abs_std = (
                                round(float(valid_stds.sum()), 6)
                                if len(valid_stds) > 0
                                else np.nan
                            )
                            community_cv = (
                                round(float(sum_abs_std / mean_total_abs), 6)
                                if len(valid_stds) > 0 and mean_total_abs > 0
                                else np.nan
                            )
                            n_taxa_clean = int(len(taxon_rows))
                        else:
                            sum_abs_std = np.nan
                            community_cv = np.nan
                            weighted_log_sd = np.nan
                            n_taxa_clean = 0

                        rel_matrix = (
                            df_core.pivot_table(
                                index="day",
                                columns="taxon",
                                values="rel_abund_clean",
                                aggfunc="sum",
                            )
                            .reindex(index=days, columns=core_taxa)
                            .fillna(0.0)
                        )
                        if len(days) >= 2 and rel_matrix.shape[1] > 0:
                            rel_stds = rel_matrix.std(axis=0, ddof=1)
                            sum_rel_std = round(float(rel_stds.sum()), 6)
                            community_cv_rel = round(
                                float(sum_rel_std / 100.0), 6)
                        else:
                            sum_rel_std = np.nan
                            community_cv_rel = np.nan

                        # 鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲
                        # Step 3b锛歵emporal_bc_mean
                        # 姣忓ぉ鐢?rel_abund锛堝瓨娲荤疆鍘熷€硷紝姝讳骸缃?锛夋瀯寤哄悜閲?                        # 鎵€鏈夊ぉ涓や袱BC璺濈鐨勫潎鍊?                        # 鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲
                        all_taxa_sorted = sorted(df_rep["taxon"].unique())
                        day_vecs = {}
                        for d in days:
                            df_d = df_rep[df_rep["day"] == d]
                            if (df_d.empty or
                                    df_d["abs_abund"].sum() < COLLAPSE_THRESHOLD):
                                day_vecs[d] = np.zeros(len(all_taxa_sorted))
                                continue

                            vec = []
                            for t in all_taxa_sorted:
                                row = df_d[df_d["taxon"] == t]
                                if len(row) == 0:
                                    vec.append(0.0)
                                else:
                                    r = row["rel_abund"].values[0]
                                    vec.append(float(r) if not pd.isna(r) else 0.0)
                            day_vecs[d] = np.array(vec)

                        bc_vals = []
                        for d1, d2 in combinations(days, 2):
                            bc = calc_bray_curtis(day_vecs[d1], day_vecs[d2])
                            if not np.isnan(bc):
                                bc_vals.append(bc)

                        temporal_bc_mean = (
                            round(float(np.mean(bc_vals)), 6)
                            if len(bc_vals) > 0 else np.nan
                        )
                        temporal_bc_max = (
                            round(float(np.max(bc_vals)), 6)
                            if len(bc_vals) > 0 else np.nan
                        )

                        # 鈹€鈹€ 姹囨€?鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
                        community_records.append({
                            **base,
                            "community_cv":     community_cv,
                            "community_cv_rel": community_cv_rel,
                            "sum_abs_std":      sum_abs_std,
                            "sum_rel_std":      sum_rel_std,
                            "weighted_log_sd":  weighted_log_sd,
                            "log_pseudocount":  log_pseudocount,
                            "total_biomass_cv": total_biomass_cv,
                            "total_biomass_std": total_biomass_std,
                            "temporal_bc_mean": temporal_bc_mean,
                            "temporal_bc_max":  temporal_bc_max,
                            "n_taxa_clean":     n_taxa_clean,
                            "mean_total_abs":   round(float(mean_total_abs), 6),
                            "n_days_valid":     n_days_valid,
                            "collapsed":        False,
                        })

        # 鈹€鈹€ 淇濆瓨锛堣矾寰勪笌鍘熺増瀹屽叏涓€鑷达級鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
        out_subdir = os.path.join(OUT_DIR, window_name)
        os.makedirs(out_subdir, exist_ok=True)

        df_taxon_out = pd.DataFrame(taxon_records)
        df_comm_out  = pd.DataFrame(community_records)

        if "collapsed" not in df_comm_out.columns:
            df_comm_out["collapsed"] = False
        df_comm_out["community_cv_threshold"] = COMMUNITY_CV_THRESHOLD
        df_comm_out["fluctuating_by_community_cv"] = (
            pd.to_numeric(df_comm_out["community_cv"], errors="coerce")
            >= COMMUNITY_CV_THRESHOLD
        ) & (~df_comm_out["collapsed"].astype(bool))
        df_comm_out["fluctuation_class"] = np.where(
            df_comm_out["fluctuating_by_community_cv"],
            "Fluctuation",
            "Stable",
        )

        # taxon_level 鍙繚鐣欏繀瑕佸垪
        taxon_cols = ["experiment", "condition", "community", "replica",
                      "n_days", "taxon", "abs_std", "abs_mean",
                      "log_abs_sd", "mean_abs_weight"]
        df_taxon_out = df_taxon_out[
            [c for c in taxon_cols if c in df_taxon_out.columns]
        ]

        df_taxon_out.to_csv(
            os.path.join(out_subdir, "taxon_level.csv"), index=False)
        df_comm_out.to_csv(
            os.path.join(out_subdir, "community_level.csv"), index=False)

        n_collapsed = int(df_comm_out["collapsed"].sum())
        n_total     = len(df_comm_out)
        n_active    = n_total - n_collapsed
        print(f"  taxon_level:     {len(df_taxon_out)} rows")
        print(
            f"  community_level: {n_total} rows "
            f"(collapsed: {n_collapsed}, active: {n_active})"
        )

        active = df_comm_out[~df_comm_out["collapsed"]]
        if len(active) > 0:
            for col in ["community_cv", "sum_abs_std", "total_biomass_cv",
                        "community_cv_rel", "sum_rel_std",
                        "total_biomass_std", "temporal_bc_mean",
                        "temporal_bc_max",
                        "weighted_log_sd"]:
                if col in active.columns:
                    print(f"\n  [{col}] summary:")
                    print(active[col].describe().round(4).to_string())

    print(f"\nDone. Output directory: {OUT_DIR}")
    print("  Taxa included:        all retained synthetic-community taxa")
    print("  Renormalization:      retained taxa are renormalized to 100% each day")
    print(f"  Collapse threshold:   mean_total_abs < {COLLAPSE_THRESHOLD}")
    print(f"  Fluctuation class:    community_cv >= {COMMUNITY_CV_THRESHOLD}")
    print("  community_cv:         sum_sigma_i_abs / mean_total_abs")
    print("  community_cv_rel:     sum_i sd(relative_abundance_i_clean) / 100")
    print("  sum_abs_std:          sum_sigma_i_abs")
    print("  sum_rel_std:          sum_i sd(relative_abundance_i_clean)")
    print("  weighted_log_sd:      sum_i w_i * sd(log(abs_i + eps))")
    print(f"  log pseudocount:      max(mean_total_abs * "
          f"{LOG_PSEUDOCOUNT_FRAC}, {LOG_PSEUDOCOUNT_FLOOR})")
    print("  total_biomass_cv:     std(OD_t) / mean(OD_t)")
    print("  total_biomass_std:    std(OD_t)")
    print("  temporal_bc_mean:     pairwise Bray-Curtis on all retained species")
    print("  temporal_bc_max:      max pairwise Bray-Curtis within the time window")

if __name__ == "__main__":
    main()
