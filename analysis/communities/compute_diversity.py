"""
compute_diversity.py
====================
计算每个 experiment × condition × community × replica 的多样性指标

taxon存活判定：rel_abund > REL_THRESHOLD (1%)
群落存活判定：时间窗内每天 total abs_abund 均值 >= COLLAPSE_THRESHOLD (0.05)

保留指标（两套，last_day 和 window_pooled）：
  1. richness           — 存活taxon数（rel_abund > 1% 的物种数）
  2. survival_fraction  — richness / theoretical_n（从Excel header读取）
  3. shannon            — Shannon H，基于真实 rel_abund 分布（不重归一化）
  4. simpson            — 1 - Σpᵢ²
  5. evenness           — H / ln(S)，S 为 richness（>1%的存活物种数）

last_day（→ alpha_diversity.csv）：
  基于 last_day 单天
  - richness: 该天 rel_abund > 1% 的 taxon 数
  - shannon/simpson: 该天**所有 taxon** 的 rel_abund 归一化后计算
    （不丢 < 1% 的，不在存活子集内部重新归一化）
  - evenness: shannon / ln(richness)

window_pooled（→ gamma_diversity.csv）：
  abs pooling：窗内各天 abs_abund 按 taxon 加总
  - richness: 窗内任意一天 rel_abund > 1% 的 taxon 数
  - shannon/simpson: **所有 taxon** 的 pooled abs 归一化后计算
    （不丢窗内从未 >1% 的 taxon）
  - evenness: gamma_shannon / ln(gamma_richness)
  列名加 gamma_ 前缀

设计理由（核心变更）：
  原版把 "richness 阈值过滤" 和 "概率分布构造" 耦合在一起，
  导致 Shannon/Simpson 仅基于存活子集内部重归一化，
  系统性低估了真实多样性（"差点存活"的 taxon 完全消失）。
  这版解耦：阈值只用于 richness/evenness 分母，
  shannon/simpson 用真实完整分布。

输出：
  processed/diversity/
  ├── full/
  │   ├── alpha_diversity.csv
  │   └── gamma_diversity.csv
  └── early/
      ├── alpha_diversity.csv
      └── gamma_diversity.csv
"""

import os
import re
import numpy as np
import pandas as pd
import openpyxl

BASE_DIR = os.environ.get("COMMUNITY_WORKSPACE", os.path.dirname(os.path.abspath(__file__)))
PROC_DIR = os.path.join(BASE_DIR, "processed")
OUT_DIR  = os.path.join(PROC_DIR, "diversity")

# ── 参数 ─────────────────────────────────────────────────────
REL_THRESHOLD      = 1      # taxon存活：rel_abund > 1%
COLLAPSE_THRESHOLD = 0.05   # 群落存活：mean total abs_abund >= 此值

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


# ── 理论物种数（从Excel Community header中读取 n=XX）────────
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


# ── 群落存活判定 ──────────────────────────────────────────────
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


# ── 多样性核心计算（解耦 richness 阈值 与 分布构造）──────────
def _shannon_simpson(p):
    """
    给定一个已归一化的概率向量 p（Σp=1, 元素含 0 合法），
    返回 shannon 和 simpson。
    """
    p = np.asarray(p, dtype=float)
    p_nz = p[p > 0]                   # ln(0) 按惯例 0·ln 0 = 0
    shannon = float(-np.sum(p_nz * np.log(p_nz))) if len(p_nz) > 0 else 0.0
    simpson = float(1.0 - np.sum(p ** 2))
    return shannon, simpson


def calc_diversity(weights_all, alive_mask, theoretical_n):
    """
    通用 diversity 计算（alpha / gamma 共用）。

    参数
    ----
    weights_all : 1D array
        所有 taxon 的丰度权重（alpha: 当天 rel_abund，
        gamma: pooled abs_abund）。可含 0。
    alive_mask : 1D bool array, 同长度
        哪些 taxon 算"存活"（用于 richness 计数）。
    theoretical_n : int / float
        初始理论物种数，用于 survival_fraction。

    返回
    ----
    dict with: richness, survival_fraction, shannon, simpson, evenness

    设计要点
    --------
    - richness  = sum(alive_mask)                       (阈值过滤)
    - shannon/simpson 用全部 weights_all 归一化后的分布   (不阈值过滤)
    - evenness  = shannon / ln(richness)                (S = richness)
    """
    weights_all = np.asarray(weights_all, dtype=float)
    alive_mask  = np.asarray(alive_mask,  dtype=bool)

    s = int(alive_mask.sum())
    total = float(weights_all.sum())

    if s == 0 or total <= 0:
        return dict(richness=0, survival_fraction=0.0,
                    shannon=0.0, effective_shannon=0.0,
                    simpson=0.0, evenness=0.0)

    # 真实分布，不丢任何 taxon
    p = weights_all / total
    shannon, simpson = _shannon_simpson(p)

    # evenness 分母用 richness（存活物种数）
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


# ── collapsed 群落的空记录模板 ────────────────────────────────
DIVERSITY_COLS = ["richness", "survival_fraction",
                  "shannon", "effective_shannon", "simpson", "evenness"]

def collapsed_record(base, theoretical_n, mean_total_abs, prefix=""):
    return {**base,
            **{f"{prefix}{k}": 0 if k == "richness" else 0.0
               for k in DIVERSITY_COLS},
            "theoretical_n":  theoretical_n,
            "mean_total_abs": mean_total_abs,
            "collapsed":      True}


# ── 主流程 ────────────────────────────────────────────────────
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

                        # 群落存活用时间窗判定
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

                        # ══ Alpha：last_day 单天 ══════════════
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
                            # 所有 taxon 的 rel_abund，存活 mask 用 >1%
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
                            # 每个 taxon 在窗内各天 abs_abund 加总（所有 taxon，不筛）
                            pooled_abs = (
                                df_window_valid.groupby("taxon")["abs_abund"]
                                .sum()
                                .sort_index()
                            )

                            # 存活判定：窗内**任意一天** rel_abund > 1% 的 taxon
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

        # ── 保存 ─────────────────────────────────────────────
        out_subdir = os.path.join(OUT_DIR, window_name)
        os.makedirs(out_subdir, exist_ok=True)

        pd.DataFrame(alpha_records).to_csv(
            os.path.join(out_subdir, "alpha_diversity.csv"), index=False)
        pd.DataFrame(gamma_records).to_csv(
            os.path.join(out_subdir, "gamma_diversity.csv"), index=False)
        pd.DataFrame(mean_daily_records).to_csv(
            os.path.join(out_subdir, "mean_daily_diversity.csv"), index=False)

        print(f"  alpha (last_day={last_day}):    {len(alpha_records)} 行")
        print(f"  gamma (pooled {gamma_days}): {len(gamma_records)} 行")
        print(f"  mean daily ({gamma_days}):   {len(mean_daily_records)} 行")

    print(f"\n完成！输出目录: {OUT_DIR}")
    print(f"  taxon 存活:     rel_abund > {REL_THRESHOLD}%")
    print(f"  群落存活:       mean_total_abs >= {COLLAPSE_THRESHOLD}")
    print(f"  Richness:       存活 taxon 数（阈值 {REL_THRESHOLD}%）")
    print(f"  Shannon/Simpson: 全部 taxon 真实分布（不重归一化）")
    print(f"  Evenness 分母:  ln(richness)")
    print(f"  Gamma pooling:  abs_abund 加和（biomass-weighted）")


if __name__ == "__main__":
    main()
