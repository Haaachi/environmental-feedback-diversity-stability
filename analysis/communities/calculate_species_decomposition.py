#!/usr/bin/env python3
# ===========================================================
# calculate_species_decomposition.py
# 主目录: /home/hachi/Tem_mortality_workspace
#
# 目的
# ----
# 在两个窗口内计算三套互补指标：
#
# A. species correlation / variance decomposition
#    基于连续绝对丰度 abs_abund，不设 presence 阈值。
#    collapse day 保留，因为它是真实 total biomass fluctuation。
#
#    B(t) = Σ_i X_i(t)
#
#    Var(B) = Σ_i Var(X_i) + 2Σ_{i<j} Cov(X_i, X_j)
#
#    synchrony_phi = Var(B) / (Σ_i sd(X_i))^2
#                  = (total_biomass_cv / community_cv)^2
#
#
# B. within-window turnover analysis
#    基于 rel_abund >= 1% 定义 presence。
#
#    但为了避免 low-biomass / collapsed day 的测序噪音被相对丰度
#    放大，presence 判据改为：
#
#      I_i(t) = 1 if rel_abund_i(t) >= 1%
#                    and total_biomass(t) >= MIN_TOTAL_BIOMASS_FOR_PRESENCE
#               else 0
#
#    注意：
#      - 不删除 collapsed day。
#      - collapsed day 保留在时间序列中。
#      - collapsed day 的 presence vector 全部设为 0。
#      - 因此 collapse / recovery 仍然会贡献 disappearance / appearance。
#
#
# C. biomass-based Shannon cumulative diversity
#    基于绝对 biomass，不使用 rel_abund。
#    collapse day 保留。
#
#    每日 Shannon:
#      p_i(t) = X_i(t) / Σ_j X_j(t)
#      H(t) = -Σ_i p_i(t) log(p_i(t))
#
#    累积 Shannon:
#      A_i = Σ_t X_i(t)
#      p_i_cum = A_i / Σ_j A_j
#      H_cum = -Σ_i p_i_cum log(p_i_cum)
#
#    cumulative_excess_shannon = H_cum - mean_t H(t)
#
#
# 特殊规则
# --------
# temperature W5:
#   只保留 community 3, 7, 12。
#   其他 community 整体排除，因为在 50°C 下已 collapse，
#   rel_abund 主要反映测序噪音。
#
# 但对于 community 3, 7, 12:
#   所有天都保留，包括中间真实 collapse 的 day。
#
#
# 时间窗口
# --------
# 1) last3
#    temperature: day 8-10, R1/R2/R3
#    mortality:   day 4-6,  R1 only
#
# 2) long
#    temperature: day 3-10, R1 only
#    mortality:   day 3-6,  R1 only
#
#
# 输入
# ----
# processed/{experiment}/W*/abs_abundance.csv
#
#
# 输出
# ----
# processed/species_decomposition/
#   last3/
#     community_decomposition.csv
#     pairwise_correlations.csv
#   long/
#     community_decomposition.csv
#     pairwise_correlations.csv
# ===========================================================

import os
import glob
import itertools
import numpy as np
import pandas as pd


BASE_DIR = os.environ.get("COMMUNITY_WORKSPACE", os.path.dirname(os.path.abspath(__file__)))
PROC_DIR = os.path.join(BASE_DIR, "processed")
OUT_DIR  = os.path.join(PROC_DIR, "species_decomposition")

# rel_abund 是百分比，因此 1.0 表示 1%
PRESENCE_THRESHOLD_REL = 1.0

# 只用于 presence/absence turnover。
# 不用于过滤 abs_abund decomposition。
# 如果某天 total biomass < 该阈值，则该天 presence vector 全部设为 0。
MIN_TOTAL_BIOMASS_FOR_PRESENCE = 0.01
MIN_TOTAL_BIOMASS_FOR_DIVERSITY = 0.05
COMMUNITY_CV_FLUCTUATING_THRESHOLD = 0.25
PHI_MIN_COMMUNITY_CV_FOR_ROBUST = 0.05
PHI_MIN_SUM_SPECIES_SD = 1e-9

# temperature W5 只分析这些 community
VALID_COMMUNITIES_BY_CONDITION = {
    ("temperature", "W5"): {3, 7, 12},
}


# ============================================================
# 时间窗口定义
# ============================================================

WINDOW_CONFIG = {
    "last3": {
        "temperature": {
            "days": [8, 9, 10],
            "replicas": [1, 2, 3],
        },
        "mortality": {
            "days": [4, 5, 6],
            "replicas": [1],
        },
    },
    "last4": {
        "temperature": {
            "days": [7, 8, 9, 10],
            "replicas": [1],
        },
        "mortality": {
            "days": [3, 4, 5, 6],
            "replicas": [1],
        },
    },
    "long": {
        "temperature": {
            "days": list(range(3, 11)),
            "replicas": [1],
        },
        "mortality": {
            "days": [3, 4, 5, 6],
            "replicas": [1],
        },
    },
}


# ============================================================
# 工具函数
# ============================================================

def safe_divide(a, b):
    if b is None or pd.isna(b) or abs(b) < 1e-15:
        return np.nan
    return a / b


def apply_special_community_filter(df):
    """
    应用特殊 community-level 规则。

    当前规则：
      temperature W5 只保留 community 3, 7, 12。
    """
    if df.empty:
        return df

    keep_mask = pd.Series(True, index=df.index)

    for (experiment, condition), valid_comms in VALID_COMMUNITIES_BY_CONDITION.items():
        special_mask = (
            (df["experiment"] == experiment)
            & (df["condition"] == condition)
        )

        keep_mask.loc[special_mask] = df.loc[special_mask, "community"].isin(valid_comms)

    return df.loc[keep_mask].copy()


def load_abs_abundance():
    """
    读取 processed/{experiment}/W*/abs_abundance.csv
    并应用 community-level 特殊过滤。
    """
    all_dfs = []

    for experiment in ["mortality", "temperature"]:
        pattern = os.path.join(PROC_DIR, experiment, "W*", "abs_abundance.csv")
        files = sorted(glob.glob(pattern))

        if not files:
            print(f"[warn] No files found for {experiment}: {pattern}")
            continue

        for path in files:
            df = pd.read_csv(path)

            required = {
                "experiment", "condition", "community",
                "taxon", "replica", "day", "abs_abund", "rel_abund"
            }
            missing = required - set(df.columns)
            if missing:
                raise ValueError(
                    f"{path} is missing required columns: {sorted(missing)}"
                )

            all_dfs.append(df)

    if not all_dfs:
        raise FileNotFoundError(
            f"No abs_abundance.csv files found under {PROC_DIR}/mortality "
            f"or {PROC_DIR}/temperature."
        )

    df_all = pd.concat(all_dfs, ignore_index=True)

    df_all["experiment"] = df_all["experiment"].astype(str)
    df_all["condition"]  = df_all["condition"].astype(str)
    df_all["community"]  = df_all["community"].astype(int)
    df_all["replica"]    = df_all["replica"].astype(int)
    df_all["day"]        = df_all["day"].astype(int)
    df_all["taxon"]      = df_all["taxon"].astype(str)

    df_all["abs_abund"] = pd.to_numeric(df_all["abs_abund"], errors="coerce")
    df_all["rel_abund"] = pd.to_numeric(df_all["rel_abund"], errors="coerce")

    n_before = len(df_all)
    df_all = apply_special_community_filter(df_all)
    n_after = len(df_all)

    print(
        f"  -> special community filter applied: "
        f"{n_before} rows -> {n_after} rows"
    )
    print("     temperature W5 kept communities: 3, 7, 12")

    return df_all


def make_abundance_matrix(df_group, expected_days, value_col):
    """
    将一个 experiment × condition × community × replica 的数据转成:
        rows    = day
        columns = taxon
        values  = value_col

    对 abs_abund:
      - 用于 variance/covariance decomposition 和 Shannon。
      - 有数据的 day 中，taxon NA 视为 0。
      - 整个窗口内总和为 0 的 taxon 会去掉。
      - 不根据 total biomass 删除 day。

    对 rel_abund:
      - 用于 presence/absence turnover。
      - 有数据的 day 中，taxon NA 视为 0。
      - 不主动去掉低丰度 taxon，后续用 1% 阈值判断 presence。
    """
    df = df_group[df_group["day"].isin(expected_days)].copy()

    if df.empty:
        return None

    mat = (
        df.pivot_table(
            index="day",
            columns="taxon",
            values=value_col,
            aggfunc="sum",
            dropna=False,
        )
        .reindex(expected_days)
    )

    # 去掉完全没有任何记录的 day。
    # 注意：这不是 low-biomass filtering。
    # 如果某天有记录但 total biomass 接近 0，该 day 会保留。
    valid_day_mask = mat.notna().any(axis=1)
    mat = mat.loc[valid_day_mask]

    if mat.shape[0] < 2:
        return None

    # 在有效 day 中，缺失 taxon 视为 0
    mat = mat.fillna(0.0)

    if value_col == "abs_abund":
        nonzero_species = mat.sum(axis=0) > 0
        mat = mat.loc[:, nonzero_species]

        if mat.shape[1] == 0:
            return None

    return mat


def calc_bray_curtis(vec_a, vec_b):
    a = np.array(vec_a, dtype=float)
    b = np.array(vec_b, dtype=float)
    denom = a.sum() + b.sum()
    if denom <= 0:
        return 0.0
    return float(np.abs(a - b).sum() / denom)


def compute_temporal_bc_mean(rel_mat, abs_mat=None,
                             biomass_threshold=MIN_TOTAL_BIOMASS_FOR_DIVERSITY):
    if rel_mat is None or rel_mat.shape[0] < 2:
        return np.nan

    rel_use = rel_mat.copy()
    if abs_mat is not None:
        common_days = [d for d in rel_use.index if d in abs_mat.index]
        rel_use = rel_use.loc[common_days]
        abs_use = abs_mat.reindex(index=common_days).fillna(0.0)
        total_biomass_by_day = abs_use.sum(axis=1)
        rel_use.loc[total_biomass_by_day < biomass_threshold, :] = 0.0

    vals = []
    for day_a, day_b in itertools.combinations(rel_use.index, 2):
        vals.append(calc_bray_curtis(
            rel_use.loc[day_a].to_numpy(dtype=float),
            rel_use.loc[day_b].to_numpy(dtype=float),
        ))
    if len(vals) == 0:
        return np.nan
    return float(np.mean(vals))


# ============================================================
# A. species correlation / variance decomposition
# ============================================================

def compute_decomposition(abs_mat):
    """
    输入:
        abs_mat: day × species 的绝对丰度矩阵

    返回:
        summary dict, pairwise list
    """
    X = abs_mat.to_numpy(dtype=float)
    days = list(abs_mat.index)
    taxa = list(abs_mat.columns)

    n_days, n_species = X.shape

    if n_days < 2 or n_species < 1:
        return None, []

    B = X.sum(axis=1)

    mean_total_biomass = float(np.mean(B))
    total_biomass_variance = float(np.var(B, ddof=1))
    total_biomass_sd = float(np.sqrt(total_biomass_variance))

    species_var = np.var(X, axis=0, ddof=1)
    species_sd = np.sqrt(species_var)

    sum_species_sd = float(np.sum(species_sd))
    species_variance_term = float(np.sum(species_var))

    if n_species == 1:
        cov_matrix = np.array([[species_var[0]]], dtype=float)
    else:
        cov_matrix = np.cov(X, rowvar=False, ddof=1)

    covariance_term = 0.0
    positive_covariance_sum = 0.0
    negative_covariance_sum = 0.0

    pair_rows = []

    correlations = []
    weighted_corr_numer = 0.0
    weighted_corr_denom = 0.0
    mean_pairwise_cov_values = []

    for i in range(n_species):
        for j in range(i + 1, n_species):
            cov_ij = float(cov_matrix[i, j])
            contribution_2cov = 2.0 * cov_ij

            covariance_term += contribution_2cov

            if contribution_2cov > 0:
                positive_covariance_sum += contribution_2cov
            elif contribution_2cov < 0:
                negative_covariance_sum += contribution_2cov

            sd_i = float(species_sd[i])
            sd_j = float(species_sd[j])

            if sd_i > 0 and sd_j > 0:
                corr_ij = cov_ij / (sd_i * sd_j)
                corr_ij = max(-1.0, min(1.0, corr_ij))

                correlations.append(corr_ij)

                weight = sd_i * sd_j
                weighted_corr_numer += corr_ij * weight
                weighted_corr_denom += weight
            else:
                corr_ij = np.nan
                weight = np.nan

            mean_pairwise_cov_values.append(cov_ij)

            pair_rows.append({
                "taxon_i": taxa[i],
                "taxon_j": taxa[j],
                "sd_i": sd_i,
                "sd_j": sd_j,
                "variance_i": float(species_var[i]),
                "variance_j": float(species_var[j]),
                "covariance": cov_ij,
                "correlation": corr_ij,
                "sd_product": weight,
                "contribution_2cov": contribution_2cov,
            })

    total_biomass_variance_from_terms = species_variance_term + covariance_term
    decomposition_error = (
        total_biomass_variance - total_biomass_variance_from_terms
    )

    community_cv = safe_divide(sum_species_sd, mean_total_biomass)
    total_biomass_cv = safe_divide(total_biomass_sd, mean_total_biomass)

    if sum_species_sd > 0:
        synchrony_phi = total_biomass_variance / (sum_species_sd ** 2)
    else:
        synchrony_phi = np.nan

    if not pd.isna(community_cv) and community_cv > 0:
        synchrony_phi_from_cv = (total_biomass_cv / community_cv) ** 2
    else:
        synchrony_phi_from_cv = np.nan

    phi_denominator = sum_species_sd ** 2
    phi_denominator_valid_robust = (
        n_species >= 2
        and sum_species_sd > PHI_MIN_SUM_SPECIES_SD
        and not pd.isna(community_cv)
        and community_cv >= PHI_MIN_COMMUNITY_CV_FOR_ROBUST
    )
    if phi_denominator_valid_robust:
        synchrony_phi_robust = synchrony_phi
        synchrony_phi_from_cv_robust = synchrony_phi_from_cv
    else:
        synchrony_phi_robust = np.nan
        synchrony_phi_from_cv_robust = np.nan

    if total_biomass_variance > 0:
        species_variance_fraction_of_total = (
            species_variance_term / total_biomass_variance
        )
        covariance_fraction_of_total = (
            covariance_term / total_biomass_variance
        )
    else:
        species_variance_fraction_of_total = np.nan
        covariance_fraction_of_total = np.nan

    if len(correlations) > 0:
        mean_pairwise_correlation = float(np.nanmean(correlations))
        median_pairwise_correlation = float(np.nanmedian(correlations))
    else:
        mean_pairwise_correlation = np.nan
        median_pairwise_correlation = np.nan

    if weighted_corr_denom > 0:
        weighted_mean_pairwise_correlation = (
            weighted_corr_numer / weighted_corr_denom
        )
    else:
        weighted_mean_pairwise_correlation = np.nan

    if len(mean_pairwise_cov_values) > 0:
        mean_pairwise_covariance = float(np.nanmean(mean_pairwise_cov_values))
    else:
        mean_pairwise_covariance = np.nan

    if species_variance_term + abs(covariance_term) > 0:
        relative_covariance_contribution = (
            covariance_term / (species_variance_term + abs(covariance_term))
        )
    else:
        relative_covariance_contribution = np.nan

    summary = {
        # basic
        "n_days_decomposition": n_days,
        "days_used_decomposition": ",".join(str(x) for x in days),
        "n_species_decomposition": n_species,

        # biomass
        "mean_total_biomass": mean_total_biomass,
        "total_biomass_sd": total_biomass_sd,
        "total_biomass_variance": total_biomass_variance,

        # variance/covariance decomposition
        "sum_species_sd": sum_species_sd,
        "species_variance_term": species_variance_term,
        "covariance_term": covariance_term,
        "positive_covariance_sum": positive_covariance_sum,
        "negative_covariance_sum": negative_covariance_sum,
        "total_biomass_variance_from_terms": total_biomass_variance_from_terms,
        "decomposition_error": decomposition_error,

        # CV and synchrony
        "community_cv": community_cv,
        "total_biomass_cv": total_biomass_cv,
        "synchrony_phi": synchrony_phi,
        "synchrony_phi_from_cv": synchrony_phi_from_cv,
        "phi_denominator": phi_denominator,
        "phi_min_community_cv_for_robust": PHI_MIN_COMMUNITY_CV_FOR_ROBUST,
        "phi_denominator_valid_robust": phi_denominator_valid_robust,
        "synchrony_phi_robust": synchrony_phi_robust,
        "synchrony_phi_from_cv_robust": synchrony_phi_from_cv_robust,

        # fractions
        "species_variance_fraction_of_total": species_variance_fraction_of_total,
        "covariance_fraction_of_total": covariance_fraction_of_total,
        "relative_covariance_contribution": relative_covariance_contribution,

        # species correlation
        "mean_pairwise_correlation": mean_pairwise_correlation,
        "median_pairwise_correlation": median_pairwise_correlation,
        "weighted_mean_pairwise_correlation": weighted_mean_pairwise_correlation,
        "mean_pairwise_covariance": mean_pairwise_covariance,
        "n_pairwise": len(pair_rows),
    }

    return summary, pair_rows


# ============================================================
# B. within-window turnover analysis
#    collapse-aware presence, but no day deletion
# ============================================================

def compute_turnover(rel_mat,
                     abs_mat,
                     threshold_rel=PRESENCE_THRESHOLD_REL,
                     biomass_threshold=MIN_TOTAL_BIOMASS_FOR_PRESENCE):
    """
    输入:
        rel_mat: day × species 相对丰度矩阵，单位 percent
        abs_mat: day × species 绝对丰度矩阵

    presence 定义:
        I_i(t) = 1 if rel_abund_i(t) >= threshold_rel
                      and total_biomass(t) >= biomass_threshold
                 else 0

    关键：
      - 不删除 low biomass / collapse day。
      - low biomass / collapse day 的 presence vector 设为全 0。
      - 因此 collapse / recovery 会反映为 disappearance / appearance。
    """
    if rel_mat is None or abs_mat is None:
        return {}

    # 对齐 day 和 species
    common_days = [d for d in rel_mat.index if d in abs_mat.index]
    if len(common_days) < 2:
        return {}

    # rel_mat 可能有所有 taxon，abs_mat 只保留非零 taxon。
    # 为 presence 计算，使用 rel_mat 的 columns。
    rel_use = rel_mat.loc[common_days].copy()
    abs_use = abs_mat.reindex(index=common_days).fillna(0.0)

    days = list(rel_use.index)
    n_days = len(days)

    if n_days < 2:
        return {}

    total_biomass_by_day = abs_use.sum(axis=1)
    biomass_valid_day = total_biomass_by_day >= biomass_threshold

    presence = (rel_use >= threshold_rel).astype(int)

    # collapse / low-biomass day: presence 全部设为 0
    presence.loc[~biomass_valid_day, :] = 0

    n_species_pool = presence.shape[1]

    daily_richness = presence.sum(axis=1).to_numpy(dtype=float)

    species_present_any = presence.sum(axis=0) > 0
    cumulative_richness = int(species_present_any.sum())

    mean_daily_richness = float(np.mean(daily_richness))
    min_daily_richness = float(np.min(daily_richness))
    max_daily_richness = float(np.max(daily_richness))

    if n_days >= 2:
        sd_daily_richness = float(np.std(daily_richness, ddof=1))
    else:
        sd_daily_richness = np.nan

    cumulative_excess_richness = cumulative_richness - mean_daily_richness

    if mean_daily_richness > 0:
        cumulative_richness_ratio = cumulative_richness / mean_daily_richness
    else:
        cumulative_richness_ratio = np.nan

    P = presence.to_numpy(dtype=int)

    appearance_events = 0
    disappearance_events = 0

    for t in range(1, n_days):
        prev = P[t - 1, :]
        curr = P[t, :]

        appearance_events += int(np.sum((prev == 0) & (curr == 1)))
        disappearance_events += int(np.sum((prev == 1) & (curr == 0)))

    turnover_events = appearance_events + disappearance_events

    # 保留归一化版本用于检查，但主图可直接用 turnover_events
    denom_cum = (n_days - 1) * cumulative_richness
    if denom_cum > 0:
        turnover_event_rate = turnover_events / denom_cum
        appearance_event_rate = appearance_events / denom_cum
        disappearance_event_rate = disappearance_events / denom_cum
    else:
        turnover_event_rate = np.nan
        appearance_event_rate = np.nan
        disappearance_event_rate = np.nan

    denom_daily = (n_days - 1) * mean_daily_richness
    if denom_daily > 0:
        turnover_events_per_daily_richness = turnover_events / denom_daily
        appearance_events_per_daily_richness = appearance_events / denom_daily
        disappearance_events_per_daily_richness = disappearance_events / denom_daily
    else:
        turnover_events_per_daily_richness = np.nan
        appearance_events_per_daily_richness = np.nan
        disappearance_events_per_daily_richness = np.nan

    # temporal Jaccard distance across all day pairs
    jaccard_values = []

    for a, b in itertools.combinations(range(n_days), 2):
        A = P[a, :].astype(bool)
        B = P[b, :].astype(bool)

        union = np.sum(A | B)
        intersection = np.sum(A & B)

        if union > 0:
            jaccard = 1.0 - intersection / union
            jaccard_values.append(jaccard)

    if len(jaccard_values) > 0:
        temporal_jaccard_mean = float(np.mean(jaccard_values))
        temporal_jaccard_median = float(np.median(jaccard_values))
        temporal_jaccard_max = float(np.max(jaccard_values))
    else:
        temporal_jaccard_mean = np.nan
        temporal_jaccard_median = np.nan
        temporal_jaccard_max = np.nan

    collapsed_days_for_presence = [
        int(day)
        for day, ok in zip(days, biomass_valid_day.to_numpy())
        if not ok
    ]

    summary = {
        "presence_threshold_rel": threshold_rel,
        "presence_biomass_threshold": biomass_threshold,

        "n_days_turnover": n_days,
        "days_used_turnover": ",".join(str(x) for x in days),
        "n_species_pool_turnover": n_species_pool,

        "n_collapsed_days_for_presence": len(collapsed_days_for_presence),
        "collapsed_days_for_presence": ",".join(str(x) for x in collapsed_days_for_presence),

        "min_total_biomass_for_turnover": float(np.min(total_biomass_by_day)),
        "mean_total_biomass_for_turnover": float(np.mean(total_biomass_by_day)),
        "max_total_biomass_for_turnover": float(np.max(total_biomass_by_day)),

        "mean_daily_richness": mean_daily_richness,
        "min_daily_richness": min_daily_richness,
        "max_daily_richness": max_daily_richness,
        "sd_daily_richness": sd_daily_richness,

        "cumulative_richness": cumulative_richness,
        "cumulative_excess_richness": cumulative_excess_richness,
        "cumulative_richness_ratio": cumulative_richness_ratio,

        "appearance_events": appearance_events,
        "disappearance_events": disappearance_events,
        "turnover_events": turnover_events,

        "appearance_event_rate": appearance_event_rate,
        "disappearance_event_rate": disappearance_event_rate,
        "turnover_event_rate": turnover_event_rate,

        "appearance_events_per_daily_richness": appearance_events_per_daily_richness,
        "disappearance_events_per_daily_richness": disappearance_events_per_daily_richness,
        "turnover_events_per_daily_richness": turnover_events_per_daily_richness,

        "temporal_jaccard_mean": temporal_jaccard_mean,
        "temporal_jaccard_median": temporal_jaccard_median,
        "temporal_jaccard_max": temporal_jaccard_max,
    }

    return summary


# ============================================================
# C. biomass-based Shannon cumulative diversity
# ============================================================

def shannon_from_vector(v):
    """
    v: abundance vector
    returns Shannon diversity based on normalized v.
    """
    v = np.asarray(v, dtype=float)
    v = v[np.isfinite(v)]
    total = np.sum(v)

    if total <= 0:
        return np.nan

    p = v / total
    p = p[p > 0]

    if len(p) == 0:
        return np.nan

    return float(-np.sum(p * np.log(p)))


def compute_biomass_shannon(abs_mat):
    """
    基于绝对 biomass 计算 Shannon 多样性。

    Daily Shannon:
        p_i(t) = X_i(t) / sum_j X_j(t)
        H(t) = -sum_i p_i(t) log(p_i(t))

    Cumulative Shannon:
        A_i = sum_t X_i(t)
        p_i_cum = A_i / sum_j A_j
        H_cum = -sum_i p_i_cum log(p_i_cum)

    cumulative_excess_shannon:
        H_cum - mean_t H(t)

    collapse day 保留。
    若某天 total biomass = 0，则该天 daily Shannon 为 NA，
    计算 mean_daily_shannon 时自动忽略 NA。
    """
    empty_result = {
        "mean_daily_shannon": np.nan,
        "min_daily_shannon": np.nan,
        "max_daily_shannon": np.nan,
        "sd_daily_shannon": np.nan,
        "cumulative_shannon": np.nan,
        "cumulative_excess_shannon": np.nan,
    }

    if abs_mat is None:
        return empty_result

    if abs_mat.shape[0] < 1 or abs_mat.shape[1] < 1:
        return empty_result

    X = abs_mat.to_numpy(dtype=float)
    total_biomass_by_day = np.nansum(X, axis=1)
    biomass_valid_day = total_biomass_by_day >= MIN_TOTAL_BIOMASS_FOR_DIVERSITY

    daily_shannon = np.array([
        shannon_from_vector(X[t, :]) if biomass_valid_day[t] else 0.0
        for t in range(X.shape[0])
    ])

    if np.sum(np.isfinite(daily_shannon)) == 0:
        mean_daily_shannon = np.nan
        min_daily_shannon = np.nan
        max_daily_shannon = np.nan
        sd_daily_shannon = np.nan
    else:
        mean_daily_shannon = float(np.nanmean(daily_shannon))
        min_daily_shannon = float(np.nanmin(daily_shannon))
        max_daily_shannon = float(np.nanmax(daily_shannon))

        if np.sum(np.isfinite(daily_shannon)) > 1:
            sd_daily_shannon = float(np.nanstd(daily_shannon, ddof=1))
        else:
            sd_daily_shannon = np.nan

    if np.any(biomass_valid_day):
        cumulative_biomass = np.nansum(X[biomass_valid_day, :], axis=0)
        cumulative_shannon = shannon_from_vector(cumulative_biomass)
    else:
        cumulative_shannon = 0.0

    if pd.isna(cumulative_shannon) or pd.isna(mean_daily_shannon):
        cumulative_excess_shannon = np.nan
    else:
        cumulative_excess_shannon = cumulative_shannon - mean_daily_shannon

    return {
        "mean_daily_shannon": mean_daily_shannon,
        "min_daily_shannon": min_daily_shannon,
        "max_daily_shannon": max_daily_shannon,
        "sd_daily_shannon": sd_daily_shannon,
        "cumulative_shannon": cumulative_shannon,
        "cumulative_excess_shannon": cumulative_excess_shannon,
    }


# ============================================================
# 主流程
# ============================================================

def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    print("读取 abs_abundance.csv ...")
    df_all = load_abs_abundance()

    print(f"  -> total rows after filtering: {len(df_all)}")
    print(f"  -> experiments: {sorted(df_all['experiment'].unique())}")
    print(f"  -> conditions: {sorted(df_all['condition'].unique())}")
    print(f"  -> presence threshold for turnover: {PRESENCE_THRESHOLD_REL:.2f}%")
    print(
        f"  -> biomass threshold for presence: "
        f"{MIN_TOTAL_BIOMASS_FOR_PRESENCE:.4f}"
    )
    print(
        f"  -> biomass threshold for diversity-like Shannon: "
        f"{MIN_TOTAL_BIOMASS_FOR_DIVERSITY:.4f}"
    )

    for window_name, exp_cfg in WINDOW_CONFIG.items():
        print(f"\n========== window: {window_name} ==========")

        out_window_dir = os.path.join(OUT_DIR, window_name)
        os.makedirs(out_window_dir, exist_ok=True)

        community_rows = []
        pairwise_rows = []

        for experiment in ["mortality", "temperature"]:
            if experiment not in exp_cfg:
                continue

            days = exp_cfg[experiment]["days"]
            replicas = exp_cfg[experiment]["replicas"]

            print(f"\n  {experiment}: days={days}, replicas={replicas}")

            df_exp = df_all[
                (df_all["experiment"] == experiment)
                & (df_all["day"].isin(days))
                & (df_all["replica"].isin(replicas))
            ].copy()

            if df_exp.empty:
                print("    [warn] no data")
                continue

            group_cols = ["experiment", "condition", "community", "replica"]

            n_group_done = 0
            n_group_skipped = 0

            for key, df_group in df_exp.groupby(group_cols, sort=True):
                exp_i, condition, community, replica = key

                base_info = {
                    "window": window_name,
                    "experiment": exp_i,
                    "condition": condition,
                    "community": community,
                    "replica": replica,
                    "configured_days": ",".join(str(x) for x in days),
                    "configured_replicas": ",".join(str(x) for x in replicas),
                }

                # --------------------------------------------------
                # abs_abund matrix
                # 用于 decomposition 和 biomass-based Shannon。
                # collapse day 保留。
                # --------------------------------------------------
                abs_mat = make_abundance_matrix(
                    df_group,
                    expected_days=days,
                    value_col="abs_abund",
                )

                if abs_mat is not None:
                    decomp_summary, pairs = compute_decomposition(abs_mat)
                    shannon_summary = compute_biomass_shannon(abs_mat)
                else:
                    decomp_summary, pairs = None, []
                    shannon_summary = compute_biomass_shannon(None)

                # --------------------------------------------------
                # rel_abund matrix
                # 用于 1% presence/absence turnover。
                # low-biomass day 不删除，而是 presence 全 0。
                # --------------------------------------------------
                rel_mat = make_abundance_matrix(
                    df_group,
                    expected_days=days,
                    value_col="rel_abund",
                )
                temporal_bc_mean = compute_temporal_bc_mean(
                    rel_mat, abs_mat, MIN_TOTAL_BIOMASS_FOR_DIVERSITY
                )

                turnover_summary = compute_turnover(
                    rel_mat=rel_mat,
                    abs_mat=abs_mat,
                    threshold_rel=PRESENCE_THRESHOLD_REL,
                    biomass_threshold=MIN_TOTAL_BIOMASS_FOR_PRESENCE,
                )

                # 至少有一套方法成功，就输出 community-level row
                if decomp_summary is None and not turnover_summary:
                    n_group_skipped += 1
                    print(
                        f"    [skip] {exp_i} {condition} "
                        f"community={community} R{replica}: insufficient data"
                    )
                    continue

                row = dict(base_info)

                if decomp_summary is not None:
                    row.update(decomp_summary)
                else:
                    row.update({
                        "n_days_decomposition": np.nan,
                        "days_used_decomposition": "",
                        "n_species_decomposition": np.nan,
                        "mean_total_biomass": np.nan,
                        "total_biomass_sd": np.nan,
                        "total_biomass_variance": np.nan,
                        "sum_species_sd": np.nan,
                        "species_variance_term": np.nan,
                        "covariance_term": np.nan,
                        "positive_covariance_sum": np.nan,
                        "negative_covariance_sum": np.nan,
                        "total_biomass_variance_from_terms": np.nan,
                        "decomposition_error": np.nan,
                        "community_cv": np.nan,
                        "total_biomass_cv": np.nan,
                        "synchrony_phi": np.nan,
                        "synchrony_phi_from_cv": np.nan,
                        "phi_denominator": np.nan,
                        "phi_min_community_cv_for_robust":
                            PHI_MIN_COMMUNITY_CV_FOR_ROBUST,
                        "phi_denominator_valid_robust": False,
                        "synchrony_phi_robust": np.nan,
                        "synchrony_phi_from_cv_robust": np.nan,
                        "species_variance_fraction_of_total": np.nan,
                        "covariance_fraction_of_total": np.nan,
                        "relative_covariance_contribution": np.nan,
                        "mean_pairwise_correlation": np.nan,
                        "median_pairwise_correlation": np.nan,
                        "weighted_mean_pairwise_correlation": np.nan,
                        "mean_pairwise_covariance": np.nan,
                        "n_pairwise": np.nan,
                    })

                if turnover_summary:
                    row.update(turnover_summary)
                else:
                    row.update({
                        "presence_threshold_rel": PRESENCE_THRESHOLD_REL,
                        "presence_biomass_threshold": MIN_TOTAL_BIOMASS_FOR_PRESENCE,
                        "n_days_turnover": np.nan,
                        "days_used_turnover": "",
                        "n_species_pool_turnover": np.nan,
                        "n_collapsed_days_for_presence": np.nan,
                        "collapsed_days_for_presence": "",
                        "min_total_biomass_for_turnover": np.nan,
                        "mean_total_biomass_for_turnover": np.nan,
                        "max_total_biomass_for_turnover": np.nan,
                        "mean_daily_richness": np.nan,
                        "min_daily_richness": np.nan,
                        "max_daily_richness": np.nan,
                        "sd_daily_richness": np.nan,
                        "cumulative_richness": np.nan,
                        "cumulative_excess_richness": np.nan,
                        "cumulative_richness_ratio": np.nan,
                        "appearance_events": np.nan,
                        "disappearance_events": np.nan,
                        "turnover_events": np.nan,
                        "appearance_event_rate": np.nan,
                        "disappearance_event_rate": np.nan,
                        "turnover_event_rate": np.nan,
                        "appearance_events_per_daily_richness": np.nan,
                        "disappearance_events_per_daily_richness": np.nan,
                        "turnover_events_per_daily_richness": np.nan,
                        "temporal_jaccard_mean": np.nan,
                        "temporal_jaccard_median": np.nan,
                        "temporal_jaccard_max": np.nan,
                    })

                row.update(shannon_summary)
                row["temporal_bc_mean"] = temporal_bc_mean

                community_rows.append(row)

                for p in pairs:
                    pairwise_rows.append({
                        **base_info,
                        **p,
                    })

                n_group_done += 1

            print(
                f"    -> groups done: {n_group_done}, "
                f"skipped: {n_group_skipped}"
            )

        community_df = pd.DataFrame(community_rows)
        pairwise_df = pd.DataFrame(pairwise_rows)

        community_path = os.path.join(
            out_window_dir, "community_decomposition.csv"
        )
        pairwise_path = os.path.join(
            out_window_dir, "pairwise_correlations.csv"
        )
        fluctuating_path = os.path.join(
            out_window_dir, "community_decomposition_fluctuating_cv.csv"
        )
        fluctuating_phi_path = os.path.join(
            out_window_dir, "community_decomposition_phi_fluctuating_cv.csv"
        )
        phi_robust_path = os.path.join(
            out_window_dir, "community_decomposition_phi_robust.csv"
        )
        fluctuating_phi_robust_path = os.path.join(
            out_window_dir,
            "community_decomposition_phi_robust_fluctuating_cv.csv",
        )

        if not community_df.empty:
            community_df = community_df.sort_values([
                "experiment", "condition", "community", "replica"
            ])
            community_df["community_cv_threshold"] = (
                COMMUNITY_CV_FLUCTUATING_THRESHOLD
            )
            community_df["fluctuating_by_community_cv"] = (
                pd.to_numeric(community_df["community_cv"], errors="coerce")
                >= community_df["community_cv_threshold"]
            )
            community_df["fluctuating_composite"] = (
                community_df["fluctuating_by_community_cv"]
            )
            community_df.to_csv(community_path, index=False)

            print(f"\n  -> saved: {community_path}")
            print(f"     community rows: {len(community_df)}")

            fluctuating_df = community_df[
                community_df["fluctuating_composite"]
            ].copy()
            fluctuating_df.to_csv(fluctuating_path, index=False)
            print(f"  -> saved: {fluctuating_path}")

            fluctuating_phi_df = fluctuating_df[
                pd.to_numeric(
                    fluctuating_df["synchrony_phi"], errors="coerce"
                ).notna()
            ].copy()
            fluctuating_phi_df.to_csv(fluctuating_phi_path, index=False)
            print(f"  -> saved: {fluctuating_phi_path}")

            phi_robust_df = community_df[
                pd.to_numeric(
                    community_df["synchrony_phi_robust"], errors="coerce"
                ).notna()
            ].copy()
            phi_robust_df.to_csv(phi_robust_path, index=False)
            print(f"  -> saved: {phi_robust_path}")

            fluctuating_phi_robust_df = fluctuating_df[
                pd.to_numeric(
                    fluctuating_df["synchrony_phi_robust"], errors="coerce"
                ).notna()
            ].copy()
            fluctuating_phi_robust_df.to_csv(
                fluctuating_phi_robust_path, index=False
            )
            print(f"  -> saved: {fluctuating_phi_robust_path}")
            print(
                "     community rows with "
                f"community_cv >= {COMMUNITY_CV_FLUCTUATING_THRESHOLD}: "
                f"{len(fluctuating_df)}"
            )
            print(
                "     robust phi rows with "
                f"community_cv >= {PHI_MIN_COMMUNITY_CV_FOR_ROBUST}: "
                f"{len(phi_robust_df)}"
            )
        else:
            print(f"\n  [warn] no community rows for {window_name}")

        if not pairwise_df.empty:
            pairwise_df = pairwise_df.sort_values([
                "experiment", "condition", "community",
                "replica", "taxon_i", "taxon_j"
            ])
            pairwise_df.to_csv(pairwise_path, index=False)

            print(f"  -> saved: {pairwise_path}")
            print(f"     pairwise rows: {len(pairwise_df)}")
        else:
            print(f"  [warn] no pairwise rows for {window_name}")

    print("\n完成！")
    print(f"输出目录: {OUT_DIR}")
    print("\ncommunity_decomposition.csv 中包含：")
    print("  A. abs_abund-based variance/covariance/synchrony decomposition")
    print("  B. collapse-aware rel_abund >= 1% within-window turnover analysis")
    print("  C. abs_abund biomass-based cumulative Shannon diversity")
    print(
        "  Extra: community_decomposition_fluctuating_cv.csv keeps only "
        f"community_cv >= {COMMUNITY_CV_FLUCTUATING_THRESHOLD}"
    )
    print("  Extra: community_decomposition_phi_fluctuating_cv.csv is kept for compatibility")
    print(
        "  Extra: community_decomposition_phi_robust.csv excludes near-zero "
        f"phi denominators (community_cv < {PHI_MIN_COMMUNITY_CV_FOR_ROBUST})"
    )
    print("\n特殊规则：")
    print("  temperature W5 only community 3, 7, 12")
    print("  low-biomass days are NOT deleted")
    print("  low-biomass days only get all-zero presence vectors in turnover analysis")


if __name__ == "__main__":
    main()
