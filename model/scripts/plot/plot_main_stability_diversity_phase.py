"""Phase maps for stability-diversity relationships.

Default input uses the manuscript main CR-pH model dataset:

    data/runs/rh_gnormal_sd3_R1_sparse6_tri_S12_h015_g100_e5000

The script recomputes oscillation labels with CV > 0.1 and diversity with an
absolute B_i > 1e-3 activity threshold, then plots:

1. fluctuation fraction
2. diversity difference: mean(diversity | oscillating) - mean(diversity | stable)
3. distribution-free effect size: rank-biserial / Cliff's delta
4. Spearman rho between CV and diversity among non-collapsed communities
5. standard-deviation participation ratio among oscillating communities
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
from matplotlib.colors import LinearSegmentedColormap, ListedColormap, Normalize, TwoSlopeNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy import ndimage, stats


matplotlib.rcParams.update(
    {
        "font.family": "Arial",
        "font.size": 16,
        "axes.titlesize": 16,
        "axes.labelsize": 16,
        "xtick.labelsize": 14,
        "ytick.labelsize": 14,
        "axes.linewidth": 1.2,
        "xtick.major.width": 1.2,
        "ytick.major.width": 1.2,
        "xtick.major.size": 5,
        "ytick.major.size": 5,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
        "figure.dpi": 150,
        "savefig.dpi": 600,
        "savefig.bbox": "tight",
        "savefig.facecolor": "white",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_RUN_DIR = ROOT / "data" / "runs" / "rh_gnormal_sd3_R1_sparse6_tri_S12_h015_g100_e5000"
DEFAULT_OUT_DIR = ROOT / "figures" / "main_stability_diversity" / "rh_gnormal_sd3_R1_sparse6_tri_S12_h015_g100_e5000_cv010_oscmin010"


def cv_suffix(cv_thresh: float) -> str:
    return f"cv{int(round(cv_thresh * 100)):03d}"


def white_to_colour(hex_end: str, name: str) -> LinearSegmentedColormap:
    return LinearSegmentedColormap.from_list(
        name, [(1.0, 1.0, 1.0), matplotlib.colors.to_rgb(hex_end)], N=256
    )


CMAP_OSC = white_to_colour("#E84A1D", "w_orangered")
CMAP_PR = white_to_colour("#6A3D9A", "w_purple")
CMAP_GREENS = white_to_colour("#117733", "w_green")
CMAP_GREY = ListedColormap(["0.84"])


def make_edges(values: np.ndarray) -> np.ndarray:
    values = values.astype(float)
    delta = np.diff(values)
    return np.concatenate(
        [[values[0] - delta[0] / 2], (values[:-1] + values[1:]) / 2, [values[-1] + delta[-1] / 2]]
    )


def shannon(values: np.ndarray, axis: int = -1) -> np.ndarray:
    total = values.sum(axis=axis, keepdims=True)
    p = np.divide(values, total, out=np.zeros_like(values, dtype=float), where=total > 0)
    with np.errstate(divide="ignore", invalid="ignore"):
        log_p = np.where(p > 0, np.log(p), 0.0)
    return -(p * log_p).sum(axis=axis)


def std_participation_ratio(x: np.ndarray) -> np.ndarray:
    sigma = np.std(x, axis=-1, ddof=1)
    sum_sigma = sigma.sum(axis=-1)
    sum_sigma_sq = (sigma**2).sum(axis=-1)
    return np.divide(sum_sigma**2, sum_sigma_sq, out=np.full_like(sum_sigma, np.nan), where=sum_sigma_sq > 0)


def one_minus_max_std_share(x: np.ndarray) -> np.ndarray:
    sigma = np.std(x, axis=-1, ddof=1)
    sum_sigma = sigma.sum(axis=-1)
    max_share = np.divide(
        sigma.max(axis=-1),
        sum_sigma,
        out=np.full_like(sum_sigma, np.nan),
        where=sum_sigma > 0,
    )
    return 1.0 - max_share


def symmetric_norm(values: np.ndarray, min_abs: float = 0.1, percentile: float = 95.0) -> TwoSlopeNorm:
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        vmax = 1.0
    else:
        vmax = float(np.nanpercentile(np.abs(finite), percentile))
        vmax = max(vmax, min_abs)
    return TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax)


def positive_norm(values: np.ndarray, min_max: float = 1.0, percentile: float = 95.0) -> Normalize:
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        vmax = min_max
    else:
        vmax = float(np.nanpercentile(finite, percentile))
        vmax = max(vmax, min_max)
    return Normalize(0.0, vmax)


def fraction_norm(values: np.ndarray, fallback_max: float = 1.0, percentile: float = 98.0) -> Normalize:
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        vmax = fallback_max
    else:
        vmax = float(np.nanpercentile(finite, percentile))
        vmax = min(max(vmax, 0.05), fallback_max)
    return Normalize(0.0, vmax)


def smooth_nan(values: np.ndarray, sigma: float) -> np.ndarray:
    if sigma <= 0:
        return values
    finite = np.isfinite(values)
    if not finite.any():
        return values
    filled = np.where(finite, values, 0.0)
    weights = finite.astype(float)
    smoothed = ndimage.gaussian_filter(filled, sigma=sigma, mode="nearest")
    norm = ndimage.gaussian_filter(weights, sigma=sigma, mode="nearest")
    out = np.divide(smoothed, norm, out=np.full_like(values, np.nan, dtype=float), where=norm > 1e-12)
    out[~finite & (norm <= 1e-12)] = np.nan
    return out


def plot_phase(
    z: np.ndarray,
    r_axis: np.ndarray,
    h_axis: np.ndarray,
    out_path: Path,
    cmap,
    norm,
    grey_mask: np.ndarray | None = None,
) -> None:
    fig, ax = plt.subplots(figsize=(3.5, 3.2), constrained_layout=True)
    xe = make_edges(h_axis)
    ye = make_edges(r_axis)

    if grey_mask is not None:
        grey = np.where(grey_mask, 0.0, np.nan)
        ax.pcolormesh(xe, ye, grey, cmap=CMAP_GREY, vmin=-1, vmax=1, shading="flat", rasterized=True)

    pcm = ax.pcolormesh(
        xe,
        ye,
        np.ma.masked_invalid(z),
        cmap=cmap,
        norm=norm,
        shading="flat",
        rasterized=True,
    )
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="4.5%", pad=0.08)
    cb = fig.colorbar(pcm, cax=cax)
    cb.outline.set_linewidth(1.0)
    cb.ax.tick_params(labelsize=12, length=3, direction="in", width=1.0)

    ax.set_xlim(xe[0], xe[-1])
    ax.set_ylim(ye[0], ye[-1])
    ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(4))
    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator(4))
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_title("")
    for spine in ax.spines.values():
        spine.set_linewidth(1.2)
        spine.set_color("0.15")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(out_path) + ".pdf")
    fig.savefig(str(out_path) + ".png")
    plt.close(fig)


def compare_groups(
    metric: np.ndarray,
    cv: np.ndarray,
    is_osc: np.ndarray,
    is_stable: np.ndarray,
    valid: np.ndarray,
    min_group: int,
) -> dict[str, np.ndarray]:
    n_r, n_h, _ = metric.shape
    delta_mean = np.full((n_r, n_h), np.nan, dtype=float)
    cohens_d = np.full((n_r, n_h), np.nan, dtype=float)
    rank_biserial = np.full((n_r, n_h), np.nan, dtype=float)
    spearman_rho = np.full((n_r, n_h), np.nan, dtype=float)
    n_osc = is_osc.sum(axis=-1).astype(float)
    n_stable = is_stable.sum(axis=-1).astype(float)

    for ir in range(n_r):
        for ih in range(n_h):
            osc_mask = is_osc[ir, ih] & np.isfinite(metric[ir, ih])
            stable_mask = is_stable[ir, ih] & np.isfinite(metric[ir, ih])
            no = int(osc_mask.sum())
            ns = int(stable_mask.sum())
            if no >= min_group and ns >= min_group:
                xo = metric[ir, ih, osc_mask]
                xs = metric[ir, ih, stable_mask]
                delta_mean[ir, ih] = xo.mean() - xs.mean()

                pooled_var = ((no - 1) * xo.var(ddof=1) + (ns - 1) * xs.var(ddof=1)) / (no + ns - 2)
                if pooled_var > 0:
                    cohens_d[ir, ih] = delta_mean[ir, ih] / np.sqrt(pooled_var)

                u_stat = stats.mannwhitneyu(xo, xs, alternative="two-sided").statistic
                rank_biserial[ir, ih] = 2.0 * u_stat / (no * ns) - 1.0

            rho_mask = valid[ir, ih] & np.isfinite(metric[ir, ih]) & np.isfinite(cv[ir, ih])
            if int(rho_mask.sum()) >= max(20, min_group):
                x = cv[ir, ih, rho_mask]
                y = metric[ir, ih, rho_mask]
                if x.std() > 0 and y.std() > 0:
                    rho, _ = stats.spearmanr(x, y)
                    spearman_rho[ir, ih] = rho

    return {
        "delta_mean": delta_mean,
        "cohens_d": cohens_d,
        "rank_biserial": rank_biserial,
        "spearman_rho": spearman_rho,
        "n_osc": n_osc,
        "n_stable": n_stable,
    }


def find_scan_name(run_dir: Path) -> str:
    if (run_dir / "phase_r_h_meta.h5").exists():
        return "phase_r_h"
    if (run_dir / "phase_h_f_meta.h5").exists():
        return "phase_h_f"
    meta_files = sorted(run_dir.glob("*_meta.h5"))
    if meta_files:
        return meta_files[0].name.removesuffix("_meta.h5")
    raise FileNotFoundError(f"No *_meta.h5 file found in {run_dir}")


def load_meta(run_dir: Path, scan_name: str) -> tuple[np.ndarray, np.ndarray, int, int, int]:
    meta_path = run_dir / f"{scan_name}_meta.h5"
    with h5py.File(meta_path, "r") as h5:
        axis1 = h5["axis1"][...]
        axis2 = h5["axis2"][...]
        n_rows = len(axis1)
        n_axis2 = len(axis2)
        if "config_json" in h5:
            config_text = h5["config_json"][()]
        elif "config_json" in h5.attrs:
            config_text = h5.attrs["config_json"]
        else:
            config_text = None
        if isinstance(config_text, bytes):
            config_text = config_text.decode("utf-8")
        if config_text is not None:
            cfg = json.loads(str(config_text))
            species_count = int(cfg["model"]["S"])
        elif "S" in h5.attrs:
            species_count = int(h5.attrs["S"])
        else:
            species_count = -1
    return axis1, axis2, n_rows, n_axis2, species_count


def row_path(run_dir: Path, scan_name: str, row_index_1based: int) -> Path:
    return run_dir / "raw" / f"{scan_name}_rows" / f"{scan_name}_row_{row_index_1based:04d}.h5"


def legacy_phase2_row_path(run_dir: Path, row_index_1based: int) -> Path:
    return run_dir / "phase2_rows" / f"phase2_row_{row_index_1based:03d}.h5"


def compute_maps(
    run_dir: Path,
    b_thresh: float,
    cv_thresh: float,
    collapse_thresh: float,
    min_group: int,
) -> dict[str, np.ndarray]:
    scan_name = find_scan_name(run_dir)
    r_axis, h_axis, n_rows, n_axis2, species_count = load_meta(run_dir, scan_name)

    first = None
    for ir in range(n_rows):
        path = row_path(run_dir, scan_name, ir + 1)
        if path.exists():
            with h5py.File(path, "r") as h5:
                first = h5["B_last10"].shape
            break
    if first is None:
        raise FileNotFoundError(f"No raw rows found under {run_dir}")
    _, n_ens, inferred_species_count, _ = first
    if species_count <= 0:
        species_count = inferred_species_count

    shape = (n_rows, n_axis2, n_ens)
    community_cv = np.full(shape, np.nan, dtype=np.float32)
    is_collapsed = np.zeros(shape, dtype=bool)
    is_osc = np.zeros(shape, dtype=bool)
    pr_sigma = np.full(shape, np.nan, dtype=np.float32)
    one_minus_max_pi_B = np.full(shape, np.nan, dtype=np.float32)
    metrics = {
        "richness_last_absB": np.full(shape, np.nan, dtype=np.float32),
        "richness_mean_last3_absB": np.full(shape, np.nan, dtype=np.float32),
        "richness_cum3_absB": np.full(shape, np.nan, dtype=np.float32),
        "shannon_last_absB": np.full(shape, np.nan, dtype=np.float32),
        "shannon_mean_last3_absB": np.full(shape, np.nan, dtype=np.float32),
        "shannon_cum3_absB": np.full(shape, np.nan, dtype=np.float32),
    }
    completed_axis2 = np.zeros((n_rows, n_axis2), dtype=bool)

    for ir in range(n_rows):
        path = row_path(run_dir, scan_name, ir + 1)
        if not path.exists():
            print(f"[missing] row {ir + 1:04d}")
            continue

        with h5py.File(path, "r") as h5:
            b10 = h5["B_last10"][...].astype(np.float64)
            complete = h5["completed_axis2"][...].astype(bool)
        completed_axis2[ir] = complete
        if not complete.any():
            continue

        b3 = b10[..., -3:]
        total3 = b3.sum(axis=2)
        total_mean = total3.mean(axis=-1)
        collapsed = total_mean < collapse_thresh
        sigma3 = b3.std(axis=-1, ddof=1)
        sum_sigma3 = sigma3.sum(axis=-1)
        cv3 = np.divide(sum_sigma3, total_mean, out=np.zeros_like(total_mean), where=total_mean > 0)
        osc = (~collapsed) & (cv3 > cv_thresh)

        b3_active = np.where(b3 > b_thresh, b3, 0.0)
        b10_active = np.where(b10 > b_thresh, b10, 0.0)
        rich_daily = (b3_active > 0).sum(axis=2)
        sh_daily = shannon(np.moveaxis(b3_active, -1, -2), axis=-1)

        community_cv[ir] = cv3.astype(np.float32)
        is_collapsed[ir] = collapsed
        is_osc[ir] = osc
        metrics["richness_last_absB"][ir] = rich_daily[..., -1].astype(np.float32)
        metrics["richness_mean_last3_absB"][ir] = rich_daily.mean(axis=-1).astype(np.float32)
        metrics["richness_cum3_absB"][ir] = (b3_active.sum(axis=-1) > 0).sum(axis=-1).astype(np.float32)
        metrics["shannon_last_absB"][ir] = shannon(b3_active[..., -1], axis=-1).astype(np.float32)
        metrics["shannon_mean_last3_absB"][ir] = sh_daily.mean(axis=-1).astype(np.float32)
        metrics["shannon_cum3_absB"][ir] = shannon(b3_active.sum(axis=-1), axis=-1).astype(np.float32)
        pr_sigma[ir] = std_participation_ratio(b10_active).astype(np.float32)
        one_minus_max_pi_B[ir] = one_minus_max_std_share(b10_active).astype(np.float32)

        if (ir + 1) % 10 == 0 or ir == 0:
            print(f"[row {ir + 1:04d}] osc={osc.mean():.3f}, collapsed={collapsed.mean():.3f}")

    noncollapsed = (~is_collapsed) & completed_axis2[..., None]
    stable = noncollapsed & (~is_osc)
    osc_frac = np.where(completed_axis2, is_osc.mean(axis=-1), np.nan)
    collapse_frac = np.where(completed_axis2, is_collapsed.mean(axis=-1), np.nan)
    grey_group_mask = (is_osc.sum(axis=-1) < min_group) | (stable.sum(axis=-1) < min_group)

    out = {
        "axis1": r_axis,
        "axis2": h_axis,
        "community_CV_B": community_cv,
        "is_collapsed_B": is_collapsed,
        "is_oscillating_B": is_osc,
        "is_stable_B": stable,
        "fluctuation_fraction": osc_frac.astype(np.float32),
        "collapse_fraction_B": collapse_frac.astype(np.float32),
        "std_participation_ratio_B": pr_sigma,
        "std_participation_ratio_B_osc_mean": np.nanmean(np.where(is_osc, pr_sigma, np.nan), axis=-1).astype(np.float32),
        "std_participation_ratio_B_osc_norm_mean": (
            np.nanmean(np.where(is_osc, pr_sigma, np.nan), axis=-1) / species_count
        ).astype(np.float32),
        "one_minus_max_std_share_B": one_minus_max_pi_B,
        "one_minus_max_std_share_B_osc_mean": (
            np.nanmean(np.where(is_osc, one_minus_max_pi_B, np.nan), axis=-1)
        ).astype(np.float32),
        "grey_group_mask": grey_group_mask,
    }
    out.update(metrics)

    for metric_name, metric_values in metrics.items():
        comp = compare_groups(metric_values, community_cv, is_osc, stable, noncollapsed, min_group)
        for key, value in comp.items():
            out[f"{key}_{metric_name}"] = value.astype(np.float32)

    return out


def compute_maps_legacy_phase2(
    run_dir: Path,
    b_thresh: float,
    cv_thresh: float,
    collapse_thresh: float,
    min_group: int,
) -> dict[str, np.ndarray]:
    meta_path = run_dir / "CR_phase2_meta.h5"
    with h5py.File(meta_path, "r") as h5:
        r_axis = h5["rmean_list"][...]
        h_axis = h5["hstress_list"][...]
        n_rows = int(h5.attrs["n_ax1"])
        n_axis2 = int(h5.attrs["n_ax2"])
        n_ens = int(h5.attrs["n_ensemble"])
        species_count = int(h5.attrs["S"])

    shape = (n_rows, n_axis2, n_ens)
    community_cv = np.full(shape, np.nan, dtype=np.float32)
    is_collapsed = np.zeros(shape, dtype=bool)
    is_osc = np.zeros(shape, dtype=bool)
    pr_sigma = np.full(shape, np.nan, dtype=np.float32)
    one_minus_max_pi_B = np.full(shape, np.nan, dtype=np.float32)
    metrics = {
        "richness_last_absB": np.full(shape, np.nan, dtype=np.float32),
        "richness_mean_last3_absB": np.full(shape, np.nan, dtype=np.float32),
        "richness_cum3_absB": np.full(shape, np.nan, dtype=np.float32),
        "shannon_last_absB": np.full(shape, np.nan, dtype=np.float32),
        "shannon_mean_last3_absB": np.full(shape, np.nan, dtype=np.float32),
        "shannon_cum3_absB": np.full(shape, np.nan, dtype=np.float32),
    }
    completed_axis2 = np.zeros((n_rows, n_axis2), dtype=bool)

    for ir in range(n_rows):
        path = legacy_phase2_row_path(run_dir, ir + 1)
        if not path.exists():
            print(f"[missing] legacy row {ir + 1:03d}")
            continue

        with h5py.File(path, "r") as h5:
            b3 = h5["B_last3"][...].astype(np.float64)
        completed_axis2[ir] = True

        total3 = b3.sum(axis=2)
        total_mean = total3.mean(axis=-1)
        collapsed = total_mean < collapse_thresh
        sigma3 = b3.std(axis=-1, ddof=1)
        sum_sigma3 = sigma3.sum(axis=-1)
        cv3 = np.divide(sum_sigma3, total_mean, out=np.zeros_like(total_mean), where=total_mean > 0)
        osc = (~collapsed) & (cv3 > cv_thresh)

        b3_active = np.where(b3 > b_thresh, b3, 0.0)
        rich_daily = (b3_active > 0).sum(axis=2)
        sh_daily = shannon(np.moveaxis(b3_active, -1, -2), axis=-1)
        b3_cum = b3.sum(axis=-1)
        b3_cum_active = np.where(b3_cum > b_thresh, b3_cum, 0.0)

        community_cv[ir] = cv3.astype(np.float32)
        is_collapsed[ir] = collapsed
        is_osc[ir] = osc
        metrics["richness_last_absB"][ir] = rich_daily[..., -1].astype(np.float32)
        metrics["richness_mean_last3_absB"][ir] = rich_daily.mean(axis=-1).astype(np.float32)
        metrics["richness_cum3_absB"][ir] = (b3_cum_active > 0).sum(axis=-1).astype(np.float32)
        metrics["shannon_last_absB"][ir] = shannon(b3_active[..., -1], axis=-1).astype(np.float32)
        metrics["shannon_mean_last3_absB"][ir] = sh_daily.mean(axis=-1).astype(np.float32)
        metrics["shannon_cum3_absB"][ir] = shannon(b3_cum_active, axis=-1).astype(np.float32)
        pr_sigma[ir] = std_participation_ratio(b3_active).astype(np.float32)
        one_minus_max_pi_B[ir] = one_minus_max_std_share(b3_active).astype(np.float32)

        if (ir + 1) % 10 == 0 or ir == 0:
            print(f"[legacy row {ir + 1:03d}] osc={osc.mean():.3f}, collapsed={collapsed.mean():.3f}")

    noncollapsed = (~is_collapsed) & completed_axis2[..., None]
    stable = noncollapsed & (~is_osc)
    osc_frac = np.where(completed_axis2, is_osc.mean(axis=-1), np.nan)
    collapse_frac = np.where(completed_axis2, is_collapsed.mean(axis=-1), np.nan)
    grey_group_mask = (is_osc.sum(axis=-1) < min_group) | (stable.sum(axis=-1) < min_group)
    pr_osc_mean = np.nanmean(np.where(is_osc, pr_sigma, np.nan), axis=-1)

    out = {
        "axis1": r_axis,
        "axis2": h_axis,
        "community_CV_B": community_cv,
        "is_collapsed_B": is_collapsed,
        "is_oscillating_B": is_osc,
        "is_stable_B": stable,
        "fluctuation_fraction": osc_frac.astype(np.float32),
        "collapse_fraction_B": collapse_frac.astype(np.float32),
        "std_participation_ratio_B": pr_sigma,
        "std_participation_ratio_B_osc_mean": pr_osc_mean.astype(np.float32),
        "std_participation_ratio_B_osc_norm_mean": (pr_osc_mean / species_count).astype(np.float32),
        "one_minus_max_std_share_B": one_minus_max_pi_B,
        "one_minus_max_std_share_B_osc_mean": (
            np.nanmean(np.where(is_osc, one_minus_max_pi_B, np.nan), axis=-1)
        ).astype(np.float32),
        "grey_group_mask": grey_group_mask,
    }
    out.update(metrics)

    for metric_name, metric_values in metrics.items():
        comp = compare_groups(metric_values, community_cv, is_osc, stable, noncollapsed, min_group)
        for key, value in comp.items():
            out[f"{key}_{metric_name}"] = value.astype(np.float32)

    return out


def save_plots(results: dict[str, np.ndarray], out_dir: Path, suffix: str) -> None:
    r_axis = results["axis1"]
    h_axis = results["axis2"]
    grey_mask = results["grey_group_mask"]
    osc_frac_min = float(results.get("osc_frac_min", 0.10))
    smooth_sigma = float(results.get("smooth_sigma", 0.0))
    osc_frac_gate = results["fluctuation_fraction"] > osc_frac_min
    osc_gate_grey_mask = grey_mask | (~osc_frac_gate)
    
    def zmap(name: str) -> np.ndarray:
        return smooth_nan(results[name], smooth_sigma)

    plot_phase(
        zmap("fluctuation_fraction"),
        r_axis,
        h_axis,
        out_dir / f"01_fluctuation_fraction_{suffix}",
        CMAP_OSC,
        fraction_norm(zmap("fluctuation_fraction")),
    )
    plot_phase(
        zmap("collapse_fraction_B"),
        r_axis,
        h_axis,
        out_dir / "01b_collapse_fraction_B",
        "Greys",
        fraction_norm(zmap("collapse_fraction_B")),
    )
    plot_phase(
        zmap("std_participation_ratio_B_osc_mean"),
        r_axis,
        h_axis,
        out_dir / "03_std_participation_ratio_B_osc_mean",
        CMAP_PR,
        positive_norm(zmap("std_participation_ratio_B_osc_mean"), min_max=1.0),
        grey_mask=osc_gate_grey_mask,
    )
    plot_phase(
        zmap("std_participation_ratio_B_osc_norm_mean"),
        r_axis,
        h_axis,
        out_dir / "03b_std_participation_ratio_B_osc_norm_mean",
        CMAP_PR,
        Normalize(0, 1),
        grey_mask=osc_gate_grey_mask,
    )
    plot_phase(
        zmap("one_minus_max_std_share_B_osc_mean"),
        r_axis,
        h_axis,
        out_dir / "03c_one_minus_max_std_share_B_osc_mean",
        CMAP_PR,
        Normalize(0, 1),
        grey_mask=osc_gate_grey_mask,
    )

    metric_names = [
        "richness_last_absB",
        "richness_mean_last3_absB",
        "richness_cum3_absB",
        "shannon_last_absB",
        "shannon_mean_last3_absB",
        "shannon_cum3_absB",
    ]
    for metric_name in metric_names:
        delta = zmap(f"delta_mean_{metric_name}")
        cliff = zmap(f"rank_biserial_{metric_name}")
        cohen = zmap(f"cohens_d_{metric_name}")
        rho = zmap(f"spearman_rho_{metric_name}")

        plot_phase(
            delta,
            r_axis,
            h_axis,
            out_dir / "02_delta_mean" / f"02_delta_mean_{metric_name}",
            "RdBu_r",
            symmetric_norm(delta),
            grey_mask=osc_gate_grey_mask,
        )
        plot_phase(
            cliff,
            r_axis,
            h_axis,
            out_dir / "02_rank_biserial" / f"02_rank_biserial_{metric_name}",
            "RdBu_r",
            symmetric_norm(cliff, min_abs=0.25),
            grey_mask=osc_gate_grey_mask,
        )
        plot_phase(
            cohen,
            r_axis,
            h_axis,
            out_dir / "02_cohens_d" / f"02_cohens_d_{metric_name}",
            "RdBu_r",
            symmetric_norm(cohen, min_abs=0.5),
            grey_mask=osc_gate_grey_mask,
        )
        plot_phase(
            rho,
            r_axis,
            h_axis,
            out_dir / "02_spearman_cv_diversity" / f"02_spearman_cv_{metric_name}",
            "RdBu_r",
            symmetric_norm(rho, min_abs=0.25),
            grey_mask=(~np.isfinite(rho)) | (~osc_frac_gate),
        )


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot stability-diversity phase maps.")
    parser.add_argument("--run-dir", default=str(DEFAULT_RUN_DIR))
    parser.add_argument("--out-dir", default=str(DEFAULT_OUT_DIR))
    parser.add_argument("--b-thresh", type=float, default=1e-3)
    parser.add_argument("--cv-thresh", type=float, default=1e-1)
    parser.add_argument("--collapse-thresh", type=float, default=1e-3)
    parser.add_argument("--min-group", type=int, default=10)
    parser.add_argument("--osc-frac-min", type=float, default=0.10)
    parser.add_argument("--smooth-sigma", type=float, default=0.0)
    args = parser.parse_args()

    run_dir = Path(args.run_dir)
    if not run_dir.is_absolute():
        run_dir = ROOT / run_dir
    out_dir = Path(args.out_dir)
    if not out_dir.is_absolute():
        out_dir = ROOT / out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    if (run_dir / "CR_phase2_meta.h5").exists() and (run_dir / "phase2_rows").exists():
        print("Detected legacy phase2 format.")
        results = compute_maps_legacy_phase2(
            run_dir=run_dir,
            b_thresh=args.b_thresh,
            cv_thresh=args.cv_thresh,
            collapse_thresh=args.collapse_thresh,
            min_group=args.min_group,
        )
    else:
        print("Detected reorganized raw format.")
        results = compute_maps(
            run_dir=run_dir,
            b_thresh=args.b_thresh,
            cv_thresh=args.cv_thresh,
            collapse_thresh=args.collapse_thresh,
            min_group=args.min_group,
        )
    results["osc_frac_min"] = float(args.osc_frac_min)
    results["smooth_sigma"] = float(args.smooth_sigma)
    suffix = cv_suffix(args.cv_thresh)
    np.savez_compressed(out_dir / f"stability_diversity_metrics_{suffix}_absB.npz", **results)
    save_plots(results, out_dir, suffix)
    print(f"Saved figures and metrics to: {out_dir}")


if __name__ == "__main__":
    main()
