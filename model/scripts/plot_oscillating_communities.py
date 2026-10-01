#!/usr/bin/env python
"""Plot all oscillating communities from the S-scaling HDF5 outputs.

The figure style follows ``archive/old_figures/serials/cr_model_plot.py``:

    fig1: B_i at end of each cycle, log y
    fig2: N_i at end of each cycle, log y
    fig3: pH at end of each cycle

Only saved end-of-cycle points are used. No ODE re-integration is performed.
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "1")

ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np


DEFAULT_RUN_DIR = ROOT / "data" / "runs" / "s_scaling_r1_h005"
DEFAULT_OUT_DIR = ROOT / "figures" / "s_scaling_oscillating_serials"

FIG_SIZE_MM = 50.0
FIG_SIZE_IN = FIG_SIZE_MM / 25.4
FONT_SIZE = 9
SPINE_LW = 1.4
TICK_LW = 1.1
TICK_LEN_MAJ = 4.0
TICK_LEN_MIN = 2.5
LW_FIG1 = 1.2
MS_FIG1 = 2.0
LW_FIG2 = 1.2
LW_FIG3 = 1.2
MS_FIG3 = 2.0
PH_COLOR = "#1f3a68"
PH_REF_COLOR = "#999999"

Y_MIN = 1e-6
Y_MAX = 1e1
PH_MIN = 3.0
PH_MAX = 11.0
PH_TICKS = [3, 7, 11]


def endpoint_metrics(x: np.ndarray, cv_thresh: float, collapse_thresh: float,
                     last_k: int) -> tuple[np.ndarray, np.ndarray]:
    xw = x[..., -last_k:]
    per_cycle_total = xw.sum(axis=-2)
    total = per_cycle_total.mean(axis=-1)
    is_collapsed = total < collapse_thresh
    sp_std = xw.std(axis=-1, ddof=1)
    safe_total = np.where(total > 0, total, 1.0)
    community_cv = np.where(is_collapsed, 0.0, sp_std.sum(axis=-1) / safe_total)
    is_osc = (~is_collapsed) & (community_cv >= cv_thresh)
    return is_osc, community_cv


def style_axis(ax):
    ax.tick_params(which="major", direction="in",
                   labelsize=FONT_SIZE, width=TICK_LW, length=TICK_LEN_MAJ)
    ax.tick_params(which="minor", direction="in",
                   labelsize=FONT_SIZE, width=TICK_LW, length=TICK_LEN_MIN)
    for sp in ax.spines.values():
        sp.set_linewidth(SPINE_LW)
        sp.set_color("black")
    ax.set_facecolor("white")


def new_canvas():
    fig, ax = plt.subplots(figsize=(FIG_SIZE_IN, FIG_SIZE_IN))
    fig.patch.set_facecolor("white")
    fig.subplots_adjust(left=0.22, right=0.96, bottom=0.20, top=0.96)
    return fig, ax


def save_figure(fig, out_dir: Path, base: str, suffix: str,
                save_png: bool, save_pdf: bool, dpi: int) -> None:
    if save_pdf:
        fig.savefig(out_dir / f"{base}{suffix}.pdf", facecolor="white")
    if save_png:
        fig.savefig(out_dir / f"{base}{suffix}.png", facecolor="white", dpi=dpi)
    plt.close(fig)


def make_figs(
    B_hist: np.ndarray,
    N_hist: np.ndarray,
    p_end: np.ndarray,
    base: str,
    out_dir: Path,
    save_png: bool,
    save_pdf: bool,
    dpi: int,
) -> None:
    plt.rcParams["font.family"] = "Arial"
    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["ps.fonttype"] = 42

    S, n_cycles = B_hist.shape
    days = np.arange(1, n_cycles + 1)
    p_env = 7.0
    cmap = plt.get_cmap("tab10", S) if S <= 10 else plt.get_cmap("tab20", S)
    colors = [cmap(i % cmap.N) for i in range(S)]

    fig1, ax1 = new_canvas()
    ax1.axhline(1e-3, color=PH_REF_COLOR, lw=0.8, ls="--", zorder=2)
    for i in range(S):
        mask = B_hist[i] > 0
        if mask.any():
            ax1.semilogy(days[mask], B_hist[i, mask],
                         "-o", ms=MS_FIG1, color=colors[i],
                         markeredgewidth=0, alpha=0.95, lw=LW_FIG1, zorder=3)
    ax1.set_xlim(0, n_cycles + 1)
    ax1.set_ylim(Y_MIN, Y_MAX)
    ax1.xaxis.set_major_locator(ticker.FixedLocator([0, n_cycles // 2, n_cycles]))
    ax1.yaxis.set_major_locator(ticker.FixedLocator([1e-6, 1e-3, 1e0]))
    ax1.yaxis.set_minor_locator(ticker.NullLocator())
    style_axis(ax1)
    save_figure(fig1, out_dir, base, "_fig1_B_end", save_png, save_pdf, dpi)

    fig2, ax2 = new_canvas()
    ax2.axhline(1e-3, color=PH_REF_COLOR, lw=0.8, ls="--", zorder=2)
    for i in range(S):
        mask = N_hist[i] > 0
        if mask.any():
            ax2.semilogy(days[mask], N_hist[i, mask],
                         "-o", ms=MS_FIG1, color=colors[i],
                         markeredgewidth=0, alpha=0.95, lw=LW_FIG2, zorder=3)
    ax2.set_xlim(0, n_cycles + 1)
    ax2.set_ylim(Y_MIN, Y_MAX)
    ax2.xaxis.set_major_locator(ticker.FixedLocator([0, n_cycles // 2, n_cycles]))
    ax2.yaxis.set_major_locator(ticker.FixedLocator([1e-6, 1e-3, 1e0]))
    ax2.yaxis.set_minor_locator(ticker.NullLocator())
    style_axis(ax2)
    save_figure(fig2, out_dir, base, "_fig2_N_end", save_png, save_pdf, dpi)

    fig3, ax3 = new_canvas()
    ax3.axhline(p_env, color=PH_REF_COLOR, lw=0.8, ls="--", zorder=2)
    ax3.plot(days, p_end, "-o", ms=MS_FIG3, color=PH_COLOR,
             markeredgewidth=0, lw=LW_FIG3, zorder=3)
    ax3.set_xlim(0, n_cycles + 1)
    ax3.set_ylim(PH_MIN, PH_MAX)
    ax3.xaxis.set_major_locator(ticker.FixedLocator([0, n_cycles // 2, n_cycles]))
    ax3.yaxis.set_major_locator(ticker.FixedLocator(PH_TICKS))
    ax3.yaxis.set_minor_locator(ticker.NullLocator())
    style_axis(ax3)
    save_figure(fig3, out_dir, base, "_fig3_pH_end", save_png, save_pdf, dpi)


def plot_file(path: Path, out_root: Path, args) -> int:
    count = 0
    with h5py.File(path, "r") as h5:
        S = int(h5.attrs["S"])
        n_exp = int(h5.attrs["n_exp"])
        B = h5["B_end"]
        N = h5["N_end"]
        p = h5["p_end"]
        seeds = h5["seed"][...]

        if args.use_saved_labels and "metrics/B/is_oscillating" in h5:
            is_osc = h5["metrics/B/is_oscillating"][...].astype(bool)
            cv = h5["metrics/B/community_CV"][...]
        else:
            is_osc, cv = endpoint_metrics(
                B[...], cv_thresh=args.cv_thresh,
                collapse_thresh=args.collapse_thresh,
                last_k=args.metric_last_k,
            )

        idxs = np.flatnonzero(is_osc)
        if args.max_per_s is not None:
            idxs = idxs[:args.max_per_s]
        if idxs.size == 0:
            print(f"[skip] S={S}: no oscillating communities")
            return 0

        out_dir = out_root / f"S_{S:04d}"
        out_dir.mkdir(parents=True, exist_ok=True)

        for ie in idxs:
            B_hist = B[ie, :, :].astype(float)
            N_hist = N[ie, :, :].astype(float)
            p_end = p[ie, :].astype(float)

            surv = int((B_hist[:, -1] > args.alive_thresh).sum())
            base = (
                f"S{S:04d}_n{n_exp}_ens{ie:03d}_seed{int(seeds[ie])}"
                f"_cv{cv[ie]:.3f}_surv{surv:03d}"
            )
            make_figs(
                B_hist, N_hist, p_end, base, out_dir,
                save_png=args.png, save_pdf=args.pdf, dpi=args.dpi,
            )
            count += 1
            print(f"[OK] S={S} ens={ie} seed={int(seeds[ie])} cv={cv[ie]:.3f}")
    return count


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Plot all oscillating communities from S-scaling HDF5 files.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--run-dir", default=str(DEFAULT_RUN_DIR))
    p.add_argument("--out-dir", default=str(DEFAULT_OUT_DIR))
    p.add_argument("--cv-thresh", type=float, default=0.01)
    p.add_argument("--collapse-thresh", type=float, default=1e-3)
    p.add_argument("--metric-last-k", type=int, default=3)
    p.add_argument("--alive-thresh", type=float, default=1e-3)
    p.add_argument("--max-per-s", type=int, default=None)
    p.add_argument("--n-index-env", default=None,
                   help="Read exponent n from env var and only plot S=2^n.")
    p.add_argument("--n-values", default=None,
                   help="Comma/range list of exponents n to plot, e.g. 5-8.")
    p.add_argument("--use-saved-labels", action="store_true",
                   help="Use metrics/B/is_oscillating if present.")
    p.add_argument("--png", action="store_true", default=True)
    p.add_argument("--no-png", dest="png", action="store_false")
    p.add_argument("--pdf", action="store_true", default=True)
    p.add_argument("--no-pdf", dest="pdf", action="store_false")
    p.add_argument("--dpi", type=int, default=300)
    return p.parse_args()


def parse_n_values(text: str | None) -> set[int] | None:
    if text is None:
        return None
    out: set[int] = set()
    for part in text.split(","):
        part = part.strip()
        if not part:
            continue
        if "-" in part:
            a, b = part.split("-", 1)
            out.update(range(int(a), int(b) + 1))
        else:
            out.add(int(part))
    return out


def main() -> None:
    args = parse_args()
    run_dir = Path(args.run_dir)
    out_dir = Path(args.out_dir)
    if not run_dir.is_absolute():
        run_dir = ROOT / run_dir
    if not out_dir.is_absolute():
        out_dir = ROOT / out_dir

    raw_dir = run_dir / "raw"
    files = sorted(raw_dir.glob("S_*.h5"))
    if not files:
        raise FileNotFoundError(f"No S_*.h5 files found in {raw_dir}")

    if args.n_index_env:
        env_val = os.environ.get(args.n_index_env)
        if env_val is None:
            raise ValueError(f"Env var not set: {args.n_index_env}")
        wanted_n = {int(env_val)}
    else:
        wanted_n = parse_n_values(args.n_values)

    if wanted_n is not None:
        wanted_s = {2 ** n for n in wanted_n}
        files = [p for p in files if int(p.stem.split("_", 1)[1]) in wanted_s]
        if not files:
            raise FileNotFoundError(
                f"No matching S files for n={sorted(wanted_n)} in {raw_dir}"
            )

    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"Input : {raw_dir}")
    print(f"Output: {out_dir}")
    print(f"Files : {len(files)}")
    if wanted_n is not None:
        print(f"n     : {sorted(wanted_n)}")

    total = 0
    for path in files:
        total += plot_file(path, out_dir, args)
    print(f"Done. Plotted {total} oscillating communities.")


if __name__ == "__main__":
    main()
