#!/usr/bin/env python
"""Plot community-CV histograms from S-scaling HDF5 outputs."""

from __future__ import annotations

import argparse
import csv
import os
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_RUN_DIR = ROOT / "data" / "runs" / "s_scaling_r1_h005"
DEFAULT_OUT_DIR = ROOT / "figures" / "s_scaling_cv_hist"

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def compute_community_cv(
    x: np.ndarray,
    collapse_thresh: float,
    last_k: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Return community_CV and collapsed flags for x=(ens,S,T)."""
    xw = x[..., -last_k:]
    per_cycle_total = xw.sum(axis=-2)
    total = per_cycle_total.mean(axis=-1)
    collapsed = total < collapse_thresh
    sp_std = xw.std(axis=-1, ddof=1)
    safe_total = np.where(total > 0, total, 1.0)
    cv = np.where(collapsed, 0.0, sp_std.sum(axis=-1) / safe_total)
    return cv, collapsed


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


def load_cv_rows(args) -> list[dict]:
    run_dir = Path(args.run_dir)
    if not run_dir.is_absolute():
        run_dir = ROOT / run_dir
    files = sorted((run_dir / "raw").glob("S_*.h5"))
    if not files:
        raise FileNotFoundError(f"No S_*.h5 files found under {run_dir / 'raw'}")

    wanted_n = parse_n_values(args.n_values)
    if wanted_n is not None:
        wanted_s = {2 ** n for n in wanted_n}
        files = [p for p in files if int(p.stem.split("_", 1)[1]) in wanted_s]

    rows = []
    for path in files:
        with h5py.File(path, "r") as h5:
            S = int(h5.attrs["S"])
            n = int(h5.attrs["n_exp"])
            if args.readout == "B":
                x = h5["B_end"][...].astype(float)
            elif args.readout == "N":
                x = h5["N_end"][...].astype(float)
            else:
                raise ValueError(f"Unknown readout: {args.readout}")
            cv, collapsed = compute_community_cv(
                x, collapse_thresh=args.collapse_thresh, last_k=args.last_k)
            rows.append({
                "path": path,
                "S": S,
                "n": n,
                "cv": cv,
                "collapsed": collapsed,
            })
    rows.sort(key=lambda r: r["S"])
    return rows


def write_quantiles(rows: list[dict], out_dir: Path, readout: str) -> Path:
    q_values = [0, 0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 1.0]
    out_path = out_dir / f"community_cv_{readout}_quantiles.csv"
    with open(out_path, "w", newline="") as f:
        fieldnames = [
            "n", "S", "n_total", "n_noncollapsed", "collapse_frac",
            "mean", "median", "max",
        ] + [f"q{int(q * 100):03d}" for q in q_values]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            cv = row["cv"]
            collapsed = row["collapsed"]
            valid = cv[~collapsed]
            if valid.size == 0:
                valid = np.array([np.nan])
            record = {
                "n": row["n"],
                "S": row["S"],
                "n_total": int(cv.size),
                "n_noncollapsed": int((~collapsed).sum()),
                "collapse_frac": float(collapsed.mean()),
                "mean": float(np.nanmean(valid)),
                "median": float(np.nanmedian(valid)),
                "max": float(np.nanmax(valid)),
            }
            for q in q_values:
                record[f"q{int(q * 100):03d}"] = float(np.nanquantile(valid, q))
            writer.writerow(record)
    return out_path


def plot_grid(rows: list[dict], args, out_dir: Path) -> Path:
    n_panels = len(rows)
    n_cols = min(3, n_panels)
    n_rows = int(np.ceil(n_panels / n_cols))
    fig, axes = plt.subplots(
        n_rows, n_cols,
        figsize=(3.2 * n_cols, 2.6 * n_rows),
        squeeze=False,
    )

    positive = np.concatenate([
        row["cv"][~row["collapsed"]][row["cv"][~row["collapsed"]] > 0]
        for row in rows
    ])
    if positive.size == 0:
        bins = np.linspace(0, 1, args.bins)
        log_x = False
    else:
        log_x = args.log_x
        if log_x:
            lo = max(float(np.nanmin(positive)), 1e-5)
            hi = max(float(np.nanmax(positive)), lo * 1.01)
            bins = np.logspace(np.log10(lo), np.log10(hi), args.bins)
        else:
            hi = float(np.nanquantile(positive, 0.995))
            hi = max(hi, args.cv_thresh * 2, 0.02)
            bins = np.linspace(0, hi, args.bins)

    for ax, row in zip(axes.ravel(), rows):
        cv = row["cv"]
        values = cv[~row["collapsed"]] if args.exclude_collapsed else cv
        values = values[np.isfinite(values)]
        if args.log_x:
            values = values[values > 0]
        ax.hist(values, bins=bins, color="#4C78A8", alpha=0.85,
                edgecolor="white", linewidth=0.4)
        ax.axvline(args.cv_thresh, color="#D62728", lw=1.2, ls="--",
                   label=f"{args.cv_thresh:g}")
        if args.alt_cv_thresh is not None:
            ax.axvline(args.alt_cv_thresh, color="#333333", lw=1.0, ls=":")
        if args.log_x:
            ax.set_xscale("log")
        ax.set_title(
            f"S={row['S']}  collapse={row['collapsed'].mean():.2f}",
            fontsize=9,
        )
        ax.tick_params(labelsize=8)
        ax.set_ylabel("count", fontsize=8)
        ax.set_xlabel(f"community_CV_{args.readout}", fontsize=8)

    for ax in axes.ravel()[n_panels:]:
        ax.axis("off")

    fig.tight_layout()
    suffix = "logx" if args.log_x else "linear"
    excl = "noncollapsed" if args.exclude_collapsed else "all"
    out_path = out_dir / f"community_cv_{args.readout}_hist_{suffix}_{excl}.png"
    fig.savefig(out_path, dpi=args.dpi)
    plt.close(fig)
    return out_path


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Plot community-CV histograms for S-scaling outputs.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--run-dir", default=str(DEFAULT_RUN_DIR))
    p.add_argument("--out-dir", default=str(DEFAULT_OUT_DIR))
    p.add_argument("--readout", choices=["B", "N"], default="B")
    p.add_argument("--n-values", default=None, help="Exponent list/range, e.g. 0,4-8")
    p.add_argument("--last-k", type=int, default=3)
    p.add_argument("--collapse-thresh", type=float, default=1e-3)
    p.add_argument("--cv-thresh", type=float, default=0.01)
    p.add_argument("--alt-cv-thresh", type=float, default=0.1)
    p.add_argument("--bins", type=int, default=40)
    p.add_argument("--log-x", action="store_true")
    p.add_argument("--include-collapsed", dest="exclude_collapsed",
                   action="store_false")
    p.set_defaults(exclude_collapsed=True)
    p.add_argument("--dpi", type=int, default=300)
    return p.parse_args()


def main() -> None:
    args = parse_args()
    out_dir = Path(args.out_dir)
    if not out_dir.is_absolute():
        out_dir = ROOT / out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    rows = load_cv_rows(args)
    hist_path = plot_grid(rows, args, out_dir)
    q_path = write_quantiles(rows, out_dir, args.readout)
    print(f"Wrote histogram: {hist_path}")
    print(f"Wrote quantiles: {q_path}")


if __name__ == "__main__":
    main()
