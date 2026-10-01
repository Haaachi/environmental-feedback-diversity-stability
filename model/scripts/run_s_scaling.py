#!/usr/bin/env python
"""Run S-scaling simulations at fixed r and h.

This script is intentionally small and HPC-friendly:

    python scripts/run_s_scaling.py --workers $SLURM_CPUS_PER_TASK

For SLURM arrays, use array indices 0-8:

    sbatch --array=0-8 ... --wrap="python scripts/run_s_scaling.py \
        --n-index-env SLURM_ARRAY_TASK_ID --workers $SLURM_CPUS_PER_TASK"

Each S value is written to a separate HDF5 file under
``data/runs/s_scaling_r1_h005/raw/``.  This avoids concurrent HDF5 writes when
using a job array.  Run the script again with ``--summary-only`` after an array
finishes to regenerate the global summary CSV.
"""

from __future__ import annotations

import argparse
import json
import os
import socket
import sys
import time
from datetime import datetime
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "1")

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = PROJECT_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

import h5py
import numpy as np
from joblib import Parallel, delayed
from scipy.integrate import solve_ivp

from crmodel.io.config import load_config
from crmodel.model.core import _make_rhs_closure
from crmodel.model.resources import build_resource_vectors
from crmodel.model.traits import draw_quenched_traits


DEFAULT_RUN_DIR = PROJECT_ROOT / "data" / "runs" / "s_scaling_r1_h005"


def _slurm_workers() -> int:
    if "SLURM_CPUS_PER_TASK" in os.environ:
        return max(1, int(os.environ["SLURM_CPUS_PER_TASK"]))
    return max(1, (os.cpu_count() or 2) - 1)


def _compression_kwargs(name: str) -> dict:
    if name == "none":
        return {}
    if name == "lzf":
        return {"compression": "lzf"}
    if name == "gzip":
        return {"compression": "gzip", "compression_opts": 4}
    raise ValueError(f"Unsupported compression: {name}")


def _parse_n_values(text: str | None, n_min: int, n_max: int) -> list[int]:
    if text is None:
        return list(range(n_min, n_max + 1))
    out: list[int] = []
    for part in text.split(","):
        part = part.strip()
        if not part:
            continue
        if "-" in part:
            a, b = part.split("-", 1)
            out.extend(range(int(a), int(b) + 1))
        else:
            out.append(int(part))
    bad = [n for n in out if n < n_min or n > n_max]
    if bad:
        raise ValueError(f"n values out of range {n_min}..{n_max}: {bad}")
    return sorted(set(out))


def _build_ph_shape(cfg: dict, S: int) -> np.ndarray:
    mode = cfg.get("ph_shaping", "all")
    if mode == "all":
        return np.ones(S, dtype=np.float64)
    if mode == "first_only":
        ph_shape = np.zeros(S, dtype=np.float64)
        ph_shape[0] = 1.0
        return ph_shape
    raise ValueError(f"Unknown ph_shaping: {mode}")


def seed_for(seed_base: int, n_exp: int, i_ens: int) -> int:
    return int(seed_base + (n_exp + 1) * 1_000_000 + (i_ens + 1))


def simulate_one_community(
    seed: int,
    S: int,
    M: int,
    r_mean: float,
    h_stress: float,
    f: float,
    cfg: dict,
    R0: np.ndarray,
    gamma: np.ndarray,
    ph_shape: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Simulate one community and return all daily endpoints.

    Returns
    -------
    N_end : ndarray, shape (S, n_cycles)
    B_end : ndarray, shape (S, n_cycles)
        End-of-cycle biomass readout. Each cycle starts from
        ``B_i(0)=N_i(0)`` and then integrates ``dB_i/dt=N_i*g_i``.
    R_end : ndarray, shape (M, n_cycles)
    p_end : ndarray, shape (n_cycles,)
    p_opt : ndarray, shape (S,)
    v_max : ndarray, shape (S, M)
    """
    model = cfg["model"]
    traits = cfg["traits"]
    sim = cfg["simulation"]

    n_cycles = int(sim["n_cycles"])
    tau = float(model["tau"])
    Km = float(model["Km"])
    p_tol = float(model["p_tol"])
    delta = float(model["delta"])
    p_env = float(model["p_env"])
    D = float(model["D"])
    suicide_thresh = float(sim.get("suicide_thresh", 1e-7))
    rtol = float(sim["rtol"])
    atol = float(sim["atol"])
    max_step = tau * float(sim["max_step_frac"])

    v_max, p_opt = draw_quenched_traits(
        seed=seed,
        n_species=S,
        n_resources=M,
        r_mean=r_mean,
        growth_cv=float(traits.get("growth_cv", 0.0)),
        resource_degree=int(traits["resource_degree"]),
        p_opt_min=float(traits["p_opt_min"]),
        p_opt_max=float(traits["p_opt_max"]),
        preference_model=str(traits.get("preference_model", "sparse_random_links")),
    )

    rhs = _make_rhs_closure(
        S, M, Km, p_tol, h_stress, delta, p_env,
        v_max, p_opt, gamma, ph_shape,
    )

    N = 1e-3 * np.ones(S, dtype=np.float64)
    N_save = np.zeros((S, n_cycles), dtype=np.float64)
    B_save = np.zeros((S, n_cycles), dtype=np.float64)
    R_save = np.zeros((M, n_cycles), dtype=np.float64)
    p_save = np.zeros(n_cycles, dtype=np.float64)

    state_len = 2 * S + M + 1
    t_span = (0.0, tau)

    for nc in range(n_cycles):
        y0 = np.empty(state_len, dtype=np.float64)
        y0[:S] = N
        y0[S:S + M] = R0
        y0[S + M] = p_env
        # B is the per-cycle biomass readout.  Initialize at the diluted
        # inoculum N(0), then add gross within-cycle production.
        y0[S + M + 1:] = N

        sol = solve_ivp(
            rhs, t_span, y0,
            method="LSODA",
            rtol=rtol, atol=atol, max_step=max_step,
            dense_output=False, vectorized=False,
        )
        if not sol.success:
            sol = solve_ivp(
                rhs, t_span, y0,
                method="BDF",
                rtol=rtol, atol=atol, max_step=max_step,
            )
            if not sol.success:
                y_end = y0.copy()
                y_end[:S] = 0.0
            else:
                y_end = sol.y[:, -1]
        else:
            y_end = sol.y[:, -1]

        N_end = y_end[:S].copy()
        R_end = y_end[S:S + M].copy()
        B_end = y_end[S + M + 1:].copy()
        p_end = float(y_end[S + M])

        N_end[N_end < suicide_thresh] = 0.0
        R_end[R_end < 0.0] = 0.0

        N_save[:, nc] = N_end
        B_save[:, nc] = B_end
        R_save[:, nc] = R_end
        p_save[nc] = p_end

        N = f * (N_end + D)
        N[N < 1e-12] = 0.0

    return N_save, B_save, R_save, p_save, p_opt, v_max


def endpoint_metrics(
    x: np.ndarray,
    cv_thresh: float,
    collapse_thresh: float,
    last_k: int = 3,
) -> dict[str, np.ndarray | float]:
    """Compute per-community last-k endpoint fluctuation metrics."""
    xw = x[..., -last_k:]
    per_cycle_total = xw.sum(axis=-2)
    total = per_cycle_total.mean(axis=-1)
    is_collapsed = total < collapse_thresh

    sp_std = xw.std(axis=-1, ddof=1)
    safe_total = np.where(total > 0, total, 1.0)
    community_cv = np.where(is_collapsed, 0.0, sp_std.sum(axis=-1) / safe_total)
    is_osc = (~is_collapsed) & (community_cv >= cv_thresh)

    return {
        "total": total,
        "is_collapsed": is_collapsed,
        "community_CV": community_cv,
        "is_oscillating": is_osc,
        "osc_frac": float(np.mean(is_osc)),
        "collapse_frac": float(np.mean(is_collapsed)),
        "mean_community_CV": float(np.mean(community_cv)),
        "median_community_CV": float(np.median(community_cv)),
    }


def write_one_s(
    n_exp: int,
    args: argparse.Namespace,
    cfg: dict,
    R0: np.ndarray,
    gamma: np.ndarray,
) -> Path:
    S = 2 ** n_exp
    M = int(cfg["model"]["M"])
    run_dir = Path(args.run_dir)
    raw_dir = run_dir / "raw"
    raw_dir.mkdir(parents=True, exist_ok=True)
    out_path = raw_dir / f"S_{S:04d}.h5"

    if out_path.exists() and not args.overwrite:
        print(f"[skip] S={S}: {out_path} exists", flush=True)
        return out_path
    if out_path.exists() and args.overwrite:
        out_path.unlink()

    dtype = np.float32 if args.dtype == "float32" else np.float64
    n_ens = int(args.n_ensemble)
    n_cycles = int(args.n_cycles)
    seed_base = int(args.seed_base)
    seeds = np.array(
        [seed_for(seed_base, n_exp, i) for i in range(n_ens)],
        dtype=np.int64,
    )
    ph_shape = _build_ph_shape(cfg, S)

    print(
        f"[run] n={n_exp}, S={S}, ensembles={n_ens}, "
        f"cycles={n_cycles}, workers={args.workers}",
        flush=True,
    )
    t0 = time.time()

    results = Parallel(n_jobs=args.workers, backend="loky", verbose=0)(
        delayed(simulate_one_community)(
            int(seed), S, M, float(args.r_mean), float(args.h_stress),
            float(args.f), cfg, R0, gamma, ph_shape,
        )
        for seed in seeds
    )

    N_end = np.zeros((n_ens, S, n_cycles), dtype=dtype)
    B_end = np.zeros((n_ens, S, n_cycles), dtype=dtype)
    R_end = np.zeros((n_ens, M, n_cycles), dtype=dtype)
    p_end = np.zeros((n_ens, n_cycles), dtype=dtype)
    p_opt = np.zeros((n_ens, S), dtype=dtype)
    v_max = np.zeros((n_ens, S, M), dtype=dtype)

    for i, (Ns, Bs, Rs, ps, po, vm) in enumerate(results):
        N_end[i] = Ns
        B_end[i] = Bs
        R_end[i] = Rs
        p_end[i] = ps
        p_opt[i] = po
        v_max[i] = vm

    metrics_B = endpoint_metrics(
        B_end, cv_thresh=float(args.cv_thresh),
        collapse_thresh=float(args.collapse_thresh),
        last_k=int(args.metric_last_k),
    )
    metrics_N = endpoint_metrics(
        N_end, cv_thresh=float(args.cv_thresh),
        collapse_thresh=float(args.collapse_thresh),
        last_k=int(args.metric_last_k),
    )
    metrics_B_alt = endpoint_metrics(
        B_end, cv_thresh=float(args.alt_cv_thresh),
        collapse_thresh=float(args.collapse_thresh),
        last_k=int(args.metric_last_k),
    )

    ck = _compression_kwargs(args.compression)
    tmp_path = raw_dir / f".S_{S:04d}.{os.getpid()}.tmp.h5"
    if tmp_path.exists():
        tmp_path.unlink()

    with h5py.File(tmp_path, "w") as h5:
        h5.attrs["created_at"] = datetime.now().isoformat(timespec="seconds")
        h5.attrs["hostname"] = socket.gethostname()
        h5.attrs["command"] = " ".join(sys.argv)
        h5.attrs["config_path"] = str(args.config)
        h5.attrs["config_json"] = json.dumps(cfg, sort_keys=True, default=str)
        h5.attrs["n_exp"] = n_exp
        h5.attrs["S"] = S
        h5.attrs["M"] = M
        h5.attrs["r_mean"] = float(args.r_mean)
        h5.attrs["h_stress"] = float(args.h_stress)
        h5.attrs["f"] = float(args.f)
        h5.attrs["n_ensemble"] = n_ens
        h5.attrs["n_cycles"] = n_cycles
        h5.attrs["metric_last_k"] = int(args.metric_last_k)
        h5.attrs["cv_thresh"] = float(args.cv_thresh)
        h5.attrs["alt_cv_thresh"] = float(args.alt_cv_thresh)
        h5.attrs["collapse_thresh"] = float(args.collapse_thresh)
        h5.attrs["B_definition"] = "B_i(0)=N_i(0); dB_i/dt=N_i*g_i"
        h5.attrs["osc_frac_B"] = metrics_B["osc_frac"]
        h5.attrs["osc_frac_N"] = metrics_N["osc_frac"]
        h5.attrs["osc_frac_B_alt"] = metrics_B_alt["osc_frac"]

        h5.create_dataset("seed", data=seeds, **ck)
        h5.create_dataset("R0", data=R0.astype(dtype), **ck)
        h5.create_dataset("gamma", data=gamma.astype(dtype), **ck)
        h5.create_dataset("ph_shape", data=ph_shape.astype(dtype), **ck)
        h5.create_dataset("p_opt", data=p_opt, **ck)
        h5.create_dataset("v_max", data=v_max, **ck)
        h5.create_dataset("N_end", data=N_end, **ck)
        h5.create_dataset("B_end", data=B_end, **ck)
        h5.create_dataset("R_end", data=R_end, **ck)
        h5.create_dataset("p_end", data=p_end, **ck)

        mg = h5.create_group("metrics")
        for prefix, metrics in (
            ("B", metrics_B),
            ("N", metrics_N),
            ("B_alt", metrics_B_alt),
        ):
            g = mg.create_group(prefix)
            for key, value in metrics.items():
                if np.isscalar(value):
                    g.attrs[key] = value
                else:
                    g.create_dataset(key, data=value, **ck)

    tmp_path.replace(out_path)
    print(
        f"[done] S={S}: osc_frac_B={metrics_B['osc_frac']:.3f}, "
        f"osc_frac_N={metrics_N['osc_frac']:.3f}, "
        f"elapsed={time.time() - t0:.1f}s",
        flush=True,
    )
    return out_path


def write_summary(run_dir: Path) -> Path:
    rows = []
    raw_dir = run_dir / "raw"
    for path in sorted(raw_dir.glob("S_*.h5")):
        with h5py.File(path, "r") as h5:
            row = {
                "n": int(h5.attrs["n_exp"]),
                "S": int(h5.attrs["S"]),
                "n_ensemble": int(h5.attrs["n_ensemble"]),
                "n_cycles": int(h5.attrs["n_cycles"]),
                "r_mean": float(h5.attrs["r_mean"]),
                "h_stress": float(h5.attrs["h_stress"]),
                "f": float(h5.attrs["f"]),
                "cv_thresh": float(h5.attrs["cv_thresh"]),
                "alt_cv_thresh": float(h5.attrs["alt_cv_thresh"]),
                "osc_frac_B": float(h5["metrics/B"].attrs["osc_frac"]),
                "osc_frac_N": float(h5["metrics/N"].attrs["osc_frac"]),
                "osc_frac_B_alt": float(h5["metrics/B_alt"].attrs["osc_frac"]),
                "collapse_frac_B": float(h5["metrics/B"].attrs["collapse_frac"]),
                "collapse_frac_N": float(h5["metrics/N"].attrs["collapse_frac"]),
                "mean_community_CV_B": float(h5["metrics/B"].attrs["mean_community_CV"]),
                "median_community_CV_B": float(h5["metrics/B"].attrs["median_community_CV"]),
            }
            rows.append(row)

    summary_path = run_dir / "summary_s_scaling.csv"
    if not rows:
        print(f"[summary] no S_*.h5 files found in {raw_dir}", flush=True)
        return summary_path

    import csv

    rows.sort(key=lambda r: r["S"])
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    with open(summary_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    print(f"[summary] wrote {summary_path}", flush=True)
    for row in rows:
        print(
            f"  S={row['S']:4d}  osc_frac_B={row['osc_frac_B']:.3f}  "
            f"collapse_frac_B={row['collapse_frac_B']:.3f}",
            flush=True,
        )
    return summary_path


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Scan S=2^n at fixed r=1, h=0.05 and save all endpoints.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--config", default=str(PROJECT_ROOT / "configs" / "base.yaml"))
    p.add_argument("--run-dir", default=str(DEFAULT_RUN_DIR))
    p.add_argument("--n-min", type=int, default=0)
    p.add_argument("--n-max", type=int, default=8)
    p.add_argument("--n-values", default=None, help="Comma/range list, e.g. 0,1,4-8")
    p.add_argument("--n-index-env", default=None,
                   help="Read one exponent n from an env var, e.g. SLURM_ARRAY_TASK_ID")
    p.add_argument("--summary-only", action="store_true")
    p.add_argument("--no-summary", action="store_true")
    p.add_argument("--workers", type=int, default=_slurm_workers())
    p.add_argument("--n-ensemble", type=int, default=100)
    p.add_argument("--n-cycles", type=int, default=30)
    p.add_argument("--seed-base", type=int, default=0)
    p.add_argument("--r-mean", type=float, default=1.0)
    p.add_argument("--h-stress", type=float, default=0.05)
    p.add_argument("--f", type=float, default=1.0 / 30.0)
    p.add_argument("--cv-thresh", type=float, default=0.01)
    p.add_argument("--alt-cv-thresh", type=float, default=0.1)
    p.add_argument("--collapse-thresh", type=float, default=1e-3)
    p.add_argument("--metric-last-k", type=int, default=3)
    p.add_argument("--dtype", choices=["float32", "float64"], default="float32")
    p.add_argument("--compression", choices=["lzf", "gzip", "none"], default="lzf")
    p.add_argument("--overwrite", action="store_true")
    p.add_argument("--dry-run", action="store_true")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    run_dir_arg = Path(args.run_dir)
    config_arg = Path(args.config)
    if not run_dir_arg.is_absolute():
        run_dir_arg = PROJECT_ROOT / run_dir_arg
    if not config_arg.is_absolute():
        config_arg = PROJECT_ROOT / config_arg
    args.run_dir = str(run_dir_arg.resolve())
    args.config = str(config_arg.resolve())

    run_dir = Path(args.run_dir)
    if args.summary_only:
        write_summary(run_dir)
        return

    cfg = load_config(args.config)
    cfg["model"]["S"] = None
    cfg["simulation"]["n_ensemble"] = int(args.n_ensemble)
    cfg["simulation"]["n_cycles"] = int(args.n_cycles)
    cfg["simulation"]["n_saved_days"] = int(args.n_cycles)
    cfg["simulation"]["seed_base"] = int(args.seed_base)

    M = int(cfg["model"]["M"])
    R0, gamma = build_resource_vectors(cfg["resources"], M)

    if args.n_index_env:
        env_val = os.environ.get(args.n_index_env)
        if env_val is None:
            raise ValueError(f"Env var not set: {args.n_index_env}")
        n_values = [int(env_val)]
    else:
        n_values = _parse_n_values(args.n_values, args.n_min, args.n_max)

    print("=" * 72)
    print("S-scaling CR-pH scan")
    print(f"  config       : {args.config}")
    print(f"  run_dir      : {run_dir}")
    print(f"  n values     : {n_values}")
    print(f"  S values     : {[2 ** n for n in n_values]}")
    print(f"  r_mean/h/f   : {args.r_mean} / {args.h_stress} / {args.f}")
    print(f"  M            : {M}")
    print(f"  gamma range  : [{gamma.min():.4g}, {gamma.max():.4g}]")
    print(f"  R0 unique    : {sorted(set(np.round(R0, 12)))}")
    print(f"  ensembles    : {args.n_ensemble}")
    print(f"  cycles saved : {args.n_cycles}")
    print(f"  cv threshold : {args.cv_thresh} (alt {args.alt_cv_thresh})")
    print(f"  workers      : {args.workers}")
    print("=" * 72)

    if args.dry_run:
        return

    run_dir.mkdir(parents=True, exist_ok=True)
    with open(run_dir / "run_config_s_scaling.json", "w") as f:
        json.dump(
            {
                "command": " ".join(sys.argv),
                "config": cfg,
                "n_values": n_values,
                "S_values": [2 ** n for n in n_values],
                "r_mean": args.r_mean,
                "h_stress": args.h_stress,
                "f": args.f,
                "cv_thresh": args.cv_thresh,
                "alt_cv_thresh": args.alt_cv_thresh,
            },
            f,
            indent=2,
            sort_keys=True,
            default=str,
        )

    for n_exp in n_values:
        write_one_s(n_exp, args, cfg, R0, gamma)

    if not args.no_summary and not args.n_index_env:
        write_summary(run_dir)
    elif args.n_index_env:
        print(
            "[summary] skipped because --n-index-env is set; "
            "run with --summary-only after the array finishes.",
            flush=True,
        )


if __name__ == "__main__":
    main()
