#!/usr/bin/env python
"""Unified entry point for CR-pH feedback model simulations.

Usage
-----
Simulate one or more rows of a phase diagram::

    python run.py simulate --config configs/phase_r_h.yaml --rows 1-10 --workers 32

Post-process raw HDF5 row files into compressed NPZ metrics::

    python run.py postprocess --config configs/phase_r_h.yaml

Check completion status of a run::

    python run.py status --config configs/phase_r_h.yaml

Simulate with a sensitivity override::

    python run.py simulate --config configs/phase_r_h.yaml \\
        --overlay configs/sensitivity/growth_cv_03.yaml --workers 32

SLURM array job (one row per task)::

    python run.py simulate --config configs/phase_r_h.yaml \\
        --row-index-env SLURM_ARRAY_TASK_ID --workers $SLURM_CPUS_PER_TASK
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

PROJECT_ROOT = Path(__file__).resolve().parent
SRC_DIR = PROJECT_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))


# ── Config & scan helpers ─────────────────────────────────────────────

def load_and_merge(args) -> dict:
    from crmodel.io.config import load_config, deep_update
    cfg = load_config(args.config)
    if hasattr(args, "overlay") and args.overlay:
        overlay = load_config(args.overlay)
        overlay.pop("_config_path", None)
        overlay.pop("extends", None)
        cfg = deep_update(cfg, overlay)
        cfg["_overlay_path"] = str(args.overlay)
    for attr, keys in [
        ("n_ensemble", ("simulation", "n_ensemble")),
        ("n_cycles", ("simulation", "n_cycles")),
        ("n_saved_days", ("simulation", "n_saved_days")),
        ("seed_base", ("simulation", "seed_base")),
        ("growth_cv", ("traits", "growth_cv")),
        ("resource_degree", ("traits", "resource_degree")),
    ]:
        val = getattr(args, attr, None)
        if val is not None:
            node = cfg
            for k in keys[:-1]:
                node = node.setdefault(k, {})
            node[keys[-1]] = type(node.get(keys[-1], val))(val)
    return cfg


def build_scan(cfg: dict, n_grid_override: int | None = None) -> dict:
    import numpy as np
    scan_cfg = cfg["scan"]
    fixed = cfg.get("fixed", {})
    name = str(scan_cfg["name"])
    n_grid = n_grid_override or int(scan_cfg.get("n_grid", 100))

    if name == "phase_r_h":
        axis1 = np.linspace(float(scan_cfg["r_min"]),
                            float(scan_cfg["r_max"]), n_grid)
        axis2 = np.linspace(float(scan_cfg["h_min"]),
                            float(scan_cfg["h_max"]), n_grid)
        return dict(name=name, axis1_name="r_mean", axis2_name="h_stress",
                    axis1=axis1, axis2=axis2,
                    fixed_params={"f": float(fixed["f"])})

    if name == "phase_r_f":
        axis1 = np.linspace(float(scan_cfg["r_min"]),
                            float(scan_cfg["r_max"]), n_grid)
        axis2 = np.logspace(float(scan_cfg["log10_f_min"]),
                            float(scan_cfg["log10_f_max"]), n_grid)
        return dict(name=name, axis1_name="r_mean", axis2_name="f",
                    axis1=axis1, axis2=axis2,
                    fixed_params={"h_stress": float(fixed["h_stress"])})

    if name == "phase_h_f":
        axis1 = np.linspace(float(scan_cfg["h_min"]),
                            float(scan_cfg["h_max"]),
                            int(scan_cfg.get("n_h", n_grid)))
        axis2 = np.array(scan_cfg["f_values"], dtype=np.float64)
        return dict(name=name, axis1_name="h_stress", axis2_name="f",
                    axis1=axis1, axis2=axis2,
                    fixed_params={"r_mean": float(fixed["r_mean"])})

    raise ValueError(f"Unknown scan: {name}")


def resolve_run_dir(args, cfg: dict) -> Path:
    if args.run_dir:
        d = Path(args.run_dir)
    elif os.environ.get("CR_RUN_DIR"):
        d = Path(os.environ["CR_RUN_DIR"])
    else:
        d = Path(cfg.get("output", {}).get("run_dir", "data/runs/default"))
    if hasattr(args, "overlay") and args.overlay:
        overlay_stem = Path(args.overlay).stem
        d = d.parent / f"{d.name}_{overlay_stem}"
    if not d.is_absolute():
        d = PROJECT_ROOT / d
    return d


# ── Seed generation ───────────────────────────────────────────────────

def seed_for(seed_base: int, i1: int, i2: int, i_ens: int) -> int:
    return int(seed_base + (i1 + 1) * 1_000_000
               + (i2 + 1) * 10_000 + (i_ens + 1))


# ── Simulate ──────────────────────────────────────────────────────────

def _build_ph_shape(cfg: dict, S: int):
    import numpy as np
    mode = cfg.get("ph_shaping", "all")
    if mode == "all":
        return np.ones(S, dtype=np.float64)
    if mode == "first_only":
        ps = np.zeros(S, dtype=np.float64)
        ps[0] = 1.0
        return ps
    raise ValueError(f"Unknown ph_shaping: {mode}")


def _ph_response_mode(cfg: dict) -> int:
    mode = str(cfg.get("ph_response",
                       cfg.get("model", {}).get("ph_response",
                                               "triangular"))).lower()
    aliases = {
        "triangular": 0,
        "triangle": 0,
        "linear": 0,
        "tent": 0,
        "top_hat": 1,
        "tophat": 1,
        "boxcar": 1,
        "window": 1,
        "binary": 1,
    }
    if mode not in aliases:
        raise ValueError(f"Unknown ph_response: {mode}")
    return aliases[mode]


def run_one(seed, r_mean, params, cfg, R0, gamma, ph_shape):
    from crmodel.model.core import run_community_endpoints
    from crmodel.model.resources import draw_community_gamma, is_gamma_per_community
    from crmodel.model.traits import draw_quenched_traits
    model = cfg["model"]
    traits = cfg["traits"]
    sim = cfg["simulation"]
    S, M = int(model["S"]), int(model["M"])
    v_max, p_opt = draw_quenched_traits(
        seed=seed, n_species=S, n_resources=M,
        r_mean=float(r_mean),
        growth_cv=float(traits.get("growth_cv", 0.0)),
        resource_degree=int(traits["resource_degree"]),
        p_opt_min=float(traits["p_opt_min"]),
        p_opt_max=float(traits["p_opt_max"]),
        preference_model=str(traits.get("preference_model",
                                        "sparse_random_links")),
    )
    gamma_used = (
        draw_community_gamma(cfg["resources"], M, seed)
        if is_gamma_per_community(cfg["resources"])
        else gamma
    )
    N_save, B_save, p_save = run_community_endpoints(
        S=S, M=M,
        Km=float(model["Km"]),
        p_tol=float(model["p_tol"]),
        h_stress=float(params["h_stress"]),
        delta=float(model["delta"]),
        p_env=float(model["p_env"]),
        D=float(model["D"]),
        tau=float(model["tau"]),
        f=float(params["f"]),
        R0=R0, gamma=gamma_used, v_max=v_max, p_opt=p_opt,
        n_cycles=int(sim["n_cycles"]),
        n_save=int(sim["n_saved_days"]),
        suicide_thresh=float(sim.get("suicide_thresh", 1e-6)),
        rtol=float(sim["rtol"]),
        atol=float(sim["atol"]),
        max_step=float(model["tau"]) * float(sim["max_step_frac"]),
        ph_shape=ph_shape,
        ph_response_mode=_ph_response_mode(cfg),
    )
    return N_save, B_save, p_save, gamma_used


def cmd_simulate(args):
    import h5py
    import numpy as np
    from joblib import Parallel, delayed
    from crmodel.model.resources import build_resource_vectors, is_gamma_per_community
    from crmodel.io.manifest import snapshot_run

    cfg = load_and_merge(args)
    scan = build_scan(cfg, args.n_grid)
    run_dir = resolve_run_dir(args, cfg)
    rows_dir = run_dir / "raw" / f"{scan['name']}_rows"
    rows_dir.mkdir(parents=True, exist_ok=True)
    (run_dir / "processed").mkdir(parents=True, exist_ok=True)

    S = int(cfg["model"]["S"])
    M = int(cfg["model"]["M"])
    R0, gamma = build_resource_vectors(cfg["resources"], M)
    gamma_per_community = is_gamma_per_community(cfg["resources"])
    ph_shape = _build_ph_shape(cfg, S)
    snapshot_run(run_dir, cfg["_config_path"], cfg, " ".join(sys.argv))
    _save_metadata(run_dir, rows_dir, cfg, scan, R0, gamma)

    rows = _resolve_rows(args, len(scan["axis1"]))

    print("=" * 72)
    print("CRmodel — simulate")
    print(f"  config     : {cfg['_config_path']}")
    if cfg.get("_overlay_path"):
        print(f"  overlay    : {cfg['_overlay_path']}")
    print(f"  run_dir    : {run_dir}")
    print(f"  scan       : {scan['name']} "
          f"({scan['axis1_name']} x {scan['axis2_name']})")
    print(f"  grid       : {len(scan['axis1'])} x {len(scan['axis2'])}")
    print(f"  ensemble   : {cfg['simulation']['n_ensemble']}")
    print(f"  cycles/save: {cfg['simulation']['n_cycles']} "
          f"/ {cfg['simulation']['n_saved_days']}")
    if gamma_per_community:
        res = cfg["resources"]
        print("  gamma      : Normal("
              f"{float(res.get('gamma_mean', 0.0))}, "
              f"{float(res['gamma_sd'])}^2) per community")
    else:
        print(f"  gamma      : {list(set(gamma))}")
    print(f"  R0         : {list(set(R0))}")
    print(f"  growth_cv  : {cfg['traits'].get('growth_cv', 0.0)}")
    print(f"  ph_shaping : {cfg.get('ph_shaping', 'all')}")
    print(f"  ph_response: {cfg.get('ph_response', cfg.get('model', {}).get('ph_response', 'triangular'))}")
    print(f"  workers    : {args.workers}")
    print(f"  rows       : {rows[:8]}{' ...' if len(rows) > 8 else ''}")
    print("=" * 72)

    if args.dry_run:
        return

    if not args.no_warmup:
        print("Warming up Numba JIT...", flush=True)
        warm_params = dict(scan["fixed_params"])
        warm_params[scan["axis1_name"]] = float(scan["axis1"][0])
        warm_params[scan["axis2_name"]] = float(scan["axis2"][0])
        if scan["axis1_name"] == "r_mean":
            warm_r_mean = float(scan["axis1"][0])
        else:
            warm_r_mean = float(scan["fixed_params"]["r_mean"])
        warm_cfg = json.loads(json.dumps(cfg))
        warm_cfg["simulation"]["n_cycles"] = 2
        warm_cfg["simulation"]["n_saved_days"] = 1
        _ = run_one(1, warm_r_mean, warm_params, warm_cfg, R0, gamma, ph_shape)
        print("Warmup done.", flush=True)

    sim = cfg["simulation"]
    n_ens = int(sim["n_ensemble"])
    n_save = int(sim["n_saved_days"])
    seed_base = int(sim.get("seed_base", 0))
    dtype = np.float32 if args.dtype == "float32" else np.float64
    ck = _compression_kwargs(args.compression)

    for row_1based in rows:
        i1 = row_1based - 1
        row_path = rows_dir / f"{scan['name']}_row_{row_1based:04d}.h5"
        if row_path.exists() and args.overwrite:
            row_path.unlink()

        t0 = time.time()
        axis1_val = float(scan["axis1"][i1])
        axis2 = scan["axis2"]
        n_axis2 = len(axis2)

        # Determine r_mean: from axis1 if it *is* r_mean, else from fixed
        if scan["axis1_name"] == "r_mean":
            r_mean = axis1_val
        else:
            r_mean = float(scan["fixed_params"]["r_mean"])

        with h5py.File(row_path, "a") as h5:
            if "completed_axis2" not in h5:
                h5.attrs["row_index"] = row_1based
                h5.attrs["axis1_name"] = scan["axis1_name"]
                h5.attrs["axis1_value"] = axis1_val
                h5.attrs["axis2_name"] = scan["axis2_name"]
                h5.attrs["created_at"] = datetime.now().isoformat(
                    timespec="seconds")
                h5.create_dataset(
                    "N_last10", shape=(n_axis2, n_ens, S, n_save),
                    dtype=dtype,
                    chunks=(1, min(64, n_ens), S, n_save), **ck)
                h5.create_dataset(
                    "B_last10", shape=(n_axis2, n_ens, S, n_save),
                    dtype=dtype,
                    chunks=(1, min(64, n_ens), S, n_save), **ck)
                h5.create_dataset(
                    "p_end_last10", shape=(n_axis2, n_ens, n_save),
                    dtype=dtype,
                    chunks=(1, min(64, n_ens), n_save), **ck)
                if gamma_per_community:
                    h5.create_dataset(
                        "gamma", shape=(n_axis2, n_ens, M),
                        dtype=dtype,
                        chunks=(1, min(64, n_ens), M), **ck)
                h5.create_dataset(
                    "seed", shape=(n_axis2, n_ens), dtype=np.int64,
                    chunks=(1, min(256, n_ens)), **ck)
                h5.create_dataset(
                    "completed_axis2", data=np.zeros(n_axis2, dtype=bool))
                h5.flush()

            completed = h5["completed_axis2"]
            if bool(np.all(completed[...])):
                print(f"[skip] row {row_1based:04d}: already complete",
                      flush=True)
                continue

            for i2, a2_val in enumerate(axis2):
                if bool(completed[i2]):
                    continue

                params = dict(scan["fixed_params"])
                params[scan["axis1_name"]] = axis1_val
                params[scan["axis2_name"]] = float(a2_val)
                seeds = np.array(
                    [seed_for(seed_base, i1, i2, ie) for ie in range(n_ens)],
                    dtype=np.int64)

                results = Parallel(
                    n_jobs=args.workers, backend="loky", verbose=0)(
                    delayed(run_one)(
                        int(s), r_mean, params, cfg, R0, gamma, ph_shape)
                    for s in seeds)

                N_px = np.zeros((n_ens, S, n_save), dtype=dtype)
                B_px = np.zeros((n_ens, S, n_save), dtype=dtype)
                p_px = np.zeros((n_ens, n_save), dtype=dtype)
                gamma_px = (
                    np.zeros((n_ens, M), dtype=dtype)
                    if gamma_per_community else None
                )
                for ie, (Ns, Bs, ps, gs) in enumerate(results):
                    N_px[ie] = Ns
                    B_px[ie] = Bs
                    p_px[ie] = ps
                    if gamma_per_community:
                        gamma_px[ie] = gs

                h5["N_last10"][i2] = N_px
                h5["B_last10"][i2] = B_px
                h5["p_end_last10"][i2] = p_px
                if gamma_per_community:
                    h5["gamma"][i2] = gamma_px
                h5["seed"][i2] = seeds
                completed[i2] = True
                h5.attrs["updated_at"] = datetime.now().isoformat(
                    timespec="seconds")
                h5.flush()

                print(f"[row {row_1based:04d}] "
                      f"{scan['axis1_name']}={axis1_val:.4g} "
                      f"{scan['axis2_name']} {i2+1}/{n_axis2}"
                      f"={float(a2_val):.4g}", flush=True)

        print(f"[done] row {row_1based:04d} ({time.time()-t0:.1f}s)",
              flush=True)

    print("All requested rows complete.")


# ── Post-process ──────────────────────────────────────────────────────

def cmd_postprocess(args):
    import h5py
    import numpy as np
    from crmodel.analysis.metrics import compute_last3_metrics, compute_last10_metrics

    cfg = load_and_merge(args)
    scan = build_scan(cfg, args.n_grid)
    run_dir = resolve_run_dir(args, cfg)

    meta = _load_meta(run_dir, scan["name"])
    axis1 = np.asarray(meta["axis1"])
    axis2 = np.asarray(meta["axis2"])
    rows_dir_str = meta.get("rows_dir", "")
    if isinstance(rows_dir_str, bytes):
        rows_dir_str = rows_dir_str.decode()
    rows_dir = Path(rows_dir_str)
    if not rows_dir.is_absolute():
        rows_dir = run_dir / "raw" / f"{scan['name']}_rows"

    n_rows, n_axis2 = len(axis1), len(axis2)

    first_shape = None
    for i in range(n_rows):
        p = rows_dir / f"{scan['name']}_row_{i+1:04d}.h5"
        if p.exists():
            with h5py.File(p, "r") as h5:
                first_shape = h5["B_last10"].shape
            break
    if first_shape is None:
        raise FileNotFoundError(f"No row files in {rows_dir}")
    _, n_ens, _, n_save = first_shape

    LAST3_FLOAT_KEYS = [
        "total", "surv_frac", "shannon", "richness",
        "community_CV", "total_CV", "bray_curtis",
    ]
    LAST3_BOOL_KEYS = ["is_collapsed", "is_oscillating"]
    LAST10_KEYS = [
        "shannon_last", "shannon_cum", "delta_shannon",
        "richness_last", "richness_cum", "delta_richness",
        "synchrony_phi", "std_participation_ratio", "eigenvalue_dominance",
    ]

    fs = (n_rows, n_axis2, n_ens)
    last3 = {}
    for k in LAST3_FLOAT_KEYS:
        for suf in ("", "_B"):
            last3[k + suf] = np.full(fs, np.nan, dtype=np.float32)
    for k in LAST3_BOOL_KEYS:
        for suf in ("", "_B"):
            last3[k + suf] = np.zeros(fs, dtype=bool)
    last3["p_end"] = np.full(fs, np.nan, dtype=np.float32)

    turnover = {}
    for k in LAST10_KEYS:
        for suf in ("", "_B"):
            turnover[k + suf] = np.full(fs, np.nan, dtype=np.float32)

    seeds = np.zeros(fs, dtype=np.int64)
    completed = np.zeros(n_rows, dtype=bool)
    completed_ax2 = np.zeros((n_rows, n_axis2), dtype=bool)

    for i in range(n_rows):
        p = rows_dir / f"{scan['name']}_row_{i+1:04d}.h5"
        if not p.exists():
            if args.strict:
                raise FileNotFoundError(f"Missing row {i+1}: {p}")
            print(f"[missing] row {i+1:04d}")
            continue

        with h5py.File(p, "r") as h5:
            N10 = h5["N_last10"][...]
            B10 = h5["B_last10"][...]
            p10 = h5["p_end_last10"][...]
            c_ax2 = (h5["completed_axis2"][...].astype(bool)
                     if "completed_axis2" in h5
                     else np.ones(n_axis2, dtype=bool))
            if "seed" in h5:
                seeds[i] = h5["seed"][...]

        completed_ax2[i] = c_ax2
        if not bool(np.any(c_ax2)):
            continue

        N10c, B10c, p10c = N10[c_ax2], B10[c_ax2], p10[c_ax2]

        mN3 = compute_last3_metrics(N10c[..., -3:])
        mB3 = compute_last3_metrics(B10c[..., -3:])
        for k in LAST3_FLOAT_KEYS:
            last3[k][i, c_ax2] = mN3[k].astype(np.float32, copy=False)
            last3[k + "_B"][i, c_ax2] = mB3[k].astype(np.float32, copy=False)
        for k in LAST3_BOOL_KEYS:
            last3[k][i, c_ax2] = mN3[k]
            last3[k + "_B"][i, c_ax2] = mB3[k]
        last3["p_end"][i, c_ax2] = p10c[..., -1].astype(np.float32)

        mN10 = compute_last10_metrics(N10c)
        mB10 = compute_last10_metrics(B10c)
        for k in LAST10_KEYS:
            turnover[k][i, c_ax2] = mN10[k].astype(np.float32, copy=False)
            turnover[k + "_B"][i, c_ax2] = mB10[k].astype(np.float32, copy=False)

        completed[i] = bool(np.all(c_ax2))
        print(f"[done] row {i+1:04d}")

    meta_arrays = {
        "seed": seeds,
        "completed_rows": completed,
        "completed_axis2": completed_ax2,
        "axis1": axis1,
        "axis2": axis2,
        "scan_name": np.array(scan["name"]),
        "axis1_name": np.array(scan["axis1_name"]),
        "axis2_name": np.array(scan["axis2_name"]),
    }
    last3.update(meta_arrays)
    turnover.update(meta_arrays)

    proc = run_dir / "processed"
    proc.mkdir(parents=True, exist_ok=True)
    p1 = proc / f"{scan['name']}_metrics_last3.npz"
    p2 = proc / f"{scan['name']}_metrics_last10.npz"
    np.savez_compressed(p1, **last3)
    np.savez_compressed(p2, **turnover)
    print(f"Saved: {p1}")
    print(f"Saved: {p2}")


# ── Status ────────────────────────────────────────────────────────────

def cmd_status(args):
    import h5py
    import numpy as np

    cfg = load_and_merge(args)
    scan = build_scan(cfg, args.n_grid)
    run_dir = resolve_run_dir(args, cfg)
    rows_dir = run_dir / "raw" / f"{scan['name']}_rows"

    n_rows = len(scan["axis1"])
    n_axis2 = len(scan["axis2"])
    done, partial, missing = [], [], []

    for i in range(n_rows):
        p = rows_dir / f"{scan['name']}_row_{i+1:04d}.h5"
        if not p.exists():
            missing.append(i + 1)
            continue
        with h5py.File(p, "r") as h5:
            if "completed_axis2" not in h5:
                missing.append(i + 1)
                continue
            c = h5["completed_axis2"][...]
            if np.all(c):
                done.append(i + 1)
            else:
                partial.append((i + 1, int(c.sum()), n_axis2))

    total = n_rows
    print(f"Run: {run_dir}")
    print(f"Scan: {scan['name']} ({n_rows} x {n_axis2})")
    print(f"  Complete : {len(done)}/{total}")
    print(f"  Partial  : {len(partial)}/{total}")
    print(f"  Missing  : {len(missing)}/{total}")
    if partial:
        for r, ndone, ntot in partial[:10]:
            print(f"    row {r:04d}: {ndone}/{ntot}")
    if missing and len(missing) <= 20:
        print(f"  Missing rows: {missing}")


# ── Helper functions ──────────────────────────────────────────────────

def _resolve_rows(args, n_rows: int) -> list[int]:
    if args.row_index_env:
        val = os.environ.get(args.row_index_env)
        if val is None:
            raise ValueError(f"Env var not set: {args.row_index_env}")
        rows = [int(val)]
    elif args.rows:
        rows = []
        for part in args.rows.split(","):
            part = part.strip()
            if "-" in part:
                a, b = part.split("-", 1)
                rows.extend(range(int(a), int(b) + 1))
            else:
                rows.append(int(part))
    else:
        start = int(args.row_start or 1)
        end = int(args.row_end or n_rows)
        rows = list(range(start, end + 1))
    bad = [r for r in rows if r < 1 or r > n_rows]
    if bad:
        raise ValueError(f"Row indices out of range 1..{n_rows}: {bad}")
    return sorted(set(rows))


def _compression_kwargs(name: str) -> dict:
    if name == "none":
        return {}
    if name == "lzf":
        return {"compression": "lzf"}
    if name == "gzip":
        return {"compression": "gzip", "compression_opts": 4}
    raise ValueError(f"Unsupported compression: {name}")


def _save_metadata(run_dir, rows_dir, cfg, scan, R0, gamma):
    import h5py
    run_dir.mkdir(parents=True, exist_ok=True)
    meta_path = run_dir / f"{scan['name']}_meta.h5"
    if meta_path.exists():
        return
    try:
        with h5py.File(meta_path, "x") as h5:
            h5.attrs["scan_name"] = scan["name"]
            h5.attrs["axis1_name"] = scan["axis1_name"]
            h5.attrs["axis2_name"] = scan["axis2_name"]
            h5.attrs["rows_dir"] = str(rows_dir)
            h5.attrs["config_json"] = json.dumps(cfg, sort_keys=True,
                                                  default=str)
            h5.create_dataset("axis1", data=scan["axis1"])
            h5.create_dataset("axis2", data=scan["axis2"])
            h5.create_dataset("R0", data=R0)
            res = cfg.get("resources", {})
            gamma_mode = str(res.get("gamma_mode", ""))
            h5.attrs["gamma_mode"] = gamma_mode
            if gamma_mode.lower() in {
                "normal_per_community",
                "normal_community",
                "quenched_normal",
            }:
                h5.attrs["gamma_mean"] = float(res.get("gamma_mean", 0.0))
                h5.attrs["gamma_sd"] = float(res["gamma_sd"])
                h5.attrs["gamma_seed_offset"] = int(
                    res.get("gamma_seed_offset", 730_000_000)
                )
            else:
                h5.create_dataset("gamma", data=gamma)
    except FileExistsError:
        return


def _load_meta(run_dir, scan_name):
    import h5py
    p = run_dir / f"{scan_name}_meta.h5"
    with h5py.File(p, "r") as h5:
        meta = {k: h5[k][...] for k in h5.keys()}
        for k, v in h5.attrs.items():
            meta[k] = v
    return meta


def _slurm_workers() -> int:
    if "SLURM_CPUS_PER_TASK" in os.environ:
        return max(1, int(os.environ["SLURM_CPUS_PER_TASK"]))
    return max(1, (os.cpu_count() or 2) - 1)


# ── CLI ───────────────────────────────────────────────────────────────

def main():
    p = argparse.ArgumentParser(
        description="CR-pH feedback model: simulate, postprocess, status.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    sub = p.add_subparsers(dest="command", required=True)

    # ── simulate ──
    s = sub.add_parser("simulate", help="Run phase-diagram simulation rows")
    s.add_argument("--config", required=True)
    s.add_argument("--overlay", default=None,
                   help="Sensitivity overlay YAML (merged on top of config)")
    s.add_argument("--run-dir", default=None)
    s.add_argument("--rows", default=None, help="e.g. 1,3,10-20")
    s.add_argument("--row-start", type=int, default=None)
    s.add_argument("--row-end", type=int, default=None)
    s.add_argument("--row-index-env", default=None,
                   help="Read row from env var (e.g. SLURM_ARRAY_TASK_ID)")
    s.add_argument("--workers", type=int, default=_slurm_workers())
    s.add_argument("--n-grid", type=int, default=None)
    s.add_argument("--n-ensemble", type=int, default=None)
    s.add_argument("--n-cycles", type=int, default=None)
    s.add_argument("--n-saved-days", type=int, default=None)
    s.add_argument("--seed-base", type=int, default=None)
    s.add_argument("--growth-cv", type=float, default=None)
    s.add_argument("--resource-degree", type=int, default=None)
    s.add_argument("--dtype", choices=["float32", "float64"], default="float32")
    s.add_argument("--compression", choices=["lzf", "gzip", "none"],
                   default="lzf")
    s.add_argument("--overwrite", action="store_true")
    s.add_argument("--dry-run", action="store_true")
    s.add_argument("--no-warmup", action="store_true")
    s.set_defaults(func=cmd_simulate)

    # ── postprocess ──
    pp = sub.add_parser("postprocess", help="Aggregate row files into NPZ")
    pp.add_argument("--config", required=True)
    pp.add_argument("--overlay", default=None)
    pp.add_argument("--run-dir", default=None)
    pp.add_argument("--n-grid", type=int, default=None)
    pp.add_argument("--strict", action="store_true")
    # CLI defaults so load_and_merge works uniformly
    for k in ("n_ensemble", "n_cycles", "n_saved_days", "seed_base",
              "growth_cv", "resource_degree"):
        pp.add_argument(f"--{k.replace('_','-')}", default=None, type=float)
    pp.set_defaults(func=cmd_postprocess)

    # ── status ──
    st = sub.add_parser("status", help="Check run completion")
    st.add_argument("--config", required=True)
    st.add_argument("--overlay", default=None)
    st.add_argument("--run-dir", default=None)
    st.add_argument("--n-grid", type=int, default=None)
    for k in ("n_ensemble", "n_cycles", "n_saved_days", "seed_base",
              "growth_cv", "resource_degree"):
        st.add_argument(f"--{k.replace('_','-')}", default=None, type=float)
    st.set_defaults(func=cmd_status)

    args = p.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
