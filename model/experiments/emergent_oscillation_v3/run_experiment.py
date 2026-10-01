"""
Emergent oscillation experiment — formal version v2.

Protocol
--------
1. Build a SHARED single-species library (N=1000 seeds) at r=1.
   r_i is FIXED at r_mean (no growth_cv disorder) to isolate
   pH-mediated oscillation mechanism.
   For each seed, run S=1 at every target h to classify as
   stable / oscillating / collapsed.  Record all quenched traits.

2. At each h, assemble S=12 communities (matching legacy) from:
   (A) Stable-only species  → test for emergent oscillation
   (B) Oscillating-only species (where available) → test for emergent stability
   200 random draws per condition.

3. Compute per-community metrics:
   Shannon, community_CV, synchrony phi, SD participation ratio.

4. Select showcase trajectories:
   - Stable singles → oscillating community (most chaotic)
   - Oscillating singles → stable community
   Run full time series, plot ALL singles + community figures.

Parameters match legacy_data_Final:
  S_comm=12, M=24, gamma=[-8,0,+10], R0=1.0, r_i=1.0 (fixed),
  f=1/30, 30 cycles, tau=24, Km=0.3, p_env=7, p_tol=2, D=1e-6
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src"))

import numpy as np
from scipy.integrate import solve_ivp
from scipy import stats as sp_stats
from crmodel.model.core import _make_rhs_closure
from crmodel.model.traits import draw_sparse_preferences

# ── Model parameters (match legacy) ─────────────────────────────────
M = 24
Km = 0.3
p_tol = 2.0
p_env = 7.0
delta = 0.0
D = 1e-6
tau = 24.0
f = 1.0 / 30.0
n_cycles = 30
suicide_thresh = 1e-6
rtol = 1e-6
atol = 1e-9
max_step = tau * 0.02
r_mean = 1.0
growth_cv = 0.3
resource_degree = 6
p_opt_min, p_opt_max = 5.5, 8.5
CV_THRESH = 0.1
COLLAPSE_THRESH = 1e-3

R0 = np.ones(M)
gamma = np.zeros(M)
gamma[:8] = -8.0
gamma[8:12] = 0.0
gamma[12:] = 10.0

H_TARGETS = [0.0, 0.05, 0.06, 0.10, 0.15]
N_LIBRARY = 1000
N_COMM = 200
COMM_S = 12

OUT_DIR = os.path.dirname(os.path.abspath(__file__))


# ── Trait generation ─────────────────────────────────────────────────
def draw_single_species(seed):
    rng = np.random.default_rng(seed)
    v_pref = draw_sparse_preferences(rng, 1, M, resource_degree)
    p_opt = rng.uniform(p_opt_min, p_opt_max, size=1)
    r_i = r_mean
    v_max = (v_pref * r_i).astype(np.float64)
    return v_max, p_opt.astype(np.float64), float(r_i)


# ── Simulation helpers ───────────────────────────────────────────────
def run_endpoints(S, v_max, p_opt, h_stress, n_cyc=30):
    ph_shape = np.ones(S, dtype=np.float64)
    rhs = _make_rhs_closure(S, M, Km, p_tol, h_stress, delta, p_env,
                            v_max, p_opt, gamma, ph_shape)
    N = 1e-3 * np.ones(S)
    state_len = 2 * S + M + 1
    B_last10 = np.zeros((S, 10))

    for nc in range(n_cyc):
        y0 = np.empty(state_len)
        y0[:S] = N; y0[S:S+M] = R0; y0[S+M] = p_env; y0[S+M+1:] = 0.0
        sol = solve_ivp(rhs, (0, tau), y0, method="LSODA",
                        rtol=rtol, atol=atol, max_step=max_step)
        if not sol.success:
            sol = solve_ivp(rhs, (0, tau), y0, method="BDF",
                            rtol=rtol, atol=atol, max_step=max_step)
        y_end = sol.y[:, -1] if sol.success else np.zeros(state_len)
        N_end = y_end[:S].copy()
        N_end[N_end < suicide_thresh] = 0.0
        N = f * (N_end + D); N[N < 1e-12] = 0.0
        if nc >= n_cyc - 10:
            B_last10[:, nc - (n_cyc - 10)] = y_end[S+M+1:S+M+1+S]
    return B_last10


def run_full_ts(S, v_max, p_opt, h_stress, n_cyc=30, npts=200):
    """Return full time series for all cycles."""
    ph_shape = np.ones(S, dtype=np.float64)
    rhs = _make_rhs_closure(S, M, Km, p_tol, h_stress, delta, p_env,
                            v_max, p_opt, gamma, ph_shape)
    N = 1e-3 * np.ones(S)
    state_len = 2 * S + M + 1

    B_hist = np.zeros((S, n_cyc))
    pH_end = np.zeros(n_cyc)
    last_t = None
    last_N = None
    last_B_acc = None
    last_p = None

    for nc in range(n_cyc):
        is_last = (nc == n_cyc - 1)
        n_eval = npts if is_last else 30
        t_eval = np.linspace(0, tau, n_eval)

        y0 = np.empty(state_len)
        y0[:S] = N; y0[S:S+M] = R0; y0[S+M] = p_env; y0[S+M+1:] = 0.0

        sol = solve_ivp(rhs, (0, tau), y0, method="LSODA",
                        rtol=rtol, atol=atol, max_step=max_step,
                        t_eval=t_eval)
        if not sol.success:
            sol = solve_ivp(rhs, (0, tau), y0, method="BDF",
                            rtol=rtol, atol=atol, max_step=max_step,
                            t_eval=t_eval)

        if sol.success:
            y_end = sol.y[:, -1]
            if is_last:
                last_t = sol.t
                last_N = sol.y[:S, :].copy()
                last_B_acc = sol.y[S+M+1:S+M+1+S, :].copy()
                last_p = sol.y[S+M, :].copy()
        else:
            y_end = np.zeros(state_len)

        N_end = y_end[:S].copy()
        B_end = y_end[S+M+1:S+M+1+S].copy()
        N_end[N_end < suicide_thresh] = 0.0
        N = f * (N_end + D); N[N < 1e-12] = 0.0
        B_hist[:, nc] = B_end
        pH_end[nc] = float(y_end[S+M])

    return {
        "B_hist": B_hist, "pH_end": pH_end,
        "last_t": last_t, "last_N": last_N,
        "last_B_acc": last_B_acc, "last_p": last_p,
    }


# ── Metrics ──────────────────────────────────────────────────────────
def compute_metrics(B_last10):
    """From B endpoints shape (S, 10). Returns dict."""
    S = B_last10.shape[0]
    n_days = B_last10.shape[1]
    total = B_last10.sum(axis=0)
    mean_total = total.mean()
    is_collapsed = mean_total < COLLAPSE_THRESH

    if is_collapsed:
        return {"cv": 0.0, "collapsed": True, "oscillating": False,
                "shannon": 0.0, "richness": 0, "phi": np.nan, "PR": np.nan}

    sp_std = B_last10.std(axis=1, ddof=1)
    cv = sp_std.sum() / mean_total
    is_osc = cv >= CV_THRESH

    shannon_list = []
    richness_list = []
    for d in range(n_days):
        B_d = B_last10[:, d]
        tot_d = B_d.sum()
        if tot_d > 0:
            p = B_d / tot_d
            p = p[p > 0]
            shannon_list.append(-np.sum(p * np.log(p)))
            richness_list.append(int((B_d / tot_d > 0.01).sum()))
        else:
            shannon_list.append(0.0)
            richness_list.append(0)
    shannon = float(np.mean(shannon_list))
    richness = float(np.mean(richness_list))

    var_total = np.var(total, ddof=1)
    sum_std = sp_std.sum()
    denom = sum_std ** 2
    phi = var_total / denom if denom > 0 else np.nan

    sum_var = (sp_std ** 2).sum()
    PR = denom / sum_var if sum_var > 0 else np.nan

    if not is_osc:
        phi = np.nan
        PR = np.nan

    return {"cv": cv, "collapsed": False, "oscillating": is_osc,
            "shannon": shannon, "richness": richness, "phi": phi, "PR": PR}


# ══════════════════════════════════════════════════════════════════════
#  MAIN
# ══════════════════════════════════════════════════════════════════════
def main():
    print("=" * 70)
    print("EMERGENT OSCILLATION EXPERIMENT (formal)")
    print(f"r={r_mean}, M={M}, S_comm={COMM_S}, growth_cv={growth_cv}")
    print(f"h targets: {H_TARGETS}")
    print(f"Library: {N_LIBRARY} species, Communities: {N_COMM} per condition")
    print("=" * 70)

    # ── Step 1: Build shared library ─────────────────────────────────
    print(f"\n[1] Building library of {N_LIBRARY} single species...")
    library = []
    for seed in range(1, N_LIBRARY + 1):
        v_max, p_opt, r_i = draw_single_species(seed)
        entry = {"seed": seed, "v_max": v_max, "p_opt": p_opt, "r_i": r_i,
                 "status": {}}

        for h_val in H_TARGETS:
            B3 = run_endpoints(1, v_max, p_opt, h_val)
            m = compute_metrics(B3)
            if m["collapsed"]:
                entry["status"][h_val] = "collapsed"
            elif m["oscillating"]:
                entry["status"][h_val] = "oscillating"
            else:
                entry["status"][h_val] = "stable"

        library.append(entry)
        if seed % 200 == 0:
            print(f"  {seed}/{N_LIBRARY} done")

    # Print library summary
    print(f"\nLibrary summary (single-species classification):")
    print(f"  {'h':>6s} {'stable':>7s} {'osc':>5s} {'coll':>5s}")
    for h_val in H_TARGETS:
        n_s = sum(1 for e in library if e["status"][h_val] == "stable")
        n_o = sum(1 for e in library if e["status"][h_val] == "oscillating")
        n_c = sum(1 for e in library if e["status"][h_val] == "collapsed")
        print(f"  {h_val:6.2f} {n_s:7d} {n_o:5d} {n_c:5d}")

    # Save library
    lib_data = {
        "seeds": np.array([e["seed"] for e in library]),
        "r_i": np.array([e["r_i"] for e in library]),
        "p_opt": np.array([e["p_opt"][0] for e in library]),
    }
    for h_val in H_TARGETS:
        key = f"status_h{h_val:.2f}".replace(".", "")
        lib_data[key] = np.array([e["status"][h_val] for e in library])
    np.savez_compressed(os.path.join(OUT_DIR, "library.npz"), **lib_data)

    # ── Step 2: Assemble and test communities ────────────────────────
    rng = np.random.default_rng(42)
    all_results = {}

    for h_val in H_TARGETS:
        print(f"\n{'='*60}")
        print(f"h = {h_val}")
        print(f"{'='*60}")

        stable_lib = [e for e in library if e["status"][h_val] == "stable"]
        osc_lib = [e for e in library if e["status"][h_val] == "oscillating"]

        conditions = {}

        # Condition A: stable-only
        if len(stable_lib) >= COMM_S:
            print(f"\n  [A] Stable-only assembly ({len(stable_lib)} available)...")
            cond_results = []
            for ic in range(N_COMM):
                idx = rng.choice(len(stable_lib), size=COMM_S, replace=False)
                spp = [stable_lib[k] for k in idx]
                vm = np.vstack([s["v_max"] for s in spp])
                po = np.concatenate([s["p_opt"] for s in spp])
                B3 = run_endpoints(COMM_S, vm, po, h_val)
                m = compute_metrics(B3)
                m["seeds"] = [s["seed"] for s in spp]
                m["v_max"] = vm
                m["p_opt"] = po
                cond_results.append(m)
            n_osc = sum(1 for r in cond_results if r["oscillating"])
            print(f"      Oscillating: {n_osc}/{N_COMM} ({100*n_osc/N_COMM:.1f}%)")
            conditions["stable_only"] = cond_results
        else:
            print(f"  [A] Not enough stable species ({len(stable_lib)})")

        # Condition B: oscillating-only
        if len(osc_lib) >= COMM_S:
            print(f"\n  [B] Oscillating-only assembly ({len(osc_lib)} available)...")
            cond_results = []
            for ic in range(N_COMM):
                idx = rng.choice(len(osc_lib), size=COMM_S, replace=False)
                spp = [osc_lib[k] for k in idx]
                vm = np.vstack([s["v_max"] for s in spp])
                po = np.concatenate([s["p_opt"] for s in spp])
                B3 = run_endpoints(COMM_S, vm, po, h_val)
                m = compute_metrics(B3)
                m["seeds"] = [s["seed"] for s in spp]
                m["v_max"] = vm
                m["p_opt"] = po
                cond_results.append(m)
            n_osc = sum(1 for r in cond_results if r["oscillating"])
            n_sta = sum(1 for r in cond_results
                        if not r["oscillating"] and not r["collapsed"])
            print(f"      Oscillating: {n_osc}/{N_COMM} ({100*n_osc/N_COMM:.1f}%)")
            print(f"      Stable: {n_sta}/{N_COMM} ({100*n_sta/N_COMM:.1f}%)")
            conditions["osc_only"] = cond_results
        else:
            print(f"  [B] Not enough oscillating species ({len(osc_lib)})")

        # Print metrics summary per condition
        for cond_name, results in conditions.items():
            osc_results = [r for r in results if r["oscillating"]]
            sta_results = [r for r in results
                           if not r["oscillating"] and not r["collapsed"]]
            print(f"\n  {cond_name} metrics (oscillating n={len(osc_results)}):")
            if osc_results:
                cvs = [r["cv"] for r in osc_results]
                shs = [r["shannon"] for r in osc_results]
                phis = [r["phi"] for r in osc_results if np.isfinite(r["phi"])]
                prs = [r["PR"] for r in osc_results if np.isfinite(r["PR"])]
                richs = [r["richness"] for r in osc_results]
                print(f"    CV:      {np.mean(cvs):.3f} +/- {np.std(cvs):.3f}")
                print(f"    Shannon: {np.mean(shs):.3f} +/- {np.std(shs):.3f}")
                print(f"    Rich:    {np.mean(richs):.2f} +/- {np.std(richs):.2f}")
                if phis:
                    print(f"    phi:     {np.mean(phis):.3f} +/- {np.std(phis):.3f}")
                if prs:
                    print(f"    PR:      {np.mean(prs):.2f} +/- {np.std(prs):.2f}")
            if sta_results:
                print(f"  {cond_name} metrics (stable n={len(sta_results)}):")
                shs = [r["shannon"] for r in sta_results]
                richs = [r["richness"] for r in sta_results]
                print(f"    Shannon: {np.mean(shs):.3f} +/- {np.std(shs):.3f}")
                print(f"    Rich:    {np.mean(richs):.2f} +/- {np.std(richs):.2f}")

        all_results[h_val] = conditions

    # ── Step 3: Extra sampling for emergent stability ──────────────
    print(f"\n\n{'='*70}")
    print("EXTRA SAMPLING: osc-only assembly (500 trials) for emergent stability")
    print("=" * 70)
    rng_extra = np.random.default_rng(7777)
    N_EXTRA = 500
    for h_val in [0.06, 0.10, 0.15]:
        osc_lib = [e for e in library if e["status"][h_val] == "oscillating"]
        if len(osc_lib) < COMM_S:
            continue
        print(f"\n  h={h_val}: {len(osc_lib)} osc species, assembling {N_EXTRA}...")
        extra_results = []
        for ic in range(N_EXTRA):
            idx = rng_extra.choice(len(osc_lib), size=COMM_S, replace=False)
            spp = [osc_lib[k] for k in idx]
            vm = np.vstack([s["v_max"] for s in spp])
            po = np.concatenate([s["p_opt"] for s in spp])
            B10 = run_endpoints(COMM_S, vm, po, h_val)
            m = compute_metrics(B10)
            m["seeds"] = [s["seed"] for s in spp]
            m["v_max"] = vm
            m["p_opt"] = po
            extra_results.append(m)
        n_sta = sum(1 for r in extra_results
                    if not r["oscillating"] and not r["collapsed"])
        n_osc = sum(1 for r in extra_results if r["oscillating"])
        print(f"    Stable: {n_sta}/{N_EXTRA}, Oscillating: {n_osc}/{N_EXTRA}")
        if "osc_only_extra" not in all_results.get(h_val, {}):
            if h_val not in all_results:
                all_results[h_val] = {}
            all_results[h_val]["osc_only_extra"] = extra_results

    # ── Step 4: Select showcase examples ─────────────────────────────
    print(f"\n\n{'='*70}")
    print("SELECTING SHOWCASE TRAJECTORIES")
    print("=" * 70)

    showcase = {}

    # --- Emergent oscillation: stable→osc ---
    # (a) Highest CV
    for h_pick in [0.05, 0.06]:
        conds = all_results.get(h_pick, {})
        if "stable_only" not in conds:
            continue
        osc_comms = [r for r in conds["stable_only"] if r["oscillating"]]
        if not osc_comms:
            continue
        best = max(osc_comms, key=lambda r: r["cv"])
        tag = f"emergent_highCV_h{h_pick}"
        print(f"\n  {tag}: CV={best['cv']:.3f}, Shannon={best['shannon']:.2f}, "
              f"Rich={best['richness']:.1f}")
        showcase[tag] = {
            "h": h_pick, "type": "emergent",
            "v_max": best["v_max"], "p_opt": best["p_opt"],
            "seeds": best["seeds"], "cv": best["cv"],
            "shannon": best["shannon"], "richness": best["richness"],
        }
        break

    # (b) Low-diversity oscillation (lowest Shannon among oscillating)
    for h_pick in [0.10, 0.15, 0.06]:
        conds = all_results.get(h_pick, {})
        if "stable_only" not in conds:
            continue
        osc_comms = [r for r in conds["stable_only"]
                     if r["oscillating"] and r["richness"] <= 3]
        if not osc_comms:
            osc_comms = [r for r in conds["stable_only"] if r["oscillating"]]
        if not osc_comms:
            continue
        best = min(osc_comms, key=lambda r: r["shannon"])
        tag = f"emergent_lowdiv_h{h_pick}"
        print(f"\n  {tag}: CV={best['cv']:.3f}, Shannon={best['shannon']:.2f}, "
              f"Rich={best['richness']:.1f}")
        showcase[tag] = {
            "h": h_pick, "type": "emergent_lowdiv",
            "v_max": best["v_max"], "p_opt": best["p_opt"],
            "seeds": best["seeds"], "cv": best["cv"],
            "shannon": best["shannon"], "richness": best["richness"],
        }
        break

    # (c) High-diversity high-CV at h=0.15
    conds_015 = all_results.get(0.15, {})
    if "stable_only" in conds_015:
        osc_comms = [r for r in conds_015["stable_only"]
                     if r["oscillating"] and r["richness"] >= 6]
        if osc_comms:
            best = max(osc_comms, key=lambda r: r["cv"])
            tag = "emergent_highdiv_h015"
            print(f"\n  {tag}: CV={best['cv']:.3f}, Shannon={best['shannon']:.2f}, "
                  f"Rich={best['richness']:.1f}")
            showcase[tag] = {
                "h": 0.15, "type": "emergent_highdiv",
                "v_max": best["v_max"], "p_opt": best["p_opt"],
                "seeds": best["seeds"], "cv": best["cv"],
                "shannon": best["shannon"], "richness": best["richness"],
            }

    # --- Emergent stability: osc→stable ---
    for h_pick in [0.06, 0.10, 0.15]:
        for cond_key in ["osc_only_extra", "osc_only"]:
            conds = all_results.get(h_pick, {})
            if cond_key not in conds:
                continue
            sta_comms = [r for r in conds[cond_key]
                         if not r["oscillating"] and not r["collapsed"]]
            if not sta_comms:
                continue
            best = min(sta_comms, key=lambda r: r["cv"])
            tag = f"stability_h{h_pick}"
            print(f"\n  {tag}: CV={best['cv']:.3f}, Shannon={best['shannon']:.2f}, "
                  f"Rich={best['richness']:.1f}")
            showcase[tag] = {
                "h": h_pick, "type": "emergent_stability",
                "v_max": best["v_max"], "p_opt": best["p_opt"],
                "seeds": best["seeds"], "cv": best["cv"],
                "shannon": best["shannon"], "richness": best["richness"],
            }
            break
        if any(k.startswith("stability_") for k in showcase):
            break

    # Run full time series for showcase examples
    for name, sc in showcase.items():
        h_val = sc["h"]
        print(f"\n  Running full TS: {name} ...")

        ts_comm = run_full_ts(COMM_S, sc["v_max"], sc["p_opt"], h_val)
        sc["ts_comm"] = ts_comm

        single_ts_list = []
        for i, seed in enumerate(sc["seeds"]):
            entry = library[seed - 1]
            ts_1 = run_full_ts(1, entry["v_max"], entry["p_opt"], h_val)
            ts_1["seed"] = seed
            ts_1["status"] = entry["status"][h_val]
            ts_1["p_opt"] = float(entry["p_opt"][0])
            single_ts_list.append(ts_1)
        sc["singles"] = single_ts_list

    # ── Step 4: Plot ─────────────────────────────────────────────────
    print("\n\nPlotting...")
    plot_showcase(showcase)

    # ── Step 5: Write statistics report ──────────────────────────────
    lines = []
    lines.append("=" * 70)
    lines.append("EMERGENT OSCILLATION EXPERIMENT — STATISTICS REPORT")
    lines.append(f"r_mean={r_mean}, r_i=FIXED (no growth_cv disorder)")
    lines.append(f"S_comm={COMM_S}, M={M}, f={f:.4f}, n_cycles={n_cycles}")
    lines.append(f"gamma=[-8, 0, +10], R0=1.0, p_tol={p_tol}, p_env={p_env}")
    lines.append(f"Library: {N_LIBRARY} species, Communities: {N_COMM}/condition")
    lines.append("=" * 70)

    lines.append(f"\n{'h':>6s} | {'1sp_sta':>7s} {'1sp_osc':>7s} {'1sp_col':>7s} | "
                 f"{'A:sta->osc':>10s} | {'B:osc->sta':>10s}")
    lines.append("-" * 70)
    for h_val in H_TARGETS:
        n_s = sum(1 for e in library if e["status"][h_val] == "stable")
        n_o = sum(1 for e in library if e["status"][h_val] == "oscillating")
        n_c = sum(1 for e in library if e["status"][h_val] == "collapsed")
        conds = all_results.get(h_val, {})
        a_str = ""
        if "stable_only" in conds:
            n_a = sum(1 for r in conds["stable_only"] if r["oscillating"])
            a_str = f"{n_a:>4d}/200"
        b_str = ""
        if "osc_only" in conds:
            n_b = sum(1 for r in conds["osc_only"]
                      if not r["oscillating"] and not r["collapsed"])
            b_str = f"{n_b:>4d}/200"
        lines.append(f"{h_val:6.2f} | {n_s:7d} {n_o:7d} {n_c:7d} | "
                     f"{a_str:>10s} | {b_str:>10s}")

    # Per-condition detailed metrics
    for h_val in H_TARGETS:
        conds = all_results.get(h_val, {})
        if not conds:
            continue
        lines.append(f"\n--- h = {h_val} ---")
        for cond_name, results in conds.items():
            osc_r = [r for r in results if r["oscillating"]]
            sta_r = [r for r in results
                     if not r["oscillating"] and not r["collapsed"]]
            lines.append(f"  {cond_name}:")
            if osc_r:
                lines.append(f"    Oscillating (n={len(osc_r)}):")
                lines.append(f"      CV:      {np.mean([r['cv'] for r in osc_r]):.3f} +/- {np.std([r['cv'] for r in osc_r]):.3f}")
                lines.append(f"      Shannon: {np.mean([r['shannon'] for r in osc_r]):.3f} +/- {np.std([r['shannon'] for r in osc_r]):.3f}")
                lines.append(f"      Rich:    {np.mean([r['richness'] for r in osc_r]):.2f} +/- {np.std([r['richness'] for r in osc_r]):.2f}")
                phis = [r['phi'] for r in osc_r if np.isfinite(r['phi'])]
                prs = [r['PR'] for r in osc_r if np.isfinite(r['PR'])]
                if phis:
                    lines.append(f"      phi:     {np.mean(phis):.3f} +/- {np.std(phis):.3f}")
                if prs:
                    lines.append(f"      PR:      {np.mean(prs):.2f} +/- {np.std(prs):.2f}")
            if sta_r:
                lines.append(f"    Stable (n={len(sta_r)}):")
                lines.append(f"      Shannon: {np.mean([r['shannon'] for r in sta_r]):.3f} +/- {np.std([r['shannon'] for r in sta_r]):.3f}")
                lines.append(f"      Rich:    {np.mean([r['richness'] for r in sta_r]):.2f} +/- {np.std([r['richness'] for r in sta_r]):.2f}")

    # Showcase info
    for name, sc in showcase.items():
        lines.append(f"\n{'='*50}")
        lines.append(f"SHOWCASE: {name}")
        lines.append(f"  h={sc['h']}, type={sc['type']}, CV={sc['cv']:.4f}")
        lines.append(f"  Seeds: {sc['seeds']}")
        lines.append(f"  Singles (p_opt):")
        for i, ts_1 in enumerate(sc["singles"]):
            lines.append(f"    sp{i+1:02d}: seed={ts_1['seed']}, "
                         f"status={ts_1['status']}, p_opt={ts_1['p_opt']:.2f}")

    report = "\n".join(lines)
    print(f"\n\n{report}")

    with open(os.path.join(OUT_DIR, "statistics_report.txt"), "w") as fout:
        fout.write(report)
    print(f"\nReport saved: statistics_report.txt")
    print(f"All outputs in: {OUT_DIR}")


# ══════════════════════════════════════════════════════════════════════
#  PUBLICATION PLOTTING — matches cr_model_plot.py style
# ══════════════════════════════════════════════════════════════════════
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

FIG_SIZE_MM = 50.0
FIG_SIZE_IN = FIG_SIZE_MM / 25.4
FONT_SIZE = 9
SPINE_LW = 1.4
TICK_LW = 1.1
TICK_LEN_MAJ = 4.0
TICK_LEN_MIN = 2.5
LW_FIG1 = 1.2
MS_FIG1 = 2.0
LW_FIG2 = 2.0
LW_FIG3 = 1.2
MS_FIG3 = 2.0
LW_FIG4 = 2.0
PH_COLOR = '#1f3a68'
PH_REF_COLOR = '#999999'
Y_MIN = 1e-6
Y_MAX = 1e1
PH_MIN = 3.0
PH_MAX = 11.0
PH_TICKS = [3, 7, 11]


def _style_axis(ax):
    ax.tick_params(which='major', direction='in',
                   labelsize=FONT_SIZE, width=TICK_LW, length=TICK_LEN_MAJ)
    ax.tick_params(which='minor', direction='in',
                   labelsize=FONT_SIZE, width=TICK_LW, length=TICK_LEN_MIN)
    for sp in ax.spines.values():
        sp.set_linewidth(SPINE_LW)
        sp.set_color('black')
    ax.set_facecolor('white')


def _new_canvas():
    fig, ax = plt.subplots(figsize=(FIG_SIZE_IN, 40 / 25.4))
    fig.patch.set_facecolor('white')
    fig.subplots_adjust(left=0.22, right=0.96, bottom=0.20, top=0.96)
    return fig, ax


def _get_colors(S):
    cmap = plt.get_cmap('tab10', S) if S <= 10 else plt.get_cmap('tab20', S)
    return [cmap(i) for i in range(S)]


def plot_4figs(ts, S, prefix, out_dir, color_override=None):
    """Plot Fig1-4 matching cr_model_plot.py exactly.
    color_override: single color (for single-species plots matching community palette).
    """
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['ps.fonttype'] = 42

    B_hist = ts["B_hist"]
    pH_end = ts["pH_end"]
    last_t = ts["last_t"]
    last_N = ts["last_N"]
    last_p = ts["last_p"]
    if color_override is not None:
        cols = [color_override] * S
    else:
        cols = _get_colors(S)
    days = np.arange(1, n_cycles + 1)

    def _save(fig, suffix):
        fig.savefig(os.path.join(out_dir, prefix + suffix + '.pdf'),
                    facecolor='white')
        fig.savefig(os.path.join(out_dir, prefix + suffix + '.png'),
                    facecolor='white', dpi=300)
        plt.close(fig)

    # Fig1: B_i vs cycle (log-y)
    fig1, ax1 = _new_canvas()
    ax1.axhline(1e-3, color=PH_REF_COLOR, lw=0.8, ls='--', zorder=2)
    for i in range(S):
        mask = B_hist[i] > 0
        if mask.any():
            ax1.semilogy(days[mask], B_hist[i, mask],
                         '-o', ms=MS_FIG1, color=cols[i],
                         markeredgewidth=0, alpha=0.95,
                         lw=LW_FIG1, zorder=3)
    ax1.set_yscale('log')
    ax1.set_xlim(0, n_cycles + 1)
    ax1.set_ylim(Y_MIN, Y_MAX)
    ax1.xaxis.set_major_locator(ticker.FixedLocator([0, 15, 30]))
    ax1.yaxis.set_major_locator(ticker.FixedLocator([1e-6, 1e-3, 1e0]))
    ax1.yaxis.set_minor_locator(ticker.NullLocator())
    _style_axis(ax1)
    _save(fig1, '_fig1')

    # Fig2: N_i(t) last cycle (log-y)
    if last_N is not None and last_t is not None:
        fig2, ax2 = _new_canvas()
        ax2.axhline(1e-3, color=PH_REF_COLOR, lw=0.8, ls='--', zorder=2)
        for i in range(S):
            traj = last_N[i]
            if traj.max() > 1e-12:
                ax2.semilogy(last_t, traj, lw=LW_FIG2, color=cols[i],
                             alpha=0.95, zorder=3)
        ax2.set_yscale('log')
        ax2.set_xlim(0, 24)
        ax2.set_ylim(Y_MIN, Y_MAX)
        ax2.xaxis.set_major_locator(ticker.FixedLocator([0, 12, 24]))
        ax2.yaxis.set_major_locator(ticker.FixedLocator([1e-6, 1e-3, 1e0]))
        ax2.yaxis.set_minor_locator(ticker.NullLocator())
        _style_axis(ax2)
        _save(fig2, '_fig2')

    # Fig3: pH end of cycle
    fig3, ax3 = _new_canvas()
    ax3.axhline(p_env, color=PH_REF_COLOR, lw=0.8, ls='--', zorder=2)
    ax3.plot(days, pH_end, '-o', ms=MS_FIG3, color=PH_COLOR,
             markeredgewidth=0, lw=LW_FIG3, zorder=3)
    ax3.set_xlim(0, n_cycles + 1)
    ax3.set_ylim(PH_MIN, PH_MAX)
    ax3.xaxis.set_major_locator(ticker.FixedLocator([0, 15, 30]))
    ax3.yaxis.set_major_locator(ticker.FixedLocator(PH_TICKS))
    ax3.yaxis.set_minor_locator(ticker.NullLocator())
    _style_axis(ax3)
    _save(fig3, '_fig3')

    # Fig4: pH(t) last cycle
    if last_p is not None and last_t is not None:
        fig4, ax4 = _new_canvas()
        ax4.axhline(p_env, color=PH_REF_COLOR, lw=0.8, ls='--', zorder=2)
        ax4.plot(last_t, last_p, color=PH_COLOR, lw=LW_FIG4, zorder=3)
        ax4.set_xlim(0, 24)
        ax4.set_ylim(PH_MIN, PH_MAX)
        ax4.xaxis.set_major_locator(ticker.FixedLocator([0, 12, 24]))
        ax4.yaxis.set_major_locator(ticker.FixedLocator(PH_TICKS))
        ax4.yaxis.set_minor_locator(ticker.NullLocator())
        _style_axis(ax4)
        _save(fig4, '_fig4')


def plot_showcase(showcase):
    for name, sc in showcase.items():
        h_val = sc["h"]
        S_c = COMM_S
        ts_comm = sc["ts_comm"]

        case_dir = os.path.join(OUT_DIR, name)
        os.makedirs(case_dir, exist_ok=True)

        comm_colors = _get_colors(S_c)

        plot_4figs(ts_comm, S_c, "community", case_dir)
        print(f"  Plotted: {name}/community_fig{{1-4}}")

        for i, ts_1 in enumerate(sc["singles"]):
            prefix = f"single_sp{i+1:02d}_seed{ts_1['seed']}"
            plot_4figs(ts_1, 1, prefix, case_dir,
                       color_override=comm_colors[i])
        print(f"  Plotted: {name}/single_sp01-{len(sc['singles']):02d}")


if __name__ == "__main__":
    main()
