"""
ODE right-hand side and serial-dilution driver.

Implements the Consumer-Resource model with pH-mediated cross-feeding.

State vector (length 2S + M + 1):
    y[0:S]           -> N   (species abundances)
    y[S:S+M]         -> R   (resource concentrations)
    y[S+M]           -> p   (pH / environmental variable)
    y[S+M+1:2S+M+1]  -> B   (biomass accumulators, reset each cycle)

Current production version:
    - no explicit pH-buffer damping term
"""

from __future__ import annotations

import numpy as np
from numba import njit
from scipy.integrate import solve_ivp


DEFAULT_PARAMS = dict(
    S=12, M=3, Km=0.3, p_tol=2.0, h_stress=0.1,
    delta=0.0, p_env=7.0, D=1e-6, tau=24.0, f=1.0/30.0,
)


# ----------------------------------------------------------------------
#  Numba-JIT compiled RHS
# ----------------------------------------------------------------------

@njit(cache=True, fastmath=True, boundscheck=False)
def ode_rhs(t: float, y: np.ndarray,
            S: int, M: int,
            Km: float, p_tol: float, h_stress: float,
            delta: float, p_env: float,
            v_max: np.ndarray,      # (S, M)
            p_opt: np.ndarray,      # (S,)
            gamma: np.ndarray,      # (M,)
            ph_shape: np.ndarray,
            ph_response_mode: int,  # 0 triangular, 1 top-hat window
            ph_buffer_lo: float,    # lower bound of buffer zone (0 = disabled)
            ph_buffer_hi: float,    # upper bound of buffer zone
            ph_buffer_fac: float,   # damping factor inside buffer zone
            dydt: np.ndarray) -> np.ndarray:
    """In-place RHS. dydt is preallocated and overwritten."""

    p = y[S + M]
    if p < 0.0:
        p = 0.0

    # -- per-species pH response factor phi(p) --
    phi = np.empty(S)
    dp = np.empty(S)
    for i in range(S):
        dpi = p - p_opt[i]
        dp[i] = dpi
        adp = dpi if dpi >= 0.0 else -dpi
        if ph_response_mode == 1:
            val = 1.0 if adp <= p_tol else 0.0
        else:
            val = 1.0 - adp / p_tol
            if val < 0.0:
                val = 0.0
        phi[i] = val

    # -- Monod factor per resource: R/(Km+R), with R clipped at 0 --
    monod = np.empty(M)
    for j in range(M):
        Rj = y[S + j]
        if Rj < 0.0:
            Rj = 0.0
        monod[j] = Rj / (Km + Rj)

    # -- J_real, growth, dRdt, gamma_sum --
    growth    = np.zeros(S)
    gamma_sum = np.zeros(S)
    dRdt      = np.zeros(M)

    for i in range(S):
        Ni_raw = y[i]
        Ni = Ni_raw if Ni_raw >= 0.0 else 0.0
        phi_i = phi[i]
        if phi_i == 0.0:
            continue
        gsum_i  = 0.0
        ggsum_i = 0.0
        for j in range(M):
            Jij = v_max[i, j] * monod[j] * phi_i
            gsum_i  += Jij
            ggsum_i += Jij * gamma[j]
            dRdt[j] -= Ni * Jij
        growth[i]    = gsum_i
        gamma_sum[i] = ggsum_i

    # -- pH dynamics --
    dpdt = -delta * (p - p_env)
    for i in range(S):
        Ni_raw = y[i]
        Ni = Ni_raw if Ni_raw >= 0.0 else 0.0
        dpdt += Ni * gamma_sum[i] * ph_shape[i]

    # pH buffer: dampen dp/dt inside [lo, hi]
    if ph_buffer_fac > 0.0 and p >= ph_buffer_lo and p <= ph_buffer_hi:
        dpdt = dpdt / ph_buffer_fac

    # -- dN/dt and dB/dt --
    for i in range(S):
        Ni_raw = y[i]
        Ni = Ni_raw if Ni_raw >= 0.0 else 0.0
        d_stress = h_stress * dp[i] * dp[i]
        dydt[i] = Ni * (growth[i] - d_stress)
        dydt[S + M + 1 + i] = Ni * growth[i]

    for j in range(M):
        dydt[S + j] = dRdt[j]

    dydt[S + M] = dpdt

    return dydt


def _make_rhs_closure(S, M, Km, p_tol, h_stress, delta, p_env,
                      v_max, p_opt, gamma, ph_shape,
                      ph_response_mode=0,
                      ph_buffer_lo=0.0, ph_buffer_hi=0.0, ph_buffer_fac=0.0):
    dydt = np.empty(2 * S + M + 1)

    def rhs(t, y):
        return ode_rhs(t, y, S, M, Km, p_tol, h_stress, delta, p_env,
                       v_max, p_opt, gamma, ph_shape,
                       ph_response_mode,
                       ph_buffer_lo, ph_buffer_hi, ph_buffer_fac,
                       dydt).copy()
    return rhs


def run_community(S: int, M: int,
                  Km: float, p_tol: float, h_stress: float,
                  delta: float, p_env: float, D: float,
                  tau: float, f: float,
                  R0: np.ndarray, gamma: np.ndarray,
                  v_max: np.ndarray, p_opt: np.ndarray,
                  n_cycles: int, suicide_thresh: float,
                  rtol: float, atol: float, max_step: float,
                  ph_shape: np.ndarray | None = None,
                  ph_response_mode: int = 0,
                  ph_buffer_lo: float = 0.0,
                  ph_buffer_hi: float = 0.0,
                  ph_buffer_fac: float = 0.0):
    """
    Run a serial-dilution community for `n_cycles`.

    Returns:
        N_last3  (S, 3)  鈥?last 3 cycles' end-of-cycle N (ring buffer order)
        B_last3  (S, 3)  鈥?last 3 cycles' end-of-cycle B
        p_end    float   鈥?last cycle's end-of-cycle p
    """
    if ph_shape is None:
        ph_shape = np.ones(S, dtype=np.float64)
    N = 1e-3 * np.ones(S, dtype=np.float64)
    N_last3 = np.zeros((S, 3), dtype=np.float64)
    B_last3 = np.zeros((S, 3), dtype=np.float64)
    p_end = p_env

    rhs = _make_rhs_closure(S, M, Km, p_tol, h_stress, delta, p_env,
                            v_max, p_opt, gamma, ph_shape,
                            ph_response_mode,
                            ph_buffer_lo, ph_buffer_hi, ph_buffer_fac)

    t_span    = (0.0, tau)
    state_len = 2 * S + M + 1

    for nc in range(n_cycles):
        y0 = np.empty(state_len, dtype=np.float64)
        y0[:S]          = N
        y0[S:S + M]     = R0
        y0[S + M]       = p_env
        y0[S + M + 1:]  = 0.0

        sol = solve_ivp(
            rhs, t_span, y0,
            method="LSODA",
            rtol=rtol, atol=atol, max_step=max_step,
            dense_output=False, vectorized=False,
        )

        if not sol.success:
            sol = solve_ivp(
                rhs, t_span, y0,
                method="BDF", rtol=rtol, atol=atol, max_step=max_step,
            )
            if not sol.success:
                y_end = y0.copy()
                y_end[:S] = 0.0
            else:
                y_end = sol.y[:, -1]
        else:
            y_end = sol.y[:, -1]

        N_end = y_end[:S].copy()
        p_end = float(y_end[S + M])
        B_end = y_end[S + M + 1:].copy()

        # Cycle boundary
        N_end[N_end < suicide_thresh] = 0.0
        N = f * (N_end + D)
        N[N < 1e-12] = 0.0

        idx = nc % 3
        N_last3[:, idx] = N_end
        B_last3[:, idx] = B_end

    return N_last3, B_last3, p_end


def run_community_endpoints(S: int, M: int,
                            Km: float, p_tol: float, h_stress: float,
                            delta: float, p_env: float, D: float,
                            tau: float, f: float,
                            R0: np.ndarray, gamma: np.ndarray,
                            v_max: np.ndarray, p_opt: np.ndarray,
                            n_cycles: int, n_save: int,
                            suicide_thresh: float,
                            rtol: float, atol: float, max_step: float,
                            ph_shape: np.ndarray | None = None,
                            ph_response_mode: int = 0,
                            ph_buffer_lo: float = 0.0,
                            ph_buffer_hi: float = 0.0,
                            ph_buffer_fac: float = 0.0):
    """Run serial dilution and keep the final ``n_save`` daily endpoints.

    This is the production runner for large sweeps. It returns endpoints in
    chronological order, unlike the legacy ring-buffer helper above.

    Returns
    -------
    N_save : ndarray, shape (S, n_save)
        End-of-cycle abundance endpoints.
    B_save : ndarray, shape (S, n_save)
        End-of-cycle biomass/yield endpoints.
    p_save : ndarray, shape (n_save,)
        End-of-cycle pH endpoints.
    """
    if n_save < 1:
        raise ValueError("n_save must be >= 1")
    if n_save > n_cycles:
        raise ValueError("n_save cannot exceed n_cycles")

    if ph_shape is None:
        ph_shape = np.ones(S, dtype=np.float64)
    N = 1e-3 * np.ones(S, dtype=np.float64)
    N_save = np.zeros((S, n_save), dtype=np.float64)
    B_save = np.zeros((S, n_save), dtype=np.float64)
    p_save = np.zeros(n_save, dtype=np.float64)

    rhs = _make_rhs_closure(S, M, Km, p_tol, h_stress, delta, p_env,
                            v_max, p_opt, gamma, ph_shape,
                            ph_response_mode,
                            ph_buffer_lo, ph_buffer_hi, ph_buffer_fac)

    t_span = (0.0, tau)
    state_len = 2 * S + M + 1
    first_saved_cycle = n_cycles - n_save

    for nc in range(n_cycles):
        y0 = np.empty(state_len, dtype=np.float64)
        y0[:S] = N
        y0[S:S + M] = R0
        y0[S + M] = p_env
        y0[S + M + 1:] = 0.0

        sol = solve_ivp(
            rhs, t_span, y0,
            method="LSODA",
            rtol=rtol, atol=atol, max_step=max_step,
            dense_output=False, vectorized=False,
        )

        if not sol.success:
            sol = solve_ivp(
                rhs, t_span, y0,
                method="BDF", rtol=rtol, atol=atol, max_step=max_step,
            )
            if not sol.success:
                y_end = y0.copy()
                y_end[:S] = 0.0
            else:
                y_end = sol.y[:, -1]
        else:
            y_end = sol.y[:, -1]

        N_end = y_end[:S].copy()
        B_end = y_end[S + M + 1:].copy()
        p_end = float(y_end[S + M])

        N_end[N_end < suicide_thresh] = 0.0
        N = f * (N_end + D)
        N[N < 1e-12] = 0.0

        if nc >= first_saved_cycle:
            idx = nc - first_saved_cycle
            N_save[:, idx] = N_end
            B_save[:, idx] = B_end
            p_save[idx] = p_end

    return N_save, B_save, p_save
