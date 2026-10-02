# CR-pH feedback model specification

This document records the model equations and parameter conventions used by
this repository.

## Manuscript Main Run

```text
dataset id: rh_gnormal_sd3_R1_sparse6_tri_S12_h015_g100_e5000
config:     configs/phase_r_h_gnormal_sd3_R1_sparse6_tri_S12_h015_e5000.yaml
```

The main run scans the phase plane of mean growth rate `r_mean` and
pH-associated death strength `h_stress`.

```text
S = 12
M = 24
R0_j = 1.0
gamma_j ~ Normal(0, 3^2), independently drawn per community and quenched
preference model = sparse_random_links
resource_degree = 6
Phi(p) = triangular
r_mean in [0, 1.0], 100 values
h_stress in [0, 0.15], 100 values
f = 1/30
n_ensemble = 5000
n_cycles = 30
n_saved_days = 10
```

## State Variables

For species `i = 1..S` and resources `j = 1..M`, the within-cycle state is:

```text
N_i(t)  species abundance
R_j(t)  resource concentration
p(t)    shared pH-like environmental variable
B_i(t)  within-cycle cumulative biomass/yield readout
```

`B_i` is reset at the beginning of every passage cycle and is the main readout
used for manuscript phase-map classification and diversity analysis.

## Trait Draws

For each random community, species traits are drawn once and held fixed across
all passage cycles:

```text
p_opt,i ~ Uniform(5.5, 8.5)
r_i = r_mean * max(1 + growth_cv * Z_i, 0),  Z_i ~ Normal(0,1)
growth_cv = 0.3
```

Each species selects `resource_degree = 6` of the 24 resources at random.
Positive preference weights on the selected resources are drawn from a
Dirichlet distribution and then scaled by `r_i`.

For the manuscript main run, resource pH effects are also community-specific:

```text
gamma_j ~ Normal(0, 3^2)
```

The same `gamma` vector is used throughout all cycles of that community, but a
new vector is drawn for each independent community.

## Within-Cycle Dynamics

At the beginning of each cycle:

```text
R_j(0) = R0_j
p(0) = p_env = 7.0
B_i(0) = 0
```

The pH-gated uptake response is triangular:

```text
Phi_i(p) = max(0, 1 - |p - p_opt,i| / p_tol)
p_tol = 2.0
```

Resource uptake follows Monod saturation:

```text
J_ij = v_max,ij * R_j / (Km + R_j) * Phi_i(p)
Km = 0.3
g_i = sum_j J_ij
```

The ODE system is:

```text
dN_i/dt = N_i * (g_i - h_stress * (p - p_opt,i)^2)
dB_i/dt = N_i * g_i
dR_j/dt = - sum_i N_i * J_ij
dp/dt = sum_i N_i * sum_j J_ij * gamma_j - delta * (p - p_env)
```

The current manuscript model uses:

```text
delta = 0
no explicit pH-buffer damping
```

## Daily Passage Map

At the end of each cycle:

```text
N_i,end < suicide_thresh  =>  N_i,end = 0
N_i,next = f * (N_i,end + D)
N_i,next < 1e-12  =>  N_i,next = 0
```

with:

```text
suicide_thresh = 1e-7
D = 1e-6
f = 1/30
tau = 24
```

The next cycle starts with fresh resources and `p = p_env`.

## Numerical Settings

```text
solver = scipy.integrate.solve_ivp, LSODA
fallback solver = BDF
rtol = 1e-6
atol = 1e-9
max_step = 0.02 * tau
Numba is used to compile the ODE right-hand side.
```

## Output Contract

The production runner stores final daily endpoints from the last 10 cycles:

```text
N_last10
B_last10
p_end_last10
seed
gamma, for gamma-normal runs
```

Postprocessing writes:

```text
processed/phase_r_h_metrics_last3.npz
processed/phase_r_h_metrics_last10.npz
```

## Analysis Convention

The manuscript phase maps use the B-based readout.

Last-3 endpoint metrics:

```text
collapse: mean total B < 1e-3
fluctuating: non-collapsed and community_CV_B > 0.1
stable: non-collapsed and not fluctuating
```

Diversity metrics include Shannon diversity and richness. The primary
diversity-stability map uses a pairwise probability-difference effect size:

```text
Delta_P = P(H_fluctuating > H_stable) - P(H_fluctuating < H_stable)
```

Fluctuation participation is measured from temporal variation across species:

```text
PR_sigma = (sum_i sigma_i)^2 / sum_i sigma_i^2
```

where `sigma_i` is the temporal standard deviation of species `B_i` over the
final three passage endpoints, using unthresholded `B_i` and sample standard
deviation (`ddof=1`). The normalized manuscript readout is `PR = PR_sigma / S`.
The absolute biomass activity threshold used for diversity is not applied to PR.

## Non-Main Runs

The sd5/R0.5 run is retained only as an alternative/sensitivity condition:

```text
configs/phase_r_h_gnormal_sd5_R05_sparse6_tri.yaml
```

Historical comparison runs such as `gm5p5`, `gm10p10`, `gamma_zero`, `top-hat`
and S-scaling runs are included as optional analyses, not as the manuscript
main model result.
