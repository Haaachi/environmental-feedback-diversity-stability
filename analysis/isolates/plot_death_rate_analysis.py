from pathlib import Path
import sys
import os

ISOLATE_ROOT = Path(os.environ.get("ISOLATE_WORKSPACE", Path(__file__).resolve().parent))

"""
Death rate analysis: sliding-window max death rate on CFU exponent (d).
Method: for each species at each temp, scan consecutive CFU pairs
  death_rate = (d_i - d_{i+1}) / (t_{i+1} - t_i)  [orders/hour]
  max across all windows = peak instantaneous death rate.
  Species without meaningful drop (< 2 orders) get death_rate = 0.

ALL species included (death_rate=0 for no-death species).

Outputs (figures/death_rate_analysis/):
  1. death_rate_vs_pH_dev.pdf/png   — |pH-6.5| (x) vs max death rate (y), regression
  2. death_rate_by_temp.pdf/png     — jitter+mean±SEM by temperature + sig brackets
  3. death_rate_vs_pH_dev_sized     — same as 1, dot size = total drop
Exclude species10.
"""
import openpyxl
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

XLSX = str(ISOLATE_ROOT / 'Summary.xlsx')
OUTDIR = str(ISOLATE_ROOT / 'figures/death_rate_analysis')
os.makedirs(OUTDIR, exist_ok=True)

wb = openpyxl.load_workbook(XLSX, data_only=True)
temps = [20, 30, 40]
EXCLUDE = {'species10'}
CYCLE_MAX = 24
MIN_DROP_TOTAL = 1   # threshold for "real" death (>= 1 order of magnitude)
PH_CENTER = 6.5      # initial medium pH

# ── Style ─────────────────────────────────────────────────────
FONT_FAMILY = 'Arial'
FONT_AX     = 8
FONT_LABEL  = 8
FONT_ANNOT  = 6
BORDER_SIZE = 0.6
TICK_SIZE   = 0.20
TICK_LEN    = 1.0

COL_TEMP = {20: '#85C1E9', 30: '#2E6DA4', 40: '#C0392B'}
DOT_SIZE_SC = 25
DOT_ALPHA   = 0.75
JIT_SIZE    = 15
JIT_ALPHA   = 0.50
JIT_WIDTH   = 0.08
DOT_SIZE_M  = 4.0
ERR_LW      = 0.45
ERR_W       = 0.10
MEAN_STROKE = 0.7

plt.rcParams.update({
    'font.family': FONT_FAMILY,
    'font.size': FONT_AX,
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
    'axes.facecolor': 'white',
    'axes.edgecolor': 'black',
    'axes.linewidth': BORDER_SIZE,
    'axes.grid': False,
    'xtick.direction': 'in',
    'ytick.direction': 'in',
    'xtick.major.size': TICK_LEN,
    'ytick.major.size': TICK_LEN,
    'xtick.major.width': TICK_SIZE,
    'ytick.major.width': TICK_SIZE,
    'xtick.top': True,
    'ytick.right': True,
    'xtick.labeltop': False,
    'ytick.labelright': False,
    'figure.facecolor': 'white',
})

# ── Read data ─────────────────────────────────────────────────
od_times = []
for cell in wb['20OD'][1][2:]:
    if cell.value is not None and str(cell.value).strip() != '':
        try:
            h = float(cell.value)
            if h > CYCLE_MAX: break
            od_times.append(h)
        except: break

def get_cfu_times_cols(sname):
    ws = wb[sname]
    times, cols = [], []
    for idx, cell in enumerate(ws[1]):
        if cell.value is not None and str(cell.value).strip() != '':
            try:
                h = float(cell.value)
                if 0 <= h <= CYCLE_MAX:
                    times.append(h)
                    cols.append(idx + 1)
            except: pass
    return times, cols

def read_ph(sname, sid):
    ws = wb[sname]
    ph_vals = []
    for row in ws.iter_rows(min_row=2, max_row=ws.max_row, values_only=True):
        if row[0] is None: break
        if str(row[0]).strip() == sid:
            for v in row[2:2 + len(od_times)]:
                try: ph_vals.append(float(v))
                except: ph_vals.append(np.nan)
            break
    return dict(zip(od_times, ph_vals))

def read_cfu(sname, sid):
    """Return {time: log10(CFU/mL)} for sliding-window death rate calculation.
    New data format: actual colony count per 5uL spot. CFU/mL = count * 200."""
    cfu_times, cfu_cols = get_cfu_times_cols(sname)
    ws = wb[sname]
    d = {}
    for row in ws.iter_rows(min_row=2, max_row=ws.max_row):
        if row[0].value is None: break
        if str(row[0].value).strip() == sid:
            for t, col in zip(cfu_times, cfu_cols):
                v = ws.cell(row=row[0].row, column=col).value
                if v is not None and str(v).strip() != '':
                    try:
                        count = float(v)
                        if count > 0:
                            d[t] = np.log10(count * 200.0)
                    except: pass
            break
    return d

# Species list
sids = []
ws = wb['20OD']
for row in ws.iter_rows(min_row=2, max_row=ws.max_row, values_only=True):
    if row[0] is None: break
    sid = str(row[0]).strip()
    if sid not in EXCLUDE:
        sids.append(sid)

# ── Compute sliding-window max death rate (ALL species) ──────
records = []

for temp in temps:
    for sid in sids:
        cfu = read_cfu(f'{temp}CFU', sid)
        ph_lookup = read_ph(f'{temp}pH', sid)
        if not cfu:
            continue

        pts = sorted(cfu.items())  # [(t, d), ...]

        # --- Find peak and total drop ---
        peak_t, peak_d = max(pts, key=lambda x: x[1])
        pts_after_peak = [(t, d) for t, d in pts if t >= peak_t]
        min_d_after = min(d for t, d in pts_after_peak) if pts_after_peak else peak_d
        total_drop = peak_d - min_d_after

        # --- Sliding window: max instantaneous death rate ---
        best_rate = 0
        best_t_start = None
        best_t_end = None
        best_drop_seg = 0

        for i in range(len(pts) - 1):
            t1, d1 = pts[i]
            t2, d2 = pts[i + 1]
            dt = t2 - t1
            if dt <= 0:
                continue
            drop_seg = d1 - d2  # positive = death
            if drop_seg > 0:
                rate = drop_seg / dt
                if rate > best_rate:
                    best_rate = rate
                    best_t_start = t1
                    best_t_end = t2
                    best_drop_seg = drop_seg

        # Determine if this is a "real" death event
        has_death = (total_drop >= MIN_DROP_TOTAL) and (best_rate > 0)

        # pH for x-axis: use onset pH if death, else pH at 24h
        if has_death and best_t_start is not None:
            ph_ref = ph_lookup.get(best_t_start, np.nan)
            if np.isnan(ph_ref):
                for t_ph in sorted(ph_lookup.keys(),
                                   key=lambda x: abs(x - best_t_start)):
                    if not np.isnan(ph_lookup[t_ph]):
                        ph_ref = ph_lookup[t_ph]
                        break
        else:
            # No death: use pH at 24h
            ph_ref = ph_lookup.get(24.0, np.nan)
            if np.isnan(ph_ref):
                # fallback: last available pH
                for t_ph in sorted(ph_lookup.keys(), reverse=True):
                    if not np.isnan(ph_lookup[t_ph]):
                        ph_ref = ph_lookup[t_ph]
                        break

        records.append({
            'temp': temp, 'sid': sid,
            'peak_t': peak_t, 'peak_d': peak_d,
            'total_drop': total_drop,
            'max_death_rate': best_rate if has_death else 0.0,
            'has_death': has_death,
            'seg_t_start': best_t_start,
            'seg_t_end': best_t_end,
            'seg_drop': best_drop_seg if has_death else 0,
            'ph_ref': ph_ref,
            'ph_dev': abs(ph_ref - PH_CENTER) if not np.isnan(ph_ref) else np.nan,
        })

# ── Print summary ────────────────────────────────────────────
n_death = sum(1 for r in records if r['has_death'])
n_total = len(records)
print(f'Total records: {n_total} ({n_death} with death, '
      f'{n_total - n_death} without)')
print(f'Method: sliding window max death rate (consecutive 2h pairs)')
print(f'Death threshold: total peak-to-trough drop >= {MIN_DROP_TOTAL} orders')
print()
for t in temps:
    recs_d = [r for r in records if r['temp'] == t and r['has_death']]
    recs_0 = [r for r in records if r['temp'] == t and not r['has_death']]
    print(f'  {t}C: {len(recs_d)} death + {len(recs_0)} no-death = '
          f'{len(recs_d) + len(recs_0)} total')
    for r in recs_d:
        print(f'    {r["sid"]}: rate={r["max_death_rate"]:.2f} d/h '
              f'({r["seg_t_start"]:.0f}-{r["seg_t_end"]:.0f}h, '
              f'drop {r["seg_drop"]:.0f}), '
              f'pH={r["ph_ref"]:.2f}, |pH-6.5|={r["ph_dev"]:.2f}')

# ── Helpers ───────────────────────────────────────────────────
def style_ax(ax):
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color('black')
        spine.set_linewidth(BORDER_SIZE)
    ax.tick_params(axis='both', which='both', direction='in',
                   top=True, right=True,
                   labeltop=False, labelright=False,
                   length=TICK_LEN, width=TICK_SIZE, pad=3)

def add_sig_bracket(ax, label, x1, x2, y_bar, tick_h):
    ax.plot([x1, x1, x2, x2],
            [y_bar - tick_h, y_bar, y_bar, y_bar - tick_h],
            color='black', linewidth=0.3, clip_on=False)
    ax.text((x1 + x2) / 2, y_bar + tick_h * 0.5, label,
            ha='center', va='bottom', fontsize=FONT_AX,
            fontfamily=FONT_FAMILY)

def sig_label(p):
    if p > 0.05: return 'ns'
    elif p > 0.01: return '*'
    elif p > 0.001: return '**'
    else: return '***'

def format_p(p):
    if p < 0.001: return 'p < 0.001'
    elif p < 0.01: return f'p = {p:.3f}'
    else: return f'p = {p:.2f}'

# ══════════════════════════════════════════════════════════════
# Y-axis range
# ══════════════════════════════════════════════════════════════
all_rates = [r['max_death_rate'] for r in records]
y_max_data = max(all_rates) if all_rates else 1.0
y_max = np.ceil(y_max_data * 2) / 2
if y_max < y_max_data * 1.1:
    y_max += 0.5
y_ticks = np.arange(0, y_max + 0.01, 0.5)

# ══════════════════════════════════════════════════════════════
# Fig 1: Max death rate vs |pH - 6.5|  (ALL species)
# ══════════════════════════════════════════════════════════════
print('\n=== Fig 1: max death rate vs |pH - 6.5| (all species) ===')
fig, ax = plt.subplots(figsize=(55 / 25.4, 50 / 25.4))

for t in temps:
    recs = [r for r in records if r['temp'] == t and not np.isnan(r['ph_dev'])]
    x = [r['ph_dev'] for r in recs]
    y = [r['max_death_rate'] for r in recs]
    ax.scatter(x, y, s=DOT_SIZE_SC, color=COL_TEMP[t],
               alpha=DOT_ALPHA, edgecolors='white', linewidths=0.3,
               zorder=2 + temps.index(t),
               label=f'{t} °C')

# Regression line (all temps pooled, all species)
valid = [r for r in records if not np.isnan(r['ph_dev'])]
all_dev = [r['ph_dev'] for r in valid]
all_dr  = [r['max_death_rate'] for r in valid]
if len(all_dev) >= 3:
    slope, intercept, r_val, p_val_2s, _ = stats.linregress(all_dev, all_dr)
    # One-sided: H1: slope > 0
    p_val = p_val_2s / 2 if slope > 0 else 1 - p_val_2s / 2
    rho, sp_p_2s = stats.spearmanr(all_dev, all_dr)
    sp_p = sp_p_2s / 2 if rho > 0 else 1 - sp_p_2s / 2
    x_fit = np.linspace(0, max(all_dev) * 1.05, 100)
    y_fit = slope * x_fit + intercept
    ax.plot(x_fit, np.clip(y_fit, 0, None), color='grey',
            linewidth=0.6, linestyle='--', alpha=0.7, zorder=1)
    # Annotation (one-sided p)
    ann = (f'$R^2$ = {r_val**2:.2f}, {format_p(p_val)}\n'
           f'$\\rho$ = {rho:.2f}, {format_p(sp_p)}')
    ax.text(0.97, 0.97, ann, transform=ax.transAxes,
            fontsize=FONT_ANNOT, fontfamily=FONT_FAMILY,
            ha='right', va='top',
            bbox=dict(facecolor='white', edgecolor='none', alpha=0.8, pad=1.5))

ax.set_xlim(0, 4.5)
ax.set_xticks([0, 1, 2, 3, 4])
ax.set_ylim(-0.05, y_max)
ax.set_yticks(y_ticks)
ax.set_xlabel('|pH – 6.5|', fontsize=FONT_LABEL, fontfamily=FONT_FAMILY)
ax.set_ylabel('Max death rate (orders h$^{-1}$)',
              fontsize=FONT_LABEL, fontfamily=FONT_FAMILY)

# Legend
leg = ax.legend(fontsize=FONT_ANNOT, frameon=True, framealpha=0.9,
                edgecolor='none', loc='upper left',
                handletextpad=0.3, borderpad=0.3,
                labelspacing=0.2, markerscale=0.8)
leg.get_frame().set_linewidth(0)

style_ax(ax)
plt.tight_layout(pad=0.4)
for ext in ['pdf', 'png']:
    fig.savefig(os.path.join(OUTDIR, f'death_rate_vs_pH_dev.{ext}'),
                dpi=300, bbox_inches='tight', facecolor='white')
plt.close(fig)
print('Saved: death_rate_vs_pH_dev')

# ══════════════════════════════════════════════════════════════
# Fig 2: Max death rate by temperature (jitter+mean±SEM, ALL)
# ══════════════════════════════════════════════════════════════
print('\n=== Fig 2: max death rate by temperature (all species) ===')
fig, ax = plt.subplots(figsize=(40 / 25.4, 50 / 25.4))

rng = np.random.RandomState(42)

vals_per_temp = {}
for k, t in enumerate(temps):
    vals = [r['max_death_rate'] for r in records if r['temp'] == t]
    vals_per_temp[t] = vals
    x_pos = k + 1
    if vals:
        jx = x_pos + rng.uniform(-JIT_WIDTH, JIT_WIDTH, len(vals))
        ax.scatter(jx, vals, s=JIT_SIZE, color=COL_TEMP[t],
                   alpha=JIT_ALPHA, edgecolors='none', zorder=2)

for k, t in enumerate(temps):
    x_pos = k + 1
    vals = np.array(vals_per_temp[t])
    if len(vals) > 0:
        mu = np.mean(vals)
        sem = np.std(vals, ddof=1) / np.sqrt(len(vals)) if len(vals) > 1 else 0
        ax.errorbar(x_pos, mu, yerr=sem, fmt='none',
                    ecolor=COL_TEMP[t], elinewidth=ERR_LW,
                    capsize=ERR_W * 25, capthick=ERR_LW, zorder=3)
        ax.plot(x_pos, mu, 'o', markersize=DOT_SIZE_M,
                markerfacecolor='white', markeredgecolor=COL_TEMP[t],
                markeredgewidth=MEAN_STROKE, zorder=4)

# Significance tests (Wilcoxon rank-sum since many zeros, unpaired ok)
# Use Mann-Whitney U for independent samples (different species may differ)
# But same species → paired Wilcoxon
Y_MAX_FIG = y_max + 0.3  # extra room for brackets

sig_pairs = [(20, 30, 1, 2), (30, 40, 2, 3)]
for i, (t1, t2, x1, x2) in enumerate(sig_pairs):
    # Paired: same species across temperatures
    common = sorted(set(r['sid'] for r in records if r['temp'] == t1) &
                    set(r['sid'] for r in records if r['temp'] == t2))
    v1_dict = {r['sid']: r['max_death_rate'] for r in records if r['temp'] == t1}
    v2_dict = {r['sid']: r['max_death_rate'] for r in records if r['temp'] == t2}
    v1 = np.array([v1_dict[s] for s in common])
    v2 = np.array([v2_dict[s] for s in common])
    diff = v2 - v1

    # Many zeros → use Wilcoxon signed-rank on non-zero diffs
    nonzero = diff[diff != 0]
    if len(nonzero) >= 3:
        stat_res = stats.wilcoxon(v2[diff != 0], v1[diff != 0],
                                  alternative='greater')
        p = stat_res.pvalue
        method = 'Wilcoxon signed-rank (one-sided)'
    elif len(nonzero) > 0:
        # Too few non-zero: sign test
        n_pos = np.sum(nonzero > 0)
        n_neg = np.sum(nonzero < 0)
        p = stats.binom_test(n_pos, n_pos + n_neg, 0.5)
        method = 'Sign test'
    else:
        p = 1.0
        method = 'No difference'

    label = sig_label(p)
    print(f'  {t1}C vs {t2}C: n_paired={len(common)}, '
          f'n_nonzero_diff={len(nonzero)}, '
          f'{method}, p={p:.4f} -> {label}')

    y_bar = Y_MAX_FIG * (0.82 + i * 0.08)
    tick_h = Y_MAX_FIG * 0.02
    add_sig_bracket(ax, label, x1, x2, y_bar, tick_h)

ax.set_xlim(0.4, 3.6)
ax.set_xticks([1, 2, 3])
ax.set_xticklabels(['20 °C', '30 °C', '40 °C'], fontsize=FONT_AX)
ax.set_ylim(-0.05, Y_MAX_FIG)
ax.set_yticks(y_ticks)
ax.set_xlabel('Temperature', fontsize=FONT_LABEL, fontfamily=FONT_FAMILY)
ax.set_ylabel('Max death rate (orders h$^{-1}$)',
              fontsize=FONT_LABEL, fontfamily=FONT_FAMILY)

# n= annotation under each group
for k, t in enumerate(temps):
    n = len(vals_per_temp[t])
    ax.text(k + 1, -0.15, f'n={n}', ha='center', va='top',
            fontsize=FONT_ANNOT, fontfamily=FONT_FAMILY, color='grey')

style_ax(ax)
plt.tight_layout(pad=0.4)
for ext in ['pdf', 'png']:
    fig.savefig(os.path.join(OUTDIR, f'death_rate_by_temp.{ext}'),
                dpi=300, bbox_inches='tight', facecolor='white')
plt.close(fig)
print('Saved: death_rate_by_temp')

# ══════════════════════════════════════════════════════════════
# Fig 3: Max death rate vs |pH-6.5|, size = total drop (ALL)
# ══════════════════════════════════════════════════════════════
print('\n=== Fig 3: max death rate vs |pH-6.5|, size = total drop ===')
fig, ax = plt.subplots(figsize=(55 / 25.4, 50 / 25.4))

for t in temps:
    recs = [r for r in records if r['temp'] == t and not np.isnan(r['ph_dev'])]
    x = [r['ph_dev'] for r in recs]
    y = [r['max_death_rate'] for r in recs]
    # Size: min 8 for no-drop, scale with total_drop
    sizes = [max(8, r['total_drop'] * 12) for r in recs]
    ax.scatter(x, y, s=sizes, color=COL_TEMP[t],
               alpha=DOT_ALPHA, edgecolors='white', linewidths=0.3,
               zorder=2 + temps.index(t),
               label=f'{t} °C')

# Size legend
for drop_val in [2, 4, 6]:
    ax.scatter([], [], s=drop_val * 12, color='grey', alpha=0.5,
               edgecolors='white', linewidths=0.3,
               label=f'drop = {drop_val}')

ax.set_xlim(0, 4.5)
ax.set_xticks([0, 1, 2, 3, 4])
ax.set_ylim(-0.05, y_max)
ax.set_yticks(y_ticks)
ax.set_xlabel('|pH – 6.5|', fontsize=FONT_LABEL, fontfamily=FONT_FAMILY)
ax.set_ylabel('Max death rate (orders h$^{-1}$)',
              fontsize=FONT_LABEL, fontfamily=FONT_FAMILY)

leg = ax.legend(fontsize=FONT_ANNOT, frameon=True, framealpha=0.9,
                edgecolor='none', loc='upper left',
                handletextpad=0.3, borderpad=0.3,
                labelspacing=0.2, markerscale=0.8)
leg.get_frame().set_linewidth(0)

style_ax(ax)
plt.tight_layout(pad=0.4)
for ext in ['pdf', 'png']:
    fig.savefig(os.path.join(OUTDIR, f'death_rate_vs_pH_dev_sized.{ext}'),
                dpi=300, bbox_inches='tight', facecolor='white')
plt.close(fig)
print('Saved: death_rate_vs_pH_dev_sized')

# ══════════════════════════════════════════════════════════════
# Summary stats
# ══════════════════════════════════════════════════════════════
print('\n=== Summary ===')
print(f'Method: sliding window max death rate (consecutive 2h pairs)')
print(f'Death threshold: total peak-to-trough drop >= {MIN_DROP_TOTAL} orders')
print(f'Total species*temp: {len(records)} '
      f'({n_death} with death, {len(records) - n_death} without)')
for t in temps:
    vals = [r['max_death_rate'] for r in records if r['temp'] == t]
    vals_d = [v for v in vals if v > 0]
    if vals:
        print(f'  {t}C: n={len(vals)}, n_death={len(vals_d)}, '
              f'mean={np.mean(vals):.3f}+-{np.std(vals,ddof=1)/np.sqrt(len(vals)):.3f} d/h')

# ══════════════════════════════════════════════════════════════
# Multiple regression (OLS, one-sided)
#   death_rate ~ |pH-6.5| + T
#   H1: beta_pH > 0, beta_T > 0
# ══════════════════════════════════════════════════════════════
from scipy.stats import f as f_dist, t as t_dist

valid_all = [r for r in records if not np.isnan(r['ph_dev'])]
n_all = len(valid_all)
Y_rate = np.array([r['max_death_rate'] for r in valid_all])
X_ph   = np.array([r['ph_dev'] for r in valid_all])
X_T    = np.array([r['temp'] for r in valid_all], dtype=float)

print('\n' + '=' * 60)
print('Multiple Regression: death_rate ~ |pH-6.5| + T  (one-sided)')
print('=' * 60)

X_all = np.column_stack([np.ones(n_all), X_ph, X_T])
beta = np.linalg.lstsq(X_all, Y_rate, rcond=None)[0]
Y_pred = X_all @ beta
SS_res = np.sum((Y_rate - Y_pred)**2)
SS_tot = np.sum((Y_rate - np.mean(Y_rate))**2)
R2 = 1 - SS_res / SS_tot
k = 2
R2_adj = 1 - (1 - R2) * (n_all - 1) / (n_all - k - 1)

# F-test
MS_reg = (SS_tot - SS_res) / k
MS_res = SS_res / (n_all - k - 1)
F_stat = MS_reg / MS_res
p_model = 1 - f_dist.cdf(F_stat, k, n_all - k - 1)

# Coefficient SE & one-sided p
sigma2 = SS_res / (n_all - k - 1)
cov_beta = sigma2 * np.linalg.inv(X_all.T @ X_all)
se_beta = np.sqrt(np.diag(cov_beta))
t_stats = beta / se_beta
# One-sided: H1: beta > 0 → p = 1 - CDF(t)
p_one = 1 - t_dist.cdf(t_stats, n_all - k - 1)

print(f'  n = {n_all}')
print(f'  R2 = {R2:.3f}, R2_adj = {R2_adj:.3f}')
print(f'  F({k}, {n_all-k-1}) = {F_stat:.2f}, {format_p(p_model)}')
print(f'  Coefficients (one-sided H1: coef > 0):')
print(f'    Intercept:  {beta[0]:+.4f} (SE={se_beta[0]:.4f})')
print(f'    |pH-6.5|:   {beta[1]:+.4f} (SE={se_beta[1]:.4f}, '
      f'p_1sided={p_one[1]:.4f})')
print(f'    Temp (C):    {beta[2]:+.5f} (SE={se_beta[2]:.5f}, '
      f'p_1sided={p_one[2]:.4f})')

# Also Spearman for |pH-6.5| (one-sided)
rho_ph, sp_p_2s = stats.spearmanr(X_ph, Y_rate)
sp_p_1s = sp_p_2s / 2 if rho_ph > 0 else 1 - sp_p_2s / 2
print(f'  Spearman(|pH-6.5|, rate): rho={rho_ph:.3f}, p_1sided={sp_p_1s:.4f}')

# ══════════════════════════════════════════════════════════════
# Fig 4a: Bivariate FULL — with stats annotation
# ══════════════════════════════════════════════════════════════
print('\n=== Fig 4a: bivariate FULL (with stats) ===')

def p_str(p):
    if p < 0.001: return '< 0.001'
    elif p < 0.01: return f'= {p:.3f}'
    else: return f'= {p:.2f}'

fig, ax = plt.subplots(figsize=(55 / 25.4, 50 / 25.4))

for t in temps:
    recs = [r for r in records if r['temp'] == t and not np.isnan(r['ph_dev'])]
    x = np.array([r['ph_dev'] for r in recs])
    y = np.array([r['max_death_rate'] for r in recs])
    ax.scatter(x, y, s=DOT_SIZE_SC, color=COL_TEMP[t],
               alpha=DOT_ALPHA, edgecolors='white', linewidths=0.3,
               zorder=2 + temps.index(t), label=f'{t} °C')

    x_line = np.linspace(0, 4.2, 50)
    y_line = beta[0] + beta[1] * x_line + beta[2] * t
    ax.plot(x_line, np.clip(y_line, 0, None), color=COL_TEMP[t],
            linewidth=0.7, linestyle='--', alpha=0.7, zorder=1)

ann_lines = [
    f'$R^2$ = {R2:.2f} (adj. {R2_adj:.2f})',
    r'$\beta_{|\Delta pH|}$' + f' = {beta[1]:.3f},  p {p_str(p_one[1])}',
    r'$\beta_T$' + f' = {beta[2]:.4f},  p {p_str(p_one[2])}',
    f'(one-sided)',
]
ann = '\n'.join(ann_lines)
ax.text(0.97, 0.97, ann, transform=ax.transAxes,
        fontsize=FONT_ANNOT, fontfamily=FONT_FAMILY,
        ha='right', va='top', linespacing=1.4,
        bbox=dict(facecolor='white', edgecolor='grey', alpha=0.9,
                  pad=2, boxstyle='round,pad=0.3', linewidth=0.3))

ax.set_xlim(0, 4.5)
ax.set_xticks([0, 1, 2, 3, 4])
ax.set_ylim(-0.05, y_max)
ax.set_yticks(y_ticks)
ax.set_xlabel('|pH – 6.5|', fontsize=FONT_LABEL, fontfamily=FONT_FAMILY)
ax.set_ylabel('Max death rate (orders h$^{-1}$)',
              fontsize=FONT_LABEL, fontfamily=FONT_FAMILY)

leg = ax.legend(fontsize=FONT_ANNOT, frameon=True, framealpha=0.9,
                edgecolor='none', loc='upper left',
                handletextpad=0.3, borderpad=0.3,
                labelspacing=0.2, markerscale=0.8)
leg.get_frame().set_linewidth(0)

style_ax(ax)
plt.tight_layout(pad=0.4)
for ext in ['pdf', 'png']:
    fig.savefig(os.path.join(OUTDIR, f'death_rate_bivariate.{ext}'),
                dpi=600, bbox_inches='tight', facecolor='white')
plt.close(fig)
print('Saved: death_rate_bivariate')

# ══════════════════════════════════════════════════════════════
# Fig 4b: Bivariate CLEAN — no text, no axes labels
# ══════════════════════════════════════════════════════════════
print('\n=== Fig 4b: bivariate CLEAN ===')
fig, ax = plt.subplots(figsize=(55 / 25.4, 50 / 25.4))

for t in temps:
    recs = [r for r in records if r['temp'] == t and not np.isnan(r['ph_dev'])]
    x = np.array([r['ph_dev'] for r in recs])
    y = np.array([r['max_death_rate'] for r in recs])
    ax.scatter(x, y, s=DOT_SIZE_SC, color=COL_TEMP[t],
               alpha=DOT_ALPHA, edgecolors='white', linewidths=0.3,
               zorder=2 + temps.index(t))

    x_line = np.linspace(0, 4.2, 50)
    y_line = beta[0] + beta[1] * x_line + beta[2] * t
    ax.plot(x_line, np.clip(y_line, 0, None), color=COL_TEMP[t],
            linewidth=0.7, linestyle='--', alpha=0.7, zorder=1)

ax.set_xlim(0, 4.5)
ax.set_ylim(-0.05, y_max)

# No labels, no title, no tick labels
ax.set_xlabel('')
ax.set_ylabel('')
ax.set_title('')
ax.set_xticklabels([])
ax.set_yticklabels([])

for spine in ax.spines.values():
    spine.set_visible(True)
    spine.set_color('black')
    spine.set_linewidth(BORDER_SIZE)

ax.tick_params(axis='both', which='both', direction='in',
               top=True, right=True,
               labeltop=False, labelright=False,
               labelbottom=False, labelleft=False,
               length=TICK_LEN, width=TICK_SIZE, pad=2)

plt.tight_layout(pad=0.4)
for ext in ['pdf', 'png']:
    fig.savefig(os.path.join(OUTDIR, f'death_rate_bivariate_clean.{ext}'),
                dpi=600, bbox_inches='tight', facecolor='white')
plt.close(fig)
print('Saved: death_rate_bivariate_clean')

# ══════════════════════════════════════════════════════════════
# Fig 5: CFU(24h)/CFU_max ratio distribution by temperature
#   log10(CFU_24h / CFU_max):  0 = no decline, negative = death
# ══════════════════════════════════════════════════════════════
print('\n=== Fig 5: CFU(24h)/CFU_max ratio distribution ===')
from scipy.stats import gaussian_kde

ratio_data = {t: [] for t in temps}
for temp in temps:
    for sid in sids:
        cfu = read_cfu(f'{temp}CFU', sid)
        if not cfu:
            continue
        pts = sorted(cfu.items())
        # CFU at 24h (or closest available)
        t24_candidates = [t for t in cfu.keys() if abs(t - 24) <= 2]
        if not t24_candidates:
            continue
        t24 = min(t24_candidates, key=lambda x: abs(x - 24))
        log_24 = cfu[t24]
        # CFU_max (peak across all time points)
        log_max = max(cfu.values())
        log_ratio = log_24 - log_max  # always <= 0
        ratio_data[temp].append(log_ratio)
        if log_ratio < -1:
            print(f'  {sid}@{temp}C: log10(CFU24/CFUmax) = {log_ratio:+.2f}')

# --- Histogram: x = log10(CFU24/CFUmax), y = count, colors = temperature ---
TEMP_LABELS = {20: '20 °C', 30: '30 °C', 40: '40 °C'}

all_ratio_vals = [v for vals in ratio_data.values() for v in vals]
bin_edges = np.arange(np.floor(min(all_ratio_vals)) - 0.5, 1.5, 1.0)

fig, ax = plt.subplots(figsize=(18 / 25.4, 50 / 25.4))

all_vals = np.concatenate([np.array(ratio_data[t]) for t in temps])
ax.hist(all_vals, bins=bin_edges, orientation='horizontal',
        color='#555555', alpha=0.75,
        edgecolor='white', linewidth=0.4, zorder=2)

print(f'  All pooled: n={len(all_vals)}, mean={np.mean(all_vals):.2f}, '
      f'median={np.median(all_vals):.2f}')

# Threshold line at 10^-1 = ratio -1
ax.axhline(-1, color='black', linewidth=0.5, linestyle='--', alpha=0.6, zorder=3)

ax.set_ylabel('log$_{10}$(CFU$_{24h}$ / CFU$_{max}$)',
              fontsize=FONT_LABEL, fontfamily=FONT_FAMILY)
ax.set_xlabel('Count', fontsize=FONT_LABEL, fontfamily=FONT_FAMILY)
ax.set_xlim(left=0)

style_ax(ax)
plt.tight_layout(pad=0.4)
for ext in ['pdf', 'png']:
    fig.savefig(os.path.join(OUTDIR, f'cfu_ratio_24h_max.{ext}'),
                dpi=300, bbox_inches='tight', facecolor='white')
plt.close(fig)
print('Saved: cfu_ratio_24h_max')

# ══════════════════════════════════════════════════════════════
# Fig 6: CFU overview — all species grey, two mean lines
#   Group A: CFU24/CFUmax < -1 (death)   Group B: >= -1 (no death)
#   Background gradient like death example panels
# ══════════════════════════════════════════════════════════════
print('\n=== Fig 6: CFU overview (grey + two mean lines + gradient) ===')
import matplotlib.ticker as ticker

# Classify each (temp, sid) into death vs no-death based on ratio
RATIO_THRESHOLD = -1.0  # log10 scale

# Collect all CFU series with group labels
death_series = []   # list of (times, cfu_vals) for death group
nodeath_series = []

for temp in temps:
    cfu_sheet = f'{temp}CFU'
    ws_c = wb[cfu_sheet]
    cfu_times_fig, cfu_cols_fig = get_cfu_times_cols(cfu_sheet)

    for row_cells in ws_c.iter_rows(min_row=2, max_row=ws_c.max_row):
        if row_cells[0].value is None:
            break
        sid = str(row_cells[0].value).strip()
        if sid in EXCLUDE:
            continue

        # Read CFU/mL series
        vals = []
        log_vals = {}
        for t_pt, col_idx in zip(cfu_times_fig, cfu_cols_fig):
            cell = ws_c.cell(row=row_cells[0].row, column=col_idx)
            if cell.value is not None and str(cell.value).strip() != '':
                try:
                    count = float(cell.value)
                    cfu_ml = count * 200.0
                    vals.append(cfu_ml)
                    if cfu_ml > 0:
                        log_vals[t_pt] = np.log10(cfu_ml)
                except:
                    vals.append(np.nan)
            else:
                vals.append(np.nan)

        # Compute ratio for classification
        if log_vals:
            t24_c = [t for t in log_vals.keys() if abs(t - 24) <= 2]
            if t24_c:
                t24 = min(t24_c, key=lambda x: abs(x - 24))
                log_ratio = log_vals[t24] - max(log_vals.values())
            else:
                log_ratio = 0.0
        else:
            log_ratio = 0.0

        series_data = (cfu_times_fig, vals)
        if log_ratio < RATIO_THRESHOLD:
            death_series.append(series_data)
        else:
            nodeath_series.append(series_data)

print(f'  Death group (ratio < {RATIO_THRESHOLD}): {len(death_series)} curves')
print(f'  No-death group: {len(nodeath_series)} curves')

# Compute mean±SEM for each group (geometric mean in log space)
def group_mean_sem(series_list, time_ref):
    """Compute geometric mean ± SEM on log scale."""
    n_t = len(time_ref)
    log_matrix = []
    for times, vals in series_list:
        row = []
        for t_ref in time_ref:
            # Find matching time
            found = False
            for t_s, v_s in zip(times, vals):
                if abs(t_s - t_ref) < 0.01 and not np.isnan(v_s) and v_s > 0:
                    row.append(np.log10(v_s))
                    found = True
                    break
            if not found:
                row.append(np.nan)
        log_matrix.append(row)

    log_arr = np.array(log_matrix)
    means, sems = [], []
    for j in range(n_t):
        col = log_arr[:, j]
        valid = col[~np.isnan(col)]
        if len(valid) > 0:
            means.append(np.mean(valid))
            sems.append(np.std(valid, ddof=1) / np.sqrt(len(valid)) if len(valid) > 1 else 0)
        else:
            means.append(np.nan)
            sems.append(np.nan)
    return np.array(means), np.array(sems)

# Use union of all CFU time points
all_cfu_times = sorted(set(t for ts, _ in death_series + nodeath_series for t in ts
                           if 0 <= t <= CYCLE_MAX))

mean_death, sem_death = group_mean_sem(death_series, all_cfu_times)
mean_nodeath, sem_nodeath = group_mean_sem(nodeath_series, all_cfu_times)

# Plot
fig_w = 55 / 25.4
fig_h = 50 / 25.4
fig, ax = plt.subplots(figsize=(fig_w, fig_h))

# --- Gradient background (growth → white → death) ---
grad_res = 256
gradient = np.zeros((1, grad_res, 4))
c_growth = np.array([0.69, 0.85, 0.95, 1.0])
c_death  = np.array([0.95, 0.75, 0.75, 1.0])
c_white  = np.array([1.0, 1.0, 1.0, 1.0])
for i in range(grad_res):
    frac = i / (grad_res - 1)
    if frac < 0.5:
        t_local = frac / 0.5
        gradient[0, i] = c_growth * (1 - t_local) + c_white * t_local
    else:
        t_local = (frac - 0.5) / 0.5
        gradient[0, i] = c_white * (1 - t_local) + c_death * t_local
gradient[0, :, 3] = 0.35
ax.imshow(gradient, extent=[0, 24, 0, 1], aspect='auto',
          transform=ax.get_xaxis_transform(),
          zorder=0, interpolation='bicubic')

# --- Individual lines (all grey) ---
GREY = '#AAAAAA'
for times, vals in death_series + nodeath_series:
    valid = [(t, v) for t, v in zip(times, vals) if not np.isnan(v) and v > 0]
    if valid:
        x, y = zip(*valid)
        ax.plot(x, y, '-', color=GREY, alpha=0.25, linewidth=0.6, zorder=1)

# --- Mean lines ---
COL_DEATH = '#C0392B'
COL_NODEATH = '#2E6DA4'
MEAN_LW = 1.8
MARKER_SZ = 5

t_arr = np.array(all_cfu_times)

# No-death mean ± SEM (errorbar)
mask = ~np.isnan(mean_nodeath)
tp_nd = t_arr[mask]
m_nd = 10**mean_nodeath[mask]
yerr_lo_nd = m_nd - 10**(mean_nodeath[mask] - sem_nodeath[mask])
yerr_hi_nd = 10**(mean_nodeath[mask] + sem_nodeath[mask]) - m_nd
ax.errorbar(tp_nd, m_nd, yerr=[yerr_lo_nd, yerr_hi_nd], fmt='o-',
            color=COL_NODEATH, linewidth=MEAN_LW, markersize=MARKER_SZ,
            markerfacecolor='white', markeredgecolor=COL_NODEATH,
            markeredgewidth=0.8, capsize=2, capthick=0.6, elinewidth=0.6,
            zorder=5, label='No death')

# Death mean ± SEM (errorbar)
mask_d = ~np.isnan(mean_death)
tp_d = t_arr[mask_d]
m_d = 10**mean_death[mask_d]
yerr_lo_d = m_d - 10**(mean_death[mask_d] - sem_death[mask_d])
yerr_hi_d = 10**(mean_death[mask_d] + sem_death[mask_d]) - m_d
ax.errorbar(tp_d, m_d, yerr=[yerr_lo_d, yerr_hi_d], fmt='s-',
            color=COL_DEATH, linewidth=MEAN_LW, markersize=MARKER_SZ,
            markerfacecolor='white', markeredgecolor=COL_DEATH,
            markeredgewidth=0.8, capsize=2, capthick=0.6, elinewidth=0.6,
            zorder=5, label='Death')

ax.set_yscale('log')
ax.set_xlim(0, 24)
ax.set_ylim(1e3, 1e10)
ax.set_xticks([0, 6, 12, 18, 24])
ax.set_yticks([1e3, 1e5, 1e7, 1e9])

ax.set_xlabel('')
ax.set_ylabel('')
ax.set_title('')

leg = ax.legend(fontsize=FONT_ANNOT, frameon=True, framealpha=0.9,
                edgecolor='none', loc='lower right',
                handletextpad=0.3, borderpad=0.3,
                labelspacing=0.2)
leg.get_frame().set_linewidth(0)

for spine in ax.spines.values():
    spine.set_visible(True)
    spine.set_color('black')
    spine.set_linewidth(BORDER_SIZE)
ax.tick_params(axis='both', which='both', direction='in',
               top=True, right=True,
               labeltop=False, labelright=False,
               length=TICK_LEN, width=TICK_SIZE, pad=2)

plt.tight_layout(pad=0.4)
for ext in ['pdf', 'png']:
    fig.savefig(os.path.join(OUTDIR, f'cfu_overview_death_groups.{ext}'),
                dpi=600, bbox_inches='tight', facecolor='white')
plt.close(fig)
print('Saved: cfu_overview_death_groups')

# ══════════════════════════════════════════════════════════════
# Fig 6b: pH overview — grey lines only
# ══════════════════════════════════════════════════════════════
print('\n=== Fig 6b: pH overview (grey lines only) ===')

# Collect all pH series across all temps
all_ph_series = []
for temp in temps:
    ws_p = wb[f'{temp}pH']
    for row in ws_p.iter_rows(min_row=2, max_row=ws_p.max_row, values_only=True):
        if row[0] is None:
            break
        sid = str(row[0]).strip()
        if sid in EXCLUDE:
            continue
        ph_vals = []
        for v in row[2:2 + len(od_times)]:
            try:
                ph_vals.append(float(v))
            except:
                ph_vals.append(np.nan)
        all_ph_series.append(ph_vals)

fig, ax = plt.subplots(figsize=(55 / 25.4, 28 / 25.4))

for ph_vals in all_ph_series:
    valid = [(t, v) for t, v in zip(od_times, ph_vals) if not np.isnan(v)]
    if valid:
        x, y = zip(*valid)
        ax.plot(x, y, '-', color='#AAAAAA', alpha=0.35, linewidth=0.6, zorder=1)

ax.set_xlim(0, 24)
ax.set_ylim(2, 10)
ax.set_xticks([0, 6, 12, 18, 24])
ax.set_yticks([3, 5, 7, 9])
ax.set_xlabel('')
ax.set_ylabel('')
ax.set_title('')

for spine in ax.spines.values():
    spine.set_visible(True)
    spine.set_color('black')
    spine.set_linewidth(BORDER_SIZE)
ax.tick_params(axis='both', which='both', direction='in',
               top=True, right=True,
               labeltop=False, labelright=False,
               length=TICK_LEN, width=TICK_SIZE, pad=2)

plt.tight_layout(pad=0.4)
for ext in ['pdf', 'png']:
    fig.savefig(os.path.join(OUTDIR, f'ph_overview_grey.{ext}'),
                dpi=600, bbox_inches='tight', facecolor='white')
plt.close(fig)
print('Saved: ph_overview_grey')

# ══════════════════════════════════════════════════════════════
# Fig 7: pH 24h distribution (death vs no-death, vertical hist)
# ══════════════════════════════════════════════════════════════
print('\n=== Fig 7: pH 24h distribution (death vs no-death) ===')

# Collect pH 24h values, classified by death/no-death
ph24_death = []
ph24_nodeath = []

idx_24 = None
for i, t in enumerate(od_times):
    if abs(t - 24.0) < 0.01:
        idx_24 = i
        break

for temp in temps:
    cfu_sheet = f'{temp}CFU'
    ph_sheet = f'{temp}pH'
    ws_c = wb[cfu_sheet]
    ws_p = wb[ph_sheet]
    cfu_times_fig, cfu_cols_fig = get_cfu_times_cols(cfu_sheet)

    # Build sid → log_ratio map
    sid_ratio = {}
    for row_cells in ws_c.iter_rows(min_row=2, max_row=ws_c.max_row):
        if row_cells[0].value is None:
            break
        sid = str(row_cells[0].value).strip()
        if sid in EXCLUDE:
            continue
        log_v = {}
        for t_pt, col_idx in zip(cfu_times_fig, cfu_cols_fig):
            cell = ws_c.cell(row=row_cells[0].row, column=col_idx)
            if cell.value is not None and str(cell.value).strip() != '':
                try:
                    count = float(cell.value)
                    if count > 0:
                        log_v[t_pt] = np.log10(count * 200.0)
                except:
                    pass
        if log_v:
            t24_c = [t for t in log_v.keys() if abs(t - 24) <= 2]
            if t24_c:
                t24 = min(t24_c, key=lambda x: abs(x - 24))
                sid_ratio[sid] = log_v[t24] - max(log_v.values())
            else:
                sid_ratio[sid] = 0.0
        else:
            sid_ratio[sid] = 0.0

    # Read pH 24h
    for row in ws_p.iter_rows(min_row=2, max_row=ws_p.max_row, values_only=True):
        if row[0] is None:
            break
        sid = str(row[0]).strip()
        if sid in EXCLUDE:
            continue
        if idx_24 is not None and len(row) > 2 + idx_24:
            v = row[2 + idx_24]
            if v is not None:
                try:
                    ph_val = float(v)
                    ratio = sid_ratio.get(sid, 0.0)
                    if ratio < RATIO_THRESHOLD:
                        ph24_death.append(ph_val)
                    else:
                        ph24_nodeath.append(ph_val)
                except:
                    pass

print(f'  Death: n={len(ph24_death)}, No-death: n={len(ph24_nodeath)}')

COL_DEATH = '#C0392B'
COL_NODEATH = '#2E6DA4'

fig, ax = plt.subplots(figsize=(18 / 25.4, 28 / 25.4))
bin_ph = np.linspace(2, 10, 13)
ax.hist([ph24_nodeath, ph24_death], bins=bin_ph, orientation='horizontal',
        color=[COL_NODEATH, COL_DEATH], alpha=0.70,
        edgecolor='white', linewidth=0.3, stacked=True,
        label=['No death', 'Death'])
ax.set_ylim(2, 10)
ax.set_yticks([3, 5, 7, 9])
ax.set_xlim(left=0)
ax.set_xlabel('Count', fontsize=FONT_LABEL, fontfamily=FONT_FAMILY)
ax.set_ylabel('')
ax.set_yticklabels([])

for spine in ax.spines.values():
    spine.set_visible(True)
    spine.set_color('black')
    spine.set_linewidth(BORDER_SIZE)
ax.tick_params(axis='both', which='both', direction='in',
               top=True, right=True,
               labeltop=False, labelright=False,
               length=TICK_LEN, width=TICK_SIZE, pad=2)

plt.tight_layout(pad=0.3)
for ext in ['pdf', 'png']:
    fig.savefig(os.path.join(OUTDIR, f'ph_dist_death_groups.{ext}'),
                dpi=600, bbox_inches='tight', facecolor='white')
plt.close(fig)
print('Saved: ph_dist_death_groups')

print(f'\nDone. All in {OUTDIR}')
