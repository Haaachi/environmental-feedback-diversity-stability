from pathlib import Path
import sys
import os

ISOLATE_ROOT = Path(os.environ.get("ISOLATE_WORKSPACE", Path(__file__).resolve().parent))

"""
24h endpoint: OD, CFU, and OD/log10(CFU) ratio
Style: R theme_pub — jitter + mean±SEM + paired significance brackets
Exclude species10.
Output: figures/endpoint_24h/
"""
import openpyxl
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
from scipy import stats

XLSX = str(ISOLATE_ROOT / 'Summary.xlsx')
OUTDIR = str(ISOLATE_ROOT / 'figures/endpoint_24h')
os.makedirs(OUTDIR, exist_ok=True)

wb = openpyxl.load_workbook(XLSX, data_only=True)
temps = [20, 30, 40]
temp_labels = ['20°C', '30°C', '40°C']
EXCLUDE = {'species10'}

# ── R theme_pub parameters ────────────────────────────────────
FONT_FAMILY = 'Arial'
FONT_AX     = 8
BORDER_SIZE = 0.6
TICK_SIZE   = 0.20
TICK_LEN    = 1.0

DOT_SIZE    = 4.0
JIT_SIZE    = 12
JIT_ALPHA   = 0.45
JIT_WIDTH   = 0.08
ERR_LW      = 0.45
ERR_W       = 0.10
MEAN_STROKE = 0.7

COL_TEMP = {20: '#85C1E9', 30: '#2E6DA4', 40: '#C0392B'}

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
def find_24h_col(sheet_name):
    ws = wb[sheet_name]
    for idx, cell in enumerate(ws[1]):
        if cell.value is not None:
            try:
                if float(cell.value) == 24.0:
                    return idx + 1
            except:
                pass
    return None

def read_24h(sheet_name):
    col = find_24h_col(sheet_name)
    if col is None:
        return {}
    ws = wb[sheet_name]
    data = {}
    for row in ws.iter_rows(min_row=2, max_row=ws.max_row, values_only=False):
        if row[0].value is None:
            break
        sid = str(row[0].value).strip()
        if sid in EXCLUDE:
            continue
        cell = ws.cell(row=row[0].row, column=col)
        if cell.value is not None and str(cell.value).strip() != '':
            try:
                data[sid] = float(cell.value)
            except:
                pass
    return data

od_raw = {t: read_24h(f'{t}OD') for t in temps}
cfu_count_raw = {t: read_24h(f'{t}CFU') for t in temps}   # actual count per 5uL

# Derived data dicts
cfu_24h = {}      # CFU/mL for plotting
ratio_24h = {}    # OD / log10(CFU)

for t in temps:
    cfu_24h[t] = {}
    ratio_24h[t] = {}
    common = set(od_raw[t].keys()) & set(cfu_count_raw[t].keys())
    for sid in common:
        count = cfu_count_raw[t][sid]
        cfu_val = count * 200.0  # count per 5uL -> CFU/mL
        cfu_24h[t][sid] = cfu_val
        if cfu_val > 0:
            log_cfu = np.log10(cfu_val)
            if log_cfu > 0:
                ratio_24h[t][sid] = od_raw[t][sid] / log_cfu

# ── Paired significance ──────────────────────────────────────
def paired_sig(data_dict, t1, t2, use_log=False, one_sided_greater=False):
    """
    Paired test on matching species.
    one_sided_greater: H1: values at t2 > values at t1.
    Returns significance label.
    """
    common = sorted(set(data_dict[t1].keys()) & set(data_dict[t2].keys()))
    v1, v2 = [], []
    for sid in common:
        a, b = data_dict[t1].get(sid), data_dict[t2].get(sid)
        if a is None or b is None or np.isnan(a) or np.isnan(b):
            continue
        if use_log:
            if a <= 0 or b <= 0:
                continue
            v1.append(np.log10(a))
            v2.append(np.log10(b))
        else:
            v1.append(a)
            v2.append(b)

    v1, v2 = np.array(v1), np.array(v2)
    diff = v1 - v2
    n = len(diff)
    print(f'  Paired n={n}, mean_diff={diff.mean():.4f}')

    sw_p = stats.shapiro(diff).pvalue if n >= 3 else 0
    normal = sw_p > 0.05

    if one_sided_greater:
        # H1: t2 > t1  →  diff_check = v2 - v1
        if normal:
            res = stats.ttest_rel(v2, v1)
            p = res.pvalue / 2
            if (v2 - v1).mean() < 0:
                p = 1 - p
            method = 'Paired t (one-sided)'
        else:
            res = stats.wilcoxon(v2, v1, alternative='greater')
            p = res.pvalue
            method = 'Wilcoxon (one-sided)'
    else:
        if normal:
            res = stats.ttest_rel(v1, v2)
            p = res.pvalue
            method = 'Paired t'
        else:
            res = stats.wilcoxon(v1, v2)
            p = res.pvalue
            method = 'Wilcoxon'

    print(f'  {method}: p={p:.4f}')
    if p > 0.05:
        label = 'ns'
    elif p > 0.01:
        label = '*'
    elif p > 0.001:
        label = '**'
    else:
        label = '***'
    print(f'  → {label}')
    return label

# ── Drawing helpers ───────────────────────────────────────────
def add_sig_bracket(ax, label, x1, x2, y_bar, tick_h, font_size=FONT_AX):
    ax.plot([x1, x1, x2, x2],
            [y_bar - tick_h, y_bar, y_bar, y_bar - tick_h],
            color='black', linewidth=0.3, clip_on=False)
    ax.text((x1 + x2) / 2, y_bar + tick_h * 0.5, label,
            ha='center', va='bottom', fontsize=font_size,
            fontfamily=FONT_FAMILY)

def plot_stat(data_dict, filename_base, y_limits, y_ticks,
              use_log_transform=False, y_tick_labels=None,
              sig_config=None):
    """
    data_dict: {temp: {sid: value}}
    sig_config: list of dicts: t1, t2, x1, x2, use_log, one_sided_greater
    """
    W_mm, H_mm = 35, 45
    fig, ax = plt.subplots(figsize=(W_mm / 25.4, H_mm / 25.4))

    vals_per_temp = {}
    for t in temps:
        raw = [v for v in data_dict[t].values()
               if v is not None and not np.isnan(v)]
        if use_log_transform:
            vals_per_temp[t] = [np.log10(v) for v in raw if v > 0]
        else:
            vals_per_temp[t] = raw

    rng = np.random.RandomState(42)
    for k, t in enumerate(temps):
        x_pos = k + 1
        vals = vals_per_temp[t]
        jx = x_pos + rng.uniform(-JIT_WIDTH, JIT_WIDTH, len(vals))
        ax.scatter(jx, vals, s=JIT_SIZE, color=COL_TEMP[t],
                   alpha=JIT_ALPHA, edgecolors='none', zorder=2)

    for k, t in enumerate(temps):
        x_pos = k + 1
        vals = np.array(vals_per_temp[t])
        mu = np.mean(vals)
        sem = np.std(vals, ddof=1) / np.sqrt(len(vals)) if len(vals) > 1 else 0
        ax.errorbar(x_pos, mu, yerr=sem, fmt='none',
                    ecolor=COL_TEMP[t], elinewidth=ERR_LW,
                    capsize=ERR_W * 25, capthick=ERR_LW, zorder=3)
        ax.plot(x_pos, mu, 'o', markersize=DOT_SIZE,
                markerfacecolor='white', markeredgecolor=COL_TEMP[t],
                markeredgewidth=MEAN_STROKE, zorder=4)

    ax.set_xlim(0.6, 3.4)
    ax.set_xticks([1, 2, 3])
    ax.set_xticklabels(temp_labels, fontsize=FONT_AX)
    ax.set_ylim(y_limits)
    ax.set_yticks(y_ticks)
    if y_tick_labels is not None:
        ax.set_yticklabels(y_tick_labels, fontsize=FONT_AX)
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
                   length=TICK_LEN, width=TICK_SIZE, pad=3)

    if sig_config is not None:
        y_range = y_limits[1] - y_limits[0]
        for i, sc in enumerate(sig_config):
            print(f'  Sig: {sc["t1"]}°C vs {sc["t2"]}°C')
            label = paired_sig(
                data_dict, sc['t1'], sc['t2'],
                use_log=sc.get('use_log', False),
                one_sided_greater=sc.get('one_sided_greater', False),
            )
            y_bar = y_limits[0] + y_range * (0.82 + i * 0.08)
            tick_h = y_range * 0.025
            add_sig_bracket(ax, label, sc['x1'], sc['x2'], y_bar, tick_h)

    plt.tight_layout(pad=0.3)
    for ext in ['pdf', 'png']:
        outpath = os.path.join(OUTDIR, f'{filename_base}.{ext}')
        fig.savefig(outpath, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f'  Saved: {filename_base}.pdf/.png\n')

# ══════════════════════════════════════════════════════════════
# 1) OD at 24h — paired two-sided
# ══════════════════════════════════════════════════════════════
print('=' * 50)
print('1) OD 24h (paired, two-sided)')
print('=' * 50)
plot_stat(
    od_raw, 'endpoint_24h_OD',
    y_limits=(0, 1.8), y_ticks=[0, 0.5, 1.0, 1.5],
    sig_config=[
        dict(t1=20, t2=30, x1=1, x2=2),
        dict(t1=30, t2=40, x1=2, x2=3),
    ],
)

# ══════════════════════════════════════════════════════════════
# 2) CFU at 24h — log10, paired, one-sided tests
#    20v30: H1: 30>20;  30v40: H1: 30>40
# ══════════════════════════════════════════════════════════════
print('=' * 50)
print('2) CFU 24h (paired, log10, one-sided)')
print('=' * 50)

NOISE_THRESHOLD = 1.0  # |Δlog₁₀(CFU)| < 1 order = technical noise, treat as tied

def paired_sig_custom(data_dict, t_higher, t_lower, use_log=False,
                      noise_threshold=0):
    """One-sided test: H1: t_higher > t_lower (on log scale if use_log).
    noise_threshold: pairs with |diff| < threshold are excluded as ties
    (within technical error of plate counting)."""
    common = sorted(set(data_dict[t_higher].keys()) & set(data_dict[t_lower].keys()))
    vh, vl = [], []
    for sid in common:
        a = data_dict[t_higher].get(sid)
        b = data_dict[t_lower].get(sid)
        if a is None or b is None or np.isnan(a) or np.isnan(b):
            continue
        if use_log:
            if a <= 0 or b <= 0:
                continue
            vh.append(np.log10(a))
            vl.append(np.log10(b))
        else:
            vh.append(a)
            vl.append(b)
    vh, vl = np.array(vh), np.array(vl)
    diff = vh - vl
    n_total = len(diff)

    # Filter out pairs within technical noise
    if noise_threshold > 0:
        mask = np.abs(diff) >= noise_threshold
        n_tied = np.sum(~mask)
        vh_f, vl_f, diff_f = vh[mask], vl[mask], diff[mask]
        print(f'  Paired n={n_total}, noise threshold={noise_threshold} '
              f'→ {n_tied} tied, {len(diff_f)} informative')
    else:
        vh_f, vl_f, diff_f = vh, vl, diff
        n_tied = 0
        print(f'  Paired n={n_total}, mean_diff={diff.mean():.4f}')

    if len(diff_f) < 3:
        # Too few informative pairs — sign test on direction
        n_pos = np.sum(diff_f > 0)
        n_neg = np.sum(diff_f < 0)
        n_test = n_pos + n_neg
        if n_test > 0:
            p = stats.binom_test(n_pos, n_test, 0.5) / 2  # one-sided
            if n_pos < n_neg:
                p = 1 - p
        else:
            p = 1.0
        method = f'Sign test (one-sided, {n_pos}+/{n_neg}-)'
    else:
        sw_p = stats.shapiro(diff_f).pvalue
        normal = sw_p > 0.05
        if normal:
            res = stats.ttest_rel(vh_f, vl_f)
            p = res.pvalue / 2
            if diff_f.mean() < 0:
                p = 1 - p
            method = 'Paired t (one-sided)'
        else:
            res = stats.wilcoxon(vh_f, vl_f, alternative='greater')
            p = res.pvalue
            method = 'Wilcoxon (one-sided)'

    print(f'  {method}: p={p:.4f}')
    if p > 0.05: return 'ns'
    elif p > 0.01: return '*'
    elif p > 0.001: return '**'
    else: return '***'

fig, ax = plt.subplots(figsize=(35/25.4, 45/25.4))
rng = np.random.RandomState(42)

vals_per_temp = {}
for t in temps:
    raw = [v for v in cfu_24h[t].values() if not np.isnan(v) and v > 0]
    vals_per_temp[t] = [np.log10(v) for v in raw]

for k, t in enumerate(temps):
    x_pos = k + 1
    vals = vals_per_temp[t]
    jx = x_pos + rng.uniform(-JIT_WIDTH, JIT_WIDTH, len(vals))
    ax.scatter(jx, vals, s=JIT_SIZE, color=COL_TEMP[t],
               alpha=JIT_ALPHA, edgecolors='none', zorder=2)

for k, t in enumerate(temps):
    x_pos = k + 1
    vals = np.array(vals_per_temp[t])
    mu = np.mean(vals)
    sem = np.std(vals, ddof=1) / np.sqrt(len(vals)) if len(vals) > 1 else 0
    ax.errorbar(x_pos, mu, yerr=sem, fmt='none',
                ecolor=COL_TEMP[t], elinewidth=ERR_LW,
                capsize=ERR_W*25, capthick=ERR_LW, zorder=3)
    ax.plot(x_pos, mu, 'o', markersize=DOT_SIZE,
            markerfacecolor='white', markeredgecolor=COL_TEMP[t],
            markeredgewidth=MEAN_STROKE, zorder=4)

ax.set_xlim(0.6, 3.4)
ax.set_xticks([1, 2, 3])
ax.set_xticklabels(temp_labels, fontsize=FONT_AX)
ax.set_ylim(3, 10)
ax.set_yticks([3, 5, 7, 9])
ax.set_yticklabels([r'$10^{3}$', r'$10^{5}$', r'$10^{7}$', r'$10^{9}$'], fontsize=FONT_AX)
ax.set_xlabel(''); ax.set_ylabel(''); ax.set_title('')
for spine in ax.spines.values():
    spine.set_visible(True); spine.set_color('black'); spine.set_linewidth(BORDER_SIZE)
ax.tick_params(axis='both', which='both', direction='in',
               top=True, right=True, labeltop=False, labelright=False,
               length=TICK_LEN, width=TICK_SIZE, pad=3)

y_range = 7  # 10-3
# 20v30: two-sided, noise-filtered
print('  CFU 20v30 (two-sided, noise-filtered):')
sig1 = paired_sig(cfu_24h, 20, 30, use_log=True)
# Note: 20v30 keeps two-sided (no directional hypothesis)
add_sig_bracket(ax, sig1, 1, 2, 3 + y_range*0.82, y_range*0.025)

# 30v40: one-sided H1: 30°C CFU > 40°C CFU
# Pairs with |Δlog₁₀| < 1 order treated as tied (within plate counting error)
print('  CFU 30v40 (one-sided, H1: 30>40, noise-filtered):')
sig2 = paired_sig_custom(cfu_24h, 30, 40, use_log=True,
                         noise_threshold=NOISE_THRESHOLD)
add_sig_bracket(ax, sig2, 2, 3, 3 + y_range*0.90, y_range*0.025)

plt.tight_layout(pad=0.3)
for ext in ['pdf', 'png']:
    fig.savefig(os.path.join(OUTDIR, f'endpoint_24h_CFU.{ext}'),
                dpi=300, bbox_inches='tight', facecolor='white')
plt.close(fig)
print('  Saved: endpoint_24h_CFU.pdf/.png\n')

# ══════════════════════════════════════════════════════════════
# 3) OD/log10(CFU) ratio — paired, one-sided (H1: 40>others)
# ══════════════════════════════════════════════════════════════
print('=' * 50)
print('3) OD/log10(CFU) ratio (paired, one-sided H1: higher temp → higher ratio)')
print('=' * 50)

fig, ax = plt.subplots(figsize=(35/25.4, 45/25.4))
rng = np.random.RandomState(42)

vals_per_temp_r = {}
for t in temps:
    vals_per_temp_r[t] = [v for v in ratio_24h[t].values()
                          if v is not None and not np.isnan(v)]

for k, t in enumerate(temps):
    x_pos = k + 1
    vals = vals_per_temp_r[t]
    jx = x_pos + rng.uniform(-JIT_WIDTH, JIT_WIDTH, len(vals))
    ax.scatter(jx, vals, s=JIT_SIZE, color=COL_TEMP[t],
               alpha=JIT_ALPHA, edgecolors='none', zorder=2)

for k, t in enumerate(temps):
    x_pos = k + 1
    vals = np.array(vals_per_temp_r[t])
    mu = np.mean(vals)
    sem = np.std(vals, ddof=1) / np.sqrt(len(vals)) if len(vals) > 1 else 0
    ax.errorbar(x_pos, mu, yerr=sem, fmt='none',
                ecolor=COL_TEMP[t], elinewidth=ERR_LW,
                capsize=ERR_W*25, capthick=ERR_LW, zorder=3)
    ax.plot(x_pos, mu, 'o', markersize=DOT_SIZE,
            markerfacecolor='white', markeredgecolor=COL_TEMP[t],
            markeredgewidth=MEAN_STROKE, zorder=4)

ax.set_xlim(0.6, 3.4)
ax.set_xticks([1, 2, 3])
ax.set_xticklabels(temp_labels, fontsize=FONT_AX)
ax.set_ylim(0, 0.35)
ax.set_yticks([0, 0.1, 0.2, 0.3])
ax.set_xlabel(''); ax.set_ylabel(''); ax.set_title('')
for spine in ax.spines.values():
    spine.set_visible(True); spine.set_color('black'); spine.set_linewidth(BORDER_SIZE)
ax.tick_params(axis='both', which='both', direction='in',
               top=True, right=True, labeltop=False, labelright=False,
               length=TICK_LEN, width=TICK_SIZE, pad=3)

y_range_r = 0.35

# 20v30: one-sided H1: 30>20
print('  Ratio 20v30 (one-sided, H1: 30>20):')
sig_r1 = paired_sig_custom(ratio_24h, 30, 20, use_log=False)
add_sig_bracket(ax, sig_r1, 1, 2, 0 + y_range_r*0.82, y_range_r*0.025)

# 30v40: one-sided H1: 40>30
print('  Ratio 30v40 (one-sided, H1: 40>30):')
sig_r2 = paired_sig_custom(ratio_24h, 40, 30, use_log=False)
add_sig_bracket(ax, sig_r2, 2, 3, 0 + y_range_r*0.90, y_range_r*0.025)

plt.tight_layout(pad=0.3)
for ext in ['pdf', 'png']:
    fig.savefig(os.path.join(OUTDIR, f'endpoint_24h_ratio.{ext}'),
                dpi=300, bbox_inches='tight', facecolor='white')
plt.close(fig)
print('  Saved: endpoint_24h_ratio.pdf/.png\n')

print('Done. All files in figures/endpoint_24h/')
