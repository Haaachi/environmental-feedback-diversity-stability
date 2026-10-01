"""Visualize emergent oscillation experiment statistics."""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import os

matplotlib.rcParams.update({
    "font.family": "Arial", "font.size": 10,
    "axes.linewidth": 1.2, "figure.dpi": 150,
    "savefig.dpi": 300, "savefig.bbox": "tight",
    "savefig.facecolor": "white", "pdf.fonttype": 42,
})

OUT_DIR = os.path.dirname(os.path.abspath(__file__))

# ── Data from the experiment ────────────────────────────────────────
hs = np.array([0.0, 0.05, 0.06, 0.10, 0.15])

# Single-species classification
n_sta_1sp = np.array([1000, 939, 892, 731, 574])
n_osc_1sp = np.array([0, 8, 27, 102, 176])
n_col_1sp = np.array([0, 53, 81, 167, 250])

# Community assembly results (out of 200)
# A: stable-only → osc
n_A_osc = np.array([5, 10, 19, 47, 89])
n_A_total = 200
# B: osc-only → stable (using extra sampling where available)
# h=0.06: 6/500 extra, h=0.10: 0/500 extra, h=0.15: 1/500 extra
n_B_sta     = np.array([np.nan, np.nan, 6, 0, 1])
n_B_total   = np.array([np.nan, np.nan, 500, 500, 500])

# Metrics: stable_only condition
# Oscillating communities
A_osc_cv   = np.array([0.241, 0.294, 0.264, 0.292, 0.387])
A_osc_cv_e = np.array([0.107, 0.235, 0.137, 0.153, 0.211])
A_osc_sh   = np.array([1.218, 1.321, 1.382, 1.086, 1.361])
A_osc_sh_e = np.array([0.435, 0.323, 0.254, 0.328, 0.342])
A_osc_ri   = np.array([5.50, 4.91, 5.18, 3.94, 5.87])
A_osc_ri_e = np.array([2.26, 1.45, 1.00, 1.16, 1.58])
A_osc_phi  = np.array([0.664, 0.761, 0.761, 0.867, 0.672])
A_osc_phi_e= np.array([0.317, 0.262, 0.258, 0.178, 0.288])
A_osc_PR   = np.array([2.84, 3.68, 3.32, 2.77, 3.64])
A_osc_PR_e = np.array([0.87, 1.14, 1.20, 1.07, 1.41])

# Stable communities (from stable_only)
A_sta_sh   = np.array([1.260, 1.159, 1.078, 0.819, 0.744])
A_sta_sh_e = np.array([0.371, 0.456, 0.486, 0.521, 0.564])
A_sta_ri   = np.array([5.22, 4.47, 4.18, 3.29, 3.05])
A_sta_ri_e = np.array([1.60, 1.68, 1.69, 1.58, 1.81])

# osc_only condition (h >= 0.06)
B_osc_cv  = np.array([np.nan, np.nan, 0.141, 0.206, 0.408])
B_osc_phi = np.array([np.nan, np.nan, 0.817, 0.957, 0.975])
B_osc_sh  = np.array([np.nan, np.nan, 1.022, 0.944, 1.320])
B_osc_PR  = np.array([np.nan, np.nan, 2.68, 2.42, 3.53])

# ── Colors ──────────────────────────────────────────────────────────
C_STA = '#4393C3'    # blue — stable
C_OSC = '#D6604D'    # red — oscillating
C_COL = '#878787'    # grey — collapsed
C_A   = '#E66101'    # orange — emergent osc (stable→osc)
C_B   = '#5E3C99'    # purple — emergent stability (osc→stable)
C_A_osc = '#D6604D'
C_A_sta = '#4393C3'
C_B_osc = '#FC8D59'


def main():
    # ════════════════════════════════════════════════════════════════
    #  Figure 1: Single-species classification + community assembly
    # ════════════════════════════════════════════════════════════════
    fig, axes = plt.subplots(1, 3, figsize=(13, 4), constrained_layout=True)

    # Panel A: Single-species stacked bar
    ax = axes[0]
    w = 0.012
    ax.bar(hs, n_sta_1sp/1000, w, color=C_STA, label='Stable')
    ax.bar(hs, n_osc_1sp/1000, w, bottom=n_sta_1sp/1000, color=C_OSC, label='Oscillating')
    ax.bar(hs, n_col_1sp/1000, w, bottom=(n_sta_1sp+n_osc_1sp)/1000, color=C_COL, label='Collapsed')
    ax.set_xlabel(r'$h_{\mathrm{stress}}$')
    ax.set_ylabel('Fraction')
    ax.set_title('(A) Single-species classification', fontsize=11)
    ax.legend(fontsize=8, loc='center right')
    ax.set_xlim(-0.02, 0.17)
    ax.set_ylim(0, 1.05)

    # Panel B: Emergent oscillation rate
    ax = axes[1]
    frac_A = n_A_osc / n_A_total
    ax.plot(hs, frac_A, 'o-', ms=7, lw=2, color=C_A,
            markeredgecolor='white', markeredgewidth=0.8, zorder=3)
    ax.fill_between(hs, 0, frac_A, alpha=0.15, color=C_A)
    ax.set_xlabel(r'$h_{\mathrm{stress}}$')
    ax.set_ylabel('Fraction oscillating')
    ax.set_title('(B) Emergent oscillation\n(stable singles → osc community)', fontsize=11)
    ax.set_xlim(-0.02, 0.17)
    ax.set_ylim(0, 0.55)
    for i, h in enumerate(hs):
        ax.annotate(f'{n_A_osc[i]}/200', (h, frac_A[i]),
                    textcoords='offset points', xytext=(8, 5), fontsize=8, color=C_A)

    # Panel C: Emergent stability rate
    ax = axes[2]
    mask_B = np.isfinite(n_B_sta)
    hs_B = hs[mask_B]
    frac_B = n_B_sta[mask_B] / n_B_total[mask_B]
    ax.plot(hs_B, frac_B, 's-', ms=7, lw=2, color=C_B,
            markeredgecolor='white', markeredgewidth=0.8, zorder=3)
    ax.fill_between(hs_B, 0, frac_B, alpha=0.15, color=C_B)
    ax.set_xlabel(r'$h_{\mathrm{stress}}$')
    ax.set_ylabel('Fraction stable')
    ax.set_title('(C) Emergent stability\n(osc singles → stable community)', fontsize=11)
    ax.set_xlim(0.04, 0.17)
    ax.set_ylim(-0.002, 0.025)
    for i, h in enumerate(hs_B):
        n = int(n_B_sta[mask_B][i])
        t = int(n_B_total[mask_B][i])
        ax.annotate(f'{n}/{t}', (h, frac_B[i]),
                    textcoords='offset points', xytext=(8, 5), fontsize=8, color=C_B)

    fig.savefig(os.path.join(OUT_DIR, 'stat_fig1_classification.png'))
    fig.savefig(os.path.join(OUT_DIR, 'stat_fig1_classification.pdf'))
    plt.close(fig)
    print('Saved: stat_fig1_classification')

    # ════════════════════════════════════════════════════════════════
    #  Figure 2: Metrics comparison (stable_only condition)
    # ════════════════════════════════════════════════════════════════
    fig2, axes2 = plt.subplots(2, 3, figsize=(13, 7.5), constrained_layout=True)

    # Row 1: Emergent oscillation metrics (stable_only → osc vs stable)
    # Shannon
    ax = axes2[0, 0]
    ax.errorbar(hs, A_osc_sh, A_osc_sh_e, fmt='o-', ms=5, lw=1.5, capsize=3,
                color=C_A_osc, label='Oscillating', zorder=3)
    ax.errorbar(hs, A_sta_sh, A_sta_sh_e, fmt='s-', ms=5, lw=1.5, capsize=3,
                color=C_A_sta, label='Stable', zorder=3)
    ax.set_ylabel('Shannon $H$')
    ax.set_title('Shannon (stable-only assembly)', fontsize=10)
    ax.legend(fontsize=8)

    # Richness
    ax = axes2[0, 1]
    ax.errorbar(hs, A_osc_ri, A_osc_ri_e, fmt='o-', ms=5, lw=1.5, capsize=3,
                color=C_A_osc, label='Oscillating', zorder=3)
    ax.errorbar(hs, A_sta_ri, A_sta_ri_e, fmt='s-', ms=5, lw=1.5, capsize=3,
                color=C_A_sta, label='Stable', zorder=3)
    ax.set_ylabel('Richness')
    ax.set_title('Richness (stable-only assembly)', fontsize=10)
    ax.legend(fontsize=8)

    # CV
    ax = axes2[0, 2]
    ax.errorbar(hs, A_osc_cv, A_osc_cv_e, fmt='o-', ms=5, lw=1.5, capsize=3,
                color=C_A_osc, zorder=3)
    ax.axhline(0.1, color='0.7', lw=0.8, ls='--')
    ax.set_ylabel('Community CV')
    ax.set_title('CV of emergent oscillators', fontsize=10)

    # Row 2: phi, PR, and comparison A vs B
    # phi
    ax = axes2[1, 0]
    ax.errorbar(hs, A_osc_phi, A_osc_phi_e, fmt='o-', ms=5, lw=1.5, capsize=3,
                color=C_A, label='A: stable→osc', zorder=3)
    mask = np.isfinite(B_osc_phi)
    ax.errorbar(hs[mask], B_osc_phi[mask], fmt='s-', ms=5, lw=1.5,
                color=C_B, label='B: osc→osc (from osc singles)', zorder=3)
    ax.set_ylabel(r'Synchrony $\varphi$')
    ax.set_xlabel(r'$h_{\mathrm{stress}}$')
    ax.set_title(r'$\varphi$ by source condition', fontsize=10)
    ax.legend(fontsize=8)
    ax.set_ylim(0.3, 1.05)

    # PR
    ax = axes2[1, 1]
    ax.errorbar(hs, A_osc_PR, A_osc_PR_e, fmt='o-', ms=5, lw=1.5, capsize=3,
                color=C_A, label='A: stable→osc', zorder=3)
    mask = np.isfinite(B_osc_PR)
    ax.errorbar(hs[mask], B_osc_PR[mask], fmt='s-', ms=5, lw=1.5,
                color=C_B, label='B: osc→osc', zorder=3)
    ax.set_ylabel('Participation ratio')
    ax.set_xlabel(r'$h_{\mathrm{stress}}$')
    ax.set_title('PR by source condition', fontsize=10)
    ax.legend(fontsize=8)
    ax.set_ylim(1, 5.5)

    # CV comparison A vs B
    ax = axes2[1, 2]
    ax.errorbar(hs, A_osc_cv, A_osc_cv_e, fmt='o-', ms=5, lw=1.5, capsize=3,
                color=C_A, label='A: stable→osc', zorder=3)
    mask = np.isfinite(B_osc_cv)
    ax.errorbar(hs[mask], B_osc_cv[mask], fmt='s-', ms=5, lw=1.5,
                color=C_B, label='B: osc→osc', zorder=3)
    ax.axhline(0.1, color='0.7', lw=0.8, ls='--')
    ax.set_ylabel('Community CV')
    ax.set_xlabel(r'$h_{\mathrm{stress}}$')
    ax.set_title('CV: emergent vs inherited osc', fontsize=10)
    ax.legend(fontsize=8)

    fig2.savefig(os.path.join(OUT_DIR, 'stat_fig2_metrics.png'))
    fig2.savefig(os.path.join(OUT_DIR, 'stat_fig2_metrics.pdf'))
    plt.close(fig2)
    print('Saved: stat_fig2_metrics')

    # ════════════════════════════════════════════════════════════════
    #  Figure 3: Dual-axis summary — single-sp osc fraction + emergent osc fraction
    # ════════════════════════════════════════════════════════════════
    fig3, ax3 = plt.subplots(figsize=(6, 4), constrained_layout=True)

    frac_1sp_osc = n_osc_1sp / 1000
    frac_1sp_col = n_col_1sp / 1000

    l1, = ax3.plot(hs, frac_1sp_osc, 'D-', ms=6, lw=1.8, color=C_OSC,
                   markeredgecolor='white', markeredgewidth=0.6,
                   label='Single-sp oscillating', zorder=3)
    l2, = ax3.plot(hs, frac_1sp_col, '^-', ms=6, lw=1.8, color=C_COL,
                   markeredgecolor='white', markeredgewidth=0.6,
                   label='Single-sp collapsed', zorder=3)
    l3, = ax3.plot(hs, frac_A, 'o-', ms=7, lw=2.2, color=C_A,
                   markeredgecolor='white', markeredgewidth=0.8,
                   label='Emergent osc (stable→osc)', zorder=4)
    ax3.set_xlabel(r'$h_{\mathrm{stress}}$', fontsize=12)
    ax3.set_ylabel('Fraction', fontsize=12)
    ax3.set_xlim(-0.01, 0.16)
    ax3.set_ylim(0, 0.52)
    ax3.legend(handles=[l1, l2, l3], fontsize=9, loc='upper left')
    ax3.set_title('Single-species vs emergent community oscillation', fontsize=11)

    fig3.savefig(os.path.join(OUT_DIR, 'stat_fig3_summary.png'))
    fig3.savefig(os.path.join(OUT_DIR, 'stat_fig3_summary.pdf'))
    plt.close(fig3)
    print('Saved: stat_fig3_summary')

    # ════════════════════════════════════════════════════════════════
    #  Figure 4: Showcase p_opt distributions
    # ════════════════════════════════════════════════════════════════
    showcase_data = {
        'emergent_highCV\n(h=0.05)': [6.34, 7.14, 6.01, 7.34, 7.78, 7.99,
                                       7.92, 7.75, 5.83, 8.08, 6.29, 7.95],
        'emergent_lowdiv\n(h=0.10)': [7.09, 6.37, 5.81, 5.83, 7.39, 5.77,
                                       5.93, 6.51, 6.13, 7.23, 7.29, 7.71],
        'emergent_highdiv\n(h=0.15)': [7.11, 6.96, 7.70, 7.26, 6.10, 7.08,
                                        7.12, 6.53, 6.77, 8.12, 6.86, 7.23],
        'stability\n(h=0.06)': [8.33, 5.65, 5.64, 5.64, 5.68, 8.35,
                                 8.35, 8.30, 8.31, 5.59, 8.30, 8.39],
    }

    fig4, axes4 = plt.subplots(1, 4, figsize=(14, 3.5), constrained_layout=True,
                                sharey=True)

    for ax, (title, popts) in zip(axes4, showcase_data.items()):
        popts = np.array(popts)
        colors_sc = [C_A if 'emergent' in title else C_B] * len(popts)
        ax.barh(range(len(popts)), popts - 7.0, left=7.0,
                color=[C_A if p >= 7 else C_STA for p in popts],
                edgecolor='white', linewidth=0.5, height=0.7)
        ax.axvline(7.0, color='0.5', lw=0.8, ls='--')
        ax.set_xlim(5.0, 9.0)
        ax.set_xlabel(r'$p_{\mathrm{opt}}$')
        ax.set_title(title, fontsize=10)
        ax.set_yticks(range(len(popts)))
        ax.set_yticklabels([f'sp{i+1:02d}' for i in range(len(popts))], fontsize=7)
        ax.invert_yaxis()
        # Add text values
        for i, p in enumerate(popts):
            ax.text(p + 0.05 if p < 8 else p - 0.05, i, f'{p:.1f}',
                    va='center', ha='left' if p < 8 else 'right', fontsize=7)

    fig4.savefig(os.path.join(OUT_DIR, 'stat_fig4_popt_distribution.png'))
    fig4.savefig(os.path.join(OUT_DIR, 'stat_fig4_popt_distribution.pdf'))
    plt.close(fig4)
    print('Saved: stat_fig4_popt_distribution')


if __name__ == '__main__':
    main()
