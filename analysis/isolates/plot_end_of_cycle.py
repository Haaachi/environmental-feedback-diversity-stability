from pathlib import Path
import sys
import os

ISOLATE_ROOT = Path(os.environ.get("ISOLATE_WORKSPACE", Path(__file__).resolve().parent))

import openpyxl
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

XLSX = str(ISOLATE_ROOT / 'Summary.xlsx')
OUTDIR = str(ISOLATE_ROOT / 'figures')
os.makedirs(OUTDIR, exist_ok=True)

wb = openpyxl.load_workbook(XLSX, data_only=True)
temps = [20, 30, 40]

# Find end-of-cycle columns: 24, 48, 72, ... (multiples of 24)
time_header_ws = wb['20OD']
all_cols = []
for idx, cell in enumerate(time_header_ws[1]):
    if cell.value is not None and str(cell.value).strip() != '':
        try:
            h = float(cell.value)
            all_cols.append((idx + 1, h))  # 1-based column index, hour
        except (ValueError, TypeError):
            pass

cycle_cols = [(col, h) for col, h in all_cols if h > 0 and h % 24 == 0]
cycle_days = [h / 24 for _, h in cycle_cols]
print(f"End-of-cycle time points: {[h for _, h in cycle_cols]} h = days {cycle_days}")

def read_cycle_values(sheet_name):
    ws = wb[sheet_name]
    data = {}
    for row in ws.iter_rows(min_row=2, max_row=ws.max_row, values_only=False):
        if row[0].value is None:
            break
        sid = row[0].value
        species_name = row[1].value
        values = []
        for col, h in cycle_cols:
            cell = ws.cell(row=row[0].row, column=col)
            if cell.value is not None and str(cell.value).strip() != '':
                try:
                    values.append(float(cell.value))
                except (ValueError, TypeError):
                    values.append(np.nan)
            else:
                values.append(np.nan)
        data[sid] = (species_name, values)
    return data

od_data = {}
ph_data = {}
for t in temps:
    od_data[t] = read_cycle_values(f'{t}OD')
    ph_data[t] = read_cycle_values(f'{t}pH')

species_ids = list(od_data[20].keys())

# Fixed axis ranges (same as within-cycle plots)
od_min = 0.02
od_max = 3.0
ph_min = 2.5
ph_max = 9.5

OD_COLOR = '#2166AC'
PH_COLOR = '#B2182B'

plt.rcParams.update({
    'font.family': 'Arial',
    'font.size': 10,
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
    'axes.linewidth': 1.0,
    'xtick.major.width': 1.0,
    'ytick.major.width': 1.0,
    'xtick.major.size': 4,
    'ytick.major.size': 4,
    'xtick.direction': 'in',
    'ytick.direction': 'in',
})

x_max = 5.5

for sid in species_ids:
    species_name = od_data[20][sid][0]

    fig, axes = plt.subplots(1, 3, figsize=(12, 3.2), sharey=False)
    fig.suptitle(species_name, fontsize=13, fontweight='bold',
                 fontstyle='italic', fontfamily='Arial', y=1.02)

    for i, t in enumerate(temps):
        ax_od = axes[i]
        ax_ph = ax_od.twinx()

        od_vals = od_data[t][sid][1]
        ph_vals = ph_data[t][sid][1]

        # OD
        valid_od = [(d, v) for d, v in zip(cycle_days, od_vals) if not np.isnan(v) and v > 0]
        if valid_od:
            x_od, y_od = zip(*valid_od)
            ax_od.plot(x_od, y_od, 'o-', color=OD_COLOR, markersize=7,
                       linewidth=1.8, markeredgecolor='white', markeredgewidth=0.5,
                       label='OD$_{600}$', zorder=3)

        # pH
        valid_ph = [(d, v) for d, v in zip(cycle_days, ph_vals) if not np.isnan(v)]
        if valid_ph:
            x_ph, y_ph = zip(*valid_ph)
            ax_ph.plot(x_ph, y_ph, 's--', color=PH_COLOR, markersize=7,
                       linewidth=1.8, markeredgecolor='white', markeredgewidth=0.5,
                       label='pH', zorder=3)

        # OD axis (linear)
        ax_od.set_ylim(0, 2.0)
        ax_od.set_xlim(0.5, x_max)
        ax_od.set_xlabel('Cycle (day)', fontfamily='Arial', fontsize=11)
        ax_od.set_title(f'{t} °C', fontsize=11, fontfamily='Arial', pad=8)
        ax_od.tick_params(axis='y', colors=OD_COLOR)
        ax_od.set_xticks([1, 2, 3, 4, 5])

        # pH axis
        ax_ph.set_ylim(ph_min, ph_max)
        ax_ph.tick_params(axis='y', colors=PH_COLOR)
        ax_ph.set_yticks(np.arange(3, 10, 1))

        if i == 0:
            ax_od.set_ylabel('OD$_{600}$', fontfamily='Arial', fontsize=11, color=OD_COLOR)
        else:
            ax_od.set_yticklabels([])

        if i == 2:
            ax_ph.set_ylabel('pH', fontfamily='Arial', fontsize=11,
                             color=PH_COLOR, rotation=270, labelpad=15)
        else:
            ax_ph.set_yticklabels([])

        for spine in ax_od.spines.values():
            spine.set_visible(True)
            spine.set_color('black')
        for spine in ax_ph.spines.values():
            spine.set_visible(True)
            spine.set_color('black')

        if i == 0:
            lines_od, labels_od = ax_od.get_legend_handles_labels()
            lines_ph, labels_ph = ax_ph.get_legend_handles_labels()
            ax_od.legend(lines_od + lines_ph, labels_od + labels_ph,
                         loc='upper left', frameon=False, fontsize=9)

    plt.tight_layout(w_pad=2.5)
    outpath = os.path.join(OUTDIR, f'cycle_{sid}_{species_name.replace(" ", "_").replace(".", "")}.png')
    fig.savefig(outpath, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f'Saved: {outpath}')

print('\nAll end-of-cycle figures saved.')
