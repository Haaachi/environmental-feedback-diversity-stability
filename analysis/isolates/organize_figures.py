from pathlib import Path
import sys
import os

ISOLATE_ROOT = Path(os.environ.get("ISOLATE_WORKSPACE", Path(__file__).resolve().parent))

"""
Organize all figures into task-based folders and export PDF copies.
Folder structure:
  figures_organized/
    1_overview/                  — overview panels (OD, pH, CFU, with/without mean)
    2_within_cycle_OD_pH/        — per-species within-cycle OD+pH dual-axis (3 temps per fig)
    3_end_of_cycle_daily/        — per-species end-of-cycle daily OD+pH (3 temps per fig)
    4_per_species_OD_CFU_pH/     — per-species per-temp OD+CFU and pH panels
    5_death_examples/            — selected death example panels
    6_endpoint_24h/              — 24h endpoint bar charts (OD, CFU)
    7_growth_rate/               — growth rate analysis figures
"""

import shutil
import glob
from PIL import Image

SRC = str(ISOLATE_ROOT / 'figures')
DST = str(ISOLATE_ROOT / 'figures_organized')

folders = {
    '1_overview': [],
    '2_within_cycle_OD_pH': [],
    '3_end_of_cycle_daily': [],
    '4_per_species_OD_CFU_pH': [],
    '5_death_examples': [],
    '6_endpoint_24h': [],
    '7_growth_rate': [],
}

# Classify files in figures/ root
for f in sorted(os.listdir(SRC)):
    fpath = os.path.join(SRC, f)
    if not os.path.isfile(fpath) or not f.endswith('.png'):
        continue

    if f.startswith('overview_'):
        folders['1_overview'].append(fpath)
    elif f.startswith('species') and not f.startswith('species') == False:
        # within-cycle OD+pH per species (species1_Klebsiella_sp1.png)
        folders['2_within_cycle_OD_pH'].append(fpath)
    elif f.startswith('cycle_'):
        folders['3_end_of_cycle_daily'].append(fpath)
    elif f.startswith('death_'):
        folders['5_death_examples'].append(fpath)
    elif f.startswith('endpoint_'):
        folders['6_endpoint_24h'].append(fpath)
    elif f.startswith('growth_rate'):
        folders['7_growth_rate'].append(fpath)

# Fix: species* files that are within-cycle
species_files = glob.glob(os.path.join(SRC, 'species*.png'))
folders['2_within_cycle_OD_pH'] = sorted(species_files)

# per_species_OD_CFU_pH subfolder
sub_dir = os.path.join(SRC, 'per_species_OD_CFU_pH')
if os.path.isdir(sub_dir):
    for f in sorted(os.listdir(sub_dir)):
        if f.endswith('.png'):
            folders['4_per_species_OD_CFU_pH'].append(os.path.join(sub_dir, f))

# Create folders, copy PNG, convert to PDF
total = 0
for folder_name, file_list in folders.items():
    out_dir = os.path.join(DST, folder_name)
    os.makedirs(out_dir, exist_ok=True)

    for src_path in file_list:
        fname = os.path.basename(src_path)
        # Copy PNG
        dst_png = os.path.join(out_dir, fname)
        shutil.copy2(src_path, dst_png)

        # Convert to PDF
        pdf_name = fname.replace('.png', '.pdf')
        dst_pdf = os.path.join(out_dir, pdf_name)
        img = Image.open(src_path)
        # Convert RGBA to RGB for PDF
        if img.mode == 'RGBA':
            bg = Image.new('RGB', img.size, (255, 255, 255))
            bg.paste(img, mask=img.split()[3])
            img = bg
        elif img.mode != 'RGB':
            img = img.convert('RGB')
        img.save(dst_pdf, 'PDF', resolution=300)
        img.close()
        total += 1

    n = len(file_list)
    print(f'{folder_name}: {n} figures ({n} PNG + {n} PDF)')

print(f'\nTotal: {total} figure pairs (PNG+PDF)')
print(f'Output: {DST}')
