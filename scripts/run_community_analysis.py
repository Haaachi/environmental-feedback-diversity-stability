"""Run the community analysis in dependency order, without changing parameters."""
from __future__ import annotations

import argparse
from pathlib import Path
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
WORKSPACE = ROOT / "analysis" / "communities"
STEPS = (
    "process_abundance.py",
    "compute_diversity.py",
    "compute_fluctuations.py",
    "calculate_species_decomposition.py",
    "calculate_mechanism_axes.py",
    "compute_community_cv_window_robustness.py",
)
FIGURES = (
    "plot_abundance.R", "plot_biomass.R", "plot_diversity.R",
    "plot_diversity_cv_halves.R", "plot_fluctuation_proportion.R",
    "plot_community_cv_sum_std_dual_axis.R", "plot_species_decomposition.R",
    "plot_community_cv_window_robustness.R", "plot_od_ph_supplement.R",
    "plot_summary_24h_temperature.R",
)

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plots", action="store_true", help="Also generate core R figure panels.")
    args = parser.parse_args()
    rscript = shutil.which("Rscript")
    if args.plots and not rscript:
        parser.error("Rscript is required for --plots. Install R and run scripts/install_r_dependencies.R.")
    for name in STEPS:
        print(f"Running {name}", flush=True)
        subprocess.run([sys.executable, str(WORKSPACE / name)], cwd=WORKSPACE, check=True)
    if args.plots:
        for name in FIGURES:
            print(f"Running {name}", flush=True)
            subprocess.run([rscript, str(WORKSPACE / name)], cwd=WORKSPACE, check=True)

if __name__ == "__main__":
    main()
