"""Batch-generate stability-diversity phase maps for available datasets."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
PLOT_SCRIPT = ROOT / "scripts" / "plot" / "plot_main_stability_diversity_phase.py"
OUT_ROOT = ROOT / "figures" / "main_stability_diversity"
CV_THRESH = 0.1
OSC_FRAC_MIN = 0.1


DATASETS = [
    ("rh_gnormal_sd3_R1_sparse6_tri_S12_h015_g100_e5000", "data/runs/rh_gnormal_sd3_R1_sparse6_tri_S12_h015_g100_e5000"),
    ("rh_gnormal_sd3_R1_sparse6_tri_S12_g100_e1000", "data/runs/rh_gnormal_sd3_R1_sparse6_tri_S12_g100_e1000"),
    ("rh_gm5p5_tri_S12_g100_e1000_runA", "data/runs/phase_r_h"),
    ("rh_gm5p5_tri_S12_g100_e1000_runB", "data/runs/v2_main/phase_r_h"),
    ("rh_gnormal_sd5_R05_sparse6_tri_S12_g100_e1000", "data/runs/rh_gnormal_sd5_R05_sparse6_tri_S12_g100_e1000"),
    ("rh_gm10p10_tri_S12_g50_e500", "data/runs/prelim/gamma_m10_p10"),
    ("rh_g0_tri_S12_g50_e500", "data/runs/prelim/gamma_zero"),
    ("rh_gm5p10_tri_S12_g50_e500", "data/runs/sensitivity/gamma_m5_p10"),
    ("rh_gm8p10_tri_S12_g50_e500_densepref", "data/runs/sensitivity/dense_pref"),
    ("rh_gm8p10_tri_S12_g50_e500_gcv0", "data/runs/sensitivity/growth_cv_00"),
    ("rh_gm8p10_tri_S1_g50_e500", "data/runs/sensitivity/S_1"),
    ("rh_gm8p10_tri_S24_g50_e500", "data/runs/sensitivity/S_24"),
    ("rh_gm8p10_tri_S36_g50_e500", "data/runs/sensitivity/S_36"),
    ("rh_gm8p10_tri_S48_g50_e500", "data/runs/sensitivity/S_48"),
    ("hf_ggrp_8_4_12_tri_S12_h100_f5_e1000_gcv0", "data/runs/phase_h_f"),
    ("rh_gm8p10_tophat_S12_g100_e1000", "data/runs/rh_gm8p10_tophat_S12_g100_e1000"),
    ("legacy_data_Final", "data/runs/legacy_data_Final"),
    ("legacy_data_Final_v1", "data/runs/legacy_data_Final_v1"),
    ("legacy_data_Final_v3", "data/runs/legacy_data_Final_v3"),
]


def has_plot_input(run_dir: Path) -> bool:
    legacy_ok = (run_dir / "CR_phase2_meta.h5").exists() and (run_dir / "phase2_rows").exists()
    reorganized_ok = (
        (
            (run_dir / "phase_r_h_meta.h5").exists()
            and (run_dir / "raw" / "phase_r_h_rows").exists()
        )
        or (
            (run_dir / "phase_h_f_meta.h5").exists()
            and (run_dir / "raw" / "phase_h_f_rows").exists()
        )
    )
    return legacy_ok or reorganized_ok


def main() -> int:
    failures = []
    suffix = f"cv{int(round(CV_THRESH * 100)):03d}_oscmin{int(round(OSC_FRAC_MIN * 100)):03d}"
    for dataset_id, run_dir_text in DATASETS:
        run_dir = ROOT / run_dir_text
        out_dir = OUT_ROOT / f"{dataset_id}_{suffix}"
        if not has_plot_input(run_dir):
            print(f"[SKIP] {dataset_id}: missing compatible raw input at {run_dir}")
            continue

        cmd = [
            sys.executable,
            str(PLOT_SCRIPT),
            "--run-dir",
            str(run_dir),
            "--out-dir",
            str(out_dir),
            "--cv-thresh",
            str(CV_THRESH),
            "--osc-frac-min",
            str(OSC_FRAC_MIN),
        ]
        print(f"\n=== {dataset_id} ===")
        print(f"run_dir: {run_dir}")
        print(f"out_dir: {out_dir}")
        completed = subprocess.run(cmd, cwd=ROOT)
        if completed.returncode != 0:
            failures.append((dataset_id, completed.returncode))

    if failures:
        print("\nFailures:")
        for dataset_id, code in failures:
            print(f"  {dataset_id}: exit {code}")
        return 1
    print("\nAll compatible datasets completed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
