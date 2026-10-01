#!/usr/bin/env python
"""Run a small end-to-end simulation and postprocessing check."""

from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CONFIG = "configs/smoke_test.yaml"
RUN_DIR = ROOT / "data" / "runs" / "smoke_test"


def run(args: list[str]) -> None:
    subprocess.run([sys.executable, *args], cwd=ROOT, check=True)


def main() -> None:
    if RUN_DIR.exists():
        shutil.rmtree(RUN_DIR)

    run([
        "run.py",
        "simulate",
        "--config",
        CONFIG,
        "--workers",
        "1",
        "--no-warmup",
    ])
    run(["run.py", "status", "--config", CONFIG])
    run(["run.py", "postprocess", "--config", CONFIG, "--strict"])

    expected = [
        RUN_DIR / "phase_r_h_meta.h5",
        RUN_DIR / "processed" / "phase_r_h_metrics_last3.npz",
        RUN_DIR / "processed" / "phase_r_h_metrics_last10.npz",
    ]
    missing = [p for p in expected if not p.exists()]
    if missing:
        raise FileNotFoundError("Smoke test missing outputs: " + ", ".join(map(str, missing)))

    print("Smoke test completed successfully.")


if __name__ == "__main__":
    main()
