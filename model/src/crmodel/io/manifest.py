"""Run manifest and config snapshot utilities."""

from __future__ import annotations

import json
import shutil
import subprocess
from datetime import datetime
from pathlib import Path


def write_manifest(path: Path, data: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        json.dump(data, f, indent=2, sort_keys=True, default=str)


def git_info() -> dict:
    info = {}
    try:
        info["hash"] = subprocess.check_output(
            ["git", "rev-parse", "HEAD"],
            stderr=subprocess.DEVNULL,
        ).decode().strip()
        info["dirty"] = bool(
            subprocess.check_output(
                ["git", "status", "--porcelain"],
                stderr=subprocess.DEVNULL,
            ).strip()
        )
    except (subprocess.CalledProcessError, FileNotFoundError):
        info["hash"] = "unknown"
        info["dirty"] = None
    return info


def snapshot_run(run_dir: Path, config_path: str | Path, cfg: dict,
                 command: str) -> Path:
    """Save a self-contained snapshot of config and metadata into run_dir."""
    snap_dir = run_dir / "snapshot"
    snap_dir.mkdir(parents=True, exist_ok=True)

    config_path = Path(config_path)
    shutil.copy2(config_path, snap_dir / config_path.name)

    if "extends" in open(config_path).read():
        parent_name = None
        import yaml
        with open(config_path) as f:
            raw = yaml.safe_load(f)
        if "extends" in raw:
            parent_path = config_path.parent / raw["extends"]
            if parent_path.exists():
                shutil.copy2(parent_path, snap_dir / parent_path.name)

    manifest = {
        "created_at": datetime.now().isoformat(timespec="seconds"),
        "command": command,
        "git": git_info(),
        "resolved_config": cfg,
    }
    write_manifest(snap_dir / "run_manifest.json", manifest)
    return snap_dir
