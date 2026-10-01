"""Dataset registry helpers for simulation outputs.

The registry lives at ``data/runs/datasets.yaml`` and is intentionally small:
it maps stable dataset ids to roots, meta files, processed files, and raw row
directories. Plotting code should resolve data paths through this module when
possible.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml


PROJECT_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_CATALOG = PROJECT_ROOT / "data" / "runs" / "datasets.yaml"


def load_catalog(path: str | Path | None = None) -> dict[str, Any]:
    """Load the dataset catalog YAML."""
    catalog_path = Path(path) if path is not None else DEFAULT_CATALOG
    if not catalog_path.is_absolute():
        catalog_path = PROJECT_ROOT / catalog_path
    with open(catalog_path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def get_dataset(
    dataset_id: str | None = None,
    catalog: dict[str, Any] | None = None,
) -> tuple[str, dict[str, Any]]:
    """Return ``(dataset_id, dataset_entry)`` from the catalog."""
    cat = catalog or load_catalog()
    selected = dataset_id or cat.get("default_dataset")
    datasets = cat.get("datasets", {})
    if selected not in datasets:
        known = ", ".join(sorted(datasets))
        raise KeyError(f"Unknown dataset '{selected}'. Known datasets: {known}")
    return selected, datasets[selected]


def project_path(path: str | Path) -> Path:
    """Resolve a project-relative path."""
    p = Path(path)
    return p if p.is_absolute() else PROJECT_ROOT / p


def dataset_root(dataset: dict[str, Any]) -> Path:
    """Resolve a dataset's canonical root directory."""
    return project_path(dataset["canonical_root"])


def phase_paths(
    dataset_id: str | None = None,
    phase: str | None = None,
    catalog: dict[str, Any] | None = None,
) -> dict[str, Path | None]:
    """Resolve root/meta/processed/raw paths for one dataset phase."""
    _, dataset = get_dataset(dataset_id, catalog)
    phases = dataset.get("phases", {})
    selected_phase = phase or next(iter(phases))
    if selected_phase not in phases:
        known = ", ".join(sorted(phases))
        raise KeyError(f"Unknown phase '{selected_phase}'. Known phases: {known}")

    root = dataset_root(dataset)
    spec = phases[selected_phase]
    out: dict[str, Path | None] = {
        "root": root,
        "meta": root / spec["meta"] if spec.get("meta") else None,
        "processed": root / spec["processed"] if spec.get("processed") else None,
        "processed_last10": (
            root / spec["processed_last10"]
            if spec.get("processed_last10")
            else None
        ),
        "raw_dir": root / spec["raw_dir"] if spec.get("raw_dir") else None,
    }
    return out


def availability_report(catalog: dict[str, Any] | None = None) -> list[dict[str, Any]]:
    """Return one availability row for every registered dataset phase."""
    cat = catalog or load_catalog()
    rows: list[dict[str, Any]] = []
    for dataset_id, dataset in cat.get("datasets", {}).items():
        root = dataset_root(dataset)
        for phase, spec in dataset.get("phases", {}).items():
            raw_dir = root / spec["raw_dir"] if spec.get("raw_dir") else None
            row_count = None
            if raw_dir is not None and raw_dir.exists():
                row_count = len(list(raw_dir.glob("*.h5")))
            rows.append(
                {
                    "dataset": dataset_id,
                    "phase": phase,
                    "role": dataset.get("role", ""),
                    "root": root,
                    "root_exists": root.exists(),
                    "meta": root / spec["meta"] if spec.get("meta") else None,
                    "meta_exists": (root / spec["meta"]).exists()
                    if spec.get("meta")
                    else None,
                    "processed": root / spec["processed"]
                    if spec.get("processed")
                    else None,
                    "processed_exists": (root / spec["processed"]).exists()
                    if spec.get("processed")
                    else None,
                    "processed_last10": root / spec["processed_last10"]
                    if spec.get("processed_last10")
                    else None,
                    "processed_last10_exists": (
                        root / spec["processed_last10"]
                    ).exists()
                    if spec.get("processed_last10")
                    else None,
                    "raw_dir": raw_dir,
                    "raw_dir_exists": raw_dir.exists() if raw_dir is not None else None,
                    "row_count": row_count,
                }
            )
    return rows
