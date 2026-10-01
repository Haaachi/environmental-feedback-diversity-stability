#!/usr/bin/env python
"""Check registered CRmodel datasets and their meta/processed files."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = PROJECT_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from crmodel.io.datasets import availability_report, load_catalog


def _mark(value: bool | None) -> str:
    if value is None:
        return "-"
    return "yes" if value else "no"


def print_table(rows: list[dict]) -> None:
    headers = [
        "dataset",
        "phase",
        "role",
        "root",
        "meta",
        "processed",
        "last10",
        "raw",
        "rows",
    ]
    table = []
    for row in rows:
        table.append(
            [
                row["dataset"],
                row["phase"],
                row["role"],
                _mark(row["root_exists"]),
                _mark(row["meta_exists"]),
                _mark(row["processed_exists"]),
                _mark(row["processed_last10_exists"]),
                _mark(row["raw_dir_exists"]),
                "" if row["row_count"] is None else str(row["row_count"]),
            ]
        )

    widths = [
        max(len(headers[i]), *(len(str(r[i])) for r in table))
        for i in range(len(headers))
    ]
    line = "  ".join(headers[i].ljust(widths[i]) for i in range(len(headers)))
    print(line)
    print("  ".join("-" * w for w in widths))
    for r in table:
        print("  ".join(str(r[i]).ljust(widths[i]) for i in range(len(headers))))


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Check data/runs/datasets.yaml availability."
    )
    parser.add_argument(
        "--catalog",
        default=PROJECT_ROOT / "data" / "runs" / "datasets.yaml",
        help="Path to datasets.yaml.",
    )
    parser.add_argument(
        "--json",
        action="store_true",
        help="Emit machine-readable JSON.",
    )
    args = parser.parse_args()

    catalog = load_catalog(args.catalog)
    rows = availability_report(catalog)

    if args.json:
        serializable = []
        for row in rows:
            serializable.append(
                {
                    k: str(v) if isinstance(v, Path) else v
                    for k, v in row.items()
                }
            )
        print(json.dumps(serializable, indent=2, sort_keys=True))
    else:
        print_table(rows)


if __name__ == "__main__":
    main()
