#!/usr/bin/env python3
"""Build a conservative merged species list for trait measurement.

The goal is to reduce trait-measurement workload across the temperature and
mortality experiments without rewriting the original species IDs. Rows are
merged only when their library-ID candidate sets strongly overlap and their
experiment/library annotations agree at species or genus level. Family-level
or conflicting rows are kept independent.
"""

from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict
from copy import copy
from pathlib import Path
from typing import Any

import pandas as pd

from english_schema import normalize_frame


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_STRICT = (
    ROOT
    / "species_annotations"
    / "library_matching"
    / "Experimental_species_library_strict_correspondence.tsv"
)
DEFAULT_LIBRARY_TAXONOMY = ROOT / "Library_taxonomy_groups_v3.xlsx"
DEFAULT_OUT_DIR = ROOT / "species_annotations" / "library_matching"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Conservatively merge experiment species into trait-measurement units."
    )
    parser.add_argument("--strict-correspondence", type=Path, default=DEFAULT_STRICT)
    parser.add_argument("--library-taxonomy", type=Path, default=DEFAULT_LIBRARY_TAXONOMY)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument(
        "--prefix",
        default="Experimental_species_trait_measurement_merged",
    )
    parser.add_argument(
        "--min-overlap-of-smaller",
        type=float,
        default=0.80,
        help="Merge eligible rows when shared library IDs cover this fraction of the smaller set.",
    )
    return parser.parse_args()


def clean(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    text = str(value).strip()
    return "" if text in {"nan", "None", "NA", "N/A", "-"} else text


def split_ids(value: object) -> list[str]:
    return [part.strip() for part in clean(value).split(";") if part.strip()]


def natural_key(value: str) -> list[Any]:
    return [int(part) if part.isdigit() else part.lower() for part in re.split(r"(\d+)", value)]


def library_id_sort_key(value: str) -> tuple[int, list[Any]]:
    if value.startswith("YB_"):
        group = 0
    elif re.match(r"^\d", value):
        group = 1
    else:
        group = 2
    return group, natural_key(value)


def sorted_library_ids(ids: set[str] | list[str]) -> list[str]:
    return sorted(ids, key=library_id_sort_key)


def tokens(value: str) -> set[str]:
    stop = {"unclassified", "sp", "species"}
    return {
        part.lower()
        for part in re.split(r"[^A-Za-z]+", clean(value))
        if part and part.lower() not in stop
    }


def annotation_tokens(annotation: str) -> set[str]:
    return tokens(annotation.split(" sp.")[0])


def annotations_compatible(a: dict[str, object], b: dict[str, object]) -> bool:
    return bool(annotation_tokens(str(a["Experimental_annotation"])) & annotation_tokens(str(b["Experimental_annotation"])))


def read_strict(path: Path) -> pd.DataFrame:
    df = normalize_frame(pd.read_csv(path, sep="\t", encoding="utf-8-sig"))
    df["Candidate_IDs"] = df["Library_IDs"].map(split_ids)
    df["Candidate_Set"] = df["Candidate_IDs"].map(set)
    df["Auto_merge_eligible"] = df["Annotation_consistency"].isin(["species_agreement", "genus_agreement"])
    return df


def read_library_taxonomy(path: Path) -> dict[str, dict[str, str]]:
    df = normalize_frame(pd.read_excel(path, sheet_name=0, header=2))
    df = df[df["Strain_ID"].notna()].copy()
    df["Strain_ID"] = df["Strain_ID"].astype(str).str.strip()

    taxonomy: dict[str, dict[str, str]] = {}
    for _, row in df.iterrows():
        strain_id = clean(row["Strain_ID"])
        taxonomy[strain_id] = {
            "Group_ID": clean(row["Group ID"]),
            "Group_type": clean(row["Group_type"]),
            "Family": clean(row["Family"]),
            "Genus": clean(row["Genus"]),
            "Species": clean(row["Species"]),
            "Blast_annotation": clean(row["Original_BLAST_annotation"]),
        }
    return taxonomy


def should_merge(a: dict[str, object], b: dict[str, object], min_overlap: float) -> bool:
    if not bool(a["Auto_merge_eligible"]) or not bool(b["Auto_merge_eligible"]):
        return False
    if not annotations_compatible(a, b):
        return False
    set_a = set(a["Candidate_Set"])
    set_b = set(b["Candidate_Set"])
    if not set_a or not set_b:
        return False
    overlap = len(set_a & set_b) / min(len(set_a), len(set_b))
    return overlap >= min_overlap


def connected_components(rows: list[dict[str, object]], min_overlap: float) -> list[list[int]]:
    parent = list(range(len(rows)))

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a: int, b: int) -> None:
        ra = find(a)
        rb = find(b)
        if ra != rb:
            parent[rb] = ra

    for i in range(len(rows)):
        for j in range(i + 1, len(rows)):
            if should_merge(rows[i], rows[j], min_overlap):
                union(i, j)

    groups: dict[int, list[int]] = defaultdict(list)
    for i in range(len(rows)):
        groups[find(i)].append(i)
    return list(groups.values())


def member_text(rows: list[dict[str, object]], experiment: str) -> str:
    members = [
        f"{row['Experimental_species_ID']}({row['Experimental_annotation']})"
        for row in rows
        if row["Experiment"] == experiment
    ]
    return "; ".join(members)


def taxon_text(ids: list[str], taxonomy: dict[str, dict[str, str]]) -> str:
    parts: list[str] = []
    for strain_id in ids:
        tax = taxonomy.get(strain_id)
        if not tax:
            parts.append(f"{strain_id}=missing_library_full_length_annotation")
            continue
        label = tax["Species"] or tax["Blast_annotation"] or tax["Genus"]
        group_id = tax["Group_ID"]
        parts.append(f"{strain_id}={label}" + (f"({group_id})" if group_id else ""))
    return "; ".join(parts)


def library_groups(ids: list[str], taxonomy: dict[str, dict[str, str]]) -> str:
    groups = {
        taxonomy[strain_id]["Group_ID"]
        for strain_id in ids
        if strain_id in taxonomy and taxonomy[strain_id]["Group_ID"]
    }
    return ";".join(sorted(groups, key=natural_key))


def representative_annotation(rows: list[dict[str, object]]) -> str:
    annotations = list(dict.fromkeys(clean(row["Experimental_annotation"]) for row in rows))
    return "; ".join(annotations)


def merge_basis(rows: list[dict[str, object]], candidate_ids: list[str], used_intersection: bool) -> tuple[str, str]:
    if any(row["Annotation_consistency"] == "annotation_conflict" for row in rows):
        return "Annotation conflict; conservatively retained separately", "review_needed"
    if any(row["Annotation_consistency"] == "family_agreement" for row in rows):
        return "Family-level agreement only; conservatively retained separately", "low"
    if len(rows) == 1:
        return "Single experiment/species retained separately", "independent"

    candidate_sets = [set(row["Candidate_Set"]) for row in rows]
    if all(candidate_sets[0] == candidate_set for candidate_set in candidate_sets[1:]):
        return "Identical candidate library ID sets", "high"
    if used_intersection:
        return "Highly overlapping candidate IDs; use shared intersection", "medium"
    return "Highly overlapping candidate IDs with empty shared intersection; use union", "low"


def build_merged_rows(
    df: pd.DataFrame,
    taxonomy: dict[str, dict[str, str]],
    min_overlap: float,
) -> list[dict[str, object]]:
    source_rows = df.to_dict("records")
    components = connected_components(source_rows, min_overlap)

    merged: list[dict[str, object]] = []
    for component in components:
        rows = [source_rows[i] for i in component]
        candidate_sets = [set(row["Candidate_Set"]) for row in rows]
        auto_group = all(row["Auto_merge_eligible"] for row in rows)

        used_intersection = False
        if len(rows) > 1 and auto_group:
            shared = set.intersection(*candidate_sets)
            if shared:
                candidate_ids = sorted_library_ids(shared)
                used_intersection = not all(shared == candidate_set for candidate_set in candidate_sets)
            else:
                candidate_ids = sorted_library_ids(set.union(*candidate_sets))
        else:
            candidate_ids = sorted_library_ids(set.union(*candidate_sets))

        basis, confidence = merge_basis(rows, candidate_ids, used_intersection)
        experiments = sorted({row["Experiment"] for row in rows})
        if experiments == ["mortality", "temperature"]:
            coverage = "both_experiments"
        elif experiments == ["temperature"]:
            coverage = "temperature_only"
        else:
            coverage = "mortality_only"

        original_ids = "; ".join(
            f"{row['Experiment']}:{row['Experimental_species_ID']}" for row in sorted(rows, key=lambda r: (r["Experiment"], natural_key(r["Experimental_species_ID"])))
        )
        notes: list[str] = []
        if len(rows) > 1:
            per_experiment_counts = defaultdict(int)
            for row in rows:
                per_experiment_counts[row["Experiment"]] += 1
            multi = [
                f"{experiment} contains {count} species"
                for experiment, count in sorted(per_experiment_counts.items())
                if count > 1
            ]
            if multi:
                notes.append("The same measurement unit contains " + ", ".join(multi) + "; original species IDs are retained")
        if any(row["Annotation_consistency"] == "annotation_conflict" for row in rows):
            notes.append("Automatic merging with other units is not recommended")

        merged.append(
            {
                "Merged_ID": "",
                "Experiment_coverage": coverage,
                "Temperature_species": member_text(rows, "temperature"),
                "Mortality_species": member_text(rows, "mortality"),
                "Merged_experimental_annotation": representative_annotation(rows),
                "Recommended_library_IDs": ";".join(candidate_ids),
                "Preferred_library_ID": candidate_ids[0] if candidate_ids else "",
                "Library_Group_ID": library_groups(candidate_ids, taxonomy),
                "Library_full_length_annotation": taxon_text(candidate_ids, taxonomy),
                "Merge_basis": basis,
                "Merge_confidence": confidence,
                "Original_species_count": len(rows),
                "Candidate_library_ID_count": len(candidate_ids),
                "Original_species_index": original_ids,
                "Notes": "; ".join(notes),
            }
        )

    def output_sort_key(row: dict[str, object]) -> tuple[int, str, list[Any]]:
        coverage_order = {"both_experiments": 0, "temperature_only": 1, "mortality_only": 2}
        return (
            coverage_order.get(str(row["Experiment_coverage"]), 9),
            str(row["Merged_experimental_annotation"]).lower(),
            natural_key(str(row["Preferred_library_ID"])),
        )

    merged.sort(key=output_sort_key)
    for index, row in enumerate(merged, start=1):
        row["Merged_ID"] = f"merged_species_{index:03d}"
    return merged


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    fields = list(rows[0].keys())
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def write_xlsx(path: Path, rows: list[dict[str, object]]) -> None:
    df = pd.DataFrame(rows)
    with pd.ExcelWriter(path, engine="openpyxl") as writer:
        df.to_excel(writer, index=False, sheet_name="trait_measurement_units")
        ws = writer.book["trait_measurement_units"]
        ws.freeze_panes = "A2"
        widths = {
            "A": 18,
            "B": 14,
            "C": 42,
            "D": 42,
            "E": 40,
            "F": 90,
            "G": 18,
            "H": 18,
            "I": 100,
            "J": 28,
            "K": 12,
            "L": 12,
            "M": 14,
            "N": 48,
            "O": 54,
        }
        for column, width in widths.items():
            ws.column_dimensions[column].width = width
        for row in ws.iter_rows():
            for cell in row:
                alignment = copy(cell.alignment)
                alignment.wrap_text = True
                alignment.vertical = "top"
                cell.alignment = alignment
        ws.auto_filter.ref = ws.dimensions


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    strict = read_strict(args.strict_correspondence)
    taxonomy = read_library_taxonomy(args.library_taxonomy)
    merged = build_merged_rows(strict, taxonomy, args.min_overlap_of_smaller)

    tsv_path = args.out_dir / f"{args.prefix}.tsv"
    xlsx_path = args.out_dir / f"{args.prefix}.xlsx"
    write_tsv(tsv_path, merged)
    write_xlsx(xlsx_path, merged)

    summary = pd.DataFrame(merged)
    print(f"Wrote {len(merged)} merged units: {tsv_path}")
    print(f"Wrote {len(merged)} merged units: {xlsx_path}")
    print("Coverage counts:", summary["Experiment_coverage"].value_counts().to_dict())
    print("Confidence counts:", summary["Merge_confidence"].value_counts().to_dict())


if __name__ == "__main__":
    main()
