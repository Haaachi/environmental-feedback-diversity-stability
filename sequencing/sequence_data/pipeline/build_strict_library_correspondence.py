#!/usr/bin/env python3
"""Build a concise experiment-species to library-ID correspondence table.

This combines:
  1. fragment-to-full-length sequence matches from match_species_to_library.py
  2. full-length BLAST annotations from Library_taxonomy_groups_v3.xlsx

The output intentionally keeps only the correspondence result, not alignment
details. Candidate IDs are ordered with YB_ IDs first, numeric IDs second, and
other prefixes last.
"""

from __future__ import annotations

import argparse
import csv
import re
from copy import copy
from pathlib import Path
from typing import Any

import pandas as pd

from english_schema import normalize_frame


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MATCHES = (
    ROOT
    / "species_annotations"
    / "library_matching"
    / "Experimental_species_to_All_Strains_library.tsv"
)
DEFAULT_LIBRARY_TAXONOMY = ROOT / "Library_taxonomy_groups_v3.xlsx"
DEFAULT_OUT_DIR = ROOT / "species_annotations" / "library_matching"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create concise strict correspondence between experiment species and library IDs."
    )
    parser.add_argument("--matches", type=Path, default=DEFAULT_MATCHES)
    parser.add_argument("--library-taxonomy", type=Path, default=DEFAULT_LIBRARY_TAXONOMY)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument(
        "--prefix",
        default="Experimental_species_library_strict_correspondence",
        help="Output file prefix.",
    )
    return parser.parse_args()


def split_ids(value: str) -> list[str]:
    return [part.strip() for part in str(value or "").split(";") if part.strip()]


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


def clean(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    text = str(value).strip()
    return "" if text in {"nan", "None", "NA", "N/A", "-"} else text


def tokens(value: str) -> set[str]:
    stop = {"unclassified", "sp", "species"}
    return {
        part.lower()
        for part in re.split(r"[^A-Za-z]+", clean(value))
        if part and part.lower() not in stop
    }


def normalize_species(value: str) -> str:
    pieces = [part.lower() for part in re.split(r"[^A-Za-z]+", clean(value)) if part]
    return " ".join(pieces[:2])


def family_from_original_taxon(value: str) -> str:
    match = re.search(r"f__([^;]+)", clean(value))
    return match.group(1).strip() if match else ""


def read_matches(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


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


def genus_family_lookup(taxonomy: dict[str, dict[str, str]]) -> dict[str, str]:
    families_by_genus: dict[str, set[str]] = {}
    for tax_row in taxonomy.values():
        family = clean(tax_row.get("Family", ""))
        if not family:
            continue
        for genus_token in tokens(tax_row.get("Genus", "")):
            families_by_genus.setdefault(genus_token, set()).add(family)
    return {
        genus: next(iter(families))
        for genus, families in families_by_genus.items()
        if len(families) == 1
    }


def experiment_family(row: dict[str, str], genus_to_family: dict[str, str]) -> str:
    family = clean(row.get("Family", ""))
    if family:
        return family.replace("Unclassified_", "")
    genus = clean(row.get("Genus", ""))
    if genus.startswith("Unclassified_"):
        return genus.replace("Unclassified_", "")
    inferred = {
        genus_to_family[token]
        for token in tokens(genus)
        if token in genus_to_family
    }
    if len(inferred) == 1:
        return next(iter(inferred))
    return family_from_original_taxon(row.get("Original_Taxon", ""))


CONSISTENCY_LABELS = {
    4: "species_agreement",
    3: "genus_agreement",
    2: "family_agreement",
    1: "annotation_conflict",
    0: "missing_library_full_length_annotation",
}


def consistency_rank(
    match_row: dict[str, str],
    tax_row: dict[str, str] | None,
    genus_to_family: dict[str, str],
) -> int:
    if tax_row is None:
        return 0

    exp_species = normalize_species(match_row.get("Species_name", ""))
    lib_species = normalize_species(tax_row.get("Species", ""))
    if exp_species and lib_species and exp_species == lib_species:
        return 4

    exp_genus = clean(match_row.get("Genus", ""))
    lib_genus = clean(tax_row.get("Genus", ""))
    if exp_genus and not exp_genus.startswith("Unclassified_"):
        if tokens(exp_genus) & tokens(lib_genus):
            return 3

    exp_family = experiment_family(match_row, genus_to_family)
    lib_family = clean(tax_row.get("Family", ""))
    if exp_family and lib_family and exp_family.lower() == lib_family.lower():
        return 2

    return 1


def taxon_text(strain_id: str, tax_row: dict[str, str] | None) -> str:
    if tax_row is None:
        return f"{strain_id}=missing_library_full_length_annotation"
    species = tax_row["Species"] or tax_row["Blast_annotation"] or tax_row["Genus"]
    group_id = tax_row["Group_ID"]
    return f"{strain_id}={species}" + (f"({group_id})" if group_id else "")


def build_correspondence(
    matches: list[dict[str, str]],
    taxonomy: dict[str, dict[str, str]],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    genus_to_family = genus_family_lookup(taxonomy)

    for match_row in matches:
        tied_ids = split_ids(match_row["Best_tied_library_IDs"])
        candidates = []
        for strain_id in tied_ids:
            tax_row = taxonomy.get(strain_id)
            rank = consistency_rank(match_row, tax_row, genus_to_family)
            candidates.append((strain_id, tax_row, rank))

        best_rank = max((rank for _, _, rank in candidates), default=0)
        if best_rank >= 2:
            kept = [(sid, tax, rank) for sid, tax, rank in candidates if rank == best_rank]
        elif any(tax is not None for _, tax, _ in candidates):
            kept = [(sid, tax, rank) for sid, tax, rank in candidates if tax is not None]
        else:
            kept = candidates

        kept.sort(key=lambda item: library_id_sort_key(item[0]))
        kept_ids = [sid for sid, _, _ in kept]
        kept_taxa = [taxon_text(sid, tax) for sid, tax, _ in kept]
        labels = sorted({CONSISTENCY_LABELS[rank] for _, _, rank in kept})
        consistency = CONSISTENCY_LABELS[best_rank] if best_rank >= 2 else ";".join(labels)

        dropped_missing = [
            sid for sid, tax, _ in candidates if tax is None and sid not in set(kept_ids)
        ]
        note = ""
        if dropped_missing:
            note = "Excluded tied candidates lacking full-length library annotations: " + ";".join(
                sorted(dropped_missing, key=library_id_sort_key)
            )

        rows.append(
            {
                "Experiment": match_row["Experiment"],
                "Experimental_species_ID": match_row["Species"],
                "Experimental_annotation": match_row["Unified_annotation"],
                "Library_IDs": ";".join(kept_ids),
                "Library_full_length_annotation": "; ".join(kept_taxa),
                "Annotation_consistency": consistency,
                "Candidate_count": len(kept_ids),
                "Notes": note,
            }
        )

    return rows


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    fields = list(rows[0].keys())
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def write_xlsx(path: Path, rows: list[dict[str, object]]) -> None:
    df = pd.DataFrame(rows)
    with pd.ExcelWriter(path, engine="openpyxl") as writer:
        df.to_excel(writer, index=False, sheet_name="correspondence")
        ws = writer.book["correspondence"]
        ws.freeze_panes = "A2"
        widths = {
            "A": 14,
            "B": 18,
            "C": 34,
            "D": 80,
            "E": 110,
            "F": 14,
            "G": 10,
            "H": 42,
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

    matches = read_matches(args.matches)
    taxonomy = read_library_taxonomy(args.library_taxonomy)
    rows = build_correspondence(matches, taxonomy)

    tsv_path = args.out_dir / f"{args.prefix}.tsv"
    xlsx_path = args.out_dir / f"{args.prefix}.xlsx"
    write_tsv(tsv_path, rows)
    write_xlsx(xlsx_path, rows)

    consistency_counts = pd.Series([row["Annotation_consistency"] for row in rows]).value_counts().to_dict()
    print(f"Wrote {len(rows)} rows: {tsv_path}")
    print(f"Wrote {len(rows)} rows: {xlsx_path}")
    print("Consistency counts:", consistency_counts)


if __name__ == "__main__":
    main()
