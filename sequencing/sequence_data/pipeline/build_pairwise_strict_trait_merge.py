#!/usr/bin/env python3
"""Build a pairwise-fragment strict merged list for trait measurement.

This is intentionally stricter than library-candidate-set merging:

* Do not merge different species within the same experiment.
* Merge across experiments only when mortality V4a is a 100% exact substring
  of the temperature V4V5 representative sequence.
* The exact cross-experiment relationship must be unique on both sides.
* Near-exact or non-unique exact relationships are kept independent and noted.
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
DEFAULT_STRICT = (
    ROOT
    / "species_annotations"
    / "library_matching"
    / "Experimental_species_library_strict_correspondence.tsv"
)
DEFAULT_LIBRARY_TAXONOMY = ROOT / "Library_taxonomy_groups_v3.xlsx"
DEFAULT_TEMPERATURE_FASTA = (
    ROOT / "16s_project" / "species_identification_standard" / "Temperature_species.fasta"
)
DEFAULT_MORTALITY_FASTA = (
    ROOT / "Mortality_analysis" / "species_identification_standard" / "Mortality_species.fasta"
)
DEFAULT_OUT_DIR = ROOT / "species_annotations" / "library_matching"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create pairwise-fragment strict trait-measurement merge table."
    )
    parser.add_argument("--strict-correspondence", type=Path, default=DEFAULT_STRICT)
    parser.add_argument("--library-taxonomy", type=Path, default=DEFAULT_LIBRARY_TAXONOMY)
    parser.add_argument("--temperature-fasta", type=Path, default=DEFAULT_TEMPERATURE_FASTA)
    parser.add_argument("--mortality-fasta", type=Path, default=DEFAULT_MORTALITY_FASTA)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument(
        "--prefix",
        default="Experimental_species_trait_measurement_pairwise_strict",
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


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1].upper()


def iter_fasta(path: Path):
    name: str | None = None
    parts: list[str] = []
    with path.open("r", encoding="utf-8-sig") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(parts).upper()
                name = line[1:].split()[0]
                parts = []
            else:
                parts.append(line)
    if name is not None:
        yield name, "".join(parts).upper()


def read_species_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    for record_id, sequence in iter_fasta(path):
        species = record_id.split("|")[0]
        records[species] = sequence
    return records


def best_substring_identity(query: str, target: str) -> dict[str, object]:
    best = {
        "identity_pct": -1.0,
        "matches": -1,
        "orientation": "forward",
        "target_start_1based": 0,
    }
    for orientation, query_seq in [
        ("forward", query.upper()),
        ("reverse_complement", reverse_complement(query)),
    ]:
        qlen = len(query_seq)
        if len(target) < qlen:
            continue
        for start in range(len(target) - qlen + 1):
            window = target[start : start + qlen]
            matches = sum(a == b for a, b in zip(query_seq, window))
            identity = (matches / qlen) * 100.0
            if matches > int(best["matches"]):
                best = {
                    "identity_pct": identity,
                    "matches": matches,
                    "orientation": orientation,
                    "target_start_1based": start + 1,
                }
    return best


def read_strict(path: Path) -> pd.DataFrame:
    df = normalize_frame(pd.read_csv(path, sep="\t", encoding="utf-8-sig"))
    df["Candidate_IDs"] = df["Library_IDs"].map(split_ids)
    df["Candidate_Set"] = df["Candidate_IDs"].map(set)
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
            "Species": clean(row["Species"]),
            "Blast_annotation": clean(row["Original_BLAST_annotation"]),
            "Genus": clean(row["Genus"]),
        }
    return taxonomy


def taxon_text(ids: list[str], taxonomy: dict[str, dict[str, str]]) -> str:
    parts: list[str] = []
    for strain_id in ids:
        tax = taxonomy.get(strain_id)
        if tax is None:
            parts.append(f"{strain_id}=missing_library_full_length_annotation")
            continue
        label = tax["Species"] or tax["Blast_annotation"] or tax["Genus"]
        group_id = tax["Group_ID"]
        parts.append(f"{strain_id}={label}" + (f"({group_id})" if group_id else ""))
    return "; ".join(parts)


def build_pairwise_table(
    temperature_seq: dict[str, str],
    mortality_seq: dict[str, str],
) -> tuple[dict[tuple[str, str], dict[str, object]], dict[str, list[str]], dict[str, list[str]]]:
    pairwise: dict[tuple[str, str], dict[str, object]] = {}
    exact_by_temp: dict[str, list[str]] = {species: [] for species in temperature_seq}
    exact_by_mort: dict[str, list[str]] = {species: [] for species in mortality_seq}

    for temp_species, temp_sequence in temperature_seq.items():
        for mort_species, mort_sequence in mortality_seq.items():
            result = best_substring_identity(mort_sequence, temp_sequence)
            pairwise[(temp_species, mort_species)] = result
            if round(float(result["identity_pct"]), 9) == 100.0:
                exact_by_temp[temp_species].append(mort_species)
                exact_by_mort[mort_species].append(temp_species)

    for values in exact_by_temp.values():
        values.sort(key=natural_key)
    for values in exact_by_mort.values():
        values.sort(key=natural_key)
    return pairwise, exact_by_temp, exact_by_mort


def best_counterpart(
    species: str,
    experiment: str,
    pairwise: dict[tuple[str, str], dict[str, object]],
) -> tuple[str, float]:
    best_species = ""
    best_identity = -1.0
    for (temp_species, mort_species), result in pairwise.items():
        if experiment == "temperature" and temp_species != species:
            continue
        if experiment == "mortality" and mort_species != species:
            continue
        counterpart = mort_species if experiment == "temperature" else temp_species
        identity = float(result["identity_pct"])
        if identity > best_identity:
            best_species = counterpart
            best_identity = identity
    return best_species, best_identity


def note_for_independent(
    species: str,
    experiment: str,
    exact_by_temp: dict[str, list[str]],
    exact_by_mort: dict[str, list[str]],
    pairwise: dict[tuple[str, str], dict[str, object]],
) -> tuple[str, str]:
    if experiment == "temperature":
        exact = exact_by_temp.get(species, [])
        if exact:
            details = []
            for mort in exact:
                details.append(f"mortality {mort} exactly matches multiple temperature {','.join(exact_by_mort[mort])}")
            return "retained_separately_nonunique_fragment", "; ".join(details)
    else:
        exact = exact_by_mort.get(species, [])
        if exact:
            details = []
            for temp in exact:
                details.append(f"temperature {temp} exactly matches multiple mortality {','.join(exact_by_temp[temp])}")
            return "retained_separately_nonunique_fragment", "; ".join(details)

    counterpart, identity = best_counterpart(species, experiment, pairwise)
    if counterpart:
        other = "mortality" if experiment == "temperature" else "temperature"
        return "retained_separately_no_unique_exact_match", f"Best cross-experiment fragment match: {other} {counterpart}, identity={identity:.3f}%"
    return "retained_separately_no_cross_experiment_match", ""


def row_from_single(
    source_row: pd.Series,
    taxonomy: dict[str, dict[str, str]],
    status: str,
    note: str,
) -> dict[str, object]:
    candidate_ids = sorted_library_ids(set(source_row["Candidate_Set"]))
    experiment = clean(source_row["Experiment"])
    return {
        "Merged_ID": "",
        "Merge_status": status,
        "Temperature_species": f"{source_row['Experimental_species_ID']}({source_row['Experimental_annotation']})"
        if experiment == "temperature"
        else "",
        "Mortality_species": f"{source_row['Experimental_species_ID']}({source_row['Experimental_annotation']})"
        if experiment == "mortality"
        else "",
        "Temperature_annotation": clean(source_row["Experimental_annotation"]) if experiment == "temperature" else "",
        "Mortality_annotation": clean(source_row["Experimental_annotation"]) if experiment == "mortality" else "",
        "Recommended_library_IDs": ";".join(candidate_ids),
        "Preferred_library_ID": candidate_ids[0] if candidate_ids else "",
        "Library_full_length_annotation": taxon_text(candidate_ids, taxonomy),
        "Fragment_alignment_basis": "Not automatically merged",
        "Original_species_count": 1,
        "Notes": note,
    }


def row_from_pair(
    temp_row: pd.Series,
    mort_row: pd.Series,
    result: dict[str, object],
    taxonomy: dict[str, dict[str, str]],
) -> dict[str, object]:
    shared = set(temp_row["Candidate_Set"]) & set(mort_row["Candidate_Set"])
    candidate_ids = sorted_library_ids(shared)
    return {
        "Merged_ID": "",
        "Merge_status": "auto_merged_unique_reciprocal_exact_fragment",
        "Temperature_species": f"{temp_row['Experimental_species_ID']}({temp_row['Experimental_annotation']})",
        "Mortality_species": f"{mort_row['Experimental_species_ID']}({mort_row['Experimental_annotation']})",
        "Temperature_annotation": clean(temp_row["Experimental_annotation"]),
        "Mortality_annotation": clean(mort_row["Experimental_annotation"]),
        "Recommended_library_IDs": ";".join(candidate_ids),
        "Preferred_library_ID": candidate_ids[0] if candidate_ids else "",
        "Library_full_length_annotation": taxon_text(candidate_ids, taxonomy),
        "Fragment_alignment_basis": (
            "mortality V4a matches the corresponding temperature V4V5 fragment exactly (100%); "
            f"start={result['target_start_1based']}; orientation={result['orientation']}; "
            "unique on both sides"
        ),
        "Original_species_count": 2,
        "Notes": "",
    }


def build_rows(
    strict: pd.DataFrame,
    taxonomy: dict[str, dict[str, str]],
    pairwise: dict[tuple[str, str], dict[str, object]],
    exact_by_temp: dict[str, list[str]],
    exact_by_mort: dict[str, list[str]],
) -> list[dict[str, object]]:
    by_key = {(row["Experiment"], row["Experimental_species_ID"]): row for _, row in strict.iterrows()}

    auto_pairs: list[tuple[str, str]] = []
    for temp_species, mort_species_list in exact_by_temp.items():
        if len(mort_species_list) != 1:
            continue
        mort_species = mort_species_list[0]
        if exact_by_mort.get(mort_species) == [temp_species]:
            auto_pairs.append((temp_species, mort_species))

    rows: list[dict[str, object]] = []
    used: set[tuple[str, str]] = set()
    for temp_species, mort_species in sorted(auto_pairs, key=lambda pair: natural_key(pair[0])):
        temp_key = ("temperature", temp_species)
        mort_key = ("mortality", mort_species)
        rows.append(
            row_from_pair(
                by_key[temp_key],
                by_key[mort_key],
                pairwise[(temp_species, mort_species)],
                taxonomy,
            )
        )
        used.add(temp_key)
        used.add(mort_key)

    for _, source_row in strict.iterrows():
        key = (source_row["Experiment"], source_row["Experimental_species_ID"])
        if key in used:
            continue
        status, note = note_for_independent(
            source_row["Experimental_species_ID"],
            source_row["Experiment"],
            exact_by_temp,
            exact_by_mort,
            pairwise,
        )
        rows.append(row_from_single(source_row, taxonomy, status, note))

    status_order = {
        "auto_merged_unique_reciprocal_exact_fragment": 0,
        "retained_separately_nonunique_fragment": 1,
        "retained_separately_no_unique_exact_match": 2,
        "retained_separately_no_cross_experiment_match": 3,
    }
    rows.sort(
        key=lambda row: (
            status_order.get(str(row["Merge_status"]), 9),
            str(row["Temperature_annotation"] or row["Mortality_annotation"]).lower(),
            natural_key(str(row["Preferred_library_ID"])),
        )
    )
    for index, row in enumerate(rows, start=1):
        row["Merged_ID"] = f"pairwise_species_{index:03d}"
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
        df.to_excel(writer, index=False, sheet_name="pairwise_strict_units")
        ws = writer.book["pairwise_strict_units"]
        ws.freeze_panes = "A2"
        widths = {
            "A": 20,
            "B": 28,
            "C": 38,
            "D": 38,
            "E": 28,
            "F": 32,
            "G": 90,
            "H": 18,
            "I": 100,
            "J": 64,
            "K": 12,
            "L": 72,
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
    temperature_seq = read_species_fasta(args.temperature_fasta)
    mortality_seq = read_species_fasta(args.mortality_fasta)
    pairwise, exact_by_temp, exact_by_mort = build_pairwise_table(temperature_seq, mortality_seq)
    rows = build_rows(strict, taxonomy, pairwise, exact_by_temp, exact_by_mort)

    tsv_path = args.out_dir / f"{args.prefix}.tsv"
    xlsx_path = args.out_dir / f"{args.prefix}.xlsx"
    write_tsv(tsv_path, rows)
    write_xlsx(xlsx_path, rows)

    summary = pd.DataFrame(rows)
    print(f"Wrote {len(rows)} pairwise-strict units: {tsv_path}")
    print(f"Wrote {len(rows)} pairwise-strict units: {xlsx_path}")
    print("Status counts:", summary["Merge_status"].value_counts().to_dict())
    print("Original species covered:", int(summary["Original_species_count"].sum()))


if __name__ == "__main__":
    main()
