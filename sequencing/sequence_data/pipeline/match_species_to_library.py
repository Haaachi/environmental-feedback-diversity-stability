#!/usr/bin/env python3
"""Match experiment species amplicons to full-length 16S library records.

The two experiment pipelines produce representative species sequences from
short 16S regions: temperature uses V4V5, mortality uses V4a. This script
compares those fragments against every record in All_Strains_Combined.fasta
and reports the library records that best explain each experimental species.

Matching is intentionally local to the provided library and does not call any
external taxonomy database. The primary pass is an exhaustive ungapped search
for the best same-length window in each full-length sequence, in both query
orientations. If numba is available, the best candidates are refined with a
semi-global gapped alignment of the full query to a local library window.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Iterable

import numpy as np
from numpy.lib.stride_tricks import sliding_window_view

try:
    from numba import njit

    HAVE_NUMBA = True
except Exception:  # pragma: no cover - optional acceleration/refinement
    HAVE_NUMBA = False


ROOT = Path(__file__).resolve().parents[1]

DEFAULT_LIBRARY = ROOT / "All_Strains_Combined.fasta"
DEFAULT_ANNOTATIONS = (
    ROOT
    / "species_annotations"
    / "unified_standard"
    / "Unified_species_annotation_mapping.tsv"
)
DEFAULT_OUT_DIR = ROOT / "species_annotations" / "library_matching"

EXPERIMENT_FASTA = {
    "temperature": (
        ROOT
        / "16s_project"
        / "species_identification_standard"
        / "Temperature_species.fasta",
        "V4V5",
    ),
    "mortality": (
        ROOT
        / "Mortality_analysis"
        / "species_identification_standard"
        / "Mortality_species.fasta",
        "V4a",
    ),
}

RC_TABLE = str.maketrans("ACGTRYKMSWBDHVNUacgtrykmswbdhvnu", "TGCAYRMKSWVHDBNAtgcayrmkswvhdbna")


@dataclass(frozen=True)
class FastaRecord:
    record_id: str
    sequence: str
    sequence_array: np.ndarray


@dataclass(frozen=True)
class QueryRecord:
    experiment: str
    region: str
    species: str
    record_id: str
    sequence: str
    feature_id: str
    representative_asv: str


@dataclass(frozen=True)
class MatchHit:
    experiment: str
    region: str
    species: str
    query_id: str
    query_length: int
    library_id: str
    library_length: int
    orientation: str
    method: str
    identity_pct: float
    query_coverage_pct: float
    alignment_length: int
    matches: int
    mismatches: int
    gap_query: int
    gap_library: int
    library_start_1based: int
    library_end_1based: int
    score: int

    @property
    def total_gaps(self) -> int:
        return self.gap_query + self.gap_library


if HAVE_NUMBA:

    @njit(cache=True)
    def semiglobal_align(q: np.ndarray, r: np.ndarray) -> tuple[int, int, int, int, int, int, int, int]:
        """Align all query bases to the best substring of r.

        Returns:
            score, matches, mismatches, gap_query, gap_library, start, end,
            alignment_length. Coordinates are 0-based half-open on r.
        """
        match_score = 2
        mismatch_score = -3
        gap_score = -5

        n = q.shape[0]
        m = r.shape[0]
        prev = np.zeros(m + 1, dtype=np.int32)
        curr = np.empty(m + 1, dtype=np.int32)
        traceback = np.zeros((n + 1, m + 1), dtype=np.uint8)

        for i in range(1, n + 1):
            curr[0] = prev[0] + gap_score
            traceback[i, 0] = 2
            qi = q[i - 1]
            for j in range(1, m + 1):
                diag = prev[j - 1] + (match_score if qi == r[j - 1] else mismatch_score)
                up = prev[j] + gap_score
                left = curr[j - 1] + gap_score
                best = diag
                direction = 1
                if up > best:
                    best = up
                    direction = 2
                if left > best:
                    best = left
                    direction = 3
                curr[j] = best
                traceback[i, j] = direction
            tmp = prev
            prev = curr
            curr = tmp

        best_j = 0
        best_score = prev[0]
        for j in range(1, m + 1):
            if prev[j] > best_score:
                best_score = prev[j]
                best_j = j

        i = n
        j = best_j
        matches = 0
        mismatches = 0
        gap_query = 0
        gap_library = 0
        while i > 0:
            direction = traceback[i, j]
            if direction == 1 and j > 0:
                if q[i - 1] == r[j - 1]:
                    matches += 1
                else:
                    mismatches += 1
                i -= 1
                j -= 1
            elif direction == 2 or j == 0:
                gap_library += 1
                i -= 1
            else:
                gap_query += 1
                j -= 1

        align_len = matches + mismatches + gap_query + gap_library
        return best_score, matches, mismatches, gap_query, gap_library, j, best_j, align_len


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Match final experiment species to full-length 16S library records."
    )
    parser.add_argument("--library", type=Path, default=DEFAULT_LIBRARY)
    parser.add_argument("--annotations", type=Path, default=DEFAULT_ANNOTATIONS)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--top-hits", type=int, default=20)
    parser.add_argument("--gapped-candidates", type=int, default=80)
    parser.add_argument("--gapped-slack-mismatches", type=int, default=5)
    parser.add_argument("--gapped-window-pad", type=int, default=40)
    parser.add_argument(
        "--no-gapped-refinement",
        action="store_true",
        help="Use only exhaustive same-length window matching.",
    )
    return parser.parse_args()


def iter_fasta(path: Path) -> Iterable[tuple[str, str]]:
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


def to_array(sequence: str) -> np.ndarray:
    return np.frombuffer(sequence.encode("ascii"), dtype=np.uint8)


def reverse_complement(sequence: str) -> str:
    return sequence.translate(RC_TABLE)[::-1].upper()


def read_library(path: Path) -> list[FastaRecord]:
    records: list[FastaRecord] = []
    for record_id, sequence in iter_fasta(path):
        records.append(FastaRecord(record_id, sequence, to_array(sequence)))
    return records


def parse_query_header(experiment: str, region: str, record_id: str, sequence: str) -> QueryRecord:
    parts = record_id.split("|")
    species = parts[0]
    feature_id = ""
    representative_asv = ""
    for part in parts[1:]:
        if "=" not in part:
            continue
        key, value = part.split("=", 1)
        if key == "FeatureID":
            feature_id = value
        elif key == "Representative_ASV":
            representative_asv = value
    return QueryRecord(
        experiment=experiment,
        region=region,
        species=species,
        record_id=record_id,
        sequence=sequence,
        feature_id=feature_id,
        representative_asv=representative_asv,
    )


def read_queries() -> list[QueryRecord]:
    queries: list[QueryRecord] = []
    for experiment, (path, region) in EXPERIMENT_FASTA.items():
        for record_id, sequence in iter_fasta(path):
            queries.append(parse_query_header(experiment, region, record_id, sequence))
    return queries


def read_annotations(path: Path) -> dict[tuple[str, str], dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    return {(row["Experiment"], row["Species"]): row for row in rows}


def best_ungapped_hit(query: QueryRecord, qseq: str, ref: FastaRecord, orientation: str) -> MatchHit:
    qarr = to_array(qseq)
    qlen = qarr.shape[0]
    if ref.sequence_array.shape[0] < qlen:
        raise ValueError(f"Library record {ref.record_id} is shorter than query {query.record_id}.")
    windows = sliding_window_view(ref.sequence_array, qlen)
    mismatches_by_pos = (windows != qarr).sum(axis=1)
    best_start = int(mismatches_by_pos.argmin())
    mismatches = int(mismatches_by_pos[best_start])
    matches = qlen - mismatches
    score = matches * 2 - mismatches * 3
    identity = (matches / qlen) * 100.0
    return MatchHit(
        experiment=query.experiment,
        region=query.region,
        species=query.species,
        query_id=query.record_id,
        query_length=qlen,
        library_id=ref.record_id,
        library_length=len(ref.sequence),
        orientation=orientation,
        method="ungapped_window",
        identity_pct=identity,
        query_coverage_pct=100.0,
        alignment_length=qlen,
        matches=matches,
        mismatches=mismatches,
        gap_query=0,
        gap_library=0,
        library_start_1based=best_start + 1,
        library_end_1based=best_start + qlen,
        score=score,
    )


def refine_gapped(hit: MatchHit, query: QueryRecord, qseq: str, ref: FastaRecord, pad: int) -> MatchHit:
    if not HAVE_NUMBA:
        return hit
    qarr = to_array(qseq)
    window_start = max(0, hit.library_start_1based - 1 - pad)
    window_end = min(len(ref.sequence), hit.library_end_1based + pad)
    ref_window = ref.sequence_array[window_start:window_end]
    score, matches, mismatches, gap_query, gap_library, start, end, align_len = semiglobal_align(
        qarr, ref_window
    )
    identity = (matches / align_len) * 100.0 if align_len else 0.0
    query_aligned = matches + mismatches + gap_library
    coverage = (query_aligned / len(qseq)) * 100.0 if qseq else 0.0
    return replace(
        hit,
        method="semiglobal_gapped",
        identity_pct=identity,
        query_coverage_pct=coverage,
        alignment_length=align_len,
        matches=matches,
        mismatches=mismatches,
        gap_query=gap_query,
        gap_library=gap_library,
        library_start_1based=window_start + start + 1,
        library_end_1based=window_start + end,
        score=score,
    )


def sort_key(hit: MatchHit) -> tuple[float, float, int, int, int, str, str]:
    return (
        -round(hit.identity_pct, 9),
        -round(hit.query_coverage_pct, 9),
        -hit.matches,
        hit.total_gaps + hit.mismatches,
        -hit.score,
        hit.library_id,
        hit.orientation,
    )


def same_best(a: MatchHit, b: MatchHit) -> bool:
    return (
        round(a.identity_pct, 9) == round(b.identity_pct, 9)
        and round(a.query_coverage_pct, 9) == round(b.query_coverage_pct, 9)
        and a.matches == b.matches
        and a.mismatches == b.mismatches
        and a.gap_query == b.gap_query
        and a.gap_library == b.gap_library
        and a.score == b.score
    )


def status_for_hit(hit: MatchHit) -> str:
    if hit.mismatches == 0 and hit.total_gaps == 0 and hit.query_coverage_pct >= 99.9:
        return "exact_fragment_match"
    if hit.identity_pct >= 99.0 and hit.query_coverage_pct >= 99.0:
        return "near_exact_fragment_match"
    if hit.identity_pct >= 97.0 and hit.query_coverage_pct >= 95.0:
        return "strong_fragment_match"
    if hit.identity_pct >= 95.0 and hit.query_coverage_pct >= 90.0:
        return "moderate_fragment_match"
    return "weak_best_available"


def hit_to_row(hit: MatchHit, rank: int | None = None) -> dict[str, object]:
    row: dict[str, object] = {
        "Experiment": hit.experiment,
        "Species": hit.species,
        "Amplicon_region": hit.region,
        "Query_ID": hit.query_id,
        "Query_length": hit.query_length,
        "Library_ID": hit.library_id,
        "Library_length": hit.library_length,
        "Orientation": hit.orientation,
        "Method": hit.method,
        "Identity_pct": f"{hit.identity_pct:.6f}",
        "Query_coverage_pct": f"{hit.query_coverage_pct:.6f}",
        "Alignment_length": hit.alignment_length,
        "Matches": hit.matches,
        "Mismatches": hit.mismatches,
        "Gap_in_query": hit.gap_query,
        "Gap_in_library": hit.gap_library,
        "Library_start_1based": hit.library_start_1based,
        "Library_end_1based": hit.library_end_1based,
        "Alignment_score": hit.score,
        "Match_status": status_for_hit(hit),
    }
    if rank is not None:
        row = {"Rank": rank, **row}
    return row


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def match_query(
    query: QueryRecord,
    library: list[FastaRecord],
    args: argparse.Namespace,
) -> list[MatchHit]:
    hits: list[MatchHit] = []
    oriented_sequences = [
        ("forward", query.sequence),
        ("reverse_complement", reverse_complement(query.sequence)),
    ]
    ref_lookup = {record.record_id: record for record in library}
    qseq_lookup = dict(oriented_sequences)

    for orientation, qseq in oriented_sequences:
        for ref in library:
            hits.append(best_ungapped_hit(query, qseq, ref, orientation))

    if args.no_gapped_refinement or not HAVE_NUMBA:
        return sorted(hits, key=sort_key)

    ungapped_sorted = sorted(hits, key=lambda h: (h.mismatches, -h.matches, h.library_id))
    best_mismatches = ungapped_sorted[0].mismatches
    selected: set[tuple[str, str]] = set()
    for hit in ungapped_sorted[: args.gapped_candidates]:
        selected.add((hit.library_id, hit.orientation))
    for hit in ungapped_sorted:
        if hit.mismatches <= best_mismatches + args.gapped_slack_mismatches:
            selected.add((hit.library_id, hit.orientation))

    refined: dict[tuple[str, str], MatchHit] = {}
    for hit in hits:
        key = (hit.library_id, hit.orientation)
        if key not in selected:
            continue
        refined[key] = refine_gapped(
            hit,
            query,
            qseq_lookup[hit.orientation],
            ref_lookup[hit.library_id],
            args.gapped_window_pad,
        )

    merged = [refined.get((hit.library_id, hit.orientation), hit) for hit in hits]
    return sorted(merged, key=sort_key)


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    library = read_library(args.library)
    queries = read_queries()
    annotations = read_annotations(args.annotations)

    top_rows: list[dict[str, object]] = []
    ranked_rows: list[dict[str, object]] = []

    for query in queries:
        ranked = match_query(query, library, args)
        top = ranked[0]
        tied = [hit for hit in ranked if same_best(hit, top)]
        top_hit_ids = ";".join(hit.library_id for hit in ranked[: args.top_hits])
        tied_ids = ";".join(hit.library_id for hit in tied)
        tied_oriented_ids = ";".join(f"{hit.library_id}:{hit.orientation}" for hit in tied)

        annotation = annotations.get((query.experiment, query.species), {})
        row = {
            "Experiment": query.experiment,
            "Species": query.species,
            "Unified_annotation": annotation.get("Unified_annotation", ""),
            "Display_Taxon": annotation.get("Display_Taxon", ""),
            "Genus": annotation.get("Genus", ""),
            "Species_name": annotation.get("Species_name", ""),
            "Feature_ID": annotation.get("Feature_ID", query.feature_id),
            "Representative_ASV": annotation.get("Representative_ASV", query.representative_asv),
            "Family": annotation.get("Family", ""),
            "Annotation_confidence": annotation.get("Confidence", ""),
            "Original_Taxon": annotation.get("Original_Taxon", ""),
            "Amplicon_region": query.region,
            "Query_length": top.query_length,
            "Best_library_ID": top.library_id,
            "Best_orientation": top.orientation,
            "Best_identity_pct": f"{top.identity_pct:.6f}",
            "Best_query_coverage_pct": f"{top.query_coverage_pct:.6f}",
            "Best_alignment_length": top.alignment_length,
            "Best_matches": top.matches,
            "Best_mismatches": top.mismatches,
            "Best_gap_in_query": top.gap_query,
            "Best_gap_in_library": top.gap_library,
            "Best_library_start_1based": top.library_start_1based,
            "Best_library_end_1based": top.library_end_1based,
            "Best_alignment_score": top.score,
            "Best_match_status": status_for_hit(top),
            "Best_tie_count": len(tied),
            "Best_tied_library_IDs": tied_ids,
            "Best_tied_library_IDs_with_orientation": tied_oriented_ids,
            f"Top_{args.top_hits}_library_IDs": top_hit_ids,
            "Matching_method": top.method,
        }
        top_rows.append(row)

        for rank, hit in enumerate(ranked[: args.top_hits], start=1):
            ranked_row = hit_to_row(hit, rank)
            ranked_row.update(
                {
                    "Unified_annotation": annotation.get("Unified_annotation", ""),
                    "Display_Taxon": annotation.get("Display_Taxon", ""),
                    "Genus": annotation.get("Genus", ""),
                    "Species_name": annotation.get("Species_name", ""),
                    "Feature_ID": annotation.get("Feature_ID", query.feature_id),
                    "Representative_ASV": annotation.get(
                        "Representative_ASV", query.representative_asv
                    ),
                }
            )
            ranked_rows.append(ranked_row)

    summary_fields = list(top_rows[0].keys())
    ranked_fields = list(ranked_rows[0].keys())
    summary_path = args.out_dir / "Experimental_species_to_All_Strains_library.tsv"
    ranked_path = args.out_dir / "Experimental_species_to_All_Strains_library_ranked_hits.tsv"
    write_tsv(summary_path, top_rows, summary_fields)
    write_tsv(ranked_path, ranked_rows, ranked_fields)

    print(f"Wrote {len(top_rows)} top matches: {summary_path}")
    print(f"Wrote {len(ranked_rows)} ranked hits: {ranked_path}")
    print(f"Gapped refinement: {'yes' if HAVE_NUMBA and not args.no_gapped_refinement else 'no'}")


if __name__ == "__main__":
    main()
