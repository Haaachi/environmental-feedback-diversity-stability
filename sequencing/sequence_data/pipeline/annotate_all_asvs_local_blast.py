#!/usr/bin/env python3
"""Annotate all ASV representative sequences with local BLAST.

This is the standard all-ASV taxonomy step for experiments that do not already
have a QIIME2 ``taxonomy.tsv``. It uses ``exported/dna-sequences.fasta`` as the
query object and emits a QIIME2-style taxonomy table:

    Feature ID    Taxon    Confidence

The ``Confidence`` value is a BLAST-derived numeric score
``identity_fraction * coverage_fraction``. It is not a sklearn classifier
probability, but keeping it numeric and in the same column lets downstream ASV
post-processing consume temperature and mortality taxonomy files identically.
"""

from __future__ import annotations

import argparse
import csv
import os
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Callable


@dataclass
class BlastHit:
    qseqid: str
    sseqid: str
    pident: float
    length: int
    mismatch: int
    gapopen: int
    qstart: int
    qend: int
    qlen: int
    qcovs: float
    evalue: float
    bitscore: float
    staxids: str
    stitle: str
    ssciname: str
    sphylum: str = ""
    sclass: str = ""
    sorder: str = ""
    sfamily: str = ""
    sgenus: str = ""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create QIIME2-style taxonomy.tsv from all ASV representative sequences."
    )
    parser.add_argument(
        "--query",
        type=Path,
        default=Path("exported/dna-sequences.fasta"),
        help="All-ASV representative sequence FASTA, normally exported/dna-sequences.fasta.",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path("exported/local_blast_taxonomy"),
        help="Output directory for taxonomy.tsv and BLAST detail tables.",
    )
    parser.add_argument("--blast-db", default="16S_ribosomal_RNA")
    parser.add_argument(
        "--blastdb-dir",
        type=Path,
        default=None,
        help="Directory containing BLAST database files. Optional if BLASTDB is already set.",
    )
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--max-target-seqs", type=int, default=10)
    parser.add_argument("--min-identity", type=float, default=90.0)
    parser.add_argument("--min-qcov", type=float, default=70.0)
    parser.add_argument(
        "--taxdump-dir",
        type=Path,
        default=None,
        help="Optional NCBI taxdump directory containing nodes.dmp and names.dmp for full lineage ranks.",
    )
    parser.add_argument(
        "--skip-blast",
        action="store_true",
        help="Reuse an existing raw BLAST TSV in --out-dir instead of running blastn.",
    )
    return parser.parse_args()


def iter_fasta_ids(path: Path) -> list[str]:
    ids: list[str] = []
    with path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if line.startswith(">"):
                ids.append(line[1:].split()[0])
    return ids


def clean_tax(value: object) -> str:
    text = str(value or "").strip()
    return "" if text in {"N/A", "NA", "nan", "None"} else text


def split_species_name(scientific_name: str) -> tuple[str, str]:
    name = clean_tax(scientific_name)
    if not name:
        return "", ""
    pieces = name.split()
    if len(pieces) >= 2 and pieces[0][0].isupper():
        return pieces[0], " ".join(pieces[:2])
    return pieces[0] if pieces else "", name


def taxonomy_string(hit: BlastHit) -> str:
    genus_from_name, species_from_name = split_species_name(hit.ssciname)
    genus = clean_tax(hit.sgenus) or genus_from_name
    levels = [
        ("d", "Bacteria"),
        ("p", hit.sphylum),
        ("c", hit.sclass),
        ("o", hit.sorder),
        ("f", hit.sfamily),
        ("g", genus),
    ]
    taxon = "; ".join(
        f"{prefix}__{clean_tax(value)}" for prefix, value in levels if clean_tax(value)
    )
    if species_from_name and genus and species_from_name.startswith(genus + " "):
        taxon = f"{taxon}; s__{species_from_name}" if taxon else f"s__{species_from_name}"
    return taxon


def load_taxdump(taxdump_dir: Path | None) -> Callable[[str], str | None]:
    if taxdump_dir is None:
        return lambda _taxids: None

    nodes_path = taxdump_dir / "nodes.dmp"
    names_path = taxdump_dir / "names.dmp"
    if not nodes_path.exists() or not names_path.exists():
        raise FileNotFoundError(
            f"taxdump directory must contain nodes.dmp and names.dmp: {taxdump_dir}"
        )

    parents: dict[str, str] = {}
    ranks: dict[str, str] = {}
    names: dict[str, str] = {}

    with nodes_path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            parts = [part.strip() for part in raw.split("|")]
            if len(parts) < 3:
                continue
            taxid, parent, rank = parts[0], parts[1], parts[2]
            parents[taxid] = parent
            ranks[taxid] = rank

    with names_path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            parts = [part.strip() for part in raw.split("|")]
            if len(parts) < 4:
                continue
            taxid, name_txt, _unique_name, name_class = parts[:4]
            if name_class == "scientific name":
                names[taxid] = name_txt

    rank_to_prefix = {
        "superkingdom": "d",
        "kingdom": "d",
        "phylum": "p",
        "class": "c",
        "order": "o",
        "family": "f",
        "genus": "g",
        "species": "s",
    }
    ordered_prefixes = ["d", "p", "c", "o", "f", "g", "s"]

    def lookup(taxids: str) -> str | None:
        first_taxid = str(taxids or "").split(";", 1)[0].strip()
        if not first_taxid or first_taxid in {"N/A", "0"}:
            return None

        found: dict[str, str] = {}
        current = first_taxid
        seen: set[str] = set()
        while current and current not in seen:
            seen.add(current)
            rank = ranks.get(current, "")
            prefix = rank_to_prefix.get(rank)
            name = names.get(current, "")
            if prefix and name and prefix not in found:
                found[prefix] = name
            parent = parents.get(current, "")
            if not parent or parent == current:
                break
            current = parent

        if not found:
            return None
        if "d" not in found and any(prefix in found for prefix in ordered_prefixes[1:]):
            found["d"] = "Bacteria"
        return "; ".join(
            f"{prefix}__{found[prefix]}" for prefix in ordered_prefixes if found.get(prefix)
        )

    return lookup


def blast_confidence(hit: BlastHit | None) -> float:
    if hit is None:
        return 0.0
    return round((hit.pident / 100.0) * (hit.qcovs / 100.0), 6)


def confidence_label(hit: BlastHit | None) -> str:
    if hit is None:
        return "no_hit"
    if hit.pident >= 99.0 and hit.qcovs >= 90 and hit.length >= 250:
        return "high_genus_possible_species"
    if hit.pident >= 97.0 and hit.qcovs >= 85:
        return "medium_genus"
    if hit.pident >= 95.0 and hit.qcovs >= 80:
        return "low_family_or_genus"
    if hit.pident >= 90.0 and hit.qcovs >= 70:
        return "very_low"
    return "unreliable"


def run_blast(args: argparse.Namespace, raw_out: Path) -> None:
    if shutil.which("blastn") is None:
        raise RuntimeError("blastn was not found on PATH.")

    env = os.environ.copy()
    if args.blastdb_dir is not None:
        existing = env.get("BLASTDB", "")
        env["BLASTDB"] = str(args.blastdb_dir) if not existing else f"{args.blastdb_dir}{os.pathsep}{existing}"

    outfmt = (
        "6 qseqid sseqid pident length mismatch gapopen qstart qend qlen qcovs "
        "evalue bitscore staxids ssciname stitle"
    )
    cmd = [
        "blastn",
        "-query",
        str(args.query),
        "-db",
        args.blast_db,
        "-outfmt",
        outfmt,
        "-max_target_seqs",
        str(args.max_target_seqs),
        "-perc_identity",
        str(args.min_identity),
        "-qcov_hsp_perc",
        str(args.min_qcov),
        "-num_threads",
        str(args.threads),
    ]
    with raw_out.open("w", encoding="utf-8", newline="") as handle:
        subprocess.run(cmd, check=True, stdout=handle, env=env)


def parse_blast(raw_out: Path) -> dict[str, list[BlastHit]]:
    standard_cols = [
        "qseqid",
        "sseqid",
        "pident",
        "length",
        "mismatch",
        "gapopen",
        "qstart",
        "qend",
        "qlen",
        "qcovs",
        "evalue",
        "bitscore",
        "staxids",
        "ssciname",
        "stitle",
    ]
    legacy_14_cols = [
        "qseqid",
        "sseqid",
        "pident",
        "length",
        "mismatch",
        "gapopen",
        "qstart",
        "qend",
        "qlen",
        "qcovs",
        "evalue",
        "bitscore",
        "stitle",
        "ssciname",
    ]
    legacy_19_cols = [
        "qseqid",
        "sseqid",
        "pident",
        "length",
        "mismatch",
        "gapopen",
        "qstart",
        "qend",
        "qlen",
        "qcovs",
        "evalue",
        "bitscore",
        "stitle",
        "ssciname",
        "sphylum",
        "sclass",
        "sorder",
        "sfamily",
        "sgenus",
    ]
    hits: dict[str, list[BlastHit]] = {}
    if not raw_out.exists() or raw_out.stat().st_size == 0:
        return hits

    with raw_out.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for values in reader:
            if len(values) == len(standard_cols):
                row = dict(zip(standard_cols, values))
            elif len(values) == len(legacy_14_cols):
                row = dict(zip(legacy_14_cols, values))
                row["staxids"] = ""
            elif len(values) == len(legacy_19_cols):
                row = dict(zip(legacy_19_cols, values))
                row["staxids"] = ""
            else:
                raise ValueError(
                    f"Unexpected BLAST column count {len(values)} in {raw_out}; "
                    f"expected {len(standard_cols)}, {len(legacy_14_cols)}, or {len(legacy_19_cols)}"
                )
            hit = BlastHit(
                qseqid=row["qseqid"],
                sseqid=row["sseqid"],
                pident=float(row["pident"]),
                length=int(row["length"]),
                mismatch=int(row["mismatch"]),
                gapopen=int(row["gapopen"]),
                qstart=int(row["qstart"]),
                qend=int(row["qend"]),
                qlen=int(row["qlen"]),
                qcovs=float(row["qcovs"]),
                evalue=float(row["evalue"]),
                bitscore=float(row["bitscore"]),
                staxids=row.get("staxids", ""),
                ssciname=row["ssciname"],
                stitle=row["stitle"],
                sphylum=row.get("sphylum", ""),
                sclass=row.get("sclass", ""),
                sorder=row.get("sorder", ""),
                sfamily=row.get("sfamily", ""),
                sgenus=row.get("sgenus", ""),
            )
            hits.setdefault(hit.qseqid, []).append(hit)

    for feature_id in hits:
        hits[feature_id].sort(key=lambda item: (-item.bitscore, -item.pident, -item.qcovs, item.evalue))
    return hits


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def build_rows(
    feature_ids: list[str],
    hits_by_feature: dict[str, list[BlastHit]],
    lineage_lookup: Callable[[str], str | None],
):
    taxonomy_rows: list[dict[str, object]] = []
    mapping_rows: list[dict[str, object]] = []
    all_hit_rows: list[dict[str, object]] = []

    for feature_id in feature_ids:
        hits = hits_by_feature.get(feature_id, [])
        best = hits[0] if hits else None
        taxon = ""
        if best is not None:
            taxon = lineage_lookup(best.staxids) or taxonomy_string(best)
        score = blast_confidence(best)
        taxonomy_rows.append(
            {
                "Feature ID": feature_id,
                "Taxon": taxon,
                "Confidence": score,
            }
        )
        mapping_rows.append(
            {
                "Feature ID": feature_id,
                "Taxon": taxon,
                "Confidence": score,
                "Blast_confidence_label": confidence_label(best),
                "Status": "success" if best is not None else "no_hit",
                "Identity_pct": best.pident if best is not None else "",
                "Coverage_pct": best.qcovs if best is not None else "",
                "Alignment_length": best.length if best is not None else "",
                "Evalue": best.evalue if best is not None else "",
                "Bitscore": best.bitscore if best is not None else "",
                "Top_accession": best.sseqid if best is not None else "",
                "Top_taxids": best.staxids if best is not None else "",
                "Scientific_name": clean_tax(best.ssciname) if best is not None else "",
                "Top_hit_description": best.stitle if best is not None else "",
                "Annotation_source": "local_blast_all_asv",
            }
        )

        for rank, hit in enumerate(hits, 1):
            all_hit_rows.append(
                {
                    "Feature ID": feature_id,
                    "Rank": rank,
                    "Accession": hit.sseqid,
                    "Taxids": hit.staxids,
                    "Identity_pct": hit.pident,
                    "Coverage_pct": hit.qcovs,
                    "Alignment_length": hit.length,
                    "Evalue": hit.evalue,
                    "Bitscore": hit.bitscore,
                    "Scientific_name": clean_tax(hit.ssciname),
                    "Taxon": lineage_lookup(hit.staxids) or taxonomy_string(hit),
                    "Top_hit_description": hit.stitle,
                }
            )

    return taxonomy_rows, mapping_rows, all_hit_rows


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    raw_out = args.out_dir / "blast_hits_raw.tsv"

    if not args.query.exists():
        raise FileNotFoundError(f"Query FASTA not found: {args.query}")

    print(f"Query FASTA: {args.query}")
    print(f"Output dir: {args.out_dir}")
    print(f"BLAST DB: {args.blast_db}")
    if args.blastdb_dir is not None:
        print(f"BLASTDB dir: {args.blastdb_dir}")

    if not args.skip_blast:
        run_blast(args, raw_out)
    elif not raw_out.exists():
        raise FileNotFoundError(f"--skip-blast was used, but raw BLAST TSV is missing: {raw_out}")

    feature_ids = iter_fasta_ids(args.query)
    hits_by_feature = parse_blast(raw_out)
    lineage_lookup = load_taxdump(args.taxdump_dir)
    taxonomy_rows, mapping_rows, all_hit_rows = build_rows(
        feature_ids, hits_by_feature, lineage_lookup
    )

    taxonomy_tsv = args.out_dir / "taxonomy.tsv"
    mapping_tsv = args.out_dir / "taxonomy_blast_mapping.tsv"
    top_hits_tsv = args.out_dir / "all_blast_hits_ranked.tsv"

    write_tsv(taxonomy_tsv, taxonomy_rows, ["Feature ID", "Taxon", "Confidence"])
    write_tsv(
        mapping_tsv,
        mapping_rows,
        [
            "Feature ID",
            "Taxon",
            "Confidence",
            "Blast_confidence_label",
            "Status",
            "Identity_pct",
            "Coverage_pct",
            "Alignment_length",
            "Evalue",
            "Bitscore",
            "Top_accession",
            "Top_taxids",
            "Scientific_name",
            "Top_hit_description",
            "Annotation_source",
        ],
    )
    write_tsv(
        top_hits_tsv,
        all_hit_rows,
        [
            "Feature ID",
            "Rank",
            "Accession",
            "Taxids",
            "Identity_pct",
            "Coverage_pct",
            "Alignment_length",
            "Evalue",
            "Bitscore",
            "Scientific_name",
            "Taxon",
            "Top_hit_description",
        ],
    )

    n_success = sum(1 for row in mapping_rows if row["Status"] == "success")
    print(f"Annotated ASVs: {n_success}/{len(feature_ids)}")
    print(f"QIIME2-style taxonomy TSV: {taxonomy_tsv}")
    print(f"BLAST mapping TSV: {mapping_tsv}")
    print(f"Ranked hit TSV: {top_hits_tsv}")


if __name__ == "__main__":
    main()
