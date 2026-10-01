"""Normalize legacy workbook and generated-table labels to English.

Only presentation labels are translated. Numeric values, identifiers and
scientific matching/merging rules are preserved.
"""
from __future__ import annotations

import pandas as pd

# Escaped legacy labels support original inputs without non-English source text.
LEGACY_LABELS = {
    "\u83cc\u682a\u7f16\u53f7": "Strain_ID",
    "\u7c7b\u578b": "Group_type",
    "BLAST\u539f\u59cb\u6ce8\u91ca": "Original_BLAST_annotation",
    "\u5b9e\u9a8cSpecies\u7f16\u53f7": "Experimental_species_ID",
    "\u5b9e\u9a8c\u6ce8\u91ca": "Experimental_annotation",
    "\u5bf9\u5e94\u5e93\u7f16\u53f7": "Library_IDs",
    "\u5e93\u5168\u957f\u6ce8\u91ca": "Library_full_length_annotation",
    "\u6ce8\u91ca\u4e00\u81f4\u6027": "Annotation_consistency",
    "\u5019\u9009\u6570": "Candidate_count",
    "\u5b9e\u9a8c": "Experiment",
    "\u5907\u6ce8": "Notes",
    "species\u4e00\u81f4": "species_agreement",
    "genus\u4e00\u81f4": "genus_agreement",
    "family\u4e00\u81f4": "family_agreement",
    "\u6ce8\u91ca\u51b2\u7a81": "annotation_conflict",
    "\u7f3a\u5c11\u5e93\u5168\u957f\u6ce8\u91ca": "missing_library_full_length_annotation",
    "\u5df2\u6392\u9664\u7f3a\u5c11\u5e93\u5168\u957f\u6ce8\u91ca\u7684\u5e76\u5217\u5019\u9009: ": "Excluded tied candidates lacking full-length library annotations: ",
    "\u6e29\u5ea6species": "Temperature_species",
    "\u6b7b\u4ea1\u7387species": "Mortality_species",
    "\u63a8\u8350\u6d4b\u91cf\u5e93\u7f16\u53f7": "Recommended_library_IDs",
    "\u4f18\u5148\u5e93\u7f16\u53f7": "Preferred_library_ID",
    "\u539f\u59cbspecies\u6570": "Original_species_count",
    "\u5408\u5e76\u540e\u5b9e\u9a8c\u6ce8\u91ca": "Merged_experimental_annotation",
    "\u5e93Group_ID": "Library_Group_ID",
    "\u5408\u5e76\u4f9d\u636e": "Merge_basis",
    "\u5408\u5e76\u7f6e\u4fe1": "Merge_confidence",
    "\u5019\u9009\u5e93\u7f16\u53f7\u6570": "Candidate_library_ID_count",
    "\u539f\u59cbspecies\u7d22\u5f15": "Original_species_index",
    "\u5b9e\u9a8c\u8986\u76d6": "Experiment_coverage",
    "\u6ce8\u91ca\u51b2\u7a81\uff0c\u4fdd\u5b88\u72ec\u7acb": "Annotation conflict; conservatively retained separately",
    "\u4ec5family\u7ea7\u4e00\u81f4\uff0c\u4fdd\u5b88\u72ec\u7acb": "Family-level agreement only; conservatively retained separately",
    "\u5355\u5b9e\u9a8c/\u5355species\u72ec\u7acb\u4fdd\u7559": "Single experiment/species retained separately",
    "\u5019\u9009\u5e93\u7f16\u53f7\u5b8c\u5168\u4e00\u81f4": "Identical candidate library ID sets",
    "\u5019\u9009\u5e93\u7f16\u53f7\u9ad8\u91cd\u53e0\uff0c\u53d6\u5171\u540c\u4ea4\u96c6": "Highly overlapping candidate IDs; use shared intersection",
    "\u5019\u9009\u5e93\u7f16\u53f7\u9ad8\u91cd\u53e0\uff0c\u4f46\u5171\u540c\u4ea4\u96c6\u4e3a\u7a7a\uff0c\u53d6\u5e76\u96c6": "Highly overlapping candidate IDs with empty shared intersection; use union",
    "\u9700\u590d\u6838": "review_needed",
    "\u72ec\u7acb": "independent",
    "\u9ad8": "high",
    "\u4e2d": "medium",
    "\u4f4e": "low",
    "\u4e24\u4e2a\u5b9e\u9a8c\u5171\u6709": "both_experiments",
    "\u4ec5\u6e29\u5ea6\u5b9e\u9a8c": "temperature_only",
    "\u4ec5\u6b7b\u4ea1\u7387\u5b9e\u9a8c": "mortality_only",
    "\u5185{count}\u4e2aspecies": " contains {count} species",
    "\u540c\u4e00\u6d4b\u91cf\u5355\u5143\u5305\u542b": "The same measurement unit contains ",
    "\uff0c\u539fspecies\u7f16\u53f7\u5df2\u4fdd\u7559": "; original species IDs are retained",
    "\u4e0d\u5efa\u8bae\u81ea\u52a8\u5408\u5e76\u5230\u5176\u4ed6\u5355\u5143": "Automatic merging with other units is not recommended",
    "\u5408\u5e76\u72b6\u6001": "Merge_status",
    "\u6e29\u5ea6\u6ce8\u91ca": "Temperature_annotation",
    "\u6b7b\u4ea1\u7387\u6ce8\u91ca": "Mortality_annotation",
    "\u7247\u6bb5\u6bd4\u5bf9\u4f9d\u636e": "Fragment_alignment_basis",
    "\u72ec\u7acb\u4fdd\u7559-\u77ed\u7247\u6bb5\u975e\u552f\u4e00": "retained_separately_nonunique_fragment",
    "\u72ec\u7acb\u4fdd\u7559-\u65e0\u552f\u4e00exact": "retained_separately_no_unique_exact_match",
    "\u72ec\u7acb\u4fdd\u7559-\u65e0\u8de8\u5b9e\u9a8c\u5339\u914d": "retained_separately_no_cross_experiment_match",
    "\u81ea\u52a8\u5408\u5e76-\u552f\u4e00\u4e92\u60e0exact\u7247\u6bb5": "auto_merged_unique_reciprocal_exact_fragment",
    "\u540c\u65f6exact\u5339\u914d": "exactly matches multiple",
    "\u6700\u4f73\u8de8\u5b9e\u9a8c\u7247\u6bb5\u5339\u914d:": "Best cross-experiment fragment match:",
    "\u672a\u81ea\u52a8\u5408\u5e76": "Not automatically merged",
    "mortality V4a \u5bf9 temperature V4V5 \u5bf9\u5e94\u7247\u6bb5 100% exact; ": "mortality V4a matches the corresponding temperature V4V5 fragment exactly (100%); ",
    "\u4e14\u4e24\u4fa7\u5747\u552f\u4e00": "unique on both sides",
    "\uff1b": "; ",
    "\u3001": ", "
}


def translate_text(value: object) -> object:
    """Translate known legacy labels and notes, retaining other values."""
    if not isinstance(value, str):
        return value
    # Exact labels avoid replacing ordinary characters in unrelated fields.
    if value in LEGACY_LABELS:
        return LEGACY_LABELS[value]
    for original, english in sorted(LEGACY_LABELS.items(), key=lambda item: -len(item[0])):
        # Single-character confidence labels are translated only by exact match.
        if len(original) > 1 or original in ("\uFF1B", "\u3001"):
            value = value.replace(original, english)
    return value


def normalize_frame(frame: pd.DataFrame) -> pd.DataFrame:
    """Accept either English or legacy labels at pipeline input boundaries."""
    aliases = {**LEGACY_LABELS, "Type": "Group_type"}
    result = frame.rename(columns=aliases).copy()
    for column in result.select_dtypes(include=["object", "string"]).columns:
        result[column] = result[column].map(translate_text)
    return result
