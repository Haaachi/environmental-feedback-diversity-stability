# Validation report (2026-10-01)

Environment: macOS arm64, Python 3.12. Exact Python package versions are listed in
`requirements-validated.txt`.

## Checks passed during organization

- Syntax checks for 53 Python files and `bash -n` checks for all shell scripts.
- All six steps of the default community analysis completed, producing abundance, diversity, instability, species decomposition, mechanism and window robustness CSV files.
- The default workflow produced 48 CSV files; 26 matched archived outputs within rtol=1e-9 and atol=1e-12.
- Nine community Python scripts had identical ASTs to their source files after excluding BASE_DIR assignments and os imports; calculation code was preserved.
- The model smoke test passed on a 2 x 2 grid with two random communities per point, including simulation, completion checks and strict post-processing.
- Standard model wheel build and local installation passed.
- Isolate `plot_growth_rate.py` read Summary.xlsx, completed statistical tests and generated PDF/PNG outputs.
- With the original exported inputs restored, the shared ASV workflow processed both experiments, retaining 25 temperature species and 38 dilution-factor species.

## Archived outputs versus current source

Twenty-two regenerated CSV files differed from archived outputs; the complete
list is in `regression_comparison.txt`. Some archived files came from earlier
classification schemes, while current source uses community CV for classification.
For example, older mechanism_metrics files contain sum_abs_std_threshold and
temporal_bc_threshold, whereas current source produces community_cv and
community_cv_threshold. The older species-decomposition field
fluctuating_by_temporal_bc was replaced by the current community_cv_threshold.

Shared continuous numeric columns were compared in five representative tables:
mechanism metrics, last3 community decomposition, per-community window CV,
replicate CV and replicate Bray-Curtis. All matched within the stated tolerances.
The window robustness fluctuation_class labels and class-specific sample sizes
and means change with the classification definition. These are differences
between current source and archived outputs; organization did not change the
classification rules or overwrite new results with older ones.

## Checks not performed

- Rscript is unavailable on this computer; R plots and systematic R statistical tests have not been executed.
- QIIME2 denoising, external BLAST database searches, SLURM submission and the full large model scan have not been executed locally.
- Other isolate plotting scripts received syntax and path checks but were not run individually.
- Two additional R plotting scripts lack original input workbooks; see `DATA.md`.

Generated validation outputs and the local virtual environment are excluded by
`.gitignore` and are not uploaded with the repository.

## English translation checks

- All 54 Python files passed syntax checks. The original 53 files retained the same computational AST structure after excluding documentation, presentation strings and the new input-label normalization wrapper.
- All 26 R files retained identical expressions after excluding comments and text literals. R execution remains unverified because Rscript is unavailable.
- The six-step community analysis was rerun; all 48 CSV outputs matched the pre-translation baseline within rtol=1e-12 and atol=1e-12.
- The three edited workbooks retained numeric values, formulas, cell styles, merged regions, column widths, hidden columns, freeze panes and filters. Text edits covered 761 cells and two taxonomy worksheet names. Changed views were rendered and inspected.
- Sequencing regression checks retained 63 strict correspondence rows, 36 candidate-based trait-measurement units and 52 pairwise-strict units, with matching identifiers, candidate sets, counts and grouping labels after translation. Both original and English taxonomy inputs were accepted; both downstream scripts accepted the English strict table.
- Isolate `plot_growth_rate.py` successfully read the translated Summary workbook and regenerated its statistics and PDF/PNG figures using the existing project environment.
- A repository text scan found no remaining Chinese prose or corrupted private-use characters. Escaped legacy labels remain only for compatibility with original inputs; worksheet font/theme metadata is preserved.

The language update preserves calculations and numerical thresholds, including
the existing 0.25/0.265 discrepancy described in `DATA.md`.

## Manuscript PR reconciliation (2026-10-02)

The manuscript phase-map plotting script now calculates PR from unthresholded
species B over the final three passages, using sample standard deviations and
division by S for the normalized map. Both reorganized and legacy phase2 input
formats follow this definition. Optional last-ten postprocessing diagnostics
remain separate from the manuscript PR.

Three regression tests passed (`python -m unittest discover -s model/tests -v`):
changing the preceding seven passages leaves manuscript PR unchanged; species
variation below the diversity activity threshold contributes to PR with matching
results across both storage formats; and zero temporal variation has undefined
participation rather than an invented finite value. Synthetic-input comparisons
also confirmed that non-PR phase-map outputs are unchanged by this correction.

The full production ensemble and manuscript PR figures have not been regenerated
locally. Rerun the main phase-map command in `model/REPRODUCE.md` from the raw
ensemble files to update the figure and NPZ outputs. The experimental CV display
threshold discrepancy and isolate death-kinetics timing discrepancy identified
during manuscript review are not resolved by this PR correction.
