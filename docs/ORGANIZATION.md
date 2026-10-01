# Organization record (2026-10-01)

| Original code directory | Repository directory | Treatment |
| --- | --- | --- |
| `code/workplace_final` | `analysis/communities` | Retained Python/R scripts and inputs; excluded generated outputs; replaced hardcoded paths |
| `code/within cycle` | `analysis/isolates` | Retained production plots and Summary input; excluded temporary _debug/_check/_verify scripts |
| `code/sequence pipeline` | `sequencing` | Retained shared ASV pipeline, QIIME2 scripts and reference inputs; excluded older standalone species-identification scripts |
| `code/model_upload_clean` | `model` | Used the existing clean model copy; retained main and sensitivity configurations, plots, SLURM scripts and trajectory experiments |

Historical scripts, duplicate HPC packages and output archives in `code/model`
were excluded. The model copy comes from the previously organized
`model_upload_clean`; its source mapping is recorded in `source_manifest.csv`.

Changes made during organization:

- Community Python scripts locate inputs from their script directory or `COMMUNITY_WORKSPACE`; R scripts use the Rscript `--file` argument or that environment variable.
- Isolate scripts locate Excel inputs and outputs using `ISOLATE_WORKSPACE` or their script directory.
- Sequencing scripts expose server paths, QIIME2 environments and thread counts through environment variables.
- The model's setuptools build backend was corrected to support standard installation.
- Dependency files, an analysis runner, Git ignore rules, input restoration notes and validation records were added.

SHA-256 values in the source manifest describe files before organization and
identify their origin. Path edits and later English translations change file
hashes. Formulas, numerical parameters, experimental windows, species-merging
rules, statistical thresholds and plot colors retain their source definitions.

The English translation covers documentation, comments, docstrings, console
messages, generated sequencing table labels and workbook text. Numeric data,
species identifiers and sequences are preserved. Sequencing readers accept both
the English labels and the original workbook/table labels.
