# Environmental feedback, diversity and stability

Research code for environmental feedback, diversity and stability in synthetic
microbial communities under temperature and dilution-factor perturbations.
The repository includes community analyses, within-cycle isolate experiments,
16S sequencing workflows and a consumer-resource model with pH feedback.
The dilution-factor experiment retains the legacy identifier `mortality` in code and data paths.

## Repository structure

```text
analysis/communities/          Community abundance, diversity, instability, variance decomposition and R plots
  data/                       Unified experimental Excel files and species annotations
analysis/isolates/             Isolate OD, CFU, pH, growth and death analyses
  Summary.xlsx                Within-cycle experimental measurements
sequencing/Tem_scripts/        Temperature QIIME2 workflows and manifest generation
sequencing/Mor_scripts/        Dilution-factor QIIME2 workflows and manifest generation
sequencing/sequence_data/      Shared ASV pipeline, reference sequences and taxonomy
model/                        Model implementation, configurations, plots and SLURM scripts
scripts/                      Community analysis runner and R dependency installation
docs/                         Data notes, organization record and validation report
```

## Installation

Python 3.11 or 3.12 is recommended. The model requires Python >= 3.10.

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

Exact package versions used for local validation are listed in
`requirements-validated.txt`. Install from that file to use the same versions.

For R plots, install R and run:

```bash
Rscript scripts/install_r_dependencies.R
```

## Community analyses

Run from the repository root:

```bash
python scripts/run_community_analysis.py
# Also generate the core R plots after completing the calculations:
python scripts/run_community_analysis.py --plots
```

Outputs are written to `analysis/communities/processed/` and
`analysis/communities/figures/`. Individual Python and R scripts can also be run
separately. Input paths default to the script directory. To use another workspace,
set `COMMUNITY_WORKSPACE` to a directory containing a `data/` subdirectory.

The default workflow uses uncorrected abundance from `process_abundance.py`.
The optional `process_abundance_copynumber.py` performs 16S copy-number correction
and writes to the same processed directory. Use a separate `COMMUNITY_WORKSPACE`
to avoid overwriting the default analysis.

## Isolate experiments

```bash
python analysis/isolates/plot_growth.py
python analysis/isolates/plot_growth_rate.py
python analysis/isolates/plot_death_rate_analysis.py
python analysis/isolates/plot_24h_endpoint.py
```

Outputs are written to `analysis/isolates/figures/`. Additional figures can be
generated with the other `plot_*.py` scripts in that directory.
Set `ISOLATE_WORKSPACE` to use an external workspace containing `Summary.xlsx`.

## Model

```bash
cd model
python scripts/smoke_test.py
```

The main model configuration is
`configs/phase_r_h_gnormal_sd3_R1_sparse6_tri_S12_h015_e5000.yaml`:
12 species, 24 resources, a 100 x 100 parameter grid and 5,000 random communities
per grid point. The full scan requires substantial computing resources; the smoke
test uses a 2 x 2 grid. See the [model overview](model/README.md),
[reproduction guide](model/REPRODUCE.md) and [model specification](model/MODEL_SPEC.md).
SLURM scripts default to the `crph-model` environment; override this with
`CRMODEL_ENV` and `CONDA_SH`.

## Sequencing

The shared pipeline uses feature tables, representative sequences and taxonomy
exported from QIIME2. Raw FASTQ files, QIIME2 artifacts, BLAST databases and large
model outputs are excluded from Git. See the [sequencing guide](sequencing/README.md)
and [data notes](docs/DATA.md) for input locations and execution order.

## Organization and validation

The original `code/` directory is preserved; this repository is a separate
organized copy. Organization preserved formulas, thresholds and experimental
design while updating paths, entry points, documentation and packaging.
See the [organization record](docs/ORGANIZATION.md),
[source manifest](docs/source_manifest.csv) and [validation report](docs/VALIDATION.md).

No open-source license has been assigned. The project owner should choose a
license before public release or granting reuse rights.
