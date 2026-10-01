# Data and analysis conventions

## Included inputs

`analysis/communities/data/` contains unified Temperature/Mortality community,
OD and pH workbooks, `Summary.xlsx`, and unified species annotations in Excel/TSV
format. `analysis/isolates/Summary.xlsx` contains separate within-cycle isolate
measurements. Both Summary files retain their respective sources; they are not
assumed to be identical.

`sequencing/sequence_data/All_Strains_Combined.fasta` contains reference isolate
sequences. `Library_taxonomy_groups_v3.xlsx` supports strict mapping between
reference isolates and experimental species and merging for trait measurement.

Workbook notes and taxonomy labels have been translated into English without
changing numeric measurements, identifiers, sequences or workbook formatting.
The sequencing readers also accept the original labels when restoring older inputs.

## Excluded inputs and outputs

- Raw FASTQ files, QIIME2 QZA/QZV artifacts, full sequencing export tables and BLAST databases.
- Large model HDF5/NPZ outputs; locations are recorded in `model/data/runs/datasets.yaml`.
- Generated processed tables, figures, caches, historical duplicate model packages and manuscript files.

Restore exported sequencing inputs from the original project directory as follows:

| Original directory, relative to the project root | Repository destination |
| --- | --- |
| `code/sequence pipeline/sequence_data/16s_project/exported/` | `sequencing/sequence_data/16s_project/exported/` |
| `code/sequence pipeline/sequence_data/Mortality_analysis/exported/` | `sequencing/sequence_data/Mortality_analysis/exported/` |

The shared ASV pipeline regenerates species-identification outputs but does not
automatically replace the unified Excel workbooks used as experimental analysis
inputs. Before replacing them, compare grouping and naming changes and update
all matching annotation inputs.

## Missing inputs for optional scripts

`rvstemperature.R` requires `analysis/communities/data/rK_expfit_skip0.xlsx`.
`plot_species_temperature.R` requires `analysis/communities/data/Speciestemperature.xlsx`.
Neither file was present in the original data directory, so these scripts are
excluded from the default reproduction entry point.

The 16S copy-number correction script first looks for
`Unified_species_annotations_with_16S_copy_number.xlsx`. This dedicated file is
currently missing, and the script may fall back to existing species annotations.
Check the available copy-number information before running correction; a fallback
result should not be treated as fully copy-number corrected.

## Experimental and model conventions

Experimental relative abundance is expressed as a percentage; presence often
uses a 1% threshold. The experimental biomass collapse threshold is generally
0.05. The model uses the B readout, a collapse threshold of 1e-3 and CV > 0.1 for
fluctuation classification. Experimental windows, presence rules and instability
thresholds retain each script's source definitions; model thresholds should not
be applied directly to experimental data.

The manuscript treats experimental instability primarily as a continuous
quantity. Fig. 1f uses a global ranking split into high- and low-instability
halves. The CV = 0.25 and accumulated Shannon = 0.8 guides describe regions in
the experimental diversity-instability plane; they do not establish discrete
dynamical classes. Some archived scripts still use CV = 0.265 for binary displays,
while others use 0.25. This translation preserves those values; the discrepancy
requires a separate methods reconciliation.

Low-OD samples excluded during sequencing species identification and collapsed
communities retained in final condition-level summaries serve different purposes.
