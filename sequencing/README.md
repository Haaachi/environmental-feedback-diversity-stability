# 16S sequencing workflow

`sequence_data/pipeline/` contains the current shared ASV implementation for the
temperature and dilution-factor experiments. The latter retains the identifier
`mortality` in scripts and data paths. See the [shared pipeline guide](sequence_data/pipeline/README.md)
for methods and species-grouping thresholds.

## Exported inputs

First restore both exported directories as described in `docs/DATA.md`, then run
from the repository root:

```bash
python sequencing/sequence_data/pipeline/synthetic_community_pipeline.py --experiment all
python sequencing/sequence_data/pipeline/build_unified_standard_annotations.py
python sequencing/sequence_data/pipeline/match_species_to_library.py
python sequencing/sequence_data/pipeline/build_strict_library_correspondence.py
python sequencing/sequence_data/pipeline/build_trait_measurement_merge.py
python sequencing/sequence_data/pipeline/build_pairwise_strict_trait_merge.py
```

Taxonomy annotation first uses restored BLAST taxonomy. New BLAST annotations
require NCBI BLAST+ and a 16S reference database. Annotation scripts provide
`--help` and `--blastdb-dir` options. Correspondence and trait-measurement outputs
use English column names and labels; readers also accept the original labels.

## QIIME2 / HPC

The temperature workflow retains the executed truncation parameters:
forward **228**, reverse **205**. The value 230 in an original comment is not
the executed value. The dilution-factor workflow retains its separate V4 parameters.

```bash
# Run on Linux / HPC with an existing QIIME2 2023.9 environment.
export WORK_DIR=/absolute/path/to/temperature/analysis
export OLD_RAW_DIR=/absolute/path/to/Temperature_V1
export NEW_RAW_DIR=/absolute/path/to/Temperature_V2
export OUT_DIR="$WORK_DIR"
export CONDA_SH=/absolute/path/to/conda/etc/profile.d/conda.sh
export QIIME2_ENV=qiime2-amplicon-2023.9
python sequencing/Tem_scripts/generate_manifest.py
bash sequencing/Tem_scripts/run_temperature_qiime2.sh
```

Before full QIIME2 processing, import paired reads and supply `demux_old.qza`,
`demux_new.qza` and `metadata.tsv`. Manifest generation does not import reads
automatically. Dilution-factor manifest generation uses `BASE_DIR`; its QIIME2
workflow uses `WORK_DIR` and `REF_DIR`.

SLURM scripts retain their original server resource settings. Adapt `#SBATCH`
directives to your cluster before submission. Submit from the repository root
after creating the log directory:

```bash
mkdir -p slurm_logs
sbatch sequencing/Tem_scripts/submit_temperature_qiime2.sh
# Or: sbatch sequencing/Mor_scripts/submit_mortality_qiime2.sh
```

When submitting from another directory, set `PIPELINE_SCRIPT_DIR` to the absolute
path of the relevant Tem_scripts/Mor_scripts directory.
