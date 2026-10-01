# Standard Synthetic-Community Pipeline

This directory is the canonical pipeline for the temperature and mortality
synthetic-community post-processing workflows.

## Shared ASV Pipeline

Run both experiments:

```powershell
python .\pipeline\synthetic_community_pipeline.py --experiment all
```

Run one experiment:

```powershell
python .\pipeline\synthetic_community_pipeline.py --experiment temperature
python .\pipeline\synthetic_community_pipeline.py --experiment mortality
```

Try a different TopN for one experiment:

```powershell
python .\pipeline\synthetic_community_pipeline.py --experiment temperature --top-n 35 --output-dir-name species_identification_top35
```

Standard output directories:

- `16s_project/species_identification_standard`
- `Mortality_analysis/species_identification_standard`

Standard file names use the same suffixes for both experiments:

- `{Experiment}_Unified_final.xlsx`
- `{Experiment}_asv_mapping.tsv`
- `{Experiment}_species_table.tsv`
- `{Experiment}_species.fasta`
- `{Experiment}_community_species.tsv`
- `{Experiment}_species_count.tsv`
- `{Experiment}_retention_qc.tsv`

The Excel workbook sheet order is identical for both experiments:

`W1 | W2 | W3 | W4 | W5 | Community_Species | Species_Members | Species_Count | Retention_QC | ASV_Mapping_Top100 | README`

## Preserved Experiment Differences

These differences are intentional and live only in `EXPERIMENTS` inside
`synthetic_community_pipeline.py`.

| Experiment | Root | Sample columns | Analysis samples | TopN |
| --- | --- | --- | --- | --- |
| Temperature | `16s_project` | `temperature_*` | W1-W4 plus non-collapsed W5 C3/C7/C12 samples; other W5 communities and OD-collapsed C7/C12 samples excluded from model/presence decisions | 45 |
| Mortality | `Mortality_analysis` | `mortality_*` | `mortality_*_R1_D1-D6` | 45 |

Shared thresholds:

- Pearson clustering: `r > 0.95`
- Presence rule: relative abundance `> 0.01` in at least 3 samples, or maximum relative abundance `>= 0.05`
- Sequence-assisted merge: same informative genus, sequence identity `>= 0.995`,
  Pearson `>= 0.90`, presence Jaccard `>= 0.85`, and minor/major total-read
  ratio `<= 0.30`
- High-correlation profile merge is also gated by taxonomy/sequence or a
  satellite-like abundance ratio, so correlation alone no longer forces a merge.

Taxon naming in the standard Excel and support tables is experiment-local:
retained ASV-species groups are labeled as `Genus sp.1`, `Genus sp.2`, etc.
The numbering is separate for temperature and mortality. Distinct retained
ASV-species groups remain distinct units even if BLAST/QIIME gives the same
scientific species or strain name.

Each standard output directory also includes `{Experiment}_asv_pair_diagnostics.tsv`
and an `ASV_Pair_Diagnostics` sheet. These tables record sequence identity,
correlation, presence Jaccard, abundance ratio, and merge/review recommendations
for all topN ASV pairs.

Temperature uses TopN 45 as the standard interpretation point. W5 communities C3, C7, and C12 are included in TopN ranking, ASV merging, presence calling, and QC only at OD-supported sample points because OD shows these communities persisted there and their sequencing profiles are meaningful. Other W5 communities are retained in exported count sheets for inspection but excluded from species-identification decisions because they collapsed and their relative-abundance profiles are noise-dominated. Within C7/C12, W5 C7 R1/R2 D7, C7 R3 D5-D7, C12 R1 D6-D8, and C12 R3 D6/D8 are also excluded from presence decisions because OD is nearly background.

## Standard BLAST Annotation

### All-ASV Taxonomy For Temperature

For temperature, annotate all ASV representative sequences first. This makes the
temperature workflow match mortality, where ASV taxonomy is available before
the synthetic-community species grouping step.

Run this on Linux from the `16s_project` directory:

```bash
python /path/to/pipeline/annotate_all_asvs_local_blast.py \
  --query exported/dna-sequences.fasta \
  --out-dir exported/local_blast_taxonomy \
  --blast-db 16S_ribosomal_RNA \
  --blastdb-dir /path/to/16S_db \
  --threads 16
```

The main output is:

- `exported/local_blast_taxonomy/taxonomy.tsv`

It has the same three columns as the mortality QIIME2 taxonomy export:

`Feature ID | Taxon | Confidence`

After copying that directory back, rerun the shared ASV pipeline. Temperature
will read `exported/local_blast_taxonomy/taxonomy.tsv` automatically.

### Retained-Species BLAST Annotation

Run after the shared ASV pipeline has produced `{Experiment}_species.fasta`.

```powershell
python .\pipeline\annotate_species_local_blast.py --experiment all --blast-db 16S_ribosomal_RNA --blastdb-dir C:\Users\ADMIN\blast_db\16S_db
```

The current machine did not have `blastn` on PATH when this pipeline was
organized, so BLAST annotation was not executed here. Once BLAST+ and the local
database are available, the command above will create matching annotation files
for both experiments under each experiment root:

- `species_annotations/local_blast_standard/{Experiment}_species_taxonomy_mapping.tsv`
- `species_annotations/local_blast_standard/{Experiment}_taxonomy_style.tsv`
- `species_annotations/local_blast_standard/{Experiment}_all_blast_hits_ranked.tsv`
- `species_annotations/local_blast_standard/{Experiment}_species_taxonomy.xlsx`

## Unified Annotation Workbook

Build the cross-experiment workbook:

```powershell
python .\pipeline\build_unified_standard_annotations.py
```

Output:

- `species_annotations/unified_standard/Unified_species_annotation_mapping.tsv`
- `species_annotations/unified_standard/Unified_community_species_annotations.tsv`
- `species_annotations/unified_standard/Unified_species_annotations.xlsx`

Before BLAST files exist, this script falls back to the standard species tables:
temperature taxonomy is blank/unclassified, and mortality keeps its QIIME
taxonomy fallback. After BLAST runs, rerun this script and it will use the
standard BLAST mapping files automatically.

