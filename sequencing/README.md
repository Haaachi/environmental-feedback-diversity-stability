# 16S sequencing workflow

`sequence_data/pipeline/` 是温度与死亡率共享的当前 ASV 处理实现。
详细方法与物种分组阈值见 [共享流程说明](sequence_data/pipeline/README.md)。

## 已导出的输入

先按 `docs/DATA.md` 恢复两个 exported 目录，再从仓库根目录运行：

```bash
python sequencing/sequence_data/pipeline/synthetic_community_pipeline.py --experiment all
python sequencing/sequence_data/pipeline/build_unified_standard_annotations.py
python sequencing/sequence_data/pipeline/match_species_to_library.py
python sequencing/sequence_data/pipeline/build_strict_library_correspondence.py
python sequencing/sequence_data/pipeline/build_trait_measurement_merge.py
python sequencing/sequence_data/pipeline/build_pairwise_strict_trait_merge.py
```

taxonomy 注释优先使用已恢复的 BLAST taxonomy；新 BLAST 注释需要另外安装 NCBI BLAST+
和 16S 参考数据库。各注释脚本提供 `--help` 与 `--blastdb-dir` 参数。

## QIIME2 / HPC

温度流程保留实际执行参数：forward truncation **228**、reverse **205**；
脚本原注释中的 230 不是执行值。死亡率流程保留独立的 V4 参数。

```bash
# 在已有 QIIME2 2023.9 环境的 Linux / HPC 上运行
export WORK_DIR=/absolute/path/to/temperature/analysis
export OLD_RAW_DIR=/absolute/path/to/Temperature_V1
export NEW_RAW_DIR=/absolute/path/to/Temperature_V2
export OUT_DIR="$WORK_DIR"
export CONDA_SH=/absolute/path/to/conda/etc/profile.d/conda.sh
export QIIME2_ENV=qiime2-amplicon-2023.9
python sequencing/Tem_scripts/generate_manifest.py
bash sequencing/Tem_scripts/run_temperature_qiime2.sh
```

运行完整 QIIME2 处理前需先导入成对 reads，提供 `demux_old.qza`、`demux_new.qza` 与
`metadata.tsv`。manifest 生成不会自动执行 reads 导入。
死亡率 manifest 使用 `BASE_DIR`；QIIME2 运行使用 `WORK_DIR` 与 `REF_DIR`。
SLURM 脚本保留原服务器的资源配置，提交前按所在集群调整 `#SBATCH`。
从仓库根目录提交，先创建日志目录：

```bash
mkdir -p slurm_logs
sbatch sequencing/Tem_scripts/submit_temperature_qiime2.sh
# 或 sbatch sequencing/Mor_scripts/submit_mortality_qiime2.sh
```

若从其他目录提交，设置 `PIPELINE_SCRIPT_DIR` 为对应 Tem_scripts/Mor_scripts 的绝对路径。
