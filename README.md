# Environmental feedback, diversity and stability

温度与死亡率扰动下，合成微生物群落的环境反馈、多样性与稳定性研究代码。
包含群落数据分析、单菌株周期内实验分析、16S 测序处理和带 pH 反馈的消费者—资源模型。

## 项目目录

```text
analysis/communities/          群落丰度、多样性、波动、方差分解与 R 绘图
  data/                       统一实验 Excel 与物种注释输入
analysis/isolates/             单菌株 OD、CFU、pH、生长和死亡率分析
  Summary.xlsx                周期内实验输入
sequencing/Tem_scripts/        温度实验 QIIME2 与 manifest 脚本
sequencing/Mor_scripts/        死亡率实验 QIIME2 与 manifest 脚本
sequencing/sequence_data/      共享 ASV 流程、菌株参考序列与分类表
model/                        独立模型程序、配置、绘图和 SLURM 提交脚本
scripts/                      群落分析入口与 R 依赖安装
docs/                         数据说明、整理记录与验证报告
```

## 安装

推荐 Python 3.11 或 3.12。模型要求 Python >= 3.10。

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

本地验证采用的具体版本记录在 `requirements-validated.txt`，需要一致版本时可用它安装。

需要生成 R 图时，另外安装 R，并执行：

```bash
Rscript scripts/install_r_dependencies.R
```

## 群落实验分析

在仓库根目录执行：

```bash
python scripts/run_community_analysis.py
# 计算结束后同时生成核心 R 图：
python scripts/run_community_analysis.py --plots
```

结果写入 `analysis/communities/processed/` 和 `analysis/communities/figures/`。
各 Python/R 脚本均可单独执行；默认输入路径相对于脚本定位。
需要使用其他数据目录时，设置 `COMMUNITY_WORKSPACE`，其中需有 `data/` 子目录。

默认使用 `process_abundance.py` 的原始丰度口径。
`process_abundance_copynumber.py` 是可选的 16S 拷贝数校正分析，会写入同一 processed 目录；
建议设置单独的 `COMMUNITY_WORKSPACE`，避免覆盖默认分析。

## 单菌株实验

```bash
python analysis/isolates/plot_growth.py
python analysis/isolates/plot_growth_rate.py
python analysis/isolates/plot_death_rate_analysis.py
python analysis/isolates/plot_24h_endpoint.py
```

结果写入 `analysis/isolates/figures/`。其他图可运行同目录下相应的 `plot_*.py`。
`ISOLATE_WORKSPACE` 可指定含有 `Summary.xlsx` 的外部工作目录。

## 模型

```bash
cd model
python scripts/smoke_test.py
```

正式模型配置为
`configs/phase_r_h_gnormal_sd3_R1_sparse6_tri_S12_h015_e5000.yaml`：
12 个物种、24 种资源、100 × 100 参数网格、每点 5000 个随机群落。
完整扫描需要大量计算资源；快速检查只运行 2 × 2 网格。
详见 [模型说明](model/README.md)、[复现步骤](model/REPRODUCE.md) 和 [模型定义](model/MODEL_SPEC.md)。
模型 SLURM 提交脚本默认环境名为 `crph-model`，可用 `CRMODEL_ENV` 与 `CONDA_SH` 覆盖。

## 测序

共享流程使用 QIIME2 导出的 feature table、代表序列与 taxonomy 文件。
测序原始 FASTQ、QIIME2 文件、BLAST 数据库和模型大体积结果未纳入 Git。
恢复位置和运行顺序见 [测序说明](sequencing/README.md) 与 [数据说明](docs/DATA.md)。

## 整理与验证

原始 `code/` 目录保留；本项目是单独的整理副本。
整理时保留计算公式、阈值和实验设计，只调整路径、入口、文档和打包配置。
详见 [整理记录](docs/ORGANIZATION.md)、[源文件清单](docs/source_manifest.csv) 和 [验证报告](docs/VALIDATION.md)。

本仓库未指定开放源代码许可；如需公开发布或授权复用，请由项目所有者确定许可。
