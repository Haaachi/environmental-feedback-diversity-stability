# 数据与计算口径

## 已纳入的输入

`analysis/communities/data/` 包含 Temperature/Mortality 的统一群落、OD、pH Excel，
`Summary.xlsx` 以及统一物种注释 Excel/TSV。
`analysis/isolates/Summary.xlsx` 是独立单菌株周期内测量表。
两处 Summary 文件分别保留原始来源，未假定它们完全相同。

`sequencing/sequence_data/All_Strains_Combined.fasta` 是菌株参考序列；
`Library_taxonomy_groups_v3.xlsx` 用于参考菌株与实验物种的严格映射和性状合并。

## 未纳入的输入与结果

- 原始 FASTQ、QIIME2 QZA/QZV、完整测序导出表、BLAST 数据库。
- 模型 HDF5/NPZ 大规模结果，恢复位置见 `model/data/runs/datasets.yaml`。
- 自动生成的 processed、figures、缓存、历史模型重复包和论文文稿。

从原课题目录恢复测序输入时，对应关系如下：

| 原目录（相对于课题根目录） | 仓库恢复位置 |
| --- | --- |
| `code/sequence pipeline/sequence_data/16s_project/exported/` | `sequencing/sequence_data/16s_project/exported/` |
| `code/sequence pipeline/sequence_data/Mortality_analysis/exported/` | `sequencing/sequence_data/Mortality_analysis/exported/` |

共享 ASV 流程会重新生成物种识别结果，但不会自动替换仓库中作为实验分析基准的统一 Excel。
若要替换，请比较分组与命名变化，更新全部匹配的注释输入后再运行群落分析。

## 可选脚本的数据缺口

`rvstemperature.R` 还需要 `analysis/communities/data/rK_expfit_skip0.xlsx`；
`plot_species_temperature.R` 还需要 `analysis/communities/data/Speciestemperature.xlsx`。
原始 data 目录没有这两个文件，因此它们不属于默认复现入口。

16S 拷贝数校正脚本优先读 `Unified_species_annotations_with_16S_copy_number.xlsx`；
该专用文件目前缺失，脚本可能回退至已有物种注释。
运行校正分析前应检查实际拷贝数信息，不能将回退结果视为已完整校正。

## 实验与模型的口径不同

实验脚本的 relative abundance 单位为百分数，presence 常使用 1%；
实验 biomass collapse 阈值通常为 0.05。
模型使用 B readout，collapse 阈值为 1e-3，fluctuation 分类使用 CV > 0.1。
实验窗口、presence 判定和波动阈值按各脚本原始定义保留，不能直接套用模型阈值。
测序物种判定排除的低 OD 样本，与最终群落统计保留的 collapse 群落也有不同用途。
