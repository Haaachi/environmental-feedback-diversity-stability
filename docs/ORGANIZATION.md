# 整理记录（2026-10-01）

| 原始代码目录 | 整理后目录 | 处理 |
| --- | --- | --- |
| `code/workplace_final` | `analysis/communities` | 保留 Python/R 与输入；移除生成结果；替换硬编码路径 |
| `code/within cycle` | `analysis/isolates` | 保留正式绘图与 Summary；不纳入临时 _debug/_check/_verify 脚本 |
| `code/sequence pipeline` | `sequencing` | 保留共享 ASV 流程、QIIME2 脚本和参考输入；不纳入旧的独立物种识别脚本 |
| `code/model_upload_clean` | `model` | 使用已有干净模型版本；保留主配置、敏感性配置、绘图、SLURM 和轨迹实验 |

不纳入 `code/model` 中的历史脚本、重复 HPC 包和输出档案。
模型发布副本来自已经整理过的 `model_upload_clean`，其来源对应清单在 `source_manifest.csv`。

修改：

- Python 群落脚本使用脚本目录或 COMMUNITY_WORKSPACE；R 使用 Rscript 的 --file 或该变量。
- 单菌株脚本以 ISOLATE_WORKSPACE 或脚本目录定位 Excel 和结果。
- 测序脚本暴露服务器路径、QIIME2 环境和线程数的环境变量。
- 修正模型 pyproject 的 setuptools build backend，以支持标准安装。
- 添加依赖文件、运行入口、Git 忽略规则、数据恢复说明和运行验证。

源文件清单中的 SHA-256 为整理前文件的哈希，便于定位来源；经过路径调整的文件哈希会变化。
计算公式、数值参数、实验窗口、物种归并规则、统计阈值和绘图颜色保持原脚本定义。
