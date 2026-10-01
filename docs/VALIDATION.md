# 验证报告（2026-10-01）

环境：macOS arm64，Python 3.12；具体 Python 包版本记录在 `requirements-validated.txt`。

## 已通过

- 53 个 Python 文件的语法检查，以及所有 Shell 脚本的 `bash -n` 检查。
- 群落实验默认分析入口：六个步骤均完成，生成丰度、多样性、波动、物种分解、机制指标和窗口稳健性 CSV。
- 默认群落分析产出 48 个 CSV；其中 26 个与原目录已有结果一致（rtol=1e-9、atol=1e-12）。
- 9 个群落 Python 脚本在忽略 BASE_DIR 路径赋值和 os 导入后，AST 与原脚本完全一致；计算代码未改变。
- 模型小规模测试：2 × 2 网格，每点 2 个随机群落，模拟、完成状态检查与 strict 后处理均通过。
- 模型标准 wheel 构建与本地安装通过。
- 单菌株 `plot_growth_rate.py` 成功读取 Summary.xlsx、完成统计检验并输出 PDF/PNG。
- 共享测序 ASV 流程在恢复原始 exported 输入后成功处理两组实验：温度保留 25 个物种，死亡率保留 38 个物种。

## 旧输出与当前代码的区别

22 个重新生成的 CSV 与原目录已有文件有差异，完整文件清单见 `regression_comparison.txt`。
原目录部分文件是早期分类方案的输出，当前源码已经使用 community CV 进行分类。
例如 mechanism_metrics 旧文件包含 sum_abs_std_threshold 和 temporal_bc_threshold，
新文件按当前源码包含 community_cv 与 community_cv_threshold；物种分解的旧字段
fluctuating_by_temporal_bc 被当前源码中的 community_cv_threshold 替代。

对机制指标、last3 群落分解、窗口逐群落 CV、replicate CV 与 replicate Bray-Curtis
五个代表性表比较共同连续数值列，均在上述容差内一致。
窗口稳健性表的 fluctuation_class 及按类别汇总的 n、均值会随分类定义而改变。
这是当前源码与旧输出的版本差异；整理过程没有改动分类规则或用旧结果覆盖新结果。

## 未验证的部分

- 当前电脑没有 Rscript，R 图和 R 系统统计检验未实际执行。
- QIIME2 denoising、外部 BLAST 数据库查询、SLURM 提交与完整大规模模型扫描未在本机执行。
- 其他单菌株绘图脚本做了语法及路径检查，未逐一执行。
- 两个额外 R 绘图脚本缺少原始输入 Excel，详情见 `DATA.md`。

验证生成的输出和本机虚拟环境被 .gitignore 排除，不随仓库上传。
