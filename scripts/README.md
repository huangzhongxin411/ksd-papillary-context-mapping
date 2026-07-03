# Script Modules

脚本按分析模块编号组织。每个模块目录建议包含：

- `README.md`: 该模块输入、输出、命令和质控标准。
- `run_*.sh`: 可复现的顶层执行脚本。
- `*.R` 或 `*.py`: 具体分析代码。
- `logs/`: 长任务运行日志可统一写入项目根目录 `logs/`。

约定：

- 原始数据只放在 `data/raw/`，不在脚本中覆盖。
- 中间标准化结果放在 `data/processed/`。
- 参考文件放在 `data/references/`。
- 表格输出放在 `results/tables/`。
- 图形输出放在 `results/figures/`。
- 每个模块都输出一个 QC 或 summary table，方便后续 manuscript 追溯。

