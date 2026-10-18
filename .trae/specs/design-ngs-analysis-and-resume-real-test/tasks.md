# Tasks

- [ ] Task 1: 盘点并固化 `ECCsplorer` 已迁移 / 未迁移边界

  - [ ] SubTask 1.1: 核对 `map`、FASTQ/BAM、gDNA control、数据库参数、输出映射的现状
  - [ ] SubTask 1.2: 明确 `PRExer`、`clu`、`all/comparative` 与 cluster handoff 仍属未完成迁移范围
  - [ ] SubTask 1.3: 将该边界写入实施说明、任务拆分和验收范围

- [ ] Task 2: 固化 `circdna.nf -> eccdna.smk` 的二代测序交接接口

  - [ ] SubTask 1.1: 为 `eccdna` 模式梳理并实现 `handoff.tsv` 与 `samples.auto.yaml` 的字段契约
  - [ ] SubTask 1.2: 在 `circdna.nf` 中补齐交接导出模块与脚本入口
  - [ ] SubTask 1.3: 校验 `group`、`data_type`、control 配对与缺失场景的导出行为

- [ ] Task 3: 完成 `eccdna.smk` 二代测序下游分析模块骨架

  - [ ] SubTask 2.1: 将 `Snakefile` 拆分为 `score.smk`、`standardize.smk`、`distribution.smk`、`deg.smk`、`visualize.smk`
  - [ ] SubTask 2.2: 将 `config/config.yaml` 扩展为支持评分、标准化、注释、差异分析和可视化参数
  - [ ] SubTask 2.3: 将 `config/samples.yaml` 调整为模板/示例，实际运行优先消费 `samples.auto.yaml`

- [ ] Task 4: 迁移二代测序分析脚本并明确模块边界

  - [ ] SubTask 3.1: 从 `ecc_pipe-main/analysis_code/Distribution.py` 提炼标准化与单样本分布分析脚本
  - [ ] SubTask 3.2: 从 `DEG.py` 提炼 burden matrix 与差异分析入口
  - [ ] SubTask 3.3: 迁移并参数化 `deseq2.R`、`edger.R`、`limma.R`、`clusterprofile.R`、`circlize.R`
  - [ ] SubTask 3.4: 明确标准 analysis BED 的列定义、输入输出路径与脚本 CLI 参数

- [ ] Task 5: 更新 `circdna.nf` 的测试与文档契约

  - [ ] SubTask 4.1: 更新 `conf/test_local.config` 及相关测试输入，使其适配最新 NGS 主链
  - [ ] SubTask 4.2: 确认整合指南中的 NGS 模块设计与实际代码入口一致
  - [ ] SubTask 4.3: 记录 `test_local` 真实测试的推荐命令、输出路径与 `-resume` 复测策略

- [ ] Task 6: 执行 `test_local` 二代测序真实测试并建立缓存

  - [ ] SubTask 5.1: 首轮执行 `nextflow run main.nf -profile test_local,docker --mode eccdna` 真实运行
  - [ ] SubTask 5.2: 记录任务数、耗时、关键输出目录与失败信息（如有）
  - [ ] SubTask 5.3: 验证 `eccsplorer`、`circlemap`、`mosdepth`、交接导出模块的真实产出

- [ ] Task 7: 使用 `-resume` 进行快速复测并验证缓存收益

  - [ ] SubTask 6.1: 在不改变上游重计算输入的前提下执行 `nextflow run main.nf -profile test_local,docker --mode eccdna -resume`
  - [ ] SubTask 6.2: 记录样本级与参考级任务的 `CACHED` 情况
  - [ ] SubTask 6.3: 验证交接文件与下游模块接口在 `-resume` 条件下仍然正确

- [ ] Task 8: 验证 `eccdna.smk` 对自动交接文件的消费能力
  - [ ] SubTask 7.1: 使用 `samples.auto.yaml` 执行 `snakemake -n`
  - [ ] SubTask 7.2: 验证 `score -> standardize -> distribution` 的规则依赖与路径解析
  - [ ] SubTask 7.3: 如真实实现已覆盖，补充 `deg` 与 `visualize` 的最小验证

# Task Dependencies

- Task 2 依赖 Task 1（先锁定 ECCsplorer 已迁移/未迁移范围，再定义交接接口）
- Task 3 依赖 Task 2（先有交接契约，再拆下游规则）
- Task 4 依赖 Task 3（模块骨架先确定，再迁移脚本）
- Task 5 可与 Task 3、Task 4 部分并行，但最终需以已落地接口为准
- Task 6 依赖 Task 2 与 Task 5（真实测试依赖稳定交接接口与测试配置）
- Task 7 依赖 Task 6（必须先有首轮真实运行缓存）
- Task 8 依赖 Task 2、Task 3，且最好在 Task 6 或 Task 7 之后进行

# 并行化说明

- `circdna.nf` 侧的交接接口设计（Task 1）与 `eccdna.smk` 下游规则骨架（Task 2）可由不同子代理并行推进
- 脚本迁移（Task 3）与文档/测试配置同步（Task 4）可并行推进
- 真实测试（Task 5）与 `-resume` 复测（Task 6）必须串行
- `snakemake -n` 验证（Task 7）可在首轮真实运行结束后尽早执行
