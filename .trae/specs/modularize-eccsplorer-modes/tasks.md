# Tasks

- [x] Task 1: 盘点并冻结四种模式的真实职责边界

  - [x] SubTask 1.1: 逐一核对 `PRExer`、`map`、`clu`、`all/comparative` 的原始输入、输出、依赖与阶段顺序
  - [x] SubTask 1.2: 标记哪些步骤应保留原始工具语义，哪些步骤应由标准 Nextflow 模块替代
  - [x] SubTask 1.3: 形成统一的模式到模块映射表

- [x] Task 2: 给出 Nextflow / Snakemake 归属分析

  - [x] SubTask 2.1: 识别适合迁移到 Nextflow 的步骤（重计算、稳定 IO、强缓存收益、资源差异明显）
  - [x] SubTask 2.2: 识别适合保留或迁移到 Snakemake 的步骤（参数敏感、统计/筛选/可视化、探索性强）
  - [x] SubTask 2.3: 为 `PRExer`、`map`、`clu`、`all/comparative` 输出明确的归属结论

- [x] Task 3: 识别“已有模块优先迁参数”的机会点

  - [x] SubTask 3.1: 盘点当前已存在且与 ECCsplorer 步骤边界相近的模块
  - [x] SubTask 3.2: 为这些模块列出可直接迁入的 ECCsplorer 参数、默认值和开关
  - [x] SubTask 3.3: 定义迁参数前后需要对比的结果差异与性能指标

- [x] Task 4: 设计共享基础模块层

  - [x] SubTask 2.1: 设计 reads 预处理、格式转换、结果提取、BLAST 注释等共享模块
  - [x] SubTask 2.2: 定义 BAM/FASTQ 统一输入接口
  - [x] SubTask 2.3: 为共享模块规划独立容器或环境边界

- [x] Task 5: 设计 `PRExer` 的替代模块链

  - [x] SubTask 3.1: 将原始 `PRExer` 需求映射到标准 Nextflow 预处理模块
  - [x] SubTask 3.2: 定义参数映射、输入输出和缓存边界
  - [x] SubTask 3.3: 明确哪些步骤不再依赖原始 ECCsplorer 总入口

- [x] Task 6: 设计 `map` 模式的最小职责模块拆分

  - [x] SubTask 4.1: 定义输入标准化层
  - [x] SubTask 4.2: 定义核心 mapping 执行层
  - [x] SubTask 4.3: 定义候选提取、junction 提取、序列提取和统计提取层
  - [x] SubTask 4.4: 定义可选 BLAST 注释层

- [x] Task 7: 设计 `clu` 模式的标准子工作流

  - [x] SubTask 5.1: 定义 cluster 输入准备模块
  - [x] SubTask 5.2: 定义核心聚类模块
  - [x] SubTask 5.3: 定义 cluster 注释模块
  - [x] SubTask 5.4: 定义 cluster 候选 BED 提取与绘图模块

- [x] Task 8: 设计 `all/comparative` 的 DAG 编排层

  - [x] SubTask 6.1: 定义 `all` 如何组合 `PRExer + map + clu`
  - [x] SubTask 6.2: 定义 comparative 的独立输入输出契约
  - [x] SubTask 6.3: 定义 comparative 的缓存、重跑和并行策略

- [x] Task 9: 设计交接与下游消费接口

  - [x] SubTask 7.1: 扩展 `handoff.tsv` 字段以支持 `clu` / comparative 输出
  - [x] SubTask 7.2: 扩展 `samples.auto.yaml` 的可选字段
  - [x] SubTask 7.3: 约束 `eccdna.smk` 如何识别和消费这些扩展字段

- [x] Task 10: 设计性能与环境优化方案

  - [x] SubTask 8.1: 为轻量模块与重型模块分别规划资源标签
  - [x] SubTask 8.2: 为不同模块规划最小容器/环境边界
  - [x] SubTask 8.3: 设计 `-resume` 下的缓存收益验证策略

- [x] Task 11: 设计实现顺序与首期验收范围
  - [x] SubTask 9.1: 确定首期必须落地的模块
  - [x] SubTask 9.2: 确定可后置的模块和重型依赖
  - [x] SubTask 9.3: 明确 `test_local` / `-resume` 的最小可行验收方案

# Task Dependencies

- Task 2 依赖 Task 1（先冻结模式职责，再做引擎归属分析）
- Task 3 依赖 Task 2（先做归属分析，再评估已有模块是否优先迁参数）
- Task 4 依赖 Task 3（共享层设计依赖已有模块复用策略）
- Task 5、Task 6、Task 7 可在 Task 4 之后并行推进
- Task 8 依赖 Task 5、Task 6、Task 7（组合层依赖各模式边界）
- Task 9 依赖 Task 6、Task 7、Task 8（交接字段依赖实际模块输出）
- Task 10 可与 Task 5、Task 6、Task 7、Task 8 部分并行
- Task 11 依赖前面所有设计结论

# 并行化说明

- `PRExer`、`map`、`clu` 三条模式设计线可并行
- 共享模块层与环境优化策略可并行推进
- `all/comparative` 必须在三条模式线收敛后再定稿

# Implementation Tasks

- [ ] Task 12: 拆分现有 `map` 黑盒模块为标准模块链

  - [ ] SubTask 12.1: 新增 `eccsplorer_input_normalize` 模块，统一 BAM / FASTQ 入口
  - [ ] SubTask 12.2: 新增 `eccsplorer_map_core` 模块，负责执行现有 ECCsplorer map 核心检测
  - [ ] SubTask 12.3: 新增 `eccsplorer_map_extract` 模块，负责导出 candidates/junction/alignment/ecc_sequences
  - [ ] SubTask 12.4: 新增 `eccsplorer_map_pipeline` 子工作流，并在 `eccdna_mode` 中接线替代旧单体调用

- [ ] Task 13: 落地参数优先迁移策略

  - [ ] SubTask 13.1: 在 `nextflow.config` 与 `modules.config` 中集中定义 ECCsplorer 参数
  - [ ] SubTask 13.2: 为已有模块补齐 `ECCsplorer` 参数语义与默认值
  - [ ] SubTask 13.3: 保持现有输出契约兼容，避免影响 `eccdna_mode` 与下游接口

- [x] Task 14: 增加 `clu` / `all-comparative` 的标准骨架

  - [x] SubTask 14.1: 新增 `eccsplorer_clu_prepare`、`eccsplorer_clu_core`、`eccsplorer_clu_candidates_plot` 模块骨架（实际实现：clu_prepare / REPEATEXPLORER2 / clu_candidates，见 test-tidehunter-module spec）
  - [x] SubTask 14.2: 新增 `subworkflows/local/eccsplorer_cluster/main.nf`（实际实现：eccsplorer_clu_slim 子工作流）
  - [ ] SubTask 14.3: 新增 `subworkflows/local/eccsplorer_all/main.nf`，定义 `PRExer + map + clu + comparative` 的 DAG 骨架
  - [x] SubTask 14.4: 为 `clu` / `all` 增加参数与默认关闭的执行开关（circle_identifier: eccsplorer_clu_slim）

- [x] Task 15: 补充测试与最小验证
  - [x] SubTask 15.1: 更新或新增最小模块测试/工作流验证配置（tests/eccsplorer_clu/test.nf）
  - [x] SubTask 15.2: 至少完成一次静态校验与一次最小运行路径验证（clu 独立测试 + 主流水线回归）
  - [x] SubTask 15.3: 记录当前实现与后续未完成范围，确保首期范围清晰
