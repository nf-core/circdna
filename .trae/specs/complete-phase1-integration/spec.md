# Complete Phase 1 Integration Spec

## Why

阶段一架构骨架已搭建完成，但真实工具集成、输出发布（publishDir）、MultiQC 报告、重复序列注释等关键功能缺失，且无法通过 `test_local.config` 端到端测试验证。需要补全这些功能使阶段一真正可用。

## What Changes

- 将 ECCsplorer / CANDIDATE_MERGE / ECC_SCORE 模块的容器镜像从 `python:3.9-slim` 替换为 `quay.io/biocontainers/python:3.9--py39` 格式
- 为所有新模块和子工作流添加 `publishDir` 配置，使结果输出到 `--outdir`
- 修复 `CIRCLE_MAP_PIPELINE` 缺失 `bed` 输出 emit 的 bug
- 修复 Integrated Mode 中 gDNA 与 eccDNA 样本 meta.id 不同导致 channel join 失败的问题
- 添加重复序列注释加载模块（读取已有的 EDTA/RepeatMasker GFF3/BED 文件并传递给 ECC_SCORE）
- 将 mosdepth、eccDNA 候选、ECC_SCORE 统计等新内容接入 MultiQC 报告
- 生成 gDNA 测试数据（使用 wgsim 从同一参考基因组模拟普通基因组 DNA reads）
- 新建 `test_integrated.config` 测试配置，支持 Integrated Mode 端到端测试
- 更新 `test_local.config` 使其支持新 mode 参数

## Impact

- Affected specs: `add-execution-modes`（已完成，无冲突）
- Affected code:
  - `modules/local/eccsplorer/main.nf` — 容器镜像替换 + publishDir
  - `modules/local/candidate_merge/main.nf` — 容器镜像替换 + publishDir
  - `modules/local/ecc_score/main.nf` — 容器镜像替换 + publishDir
  - `modules/nf-core/mosdepth/main.nf` — publishDir
  - `subworkflows/local/circle_map_pipeline/main.nf` — 修复 bed emit
  - `subworkflows/local/reference_mode/main.nf` — 重复序列注释集成
  - `subworkflows/local/eccdna_mode/main.nf` — Circle-Map bed 输出修复
  - `subworkflows/local/integrated_mode/main.nf` — channel join 逻辑修复
  - `workflows/circdna.nf` — MultiQC 集成 + 输出收集
  - `assets/multiqc_config.yml` — 新增 section 配置
  - `conf/test_local.config` — 添加 mode 参数
  - `conf/test_integrated.config` — 新建
  - `testdatasets/` — 新增 gDNA 测试数据
  - `testdatasets/samplesheet/` — 新增带 datatype 字段的 samplesheet

## ADDED Requirements

### Requirement: quay.io Python 容器镜像

所有自定义 Python 模块（ECCsplorer、CANDIDATE_MERGE、ECC_SCORE）SHALL 使用 `quay.io/biocontainers/python` 镜像，格式为 `quay.io/biocontainers/python:<version>`。

#### Scenario: 模块容器配置

- **WHEN** 检查 ECCsplorer / CANDIDATE_MERGE / ECC_SCORE 模块的 container 配置
- **THEN** 容器镜像地址 SHALL 为 `quay.io/biocontainers/python:3.9--py39` 或更高版本

### Requirement: 输出发布（publishDir）

所有新模块和子工作流的输出文件 SHALL 发布到 `params.outdir` 下对应的子目录中。

#### Scenario: Reference Mode 输出

- **WHEN** 以 `--mode reference` 运行流程
- **THEN** mosdepth 的 `.regions.bed.gz`、`.summary.txt`、`.global.dist.txt` 文件 SHALL 发布到 `outdir/reference_mode/mosdepth/`

#### Scenario: eccDNA Mode 输出

- **WHEN** 以 `--mode eccdna` 运行流程
- **THEN** ECCsplorer 候选 BED、Circle-Map 候选 BED、合并候选 BED SHALL 发布到 `outdir/eccdna_mode/` 对应子目录

#### Scenario: Integrated Mode 输出

- **WHEN** 以 `--mode integrated` 运行流程
- **THEN** ECC_SCORE 评分 BED 文件 SHALL 发布到 `outdir/integrated_mode/ecc_score/`

### Requirement: 重复序列注释加载

Reference Mode SHALL 支持加载外部 EDTA/RepeatMasker GFF3/BED 文件，并将其传递给 Integrated Mode 用于 TE repeat penalty 计算。

#### Scenario: 提供重复序列文件

- **WHEN** 用户通过 `--repeat_gff` 参数提供 GFF3/BED 文件
- **THEN** 该文件 SHALL 被加载并传递给 ECC_SCORE 模块用于 TE overlap 计算

#### Scenario: 未提供重复序列文件

- **WHEN** 用户未提供 `--repeat_gff` 参数
- **THEN** 流程 SHALL 跳过 TE penalty 计算，TE overlap ratio 设为 0

### Requirement: gDNA 测试数据生成

SHALL 使用 wgsim 或类似工具从现有参考基因组生成普通基因组 DNA（gDNA）reads，用于 Integrated Mode 测试。

#### Scenario: gDNA 数据生成

- **WHEN** 运行 gDNA 数据生成脚本
- **THEN** SHALL 生成与现有 eccDNA 测试数据相同格式、相同参考基因组的 paired-end FastQ 文件
- **AND** 文件命名 SHALL 为 `gdna_1_R1.fastq.gz` / `gdna_1_R2.fastq.gz`

### Requirement: Integrated Mode 测试配置

SHALL 新建 `test_integrated.config` 配置文件，支持 Integrated Mode 端到端测试。

#### Scenario: 运行 Integrated Mode 测试

- **WHEN** 执行 `nextflow run main.nf -profile test_integrated,docker --outdir <OUTDIR>`
- **THEN** 流程 SHALL 使用包含 gDNA 和 eccDNA 样本的 samplesheet
- **AND** 以 `--mode integrated` 运行完整的 Reference + eccDNA + Integrated 流程

## MODIFIED Requirements

### Requirement: CIRCLE_MAP_PIPELINE 输出

CIRCLE_MAP_PIPELINE 子工作流 SHALL 暴露 `bed` 输出 channel，包含 CIRCLEMAP_REALIGN 和 CIRCLEMAP_REPEATS 的 BED 文件。

#### Scenario: Circle-Map bed 输出

- **WHEN** CIRCLE_MAP_PIPELINE 运行完成
- **THEN** `CIRCLE_MAP_PIPELINE.out.bed` SHALL 返回 `tuple val(meta), path("*.bed")` 格式的 channel

### Requirement: Integrated Mode Channel Join

Integrated Mode SHALL 正确处理 gDNA 和 eccDNA 不同样本之间的 depth 比对，不依赖 meta.id 的直接 join。

#### Scenario: 不同样本 ID 的 depth 匹配

- **WHEN** gDNA 样本 ID 为 `gdna_1`，eccDNA 样本 ID 为 `circdna_1`
- **THEN** Integrated Mode SHALL 使用 eccDNA 候选区域从两个 mosdepth 结果中分别提取深度
- **AND** 不 SHALL 因为 meta.id 不匹配而丢失数据

### Requirement: MultiQC 报告集成

MultiQC 报告 SHALL 包含 mosdepth 质量统计、eccDNA 候选数量统计、ECC_SCORE 分布等新内容。

#### Scenario: MultiQC 包含新内容

- **WHEN** 流程运行完成后查看 MultiQC 报告
- **THEN** 报告 SHALL 包含 mosdepth summary 统计 section
- **AND** 报告 SHALL 包含 eccDNA 候选数量统计 section（当运行 eccDNA/integrated mode 时）

## REMOVED Requirements

无
