# circDNA 增量/减量缓存实现 - 产品需求文档

## Overview

- **摘要**: 修复 circDNA pipeline 在修改 samplesheet（增加或减少样本）时无法利用 Nextflow 缓存的问题。当前修改 CSV 文件会导致所有下游任务缓存失效，需要重新运行所有样本的计算。
- **目的**: 实现增量缓存（增加样本时旧样本复用缓存）和减量缓存（减少样本时保留样本复用缓存），大幅提升 pipeline 迭代效率。
- **目标用户**: bioinformaticians 使用 circDNA pipeline 分析 eccDNA 数据

## Goals

- 实现增量缓存：增加样本时，已存在样本的所有处理步骤复用缓存
- 实现减量缓存：减少样本时，保留样本的所有处理步骤复用缓存
- 保持现有 pipeline 功能完全兼容
- 不改变用户接口和 samplesheet 格式
- 验证 short-read (FASTQ) 和 long-read (ONT/PacBio) 两种模式

## Non-Goals (Out of Scope)

- 不修改 SAMPLESHEET_CHECK 验证逻辑
- 不改变任何 Process 的输入/输出定义
- 不修改 workflow 内部的 channel 操作逻辑（groupTuple、branch 等保持不变）
- 不优化非缓存相关的性能问题

## Background & Context

- 当前问题链路：
  1. `INPUT_CHECK` 接收 CSV 文件
  2. 调用 `SAMPLESHEET_CHECK` 进行验证，输出验证后的 CSV
  3. 使用 `SAMPLESHEET_CHECK.out.csv.splitCsv()` 创建 channel
  4. Channel 元素依赖验证后的 CSV 文件 → CSV 一变，所有哈希失效
- 解决方案（方案3）：
  1. 保留 `SAMPLESHEET_CHECK` 作为独立验证步骤（输出 versions.yml）
  2. Channel 创建改为直接从原始 CSV 解析：`file(samplesheet).splitCsv()`
  3. 每个样本用 `file()` 包装 fastq 文件，哈希只依赖于 fastq 文件内容

## Functional Requirements

- **FR-1**: INPUT_CHECK 的 `reads` channel 必须直接从原始 CSV 创建，不依赖 SAMPLESHEET_CHECK 输出
- **FR-2**: SAMPLESHEET_CHECK 仍需执行以提供版本信息和验证
- **FR-3**: 支持 FASTQ (short-read) 和 long-read (ONT/PacBio) 两种输入模式
- **FR-4**: 保持 `create_fastq_channels`、`create_bam_channels`、`create_long_read_channels` 函数签名不变
- **FR-5**: 保持 `reads` channel 的输出格式不变：`[val(meta), [reads]]`

## Non-Functional Requirements

- **NFR-1**: 代码改动最小化，仅修改 `input_check/main.nf` 中 channel 创建逻辑
- **NFR-2**: 保持与现有 test_local.config 的完全兼容
- **NFR-3**: 验证步骤必须在 channel 创建之前完成（确保文件存在性检查）

## Constraints

- **技术**: Nextflow DSL2, Groovy
- **依赖**: 必须与现有 `workflows/circdna.nf` 中的 channel 消费者兼容
- **运行环境**: 本地 test_local 配置（3个测试样本）

## Assumptions

- 假设原始 CSV 文件路径与验证后 CSV 使用相同路径
- 假设 `file(samplesheet).splitCsv()` 在文件不存在时会报错（可以替代 SAMPLESHEET_CHECK 的部分验证功能）
- 假设 test_local.config 中的测试数据有效且路径正确

## Acceptance Criteria

### AC-1: 增量缓存 - 增加样本

- **Given**: 使用 test_local.config 的 3 个样本完成首次运行
- **When**: 在 samplesheet 中增加第 4 个样本，使用 -resume 重新运行
- **Then**: 前 3 个样本的所有处理步骤显示 CACHED，仅第 4 个样本显示 NEW
- **Verification**: `programmatic`
- **Notes**: 通过 `nextflow log` 检查任务状态

### AC-2: 减量缓存 - 减少样本

- **Given**: 使用包含 3 个样本的 samplesheet 完成首次运行
- **When**: 从 samplesheet 中移除 1 个样本（保留 2 个），使用 -resume 重新运行
- **Then**: 保留的 2 个样本的所有处理步骤显示 CACHED
- **Verification**: `programmatic`

### AC-3: 同时增减缓存

- **Given**: 使用 3 个样本完成首次运行
- **When**: 增加 1 个样本并移除 1 个样本（保留 2 个），使用 -resume 重新运行
- **Then**: 保留的 2 个样本显示 CACHED，新增的 1 个样本显示 NEW
- **Verification**: `programmatic`

### AC-4: Pipeline 功能完整性

- **Given**: 修改后的代码
- **When**: 使用 test_local.config 正常运行 pipeline
- **Then**: Pipeline 完成且所有处理结果与修改前一致
- **Verification**: `programmatic`

### AC-5: Short-read 和 Long-read 兼容性

- **Given**: 修改后的代码
- **When**: 分别使用 FASTQ (short-read) 和 long-read 模式运行
- **Then**: 两种模式下 pipeline 均正常工作
- **Verification**: `programmatic`

### AC-6: 代码质量

- **Given**: 修改后的 `input_check/main.nf`
- **When**: 代码审查
- **Then**: 改动仅限 channel 创建逻辑，不引入冗余代码，保持原有代码风格
- **Verification**: `human-judgment`

## Open Questions

- [ ] `file(samplesheet).splitCsv()` 是否能完全替代 SAMPLESHEET_CHECK 的文件存在性验证？
- [ ] 是否需要保留 `SAMPLESHEET_CHECK.out.csv` 的某些验证输出？
