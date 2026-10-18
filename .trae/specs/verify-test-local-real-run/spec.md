# Verify test_local Real Run All Modes Spec

## Why

前一轮工作仅完成了 stub 模式验证（52 tasks succeeded），但 stub 模式不执行真实工具，仅产出占位文件。用户需要确认 `test_local` 配置在**真实模式**（非 stub）下能否正常执行并生成真实结果，且需要覆盖**所有模式**（`reference` 和 `eccdna`），保证新功能（ECCSPLORER 双模式支持、BAM_PREPROCESSING bug 修复等）能正常生成结果。

`test_local.config` 默认 `mode = 'eccdna'`，需通过 `--mode` 参数切换测试 `reference` 模式。两种模式覆盖不同分析链路：

- **reference 模式**：FASTQ 预处理 → BWA 比对 → BAM 排序/索引/MarkDuplicates → CIRCEXPLORER2_PARSE → Circle-Finder → AmpliconArchitect → UNICYCLER → MULTIQC
- **eccdna 模式**：FASTQ 预处理 → BWA 比对 → BAM 排序/索引/MarkDuplicates → mosdepth → ECCSPLORER → Circle-Map realign → MULTIQC

## What Changes

- **运行 test_local 真实模式 - eccdna 模式**：执行 `nextflow run main.nf -profile test_local,docker --outdir /tmp/circdna_test_local_eccdna`（默认 mode=eccdna）
- **运行 test_local 真实模式 - reference 模式**：执行 `nextflow run main.nf -profile test_local,docker --mode reference --outdir /tmp/circdna_test_local_reference`
- **验证两种模式所有关键模块产出真实结果**（非占位文件）
- **重点验证新功能**：
  - ECCSPLORER 真实执行（ENTRYPOINT 修复后）
  - BAM_PREPROCESSING 子工作流正常工作
  - MULTIQC 汇总报告生成
- **记录执行时间、任务数、失败任务**（如有）
- **若失败，记录错误信息并定位问题**

## Impact

- Affected specs: `fix-eccsplorer-entrypoint-conflict`（前置 spec，本 spec 验证其修复后真实模式可用性）、`build-eccsplorer-from-source-and-bam-support`（验证 ECCSPLORER 双模式支持）
- Affected code: 无代码变更，仅运行验证
- Affected files: 产出位于 `/tmp/circdna_test_local_eccdna/` 和 `/tmp/circdna_test_local_reference/`

## ADDED Requirements

### Requirement: test_local 真实模式端到端验证 - eccdna 模式

系统 SHALL 在 `test_local,docker` profile 下（非 stub 模式）成功运行 circdna.nf 流程的 `eccdna` 模式，所有关键模块产出真实结果文件。

#### Scenario: eccdna 模式真实运行成功

- **WHEN** 用户执行 `nextflow run main.nf -profile test_local,docker --outdir /tmp/circdna_test_local_eccdna`
- **THEN** 流程成功完成（exit code 0）
- **AND** 所有任务成功执行（无 failed tasks）
- **AND** 产出真实结果文件（非占位文件）

#### Scenario: ECCSPLORER 产出真实检测结果

- **WHEN** ECCSPLORER process 执行完成
- **THEN** 产出 `*_candidates.bed` 文件（可能为空，但不是 stub 的 0 字节占位）
- **AND** 产出 `*_junction_reads.txt` 文件
- **AND** `versions.yml` 包含真实版本号

### Requirement: test_local 真实模式端到端验证 - reference 模式

系统 SHALL 在 `test_local,docker` profile 下（非 stub 模式）成功运行 circdna.nf 流程的 `reference` 模式，所有关键模块产出真实结果文件。

#### Scenario: reference 模式真实运行成功

- **WHEN** 用户执行 `nextflow run main.nf -profile test_local,docker --mode reference --outdir /tmp/circdna_test_local_reference`
- **THEN** 流程成功完成（exit code 0）
- **AND** 所有任务成功执行（无 failed tasks）
- **AND** 产出真实结果文件（非占位文件）

#### Scenario: MULTIQC 生成汇总报告

- **WHEN** 流程完成
- **THEN** `<outdir>/multiqc/` 目录包含 `multiqc_report.html`
- **AND** 报告包含各模块的统计信息

## MODIFIED Requirements

无（本 spec 仅验证，不修改现有需求）
