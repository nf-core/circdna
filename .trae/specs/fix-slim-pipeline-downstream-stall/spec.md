# 修复 Slim 流水线下游停滞 + BAM 输出精简 Spec

## Why

1. eccsplorer_slim 流水线的 HAARZ、CANDIDATE_EXTRACT、COVERAGE_PROFILE、NORMALIZE、VISUALIZE、HTML_REPORT 6 个进程未运行
2. `results/slim_run/eccsplorer_slim/bam/` 目录中 `.sam` 与 `.bam` 并存浪费空间，且 circdna_2 样本文件缺失

## What Changes

- SEGEMEHL_ALIGN 模块补充 `--splits --briefcigar --MEOP --accuracy 95` 参数（与 ECCsplorer 源码 config.py SEGEMEHL_CMD 一致），使 segemehl 生成 `.sngl.bed` 供 HAARZ 使用
- HAARZ 模块增加 awk 行修复（复刻 eccMapper.py run_splitread_detect 的 start>end 纠错逻辑）
- SEGEMEHL_ALIGN publishDir 排除 `.sam` 文件（保留 `.sngl.bed`/`.mult.bed` 供下游），bam 目录只保留 `.bam`
- 不清理缓存：SEGEMEHL_ALIGN 的 `ext.args` 变更使缓存自然失效，3 个样本自动重跑并重新 publish，circdna_2 缺失随之补齐

## Impact

- Affected specs: eccsplorer_slim 原子化流水线
- Affected code:
  - `circdna.nf/conf/modules.config`（SEGEMEHL_ALIGN ext.args + publishDir）
  - `circdna.nf/modules/local/haarz/main.nf`（awk 修复）
  - 重新验证：Docker 全流程运行

## ADDED Requirements

### Requirement: SEGEMEHL_ALIGN 生成 split-read BED

系统 SHALL 以 ECCsplorer 源码同等的 segemehl 参数运行比对，确保生成 `.sngl.bed` 文件。

#### Scenario: 成功案例

- **WHEN** 运行 eccsplorer_map_slim 模式
- **THEN** SEGEMEHL_ALIGN 输出 `.sngl.bed`，HAARZ 及下游 6 个进程全部执行并产出文件

## MODIFIED Requirements

### Requirement: 输出目录精简

SEGEMEHL_ALIGN publishDir SHALL 排除 `.sam` 文件；`eccsplorer_slim/bam/` 目录 SHALL 仅包含 `.bam`。

#### Scenario: 成功案例

- **WHEN** 流水线完成后检查 `results/slim_run/eccsplorer_slim/bam/`
- **THEN** 每个样本仅存在一个 `.bam` 文件，无 `.sam` 文件
