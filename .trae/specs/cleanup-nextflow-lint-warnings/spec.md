# 清理 Nextflow Lint 警告 Spec

## Why

运行 `nextflow lint .` 输出 8 条警告（4 个文件），涉及未使用变量、未使用闭包参数、已废弃 API。虽然不影响运行，但污染 lint 输出、降低代码整洁度，且部分警告（如 `Channel` 废弃）在未来 Nextflow 版本可能升级为错误。

## What Changes

共 6 处最小化修改，消除全部 8 条警告：

1. **`subworkflows/local/circle_map_pipeline/main.nf:46`** — 删除未使用变量 `ch_qname_sorted_bai`（声明后从未被引用，包括 emit 和下游 join）
2. **`subworkflows/local/eccdna_mode/main.nf:39`** — 闭包参数 `meta` → `_meta`（`.map { meta, fai -> fai }` 中 `meta` 未使用）
3. **`subworkflows/local/eccdna_mode/main.nf:50`** — 闭包参数 `meta` → `_meta`（`.map { meta, fasta -> fasta }` 中 `meta` 未使用）
4. **`subworkflows/local/input_check/main.nf:12,16,20`** — `Channel.fromPath` → `channel.fromPath`（3 处，`Channel` 大写形式已废弃，Nextflow 推荐 `channel` 小写）
5. **`subworkflows/local/reference_mode/main.nf:14`** — take 参数 `repeat_gff` → `_repeat_gff`（声明但未在 main 中使用，保留参数位用于未来 repeat annotation 功能，前缀 `_` 仅抑制 lint 警告，不影响位置参数调用）
6. **`subworkflows/local/reference_mode/main.nf:31`** — 闭包参数 `meta` → `_meta`（`.map { meta, fai -> fai }` 中 `meta` 未使用）

## Impact

- Affected specs: 无（独立 lint 清理）
- Affected code:
  - `subworkflows/local/circle_map_pipeline/main.nf` — 删除 1 行
  - `subworkflows/local/eccdna_mode/main.nf` — 修改 2 行闭包参数名
  - `subworkflows/local/input_check/main.nf` — 3 处 `Channel` → `channel`
  - `subworkflows/local/reference_mode/main.nf` — 修改 take 参数名 + 1 行闭包参数名
  - `CHANGELOG.md` — 新增 PATCH 条目
  - `nextflow.config` — 版本号 PATCH bump

## Root Cause Analysis

| #   | 警告类型                       | 根因                                                                                                                           | 修复方式          |
| --- | ------------------------------ | ------------------------------------------------------------------------------------------------------------------------------ | ----------------- |
| 1   | Variable declared but not used | `ch_qname_sorted_bai` 在 `circle_map_pipeline` 中声明后从未引用（join 操作用的是 `bam_sorted_bai` 而非 `ch_qname_sorted_bai`） | 删除该行          |
| 2-3 | Parameter was not used         | `.map { meta, fai -> fai }` / `.map { meta, fasta -> fasta }` 中 `meta` 是解构第一个元素但不需要                               | 前缀 `_`          |
| 4   | `Channel` deprecated           | Nextflow DSL2 推荐 `channel` 小写作为 channel factory 入口                                                                     | 替换为 `channel`  |
| 5   | Parameter was not used         | `repeat_gff` 在 take 中声明但 main 中未引用（repeat annotation 功能未实现）                                                    | 前缀 `_` 保留位置 |
| 6   | Parameter was not used         | 同 #2                                                                                                                          | 前缀 `_`          |

## ADDED Requirements

### Requirement: Lint 零警告

`nextflow lint .` 输出中不得包含上述 8 条警告中的任何一条。

#### Scenario: nextflow lint 无警告

- **WHEN** 在 `circdna.nf` 目录下执行 `nextflow lint .`
- **THEN** 输出中不出现 "Variable was declared but not used"、"Parameter was not used"、"The use of `Channel` to access channel factories is deprecated" 三类警告
- **AND** 最终摘要显示 `✅ 81 files had no errors`（77 + 4 修复后 = 81）

#### Scenario: test_local eccdna 模式运行不受影响

- **WHEN** 执行 `nextflow run main.nf -profile test_local,docker --mode eccdna -resume`
- **THEN** Pipeline 正常完成（`Pipeline completed successfully`）
- **AND** ECCSPLORER 仍处理 `3 of 3` 样本
- **AND** 不出现新的 WARN 或 ERROR

## MODIFIED Requirements

### Requirement: circle_map_pipeline 变量声明

`CIRCLE_MAP_PIPELINE` subworkflow 中不得声明未使用的变量。`ch_qname_sorted_bai` 因从未被引用（join 操作使用 `bam_sorted_bai` 而非 `ch_qname_sorted_bai`），应删除其声明行。

### Requirement: 闭包未使用参数命名

所有 `.map { ... }` 闭包中未使用的参数 SHALL 以 `_` 前缀命名，以符合 Groovy/Nextflow 惯例并抑制 lint 警告。

### Requirement: channel factory 使用小写形式

`INPUT_CHECK` subworkflow 中 SHALL 使用 `channel.fromPath`（小写）替代已废弃的 `Channel.fromPath`（大写）。

### Requirement: reference_mode take 参数命名

`REFERENCE_MODE` subworkflow 的 `take` 参数 `repeat_gff` 在功能实现前 SHALL 命名为 `_repeat_gff` 以抑制 lint 警告，同时保留位置参数供未来 repeat annotation 功能使用。
