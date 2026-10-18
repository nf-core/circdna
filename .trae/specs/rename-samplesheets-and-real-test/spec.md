# Rename Samplesheets and Real Test Spec

## Why

当前测试 samplesheets 命名不统一（`samplesheet_*` 前缀），无法直观区分数据来源（本地/在线）和用途。同时 stub 验证已通过但尚未进行真实模式测试，需要验证 ECCsplorer BLAST 注释和 gDNA control 在真实运行中的效果。

## What Changes

- 重命名 5 个测试 samplesheets 文件，统一为 `test_<local|online>_<tags>.csv` 命名规范
- 更新所有引用这些 samplesheets 的配置文件和脚本
- 执行 `test_local` 和 `test_local_gdna` 两个 profile 的真实模式测试（非 stub）
- 验证 ECCsplorer 真实产出（`*_candidates.bed` 非空、`*_blast.m6` 非空）
- 对比有无 gDNA control 对 eccDNA 候选数量的影响

## Impact

- Affected specs: `add-eccsplorer-database-and-control`（stub 验证已完成，真实测试为 Task 21 的延伸）
- Affected code:
  - `circdna.nf/samplesheets/` — 5 个文件重命名
  - `circdna.nf/conf/test_local.config` — 引用路径更新
  - `circdna.nf/conf/test_local_gdna.config` — 引用路径更新
  - `circdna.nf/scripts/test_incremental_cache.py` — 引用路径更新
  - `circdna.nf/AGENTS.md` — 文档引用更新
  - `circdna.nf/testdatasets/README.md` — 文档引用更新

## ADDED Requirements

### Requirement: Samplesheet 命名规范

所有测试用 samplesheets SHALL 遵循 `test_<local|online>_<tags>.csv` 命名规范，其中：

- `local` 表示数据文件在本地文件系统
- `online` 表示数据文件通过 URL 远程获取
- `<tags>` 为描述性标签（如 `eccdna`、`gdna`、`bam`、`integrated`）

#### Scenario: 重命名后的文件列表

- **WHEN** 用户查看 `samplesheets/` 目录
- **THEN** 应看到以下文件（原 `samplesheet_*` 前缀已移除）：
  - `test_online.csv`（原 `samplesheet.csv`）
  - `test_local_bam.csv`（原 `samplesheet_bam_test.csv`）
  - `test_local_eccdna.csv`（原 `samplesheet_local.csv`）
  - `test_local_gdna.csv`（原 `samplesheet_local_with_gdna.csv`）
  - `test_local_integrated.csv`（原 `samplesheet_integrated.csv`）

### Requirement: 真实模式测试

系统 SHALL 通过真实模式（非 stub）验证 ECCsplorer 模块在有无 gDNA control 时的正确性。

#### Scenario: 无 control 真实运行

- **WHEN** 执行 `nextflow run main.nf -profile test_local,docker`
- **THEN** 流程成功完成，`results/test_local/eccsplorer/` 下 `*_candidates.bed` 非空，`*_blast.m6` 非空

#### Scenario: 有 control 真实运行

- **WHEN** 执行 `nextflow run main.nf -profile test_local_gdna,docker`
- **THEN** 流程成功完成，`results/test_local_gdna/eccsplorer/` 下 `*_candidates.bed` 非空，`*_blast.m6` 非空

#### Scenario: gDNA control 效果对比

- **THEN** 有 control 的候选数量应 ≤ 无 control 的候选数量（gDNA control 去除背景噪音）

## MODIFIED Requirements

### Requirement: 配置文件引用

`conf/test_local.config` 和 `conf/test_local_gdna.config` 中的 `input` 参数 SHALL 指向重命名后的 samplesheet 文件。
