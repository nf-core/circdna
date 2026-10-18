# 子工作流简并 Spec

## Why

`reference_mode` 本质是 `eccsplorer_pipeline` 的子集（仅 BAM 预处理 + mosdepth，无检测工具），100% 冗余。`eccdna_mode` 和空目录需统一清理。

## What Changes

- `eccdna_mode` → `eccsplorer_pipeline`（重命名）
- `reference_mode` 合并入 `eccsplorer_pipeline`（mode='reference' 时运行无检测工具模式）
- 删除 5 个空/冗余子工作流目录

## Impact

- 修改: `workflows/circdna.nf` — 2 个 include 合并为 1 个，mode='reference' 路由到 ECCSPLORER_PIPELINE
- 删除: `eccdna_mode/`（重命名）、`reference_mode/`、`eccsplorer_all/`、`eccsplorer_cluster/`、`eccsplorer_mapping/`、`eccsplorer_prepare/`
- **NOT BREAKING**：功能不变，仅重构

## MODIFIED Requirements

### Requirement: ECCSPLORER_PIPELINE 统一入口

The system SHALL use `ECCSPLORER_PIPELINE` for both reference and eccDNA modes.

**Scenario: mode=reference**

- **WHEN** `params.mode == 'reference'`
- **THEN** calls ECCSPLORER_PIPELINE with all detection tools disabled (run_eccsplorer=false, run_circle_map=false, etc.)
- **AND** produces BAM + mosdepth outputs only

**Scenario: mode=eccdna**

- **WHEN** `params.mode == 'eccdna'`
- **THEN** calls ECCSPLORER_PIPELINE with detection tools enabled per circle_identifier
