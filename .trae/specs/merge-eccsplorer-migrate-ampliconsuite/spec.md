# ECCsplorer 模块合并 + ampliconsuite 迁移到官方镜像 Spec

## Why

1. `eccsplorer` 和 `eccsplorer_clu` 两个模块共享同一个 Docker 镜像和 conda 依赖，应合并为单一模块减少维护成本
2. `ampliconsuite` 当前使用自定镜像 `nf-core/prepareaa:1.0.5`，biocontainers 已有官方镜像 `1.6.0`，应迁移至官方渠道并移至 `bio.nf` 模块集合留存

## What Changes

### 合并 ECCsplorer 模块

- 将 `ECCSPLORER_CLU` process 从 `eccsplorer_clu/main.nf` 迁入 `eccsplorer/main.nf`
- 删除 `modules/local/eccsplorer_clu/` 目录
- 更新 `eccdna_mode/main.nf` 中 `ECCSPLORER_CLU` 的 include 路径

### ampliconsuite 迁移

- 将 `circdna.nf/modules/local/ampliconsuite/` 移动到 `bio.nf/modules/ampliconsuite/`
- Docker 镜像替换：`nf-core/prepareaa:1.0.5` → `quay.io/biocontainers/ampliconsuite:1.6.0--pyh109da93_0`
- Conda 依赖替换：`bioconda::ampliconsuite=1.2.1` → `bioconda::ampliconsuite=1.6.0`，移除 `mosek::mosek`（mosek 许可证独立管理）
- 删除自定 Dockerfile（不再需要）
- 更新 `ampliconarchitect_pipeline/main.nf` 中的 include 路径

## Impact

- Affected code:
  - 修改: `modules/local/eccsplorer/main.nf` — 新增 ECCSPLORER_CLU process
  - 删除: `modules/local/eccsplorer_clu/` 整个目录
  - 修改: `subworkflows/local/eccdna_mode/main.nf` — include 路径
  - 移动: `modules/local/ampliconsuite/` → `bio.nf/modules/ampliconsuite/`
  - 修改: `bio.nf/modules/ampliconsuite/main.nf` — container + env
  - 删除: `modules/local/ampliconsuite/Dockerfile`
  - 修改: `subworkflows/local/ampliconarchitect_pipeline/main.nf` — include 路径
- **NOT BREAKING**：process 名称和接口不变

## ADDED Requirements

### Requirement: ECCSPLORER_CLU 并入 eccsplorer 模块

The system SHALL have all three ECCsplorer processes (ECCSPLORER, ECCSPLORER_WITH_CONTROL, ECCSPLORER_CLU) in a single `modules/local/eccsplorer/main.nf`.

**Scenario: 合并后**

- **WHEN** eccdna_mode include ECCSPLORER_CLU
- **THEN** 路径为 `from '../../../modules/local/eccsplorer/main'`
- **AND** `modules/local/eccsplorer_clu/` 目录不存在

### Requirement: ampliconsuite 使用官方 biocontainers 镜像

The system SHALL use `quay.io/biocontainers/ampliconsuite:1.6.0--pyh109da93_0` for the AMPLICONSUITE process.

**Scenario: 镜像替换**

- **WHEN** AMPLICONSUITE process 运行
- **THEN** container 为 `quay.io/biocontainers/ampliconsuite:1.6.0--pyh109da93_0`
- **AND** conda 依赖为 `bioconda::ampliconsuite=1.6.0`

### Requirement: ampliconsuite 移至 bio.nf

The system SHALL have ampliconsuite module at `bio.nf/modules/ampliconsuite/` as the canonical source.

**Scenario: 模块留存**

- **WHEN** 查看 bio.nf/modules/ampliconsuite/
- **THEN** 包含 `main.nf`、`meta.yml`、`environment.yml`（无 Dockerfile）
- **AND** `circdna.nf/modules/local/ampliconsuite/` 不存在
