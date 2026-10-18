# 同步服务器路径变更 Spec

## Why

服务器上 circdna.nf 项目路径已从 `/data1/users/siyangming/nextflow_nf_core/circdna.nf/` 迁移至 `/data1/users/siyangming/PlanteccDNADB/circdna.nf/`，但代码仓库中仍有多个文件引用旧路径，需同步更新。

## What Changes

- **SERVER_RUN_GUIDE.md**: 16+ 处 `cd /data1/users/siyangming/nextflow_nf_core/circdna.nf/` → `cd /data1/users/siyangming/PlanteccDNADB/circdna.nf/`
- **samplesheets/data_issues.txt**: 1 处旧服务器路径
- **errors.txt**: 1 处旧服务器路径
- **circdna.nf/AGENTS.md**: 跨项目路径表更新

**不变**: 本地 Mac 路径 `/Users/siyangming/nextflow_nf_core/` 仍为当前本地工作目录，无需修改。

## Impact

- Affected code: SERVER_RUN_GUIDE.md, data_issues.txt, errors.txt, AGENTS.md
- 服务器上的流程运行命令将指向正确路径
- 本地开发路径不受影响

## ADDED Requirements

### Requirement: 服务器路径同步

系统 SHALL 将所有引用旧服务器路径 `/data1/users/siyangming/nextflow_nf_core/circdna.nf/` 的文件更新为新路径 `/data1/users/siyangming/PlanteccDNADB/circdna.nf/`。

#### Scenario: 服务器运行指南路径正确

- **WHEN** 用户按照 SERVER_RUN_GUIDE.md 连接服务器执行命令
- **THEN** `cd` 命令指向 `/data1/users/siyangming/PlanteccDNADB/circdna.nf/`

#### Scenario: 本地路径不受影响

- **WHEN** 用户在本地 Mac 上运行 test_local 或 test_integrated 配置
- **THEN** 配置文件中的 `/Users/siyangming/nextflow_nf_core/` 路径仍然有效

## MODIFIED Requirements

无

## REMOVED Requirements

无
