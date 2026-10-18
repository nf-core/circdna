# ECCsplorer clu/all 模式落地 Spec

## Why

ECCsplorer clu（RepeatExplorer2 聚类）和 all（map+clu+comparative）模式的 Nextflow 模块当前均为骨架 stub（仅生成空文件），未实际执行 ECCsplorer CLI。Docker 镜像已包含 seqclust，`circle_identifier` 体系已就绪，只需将 stub 替换为真实实现并接入流程。

## What Changes

### clu 模式

- **简化架构**：用单个 `ECCSPLORER_CLU` process 替代当前 3 个 stub 子模块（prepare/core/plot）。ECCsplorer CLI 的 `--mode clu` 内部已包含完整流程，无需在 Nextflow 侧拆开
- **输入**：FASTQ（treatment）+ FASTQ（control）+ taxon 参数
- **输出**：cluster candidates csv、cluster table、HTML 摘要
- **接入**：通过 `circle_identifier` 中的 `eccsplorer` + `eccsplorer_clu=true` 参数控制
- **数据库**：PublicDB 路径通过 `eccsplorer_database` 挂载（供 seqclust 使用）

### all 模式

- **DAG 编排**：all = map (已有) + clu (新增)，在 `ECCDNA_MODE` 中顺序执行
- 不需要独立的 all 子工作流——在 `eccdna_mode` 中通过布尔开关控制

### 清理

- 移除 3 个 stub 子模块及 `eccsplorer_cluster`/`eccsplorer_all` 骨架子工作流
- 保留 stub 代码的归档（或直接删除，因为已有 git 历史）

## Impact

- Affected code:
  - 新增: `circdna.nf/modules/local/eccsplorer_clu/main.nf`
  - 删除: `modules/local/eccsplorer_clu_prepare/`, `modules/local/eccsplorer_clu_core/`, `modules/local/eccsplorer_clu_candidates_plot/`
  - 删除: `subworkflows/local/eccsplorer_cluster/`, `subworkflows/local/eccsplorer_all/`
  - 修改: `circdna.nf/subworkflows/local/eccdna_mode/main.nf` — 接入 clu 流程
  - 修改: `circdna.nf/workflows/circdna.nf` — 解析 eccsplorer_clu 参数
  - 修改: `circdna.nf/conf/modules.config` — 添加 ECCSPLORER_CLU 配置
  - 修改: `circdna.nf/conf/test_local.config` — 添加 test 参数
- **BREAKING**：`eccsplorer_cluster` / `eccsplorer_all` 子工作流被移除（未被任何地方调用）

## ADDED Requirements

### Requirement: ECCSPLORER_CLU 单 process 实现

The system SHALL provide a single `ECCSPLORER_CLU` process that calls `ECCsplorer.py --mode clu`.

**Scenario: clu 运行**

- **WHEN** 提供 treatment FASTQ + control FASTQ + `--taxon vir`
- **THEN** ECCsplorer 执行 prepare_for_clustering → RepeatExplorer 聚类
- **AND** 输出 `*_cluster_candidates.csv`、`*_comparative_cluster_table.csv`、`*_eccCL_summary.html`

**Scenario: clu 输入无需参考基因组**

- **WHEN** 运行 clu 模式
- **THEN** 不需要参考基因组 FASTA（clu 模式只用 FASTQ 输入）

### Requirement: clu 通过 eccsplorer_clu 参数控制

The system SHALL control clu mode via `params.eccsplorer_clu` (default: false).

**Scenario: 启用 clu**

- **WHEN** `--eccsplorer_clu true` 且 `circle_identifier` 含 `eccsplorer`
- **THEN** map 模式完成后自动执行 clu 模式

**Scenario: 仅 map**

- **WHEN** `--eccsplorer_clu false`（默认）
- **THEN** 仅执行 map 模式（当前行为）

### Requirement: 数据库挂载

The system SHALL make PublicDB（Dfam/RepBase）accessible to RepeatExplorer in the container.

**Scenario: 数据库路径**

- **WHEN** ECCSPLORER_CLU 启动
- **THEN** `/Users/siyangming/Library/CloudStorage/OneDrive-个人/Project_backup/PublicDB/` 可在容器内访问（通过 Docker 卷挂载或绝对路径）

### Requirement: 清理 stub 代码

The system SHALL remove the 3 stub clu sub-modules and the 2 skeleton subworkflows.

**Scenario: 清理后**

- **WHEN** 删除完成
- **THEN** `modules/local/` 下不再有 `eccsplorer_clu_prepare/`、`eccsplorer_clu_core/`、`eccsplorer_clu_candidates_plot/`
- **AND** `subworkflows/local/` 下不再有 `eccsplorer_cluster/`、`eccsplorer_all/`

### Requirement: testdata 测试验证

The system SHALL pass clu mode test with real testdata using `-resume`.

**Scenario: clu 测试通过**

- **WHEN** 执行 `nextflow run main.nf -profile test_local,docker --eccsplorer_clu true -resume`
- **THEN** ECCSPLORER map 任务 resume（缓存命中）
- **AND** ECCSPLORER_CLU 任务执行成功
- **AND** 输出 cluster candidates 文件
- **AND** 整体流程退出码 0
