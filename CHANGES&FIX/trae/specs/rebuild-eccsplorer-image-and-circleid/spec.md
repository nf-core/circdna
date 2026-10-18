# ECCsplorer 镜像完善 + circle_identifier 接入 + 测试验证 Spec

## Why

当前 `quay.io/bioinfortools/eccsplorer:2022.01.1.1` 镜像中的 RepeatExplorer2（seqclust）下载失败（repeatexplorer.org 返回 404）。需利用 `kavonrtep/repeatexplorer:2.3.8` 完善镜像。同时 ECCsplorer 目前通过 `eccsplorer_map_core` 参数控制启停，应统一到 `circle_identifier` 体系中。samplesheet 需增加 `pair` 列显式描述 eccDNA 与 gDNA 的配对关系。

**分两阶段执行**：先在 `bio.nf`（模块集合）中完善 Docker 镜像并跑通 nf-test，确认无误后再接入 `circdna.nf`。

## What Changes

### Phase 1 — bio.nf 中完善

- **Dockerfile 多阶段构建**：第一阶段从 `kavonrtep/repeatexplorer:2.3.8` 复制 `seqclust` 及 Perl/R 依赖，第二阶段在原 `condaforge/mambaforge` 基础上构建 ECCsplorer 环境
- **重建并推送镜像**：`quay.io/bioinfortools/eccsplorer:2022.01.1.1`（名称/版本不变）
- **Conda 包声明更新**：`meta.yaml` 注释说明 RepeatExplorer2 提供方式
- **nf-test 真实验证**：在 bio.nf 中添加非 stub 的 nf-test 用例，验证 ECCSPLORER process 在新镜像下实际运行成功

### Phase 2 — circdna.nf 接入

- **模块同步**：将 bio.nf 中验证通过的 `modules/eccsplorer/` 同步到 `circdna.nf/modules/local/eccsplorer/`
- **samplesheet 增加 `pair` 列**：新增可选列 `pair`，专门描述 eccDNA↔gDNA 配对关系。`group` 列保留用于通用样本分组，两者语义独立
- **circle_identifier 接入**：通过 `circle_identifier` 中是否包含 `eccsplorer` 控制启停（替代 `eccsplorer_map_core`）
- **自动 control 检测**：根据 samplesheet 的 `data_type` 区分 eccDNA/gDNA，按 `pair` 列自动配对
- **配置清理**：移除 `eccsplorer_map_core` 冗余参数
- **testdata 端到端测试**：用 `testdatasets/testdata` + `-profile test_local,docker` 验证全流程

## Impact

- Affected code:
  - Phase 1: `ECCsplorer/Dockerfile`, `ECCsplorer/conda-recipe/meta.yaml`, `bio.nf/modules/eccsplorer/tests/main.nf.test`
  - Phase 2: `circdna.nf/modules/local/eccsplorer/*`, `circdna.nf/workflows/circdna.nf`, `circdna.nf/nextflow.config`, `circdna.nf/conf/modules.config`, `circdna.nf/conf/test_local.config`, `circdna.nf/conf/server.config`, `circdna.nf/conf/test.config`, `circdna.nf/bin/check_samplesheet.py`, `circdna.nf/subworkflows/local/input_check/main.nf`, `circdna.nf/subworkflows/local/eccdna_mode/main.nf`, `circdna.nf/assets/schema_input.json`
- **NOT BREAKING**：镜像名称/版本不变；`pair` 为可选列，不影响已有 samplesheet

## ADDED Requirements

### Phase 1

#### Requirement: Dockerfile 多阶段构建集成 RepeatExplorer2

The system SHALL build `quay.io/bioinfortools/eccsplorer:2022.01.1.1` via multi-stage Dockerfile that copies `seqclust` from `kavonrtep/repeatexplorer:2.3.8`.

**Scenario: seqclust 可执行**

- **WHEN** 新镜像构建完成
- **THEN** `docker run --rm quay.io/bioinfortools/eccsplorer:2022.01.1.1 seqclust --help` 返回正常帮助信息

**Scenario: ECCsplorer.py 可执行**

- **WHEN** 新镜像构建完成
- **THEN** `docker run --rm quay.io/bioinfortools/eccsplorer:2022.01.1.1 ECCsplorer.py --help` 返回正常帮助信息

#### Requirement: bio.nf 中 nf-test 真实验证

The system SHALL have a non-stub nf-test case in `bio.nf/modules/eccsplorer/tests/main.nf.test`.

**Scenario: nf-test 真实运行通过**

- **WHEN** 执行 `nf-test test modules/eccsplorer/tests/main.nf.test`
- **THEN** ECCSPLORER process 成功完成（退出码 0）
- **AND** 输出 `*_candidates.bed` 非空
- **AND** 输出 `*_ecc_sequences.fasta` 非空

#### Requirement: Conda 包声明更新

**Scenario**: meta.yaml 注释说明 RepeatExplorer2 通过 Docker 镜像集成，conda 环境不含 seqclust

### Phase 2

#### Requirement: Samplesheet 增加 pair 列

The system SHALL support an optional `pair` column that describes eccDNA-gDNA pairing. `group` column is kept for general sample grouping — the two serve different purposes.

**Scenario: pair 列定义**

- **WHEN** samplesheet 包含 `pair` 列
- **THEN** 同一 `pair` 值的 eccDNA 和 gDNA 样本被视为配对
- **AND** `pair` 为空或不存在的样本视为无配对

**Scenario: group 列不变**

- **WHEN** samplesheet 包含 `group` 列
- **THEN** `group` 仅用于通用分组，ECCSPLORER_WITH_CONTROL 的配对不受 `group` 值影响

#### Requirement: check_samplesheet.py 支持 pair 列

**Scenario**: `bin/check_samplesheet.py` 接受 `pair` 为可选列并传递到输出

#### Requirement: input_check 读取 pair 到 meta

**Scenario**: `subworkflows/local/input_check/main.nf` 从 `pair` 列填充 `meta.pair`

#### Requirement: ECCDNA_MODE 使用 pair 列配对

The system SHALL use only `meta.pair` (not `meta.group`) for eccDNA/gDNA pairing.

**Scenario: 按 pair 配对**

- **WHEN** `data_type=eccDNA, pair=A` 和 `data_type=gDNA, pair=A`
- **THEN** 运行 ECCSPLORER_WITH_CONTROL

**Scenario: pair 不一致不配对**

- **WHEN** eccDNA `pair=A`，gDNA `pair=B`
- **THEN** 各自单独处理（eccDNA→ECCSPLORER，gDNA→仅覆盖度统计）

#### Requirement: ECCsplorer 通过 circle_identifier 选择

**Scenario: 包含 eccsplorer** → 运行（data_type + pair 自动配对）
**Scenario: 不含 eccsplorer** → 不运行

#### Requirement: 移除 eccsplorer_map_core 参数

The system SHALL remove `eccsplorer_map_core`, use `circle_identifier` default containing `eccsplorer`.

#### Requirement: testdata 端到端测试

**Scenario**: `-profile test_local,docker -resume` 成功，ECCsplorer 输出正确
