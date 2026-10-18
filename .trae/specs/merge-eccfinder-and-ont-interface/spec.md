# ecc_finder 模块合并 + 子工作流 + 三代接口 Spec

## Why

`ecc_finder` 四个子模块（map_sr / asm_sr / map_ont / asm_ont）共享同一 Docker 镜像和 Conda 环境，应合并为单一模块文件。同时需要创建子工作流统一调度二代/三代模式，为 circdnalr 分支接入长读数据预留接口。

## What Changes

### 模块合并

- 将 4 个 ecc_finder 子模块合并到 `modules/local/ecc_finder/main.nf`
- 从 bio.nf 拷贝 map_ont 和 asm_ont（二代 SR 已在 circdna.nf）
- 删除 `modules/local/ecc_finder_map_sr/` 和 `modules/local/ecc_finder_asm_sr/`
- 统一 `environment.yml`（4 个完全一致，保留一份）

### 子工作流创建

- 新建 `subworkflows/local/ecc_finder_pipeline/main.nf`
- 输入：reads + bwa_index + fasta_meta + mode 参数
- 根据 `platform`（sr / ont）和 `method`（map / asm）分发到对应 process
- 当前仅激活 sr 路径，ont 路径代码就绪但默认关闭

### 二代/三代接口

- `platform=sr`：使用 BWA index + paired FASTQ + ref → 调用 MAP_SR / ASM_SR
- `platform=ont`：使用 minimap2 index + single FASTQ + ref → 调用 MAP_ONT / ASM_ONT（预留）
- ONT 模式：模块已存在、channel 已预留，由 `ecc_finder_platform` 参数控制

### 流程接入简化

- `eccdna_mode/main.nf` 中 include 从 2 个独立模块简化为 1 个子工作流 `ECC_FINDER_PIPELINE`
- `circle_identifier` 保持 `ecc_finder_map_sr` / `ecc_finder_asm_sr` 独立控制

## Impact

- 新增: `modules/local/ecc_finder/main.nf`（4 process）、`subworkflows/local/ecc_finder_pipeline/main.nf`
- 删除: `modules/local/ecc_finder_map_sr/`、`modules/local/ecc_finder_asm_sr/`
- 修改: `subworkflows/local/eccdna_mode/main.nf`、`conf/modules.config`
- **BREAKING**: `include` 路径变更，process 名称从 `ECC_FINDER_MAP_SR` → `ECC_FINDER_MAP_SR`（不变），但来源从独立模块变为 `ecc_finder/main`

## ADDED Requirements

### Requirement: 4 个 ecc_finder process 合并到单一 main.nf

The system SHALL have `ECC_FINDER_MAP_SR`, `ECC_FINDER_ASM_SR`, `ECC_FINDER_MAP_ONT`, `ECC_FINDER_ASM_ONT` in `modules/local/ecc_finder/main.nf`.

### Requirement: ECC_FINDER_PIPELINE 子工作流

The system SHALL provide a subworkflow that routes to the correct process based on platform and method.

**Scenario: SR map**

- **WHEN** platform=sr, method=map
- **THEN** routes to ECC_FINDER_MAP_SR with BWA index + R1 + R2 + ref

**Scenario: ONT map（预留）**

- **WHEN** platform=ont, method=map
- **THEN** routes to ECC_FINDER_MAP_ONT with minimap2 index + single query + ref

### Requirement: eccdna_mode 使用子工作流

The system SHALL use `ECC_FINDER_PIPELINE` subworkflow in `eccdna_mode` instead of individual process includes.
