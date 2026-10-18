# 迁移 bio.nf 模块到 local 统一管理 + 清理归档 Spec

## Why

当前黑盒版完整工具模块（cresil / ecc_finder / eccsplorer / flye / repeatexplorer2 / ampliconsuite）散落在 `bio.nf/modules/` 下，而 circdna.nf 的 slim 原子化模块都在 `circdna.nf/modules/local/`。跨仓库引用导致：`eccsplorer_pipeline` 引用 `modules/local/eccsplorer/main` 但该黑盒模块已被 master 删除（816918e）→ **断链**；`eccsplorer_clu_slim` 与 `eccsplorer_pipeline` 仍引用 `bio.nf/modules/...` → **跨仓耦合**。目标：把 6 个模块统一迁入 `local/`，让黑盒版与 slim 版并存、通过 circle_identifier 选择运行；同时清理 `scripts/` 冗余脚本并评估当前未提交更改。

## What Changes

- **迁移 6 个 bio.nf 模块到 `circdna.nf/modules/local/`**：
  - `cresil/` → `local/cresil/`（annotate/identify/identify_wgls/trim/visualize + testdata，仅物理移动）
  - `ecc_finder/{asm_ont,asm_sr,map_ont,map_sr}` → `local/ecc_finder/{...}`（**改为 nf-core 规范：一个 process 一个子目录**），删除旧单文件 `local/ecc_finder/main.nf`（含 4 process）
  - `eccsplorer/` → `local/eccsplorer/`（黑盒，含 `ECCSPLORER` + `ECCSPLORER_WITH_CONTROL`），修复 `eccsplorer_pipeline` 断链
  - `flye/` → `local/flye/`（仅物理移动）
  - `repeatexplorer2/` → `local/repeatexplorer2/`
  - `ampliconsuite/` → `local/ampliconsuite_ec/`（**重命名避免覆盖**：local 已有 `ampliconsuite/` 实为 prepareaa，供 ampliconarchitect 使用）
- **更新所有 include 引用**（跨仓 `bio.nf/modules/...` → `modules/local/...`）：
  - `eccsplorer_pipeline/main.nf`：AMPLICONSUITE、ECCSPLORER(\_WITH_CONTROL/\_CLU)
  - `ecc_finder_pipeline/main.nf`：4 个 ecc_finder 子模块
  - `eccsplorer_clu_slim/main.nf`、`tests/eccsplorer_clu/test.nf`：REPEATEXPLORER2
- **黑盒版与 slim 版并存，circle_identifier 选择**：黑盒（`eccsplorer`/`ecc_finder_map_sr` 等）与 slim（`eccsplorer_map_slim`/`ecc_finder_map_sr_slim` 等）标识符在 workflows 中共存，各自触发对应链
- **scripts/ 归档清理**：`prepare_eccsplorer_database.sh` → `testdatasets/eccsplorer_db/`、`update_samplesheets.py` → `samplesheets/`、删除 `test_incremental_cache.py`（引用已删除旧路径）；删除 `scripts/` 目录
- **评估并提交当前未提交更改**：PRExing 预处理链（`eccsplorer_prepare_*` + `seqkit/concat` + `ecc_preprocessing`）保留；samplesheet testdata→ngs 路径更新保留；黑盒残留（`ECCSPLORER`/`CLU_*` 旧配置块）丢弃
- **BREAKING**：`local/ecc_finder/main.nf` 单文件版删除（拆分为 4 子模块）；`bio.nf/ampliconsuite` 重命名为 `local/ampliconsuite_ec`

## Impact

- 受影响能力：ecc_finder 黑盒（4 模式）、eccsplorer 黑盒、repeatexplorer2 clu、ampliconsuite、cresil/flye（待用）
- 受影响代码：
  - `subworkflows/local/ecc_finder_pipeline/main.nf`（4 子模块 include）
  - `subworkflows/local/eccsplorer_pipeline/main.nf`（ECCSPLORER/AMPLICONSUITE include）
  - `subworkflows/local/eccsplorer_clu_slim/main.nf`、`tests/eccsplorer_clu/test.nf`（REPEATEXPLORER2 include）
  - `workflows/circdna.nf`（黑盒 circle_identifier 接入）
  - `scripts/`（归档/删除）、`modules/local/`（新迁入目录）

## ADDED Requirements

### Requirement: 6 个 bio.nf 模块迁入 local

系统 SHALL 将 cresil / ecc_finder / eccsplorer / flye / repeatexplorer2 / ampliconsuite 从 `bio.nf/modules/` 迁入 `circdna.nf/modules/local/`。

#### Scenario: ecc_finder 拆分为 nf-core 规范子模块

- **WHEN** 迁移 ecc_finder
- **THEN** `local/ecc_finder/{map_sr,asm_sr,map_ont,asm_ont}/` 各含单 process 的 main.nf；删除旧 `local/ecc_finder/main.nf`（4 process 单文件）；`ecc_finder_pipeline` 改引 4 个独立子模块

#### Scenario: 黑盒 eccsplorer 就位

- **WHEN** 迁移 eccsplorer
- **THEN** `local/eccsplorer/` 含 `ECCSPLORER` 与 `ECCSPLORER_WITH_CONTROL` process，`eccsplorer_pipeline` 引用不再断链

#### Scenario: ampliconsuite 命名冲突处理

- **WHEN** 迁移 ampliconsuite
- **THEN** bio.nf 的 AmpliconSuite（`ampliconsuite:1.6.0`）迁入 `local/ampliconsuite_ec/`，不覆盖 local 已有 prepareaa 版 `local/ampliconsuite/`

#### Scenario: 仅物理移动 cresil/flye

- **WHEN** 迁移 cresil / flye
- **THEN** 目录与 testdata 完整迁入 local，不新增 circle_identifier 运行路径

### Requirement: 更新所有跨仓 include 引用

系统 SHALL 将所有 `bio.nf/modules/...` 引用改为 `modules/local/...`。

#### Scenario: 引用更新

- **WHEN** 迁移完成后检查 include
- **THEN** `grep bio.nf/modules` 在 `*.nf` 中无残留（eccsplorer_pipeline/eccsplorer_clu_slim/tests 等全部指向 local）

### Requirement: 黑盒版与 slim 版 circle_identifier 并存

系统 SHALL 使黑盒与 slim 标识符在 workflows 中共存并按选择运行。

#### Scenario: 黑盒标识符触发

- **WHEN** `--circle_identifier 'eccsplorer,ecc_finder_map_sr'`（黑盒）
- **THEN** 触发黑盒链（eccsplorer_pipeline / ecc_finder 黑盒子模块）
- **WHEN** `--circle_identifier 'eccsplorer_map_slim,ecc_finder_map_sr_slim'`（slim）
- **THEN** 触发 slim 链（eccsplorer_slim_pipeline / ecc_finder_slim_pipeline）

### Requirement: scripts 归档清理

系统 SHALL 归档有用脚本并删除过时脚本，移除 scripts 目录。

#### Scenario: 归档与删除

- **WHEN** 处理 scripts 目录
- **THEN** `prepare_eccsplorer_database.sh` → `testdatasets/eccsplorer_db/`；`update_samplesheets.py` → `samplesheets/`；`test_incremental_cache.py` 删除；`scripts/` 目录删除

### Requirement: 未提交更改评估与提交

系统 SHALL 评估并提交当前未提交更改（保留有价值、丢弃黑盒残留）。

#### Scenario: PRExing 链保留

- **WHEN** 评估未提交更改
- **THEN** `bin/eccsplorer_prepare_*`、`modules/local/eccsplorer_prepare_*`、`modules/nf-core/seqkit/concat`、`subworkflows/local/ecc_preprocessing` 保留并提交

#### Scenario: 黑盒残留丢弃

- **WHEN** 解决冲突文件
- **THEN** `conf/modules.config` 中 stash 侧 `ECCSPLORER`/`CLU_*` 旧配置丢弃；`nextflow.config` 中旧黑盒参数丢弃，保留 `run_eccprepare`/`eccprepare_*` 参数

## MODIFIED Requirements

### Requirement: ecc_finder_pipeline 子模块引用

由单文件 `local/ecc_finder/main`（4 process）改为 4 个独立子模块 `local/ecc_finder/{map_sr,asm_sr,map_ont,asm_ont}/main`。

## REMOVED Requirements

### Requirement: local/ecc_finder 单文件黑盒

**Reason**: 不符合 nf-core「一个 process 一个子目录」规范，被 4 子模块拆分版替代
**Migration**: `ecc_finder_pipeline` 改引 4 个子模块；删除 `local/ecc_finder/{main.nf,environment.yml,meta.yml}`
