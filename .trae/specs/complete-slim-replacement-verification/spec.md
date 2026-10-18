# 制作 slim 链 test\_ 测试文件并整理分支冻结修改 Spec

## Why

[slim-full-replacement.md](file:///Users/siyangming/nextflow_nf_core/.trae/documents/slim-full-replacement.md) 的 7 个工作块已实现并通过局部验证，剩余验证需真实数据（ONT/PacBio 长读 + control 配对短读）。用户已在 master 分支构建 ONT/PacBio 测试数据（提交 `7ee79a4` + `09c7473`），并指示：**制作好 `test_` 开头的测试文件（samplesheet）即可，用户自行放服务器测试**（参考基因组不下载，走 `conf/server.config` 的 `genomes` 链接）。同时评估 ECCsplorer 分支冻结修改（80 项 M/D + 未跟踪 slim 文件）保留价值并实际合并清理。

## What Changes

- **master 提交并入当前分支（先行）**：将 master 分支提交合并到 ECCsplorer 分支——当前 80 项未提交修改先暂存（stash -u），`git merge master` 后将冲突文件**以 master 内容为准**（`git checkout --theirs`），再恢复 slim 独有修改（stash pop，冲突以 master 为准）；merge 后 master 的 ONT/PacBio 数据（`testdatasets/ont|pacbio`）、`test_ont.csv`/`test_pacbio.csv`、master 版 `conf/server.config` 自动就位
- **测试数据就位（相对路径）**：确认 merge 后 `testdatasets/ont/ont_eccdna_smoke/regular.fastq.gz`、`testdatasets/pacbio/pacbio_eccdna_smoke/regular.fastq.gz`、`samplesheets/test_ont.csv`、`test_pacbio.csv`、`testdatasets/extract_test_data.sh` 就位；samplesheet 中 fastq 路径**改写为相对路径**（如 `testdatasets/ont/ont_eccdna_smoke.fastq.gz`），供用户放服务器仓库根目录即可运行
- **参考基因组走 master server.config**：采用 master 版 `conf/server.config`（`fasta_base_path=/data1/users/siyangming/PublicDB/reference` + `genomes` 块），验证运行由用户在服务器以 `-profile server` + `--genome <species>` 执行，不下载
- **制作 test\_ 测试 samplesheet**：
  - `samplesheets/test_ont.csv`：ONT 长读（smoke/regular），列 `sample,fastq_1,fastq_2,platform,protocol`，`platform=ont`，fastq 相对路径
  - `samplesheets/test_pacbio.csv`：PacBio 长读（smoke/regular），`platform=pacbio`，fastq 相对路径
  - 新 `samplesheets/test_eccsplorer_pair.csv`：短读 control 配对（gdna/circdna × 3），在 `test_real_integrated.csv` 基础上补 `pair` 列（p1/p2/p3）
- **分支冻结修改整理**：merge 后评估 ECCsplorer 分支独有修改（slim 模块链、vendored 模块、新子工作流、params 等）保留价值，产出保留/合并/丢弃清单；冲突文件按用户决策以 master 为准，其余有价值修改实际提交（不 rollback 用户修改）
- **BREAKING**：merge 会覆盖冲突文件为 master 版（工作区相关修改丢失，属用户指定决策）

## Impact

- 受影响能力：samplesheet 规范（test\_ 前缀、platform/protocol/pair 列）、ONT/PacBio 链测试入口、control/comparative/all 测试入口
- 受影响代码：
  - `samplesheets/test_ont.csv`、`samplesheets/test_pacbio.csv`（master checkout/确认）
  - `samplesheets/test_eccsplorer_pair.csv`（新建，补 pair 列）
  - `testdatasets/ont/`、`testdatasets/pacbio/`、`testdatasets/extract_test_data.sh`（checkout）
  - git：ECCsplorer 分支冻结修改整理

## ADDED Requirements

### Requirement: master 提交并入当前分支

系统 SHALL 先将 master 提交合并到 ECCsplorer 分支，冲突以 master 为准。

#### Scenario: merge 流程

- **WHEN** 当前分支存在 80 项未提交修改（M/D/??）
- **THEN** 执行 `git stash -u` 暂存 → `git merge master`（冲突文件 `git checkout --theirs` 以 master 为准，`git add`）→ 提交 merge → `git stash pop` 恢复 slim 独有修改（pop 冲突同样以 master 为准）；merge 后 `testdatasets/ont|pacbio`、`samplesheets/test_ont.csv`/`test_pacbio.csv`、master 版 `conf/server.config` 就位

#### Scenario: merge 后确认

- **WHEN** merge 完成
- **THEN** `git log` 显示 master 提交（7ee79a4/09c7473/816918e 等）已并入；数据文件与 master 版 server.config 存在；slim 独有未跟踪文件保留

#### Scenario: test_ont.csv / test_pacbio.csv 相对路径

- **WHEN** 检查 `samplesheets/test_ont.csv` 与 `test_pacbio.csv`
- **THEN** 列含 `sample,fastq_1,fastq_2,platform,protocol`；`platform` 为 `ont`/`pacbio`，`protocol=long_read`；fastq_1 路径为**相对路径**（`testdatasets/ont/ont_eccdna_smoke.fastq.gz`），用户将仓库放到服务器任意路径即可运行

#### Scenario: control 配对 samplesheet

- **WHEN** 生成 `samplesheets/test_eccsplorer_pair.csv`
- **THEN** 基于 `test_real_integrated.csv`（6 样本 gdna/circdna × 3）增加 `pair` 列（gdna_1/circdna_1=p1，以此类推），满足 clu_slim/control 按 pair join 的要求

### Requirement: 参考基因组走 master server.config

系统 SHALL 采用 master 版 `conf/server.config` 的参考基因组链接，不下载。

#### Scenario: 冲突文件以 master 为准

- **WHEN** 当前分支与 master 均修改过的文件（server.config/modules.config/nextflow.config/check_samplesheet.py 等）
- **THEN** 以 master 分支内容为准（checkout master 版覆盖工作区）

#### Scenario: 服务器运行入口

- **WHEN** 用户在服务器上执行验证
- **THEN** 以 `-profile server` + `--genome <species>`（master server.config `genomes` 块，`fasta_base_path=/data1/users/siyangming/PublicDB/reference`，如 `Arabidopsis_thaliana`→TAIR10）运行，参考基因组由 server.config 提供

### Requirement: ECCsplorer 分支冻结修改保留价值评估与实际整理

系统 SHALL 评估 ECCsplorer 分支独有冻结修改并实际合并清理。

#### Scenario: 保留价值评估

- **WHEN** 分析 `git status`（80 项 M/D/??）与 master 差异，区分：slim 链新增文件（ECCsplorer 独有）/ 冲突文件（用 master）/ 冗余废弃文件
- **THEN** 产出保留/合并/丢弃清单：slim 模块链（eccsplorer_slim/ecc_finder_slim/vendored 模块/新子工作流）为保留价值核心；冲突文件按 master 版；冗余文件标记丢弃

#### Scenario: 实际提交与清理

- **WHEN** 清单经用户确认
- **THEN** 有价值修改 git add + commit 到 ECCsplorer 分支；无价值项删除或保留待确认；不 rollback 用户修改

## MODIFIED Requirements

无。

## REMOVED Requirements

无。
