# circdnalr 同步 master 并进行长读长流程测试 Spec

## Why

`circdnalr` 分支承载长读长（PacBio/ONT）eccDNA 分析管线（CReSIL / FLED / FLYE / ECCFINDER 四引擎），但自 6401fa4 起与 `master` 分叉达 260 个文件（约 31 处冲突）。`master` 在此期间已大幅演进（ECC_SCORE 迁入 Snakemake、ONV/PacBio 内置测试数据集、samplesheet 重组、大基因组 CSI 索引等）。需要先将最新 `master` 同步进 `circdnalr`，再在服务器上用内置测试数据验证长读长管线可正常运行。

## What Changes

- 在 `circdnalr` 分支执行 `git merge master`，逐一解决约 31 处冲突（核心工作流、配置、samplesheet、校验脚本等）
- 冲突解决策略：**同时保留** master 的短读长/ECC_SCORE/Snakemake 逻辑 与 circdnalr 的长读长逻辑；二进制测试数据保留 master 版本
- 复用 master 内置 ONT/PacBio 冒烟测试数据（`testdatasets/ont/ont_eccdna_smoke.fastq.gz` 6.5MB、`testdatasets/pacbio/pacbio_eccdna_smoke.fastq.gz` 25MB）及 `testdatasets/reference/genome.fa`
- 更新/新建长读长测试配置：`protocol=ont/pacbio`、`long_read_identifier=cresil,fled,flye,eccfinder`（全部四引擎）、`entrypoint=cleaned_fastq`
- 合并结果推送到 GitHub（`origin/circdnalr`），同步到服务器 `/data1/users/siyangming/PlanteccDNADB/circdna.nf`
- 通过 SSH 在服务器上分别运行 ONT 与 PacBio 长读长测试（各含四引擎），验证流程完整执行并产出结果
- **注**：用户提到的 `metadata.csv`（`circdna.nf/circdna.nf/samplesheets/metadata.csv`）在本地与服务器均不存在（路径含重复 `circdna.nf` 段）；实际元数据源为 `SraRunInfo_eccDNA_all2.csv`。按用户确认，测试数据改用内置测试数据集，不再依赖 metadata.csv

## Impact

- **Affected specs**:
  - 长读长分析能力（CReSIL / FLED / FLYE / ECCFINDER）
  - 短读长分析能力（master 侧 ECC_SCORE / Snakemake 交接）
  - 样本表校验与解析（check_samplesheet.py / INPUT_CHECK）
- **Affected code**:
  - `workflows/circdna.nf`（长读长分支与短读长分支并存）
  - `nextflow.config` / `nextflow_schema.json`（长读长参数与 master 参数合并）
  - `conf/modules.config`（长读长模块 + master 短读长模块配置）
  - `conf/test_nanopore_lr.config` / `conf/test_pacbio_lr.config`（测试配置）
  - `bin/check_samplesheet.py` / `subworkflows/local/input_check/main.nf`
  - `samplesheets/*`（circdna*\* 与 circdnalr*\* 系列）
  - 服务器 `/data1/users/siyangming/PlanteccDNADB/circdna.nf`（同步 + 测试执行）

## ADDED Requirements

### Requirement: master 同步到 circdnalr

系统 SHALL 将 `master` 分支合并进 `circdnalr` 分支，并解决所有冲突，使合并后分支同时具备短读长与长读长能力。

#### Scenario: 成功合并

- **WHEN** 在 `circdnalr` 分支执行 `git merge master`
- **THEN** 合并完成且无未解决冲突，`git status` 干净，`nextflow config` 可正常解析，Python 校验脚本语法正确

### Requirement: 长读长四引擎测试配置

系统 SHALL 提供使用内置测试数据的 ONT 与 PacBio 长读长测试配置，覆盖全部四个检测引擎（cresil、fled、flye、eccfinder）。

#### Scenario: 测试配置可用

- **WHEN** 加载 ONT/PacBio 长读长测试配置
- **THEN** `protocol` 分别等于 `ont`/`pacbio`，`long_read_identifier` 包含 `cresil,fled,flye,eccfinder`，输入指向服务器上存在的内置 smoke 测试 FASTQ 与参考基因组

### Requirement: 服务器测试执行

系统 SHALL 通过 SSH 在服务器（192.168.16.65，代码与数据位于 `/data1/users/siyangming/PlanteccDNADB`）上运行长读长测试。

#### Scenario: ONT 测试成功

- **WHEN** 在服务器运行 ONT 长读长流程（四引擎，内置 smoke 数据）
- **THEN** 流程以成功状态结束，四个引擎均产生对应输出文件，无致命错误

#### Scenario: PacBio 测试成功

- **WHEN** 在服务器运行 PacBio 长读长流程（四引擎，内置 smoke 数据）
- **THEN** 流程以成功状态结束，四个引擎均产生对应输出文件，无致命错误

### Requirement: 测试结果记录

系统 SHALL 记录测试执行结果（命令、成功/失败状态、输出路径、主要输出文件），便于后续比对。

#### Scenario: 结果可追溯

- **WHEN** 测试完成
- **THEN** 结果记录包含运行命令、退出状态、输出目录及关键产出文件清单

## MODIFIED Requirements

### Requirement: 短读长能力保持

合并后 `circdnalr` 分支 SHALL 保持 master 的短读长能力（BWA 主链、ECC_SCORE/Snakemake 交接、samplesheet `data_type` 约定）不回退。

#### Scenario: 短读长兼容

- **WHEN** 在合并后分支上运行 `-profile test_local` 或 `nextflow config` 校验
- **THEN** 短读长相关参数与模块引用完整存在，无缺失报错

## REMOVED Requirements

无
