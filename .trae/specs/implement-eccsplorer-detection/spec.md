# 实现 ECCsplorer 真实检测逻辑 Spec

## Why

当前 [circdna.nf/modules/local/eccsplorer/main.nf](file:///Users/siyangming/nextflow_nf_core/circdna.nf/modules/local/eccsplorer/main.nf) 是占位模块，仅 `cat` 硬编码假数据，未调用真实 ECCsplorer 软件。需要参照 [crimBubble/ECCsplorer](https://github.com/crimBubble/ECCsplorer) 实现真正的 eccDNA 检测能力，使 eccdna 模式的 `eccsplorer_bed` 产物具备生物学意义。

## What Changes

### 阶段一：构建 ECCsplorer Conda 包与 Docker 镜像

- 调研 [quay.io/biocontainers](https://quay.io/organization/biocontainers) 和 [bioconda](https://bioconda.github.io/) 是否已有 ECCsplorer 镜像/包
- **若 conda 包不存在**：
  - 编写 `meta.yaml` + `build.sh` 构建 ECCsplorer conda 包，包含所有依赖（Python 3.7、numpy、biopython、scipy、pyRserve、R 包）
  - 系统级依赖（blast+、segemehl、samtools、bedtools、RepeatExplorer2、Trimmomatic）作为依赖声明或单独打包
  - 本地构建验证：`conda build` + `conda create -n test --use-local`
  - 推送 conda 包至用户 anaconda.org 频道（`https://anaconda.org/siyangming/eccsplorer`）
- **若 Docker 镜像不存在**：
  - 编写 Dockerfile 构建 ECCsplorer 镜像，基于 conda 包或直接安装所有依赖（Python 3.7、blast+、segemehl、samtools、bedtools、RepeatExplorer2、R 包、Trimmomatic）
  - 推送镜像至用户 quay.io 频道（`quay.io/siyangming/eccsplorer`）
- **若都存在**：直接使用官方镜像和包

### 阶段二：在 bio.nf 构建 eccsplorer 模块

- 在 `bio.nf/` 创建新分支 `eccsplorer`
- 在 `bio.nf/modules/eccsplorer/` 构建符合 nf-core 标准的模块：
  - `main.nf` — 调用真实 ECCsplorer.py
  - `meta.yml` — 模块元数据
  - `environment.yml` — conda 环境
  - `tests/main.nf.test` — nf-test 测试用例
  - `tests/nextflow.config` — 测试配置（含 `--platform linux/amd64`）
- 输入接口：`tuple val(meta), path(reads_r1), path(reads_r2), path(fasta)` + 可选 control reads
- 输出：`candidates.bed`（候选 eccDNA）、`junction_reads.txt`、`versions.yml`
- 使用 stub 模式支持测试

### 阶段三：拷贝模块至 circdna.nf 并接入

- 将 `bio.nf/modules/eccsplorer/` 拷贝至 `circdna.nf/modules/local/eccsplorer/`
- **BREAKING**：调整 [eccdna_mode/main.nf](file:///Users/siyangming/nextflow_nf_core/circdna.nf/subworkflows/local/eccdna_mode/main.nf) 的 ECCSPLORER 调用，输入从 `BAM + BAI` 改为 `FASTQ R1 + R2 + FASTA`（ECCsplorer 内部使用 segemehl 自行比对，不接受 BAM）
- 更新 [conf/modules.config](file:///Users/siyangming/nextflow_nf_core/circdna.nf/conf/modules.config) 添加 ECCSPLORER 资源配置
- 更新 [nextflow.config](file:///Users/siyangming/nextflow_nf_core/circdna.nf/nextflow.config) 添加 ECCsplorer 相关参数（如 `--eccsplorer_database`、`--eccsplorer_trim_reads` 等）
- 版本 bump：v4.0.0 → v4.1.0（MINOR，新功能）

### 阶段四：更新 AGENTS.md 工作流程规则

- 在 [AGENTS.md](file:///Users/siyangming/nextflow_nf_core/AGENTS.md) 新增章节 "12. 第三方模块构建工作流程"
- 规定：对于 nf-core 未提供的软件，先在 bio.nf 构建模块，再拷贝至目标流程接入
- 规定：对于 bioconda 不可用的软件，构建自定义 conda 包并推送至用户 anaconda.org 频道（`https://anaconda.org/siyangming/<software>`）
- 规定：对于 biocontainers 不可用的软件，构建自定义 Docker 镜像并推送至用户 quay.io 频道（`quay.io/siyangming/<software>`）
- 规定：conda 包与 Docker 镜像版本号保持一致，便于追溯

## Impact

- **Affected specs**: `migrate-exploration-to-snakemake`（eccdna_mode 子工作流接口变更）
- **Affected code**:
  - [circdna.nf/modules/local/eccsplorer/main.nf](file:///Users/siyangming/nextflow_nf_core/circdna.nf/modules/local/eccsplorer/main.nf) — 替换占位实现
  - [circdna.nf/subworkflows/local/eccdna_mode/main.nf](file:///Users/siyangming/nextflow_nf_core/circdna.nf/subworkflows/local/eccdna_mode/main.nf) — 输入接口调整
  - [circdna.nf/conf/modules.config](file:///Users/siyangming/nextflow_nf_core/circdna.nf/conf/modules.config) — 资源配置
  - [circdna.nf/nextflow.config](file:///Users/siyangming/nextflow_nf_core/circdna.nf/nextflow.config) — 新参数
  - [bio.nf/modules/eccsplorer/](file:///Users/siyangming/nextflow_nf_core/bio.nf/modules/eccsplorer/) — 新模块
  - [AGENTS.md](file:///Users/siyangming/nextflow_nf_core/AGENTS.md) — 新增章节 12

## ADDED Requirements

### Requirement: ECCsplorer 真实检测

系统 SHALL 在 eccdna 模式下调用真实 ECCsplorer 软件（v2022.01.1.1），从 paired-end FASTQ 数据检测 eccDNA 候选。

#### Scenario: mapping 模式（仅参考基因组）

- **WHEN** 用户提供 eccDNA FASTQ R1/R2 + 参考基因组 FASTA，且 `--run_eccsplorer true`
- **THEN** ECCSPLORER 进程调用 `python ECCsplorer.py reads_R1 reads_R2 -ref sequence.fa`，产出 `*_candidates.bed`

#### Scenario: comparative 模式（含 control reads）

- **WHEN** 用户提供 eccDNA FASTQ + control FASTQ + 参考基因组 FASTA
- **THEN** ECCSPLORER 调用 4-read 参数模式，同时执行 clustering 和 mapping 模块

#### Scenario: stub 模式

- **WHEN** 使用 `-stub-run` 运行
- **THEN** 生成空的 `*_candidates.bed` 和 `*_junction_reads.txt`，不执行真实检测

### Requirement: Conda 包可用性

系统 SHALL 提供 ECCsplorer conda 包，包含所有 Python 和 R 依赖（Python 3.7、numpy、biopython、scipy、pyRserve、R ggplot2/ggrepel/gridExtra/dplyr）。

#### Scenario: conda 包不存在于 bioconda

- **WHEN** bioconda 无 ECCsplorer 包
- **THEN** 编写 `meta.yaml` + `build.sh`，本地 `conda build` 验证
- **AND** 推送至 `https://anaconda.org/siyangming/eccsplorer`
- **AND** environment.yml 指向 `- siyangming::eccsplorer=2022.01.1.1`

### Requirement: Docker 镜像可用性

系统 SHALL 提供 ECCsplorer Docker 镜像，包含所有依赖（Python 3.7、blast+、segemehl、samtools、bedtools、RepeatExplorer2、R 包、Trimmomatic）。

#### Scenario: 镜像不存在于 biocontainers

- **WHEN** quay.io/biocontainers 无 ECCsplorer 镜像
- **THEN** 构建自定义 Dockerfile，推送至 `quay.io/siyangming/eccsplorer:<version>`
- **AND** main.nf 的 container 指向该镜像

### Requirement: bio.nf 模块构建规范

系统 SHALL 在 bio.nf 项目中构建 eccsplorer 模块，遵循 nf-core 模块标准。

#### Scenario: 模块构建

- **WHEN** 在 bio.nf 创建 eccsplorer 分支
- **THEN** 模块目录包含 main.nf、meta.yml、environment.yml、tests/main.nf.test、tests/nextflow.config
- **AND** 测试通过 stub 模式验证

### Requirement: AGENTS.md 工作流程规则

系统 SHALL 在 AGENTS.md 新增章节 12，规定第三方模块构建工作流程。

#### Scenario: nf-core 未提供的软件

- **WHEN** 需要集成 nf-core 未提供的软件工具
- **THEN** 先在 bio.nf 创建分支构建模块，测试通过后拷贝至目标流程的 modules/local/ 目录

#### Scenario: bioconda 不可用

- **WHEN** 目标软件不在 bioconda
- **THEN** 构建自定义 conda 包，推送至用户 anaconda.org 频道（`https://anaconda.org/siyangming/<software>`）

#### Scenario: biocontainers 不可用

- **WHEN** 目标软件不在 quay.io/biocontainers
- **THEN** 构建自定义 Docker 镜像，推送至用户 quay.io 频道（`quay.io/siyangming/<software>`)

## MODIFIED Requirements

### Requirement: eccdna_mode 子工作流输入接口

**修改前**：ECCSPLORER 接收 `BAM + BAI + FASTA`

**修改后**：ECCSPLORER 接收 `FASTQ R1 + R2 + FASTA`（+ 可选 control FASTQ R1/R2），因为 ECCsplorer 内部使用 segemehl 自行比对，不接受预比对的 BAM

## REMOVED Requirements

### Requirement: 占位 ECCSPLORER 模块

**Reason**: 占位模块仅产出硬编码假数据，无生物学意义
**Migration**: 替换为真实 ECCsplorer 调用，输入接口从 BAM 调整为 FASTQ

## 技术参考

- **软件仓库**: https://github.com/crimBubble/ECCsplorer
- **安装说明**: https://github.com/crimBubble/ECCsplorer/blob/master/tutorials/Installation_instructions.md
- **引用**: Mann, L., et al. BMC Bioinformatics 23, 40 (2022). https://doi.org/10.1186/s12859-021-04545-2
- **最新版本**: v2022.01.1.1 (2024-04-15)
- **依赖**: Python 3.7.10, numpy, biopython, scipy, pyRserve, R (ggplot2/ggrepel/gridExtra/dplyr), blast+, segemehl, samtools≥1.9, bedtools≥2.28.0, RepeatExplorer2, Trimmomatic (可选), seqtk (可选)
- **用户频道**:
  - conda: `https://anaconda.org/siyangming` (用于推送 eccsplorer conda 包)
  - Docker: `quay.io/siyangming` (用于推送 eccsplorer Docker 镜像)

## 模块输入输出契约

### 输入

| 通道       | 类型        | 说明                          |
| ---------- | ----------- | ----------------------------- |
| meta       | val(map)    | 样本元数据                    |
| reads_r1   | path        | eccDNA reads R1 (FASTQ/FASTA) |
| reads_r2   | path        | eccDNA reads R2 (FASTQ/FASTA) |
| fasta      | path        | 参考基因组 FASTA              |
| control_r1 | path (可选) | control reads R1              |
| control_r2 | path (可选) | control reads R2              |

### 输出

| 通道           | 文件模式               | 说明               |
| -------------- | ---------------------- | ------------------ |
| candidates_bed | `*_candidates.bed`     | eccDNA 候选区域    |
| junction_reads | `*_junction_reads.txt` | junction read 信息 |
| versions       | `versions.yml`         | 版本信息           |
