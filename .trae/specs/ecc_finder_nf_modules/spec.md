# ecc_finder Nextflow 模块构建 - 产品需求文档

## Overview
- **Summary**: 构建 ecc_finder 的 conda 和 Docker 环境，并将其封装为四个独立的 Nextflow 运行模块（map-ont、map-sr、asm-ont、asm-sr），完成模块测试和验证。
- **Purpose**: 将 ecc_finder 工具集成到 bio.nf 工作流框架中，支持在 Nextflow 管道中调用 eccDNA 检测功能，实现标准化的数据分析流程。
- **Target Users**: 生物信息学研究人员，使用 Nextflow 进行 eccDNA 数据分析的用户。

## Goals
- 构建 ecc_finder 的 conda 环境配置文件，支持 Linux 平台
- 构建 ecc_finder 的 conda 包并上传到 anaconda.org
- 构建 ecc_finder 的 Docker 镜像并推送到 quay.io
- 创建四个 Nextflow 模块：map-ont、map-sr、asm-ont、asm-sr
- 为每个模块编写测试用例，验证功能正确性
- 确保所有模块支持 task.ext.args、task.ext.prefix 和 stub 测试模式

## Non-Goals (Out of Scope)
- 修改 ecc_finder 源代码逻辑
- 创建完整的 workflow 或 subworkflow

## Background & Context
ecc_finder 是一个用于检测 eccDNA（染色体外环状 DNA）的工具，支持 Illumina 短读段和 ONT 长读段数据。项目包含四个主要命令：
- `map-ont`: 从 Nanopore 长读段识别 eccDNA 位点
- `map-sr`: 从 Illumina 短读段识别 eccDNA 位点  
- `asm-ont`: 从 Nanopore 长读段组装 eccDNA
- `asm-sr`: 从 Illumina 短读段组装 eccDNA

项目已提供 `ecc_finder.yaml` 和 `ecc_finder_MAC.yaml` 作为 conda 环境配置参考。

## Functional Requirements
- **FR-1**: 创建精简的 conda environment.yml，包含 ecc_finder 运行所需的所有依赖
- **FR-2**: 创建 Dockerfile，构建包含 ecc_finder 的容器镜像
- **FR-3**: 创建 map-ont 模块，支持从 Nanopore 长读段识别 eccDNA
- **FR-4**: 创建 map-sr 模块，支持从 Illumina 短读段识别 eccDNA
- **FR-5**: 创建 asm-ont 模块，支持从 Nanopore 长读段组装 eccDNA
- **FR-6**: 创建 asm-sr 模块，支持从 Illumina 短读段组装 eccDNA
- **FR-7**: 为每个模块编写 meta.yml，描述模块输入输出
- **FR-8**: 为每个模块编写测试用例和配置文件

## Non-Functional Requirements
- **NFR-1**: 所有模块支持 task.ext.args 传递额外参数
- **NFR-2**: 所有模块支持 task.ext.prefix 自定义输出前缀
- **NFR-3**: 所有模块支持 stub 测试模式
- **NFR-4**: 模块命名遵循 ECC_FINDER_<SUBCOMMAND> 大写命名规范
- **NFR-5**: 模块结构与现有 cresil 模块保持一致

## Constraints
- **Technical**: 使用 Python 3.8+，依赖 minimap2、bwa、samtools、bedtools、TideHunter、Genrich、fastp 等工具
- **Dependencies**: ecc_finder 源代码位于 `/Users/siyangming/nextflow_nf_core/bio.nf/ecc_finder/`
- **Module Location**: 模块文件存储在 `/Users/siyangming/nextflow_nf_core/bio.nf/modules/ecc_finder/<module_name>/`

## Assumptions
- 测试数据使用 ecc_finder 自带的 test_samples 目录中的示例数据
- 参考基因组索引文件由用户预先准备
- Docker 构建在本地完成，不涉及推送到远程仓库

## Acceptance Criteria

### AC-1: Conda 环境配置
- **Given**: 存在 ecc_finder 源代码和原始 yaml 文件
- **When**: 创建 environment.yml 并运行 `conda env create`
- **Then**: conda 环境成功创建，所有依赖安装完成
- **Verification**: `programmatic`

### AC-2: Docker 镜像构建
- **Given**: 存在 Dockerfile 和 ecc_finder 源代码
- **When**: 运行 `docker build`
- **Then**: Docker 镜像成功构建，包含所有依赖工具
- **Verification**: `programmatic`

### AC-3: map-ont 模块
- **Given**: 存在 reference.idx、query.fq、reference.fa
- **When**: 运行 ECC_FINDER_MAP_ONT 模块
- **Then**: 输出 ecc.ont.csv 和 ecc.ont.fasta 文件
- **Verification**: `programmatic`

### AC-4: map-sr 模块
- **Given**: 存在 reference.idx、query.fq1、query.fq2、reference.fa
- **When**: 运行 ECC_FINDER_MAP_SR 模块
- **Then**: 输出 ecc.sr.csv 和 ecc.sr.fasta 文件
- **Verification**: `programmatic`

### AC-5: asm-ont 模块
- **Given**: 存在 query.fq 文件
- **When**: 运行 ECC_FINDER_ASM_ONT 模块
- **Then**: 输出 ecc.asm.ont.fasta 文件
- **Verification**: `programmatic`

### AC-6: asm-sr 模块
- **Given**: 存在 query.fq1、query.fq2 文件
- **When**: 运行 ECC_FINDER_ASM_SR 模块
- **Then**: 输出 ecc.asm.sr.fasta 文件
- **Verification**: `programmatic`

### AC-7: Stub 测试模式
- **Given**: 使用 -stub 选项运行模块
- **When**: 执行测试
- **Then**: 模块成功运行，输出空的占位文件
- **Verification**: `programmatic`

## Open Questions
- [ ] 是否需要支持 macOS 平台的 conda 环境？
- [ ] 是否需要创建综合测试数据（参考基因组等）？