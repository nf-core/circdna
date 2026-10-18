# CReSIL Nextflow 模块 - Product Requirement Document

## Overview
- **Summary**: 在 bio.nf 项目的 modules 目录下构建 CReSIL 工具的 Nextflow 流程模块，包含 trim、identify、identify_wgls、annotate、visualize 五个子命令模块，使用已构建好的 Docker 和 Conda 环境，并用 example 数据进行测试。
- **Purpose**: 将 CReSIL eccDNA 检测工具集成到 Nextflow 模块化工作流中，遵循 nf-core / bio.nf 模块规范，方便在 Pipeline 中复用。
- **Target Users**: 生物信息学分析人员，使用 Nextflow 构建 eccDNA 分析流程的研究者。

## Goals
- 构建 5 个 CReSIL 子命令的独立 Nextflow 模块（trim, identify, identify_wgls, annotate, visualize）
- 每个模块遵循 bio.nf 现有模块规范（main.nf, meta.yml, environment.yml, tests/）
- 支持 Docker (quay.io/bioinfortools/cresil:1.2.0) 和 Conda (yangmingsi/cresil) 两种运行环境
- 使用 example 数据完成模块功能测试
- 提供 testdata 测试数据目录（含参考基因组等必需文件）
- 提供 stub 运行模式

## Non-Goals (Out of Scope)
- 不构建完整的端到端 subworkflow（仅构建独立的 process 模块）
- 不修改 CReSIL 源代码
- 不构建 arm64 架构的镜像（CReSIL 依赖的生物信息学工具无 arm64 版本）
- 不提供 GUI 或 Web 界面
- 不包含 Circos 可视化执行（仅生成配置文件）

## Background & Context
- CReSIL 是从 Nanopore 长读长测序数据中检测 eccDNA 的工具，版本 1.2.0
- Docker 镜像已构建并推送到 quay.io/bioinfortools/cresil:1.2.0（linux/amd64）
- Conda 包已发布到 anaconda.org/yangmingsi/cresil（linux-64）
- bio.nf 项目已有 flair、longshot、minimap2、samtools、kingfisher 等模块，遵循统一的目录结构和命名规范
- 现有模块使用 nf-core 风格：process 命名全大写 + 下划线，meta map 传递样本信息，versions 输出

## Functional Requirements

### FR-1: CReSIL TRIM 模块
- 输入：样本 meta、FASTQ/FASTA 读长文件、Minimap2 索引 (.mmi)
- 输出：trim.txt（修剪结果表）、versions.yml
- 支持参数：线程数、mapq、gap、overlap、exclude chromosome、输出目录
- 通过 `task.ext.args` 支持额外参数传递

### FR-2: CReSIL IDENTIFY 模块
- 输入：样本 meta、参考基因组 FASTA、FAI 索引、FASTQ/FASTA 读长、trim 结果表（可选）、medaka 共识模型
- 输出：eccDNA_final.txt（最终 eccDNA 列表）、versions.yml
- 支持参数：线程数、skip-variant、skip-gfa、最小区域大小、overlap 大小、断点深度、平均深度
- 通过 `task.ext.args` 支持额外参数传递

### FR-3: CReSIL IDENTIFY_WGLS 模块
- 输入：样本 meta、Minimap2 索引 (.mmi)、参考基因组 FASTA、FAI 索引、FASTQ/FASTA 读长、trim 结果表（可选）
- 输出：eccDNA_final.txt、versions.yml
- 支持参数：线程数、模式 (linkage/depth)、最小区域大小、overlap 大小、断点深度、平均深度、窗口合并参数等
- 通过 `task.ext.args` 支持额外参数传递

### FR-4: CReSIL ANNOTATE 模块
- 输入：样本 meta、identify 结果表、repeat masker BED（可选）、CpG islands BED（可选）、gene annotation BED（可选）
- 输出：注释结果（gene.annotate.txt, CpG.annotate.txt, repeat.annotate.txt, variant.annotate.txt）、versions.yml
- 支持参数：线程数
- 通过 `task.ext.args` 支持额外参数传递

### FR-5: CReSIL VISUALIZE 模块
- 输入：样本 meta、identify 结果表、eccDNA ID
- 输出：Circos 配置文件目录（for_Circos/）、versions.yml
- 支持参数：线程数、unit circos、mode circos
- 通过 `task.ext.args` 支持额外参数传递

### FR-6: 测试数据准备
- 在 modules/cresil/testdata/ 下提供测试用数据
- 包括：example FASTQ、参考基因组 FASTA + FAI、Minimap2 索引 .mmi
- 注释 BED 文件（可选，用于 annotate 测试）

### FR-7: 单元测试
- 每个模块有 tests/main.nf.test
- 包含正常运行测试和 stub 测试
- 使用 snapshot 验证输出

## Non-Functional Requirements

### NFR-1: 架构兼容性
- Docker 镜像为 linux/amd64，在 Apple Silicon Mac 上需通过 `--platform linux/amd64` 运行
- Nextflow 配置中需指定 `docker.runOptions = '--platform linux/amd64'`

### NFR-2: 遵循 bio.nf 模块规范
- 目录结构：modules/cresil/<subcommand>/
- 命名规范：process 名 `CRESIL_<SUBCOMMAND>`（全大写）
- 输入输出使用 meta map 传递样本信息
- 统一的 versions 输出格式
- conda 和 container 双配置

### NFR-3: 可复现性
- 锁定软件版本（cresil 1.2.0）
- 测试使用 snapshot 机制确保输出一致性
- 环境配置可复现（environment.yml + 固定镜像 tag）

## Constraints
- **技术**: CReSIL 仅支持 linux/amd64 架构（依赖 medaka、minimap2 等工具）
- **技术**: 必须使用 Nextflow DSL2
- **技术**: 遵循 bio.nf 现有模块的代码风格和命名规范
- **依赖**: Docker 镜像 quay.io/bioinfortools/cresil:1.2.0
- **依赖**: Conda 包 anaconda.org/yangmingsi/cresil
- **依赖**: 测试数据需要参考基因组 .mmi 索引文件

## Assumptions
- 用户在 Linux/amd64 或通过 Rosetta 2 的 macOS 上运行 Docker
- example 目录中的 exp_reads.fastq 可用于 trim 模块测试
- 需要额外生成小的参考基因组用于完整流程测试（或使用已有的 testdata）
- Nextflow 版本支持 nf-test 框架

## Acceptance Criteria

### AC-1: TRIM 模块正常运行
- **Given**: 有输入 FASTQ 文件和 Minimap2 索引
- **When**: 运行 CRESIL_TRIM 模块
- **Then**: 输出 trim.txt 文件，文件非空，包含预期的列数
- **Verification**: `programmatic`

### AC-2: IDENTIFY 模块正常运行
- **Given**: 有输入 FASTQ、参考基因组 FASTA+FAI、trim 结果
- **When**: 运行 CRESIL_IDENTIFY 模块
- **Then**: 输出 eccDNA_final.txt 文件
- **Verification**: `programmatic`

### AC-3: IDENTIFY_WGLS 模块正常运行
- **Given**: 有输入 FASTQ、Minimap2 索引、参考基因组 FASTA+FAI、trim 结果
- **When**: 运行 CRESIL_IDENTIFY_WGLS 模块
- **Then**: 输出 eccDNA_final.txt 文件
- **Verification**: `programmatic`

### AC-4: ANNOTATE 模块正常运行
- **Given**: 有 identify 结果表和注释 BED 文件
- **When**: 运行 CRESIL_ANNOTATE 模块
- **Then**: 输出注释结果文件
- **Verification**: `programmatic`

### AC-5: VISUALIZE 模块正常运行
- **Given**: 有 identify 结果表和指定的 eccDNA ID
- **When**: 运行 CRESIL_VISUALIZE 模块
- **Then**: 输出 Circos 配置文件目录
- **Verification**: `programmatic`

### AC-6: 模块遵循 bio.nf 规范
- **Given**: 所有模块已创建
- **When**: 检查目录结构、命名规范、元数据文件
- **Then**: 与现有 flair/longshot 等模块结构一致，包含 meta.yml、environment.yml、main.nf、tests/
- **Verification**: `human-judgment`

### AC-7: Docker 和 Conda 环境配置正确
- **Given**: 模块的 main.nf
- **When**: 检查 conda 和 container 指令
- **Then**: conda 指向 environment.yml 或 yangmingsi/cresil，container 指向 quay.io/bioinfortools/cresil:1.2.0
- **Verification**: `programmatic`

### AC-8: 测试通过
- **Given**: 测试数据已准备
- **When**: 运行 nf-test
- **Then**: 所有模块的测试用例（含 stub）通过
- **Verification**: `programmatic`

### AC-9: stub 模式可用
- **Given**: 使用 `-stub` 参数运行
- **When**: 运行模块
- **Then**: 快速生成占位输出文件，流程可正常通过
- **Verification**: `programmatic`

## Open Questions
- [ ] example 数据中是否有配套的参考基因组文件？如果没有，需要构建小的测试参考基因组
- [ ] annotate 模块的三个 BED 输入是否都作为可选输入？
- [ ] 是否需要同时支持 skip-variant 模式的 identify（即不依赖 medaka 的轻量模式）？
