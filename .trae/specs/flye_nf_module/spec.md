# Flye Nextflow模块 - 产品需求文档

## Overview
- **Summary**: 构建一个符合nf-core标准格式的Flye Nextflow模块，用于从头组装单分子长读长测序数据（PacBio和Oxford Nanopore）。模块接受FASTA/FASTQ格式的长读长序列作为输入，输出组装后的contigs序列文件及相关组装信息。
- **Purpose**: 为bio.nf项目提供标准化的长读长从头组装模块，支持多种长读长数据类型（PacBio CLR/HiFi、ONT raw/HQ），与现有nf-core模块格式保持一致，便于在各种基因组组装流程中复用。
- **Target Users**: 生物信息学研究人员、流程开发人员，需要进行长读长基因组组装并整合到Nextflow分析流程中。

## Goals
- 构建符合nf-core标准的flye模块（main.nf, meta.yml, environment.yml, tests/）
- 支持多种长读长数据类型输入（pacbio-raw, pacbio-hifi, nano-raw, nano-hq等）
- 输出assembly.fasta、assembly_graph.gfa、assembly_info.txt等关键文件
- 支持conda和docker两种运行环境
- 支持task.ext.args自定义参数传递
- 支持task.ext.prefix自定义输出前缀
- 支持stub测试模式
- 提供完整的nf-test测试用例
- 测试数据从GitHub Flye仓库的示例数据获取
- Docker测试配置--platform linux/amd64
- 在flye分支上开发

## Non-Goals (Out of Scope)
- 不实现多种输入类型的自动检测逻辑（由流程层控制）
- 不实现组装质量评估（如QUAST分析）
- 不构建Flye软件本身，仅封装现有工具
- 不提供GUI界面，仅提供Nextflow流程模块
- 不处理短读长数据的组装

## Background & Context
- bio.nf项目已维护多个nf-core格式的模块（cresil, ecc_finder, fastqdl, kingfisher, minimap2, samtools等）
- 现有的minimap2模块是类似的长读长比对工具，可作为结构参考
- Flye v2.9.6是最新稳定版，已在bioconda和biocontainers发布
- 注意：用户指定的conda版本是4.0.1，但Flye当前最新版本是2.9.6，Docker镜像也是2.9.6版本，以Docker镜像版本为准（2.9.6）
- Flye输出文件包括：assembly.fasta、assembly_graph.gfa、assembly_graph.gv、assembly_info.txt、flye.log
- Docker测试需要配置--platform linux/amd64以兼容Apple Silicon

## Functional Requirements
- **FR-1**: 模块接受Groovy Map meta和reads文件列表作为输入
- **FR-2**: 模块调用flye命令进行长读长从头组装
- **FR-3**: 输出组装的contig序列文件（assembly.fasta）
- **FR-4**: 输出组装图文件（assembly_graph.gfa）
- **FR-5**: 输出组装信息文件（assembly_info.txt）
- **FR-6**: 输出版本信息文件（versions.yml）
- **FR-7**: 支持通过task.ext.args传递flye额外参数
- **FR-8**: 支持通过task.ext.prefix自定义输出文件名前缀
- **FR-9**: 支持stub运行模式，生成模拟输出文件
- **FR-10**: 同时支持conda环境和Docker容器运行
- **FR-11**: 支持指定长读长数据类型（通过args参数控制，如--pacbio-hifi, --nano-raw等）

## Non-Functional Requirements
- **NFR-1**: 模块命名规范遵循nf-core标准，process名为FLYE
- **NFR-2**: 模块目录结构与现有nf-core模块保持一致（modules/flye/）
- **NFR-3**: meta.yml文档完整描述输入输出参数
- **NFR-4**: 测试用例覆盖正常运行和stub模式
- **NFR-5**: Docker镜像使用quay.io/biocontainers/flye:2.9.6--py313h7fbb527_1
- **NFR-6**: conda环境使用bioconda::flye=2.9.6（与Docker版本一致）

## Constraints
- **Technical**: 必须使用Nextflow DSL2语法，遵循nf-core模块规范
- **Business**: 模块需与bio.nf项目中现有模块风格一致
- **Dependencies**:
  - Flye v2.9.6 (bioconda::flye=2.9.6)
  - Docker镜像: quay.io/biocontainers/flye:2.9.6--py313h7fbb527_1
  - 输入数据：长读长FASTA/FASTQ文件

## Assumptions
- 用户提供的reads文件格式正确（FASTA或FASTQ，可压缩）
- 用户通过task.ext.args参数指定读长类型（--pacbio-raw, --nano-raw, --pacbio-hifi, --nano-hq等），与nf-core风格一致
- stub模式下测试不需要真实组装数据
- 测试数据使用E.coli小型长读长示例数据集（从Zenodo/GitHub获取），确保测试快速完成
- Flye的--out-dir参数设为当前工作目录
- Docker测试需使用--platform linux/amd64以确保在Apple Silicon Mac上正常运行
- 版本以Docker镜像版本2.9.6为准（与bioconda保持一致）

## Acceptance Criteria

### AC-1: 模块文件结构完整
- **Given**: 模块目录modules/flye/已创建
- **When**: 检查目录内容
- **Then**: 目录下包含main.nf, meta.yml, environment.yml, tests/main.nf.test, tests/nextflow.config文件
- **Verification**: `programmatic`
- **Notes**: 与现有模块结构一致

### AC-2: process命名和标签正确
- **Given**: main.nf文件已创建
- **When**: 检查process定义
- **Then**: process名为FLYE，tag为$meta.id，label为process_high
- **Verification**: `programmatic`

### AC-3: conda和docker环境配置正确
- **Given**: main.nf和environment.yml已创建
- **When**: 检查环境配置
- **Then**: conda使用environment.yml文件（flye=2.9.6），docker使用quay.io/biocontainers/flye:2.9.6--py313h7fbb527_1，singularity使用galaxy镜像
- **Verification**: `programmatic`

### AC-4: 输入通道定义正确
- **Given**: main.nf文件已创建
- **When**: 检查input块
- **Then**: 接受tuple val(meta), path(reads)作为输入，reads为文件或文件列表
- **Verification**: `programmatic`

### AC-5: 输出通道定义正确
- **Given**: main.nf文件已创建
- **When**: 检查output块
- **Then**: 输出assembly通道（assembly.fasta）、gfa通道（assembly_graph.gfa）、info通道（assembly_info.txt）、versions通道（versions.yml）
- **Verification**: `programmatic`

### AC-6: 支持task.ext.args和task.ext.prefix
- **Given**: main.nf文件已创建
- **When**: 检查script块
- **Then**: 使用task.ext.args传递额外参数，使用task.ext.prefix定义输出前缀（默认使用meta.id）
- **Verification**: `programmatic`

### AC-7: stub模式正常工作
- **Given**: main.nf文件已创建
- **When**: 使用-stub参数运行
- **Then**: 生成模拟的assembly.fasta、assembly_graph.gfa、assembly_info.txt和versions.yml，不执行真实组装
- **Verification**: `programmatic`

### AC-8: meta.yml文档完整
- **Given**: meta.yml文件已创建
- **When**: 检查内容
- **Then**: 包含name, description, keywords, tools, input, output, authors等完整字段
- **Verification**: `human-judgment`

### AC-9: nf-test测试用例完整
- **Given**: tests/main.nf.test已创建
- **When**: 检查测试用例
- **Then**: 至少包含正常模式和stub模式两个测试用例
- **Verification**: `programmatic`

### AC-10: stub模式测试通过
- **Given**: 模块和测试用例已创建
- **When**: 运行nf-test stub模式测试
- **Then**: 测试成功通过，snapshot匹配
- **Verification**: `programmatic`

### AC-11: Docker测试配置正确
- **Given**: tests/nextflow.config已创建
- **When**: 检查Docker配置
- **Then**: 配置docker.enabled=true和runOptions='--platform linux/amd64'
- **Verification**: `programmatic`

## Open Questions
- 无（已全部确认：通过args控制读长类型，使用E.coli小型FASTA示例数据）
