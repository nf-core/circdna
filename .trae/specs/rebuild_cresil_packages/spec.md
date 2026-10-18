# Rebuild and Publish CReSIL Docker/Conda Packages - Product Requirement Document

## Overview
- **Summary**: 将修复了 `KeyError: 'name'` bug 的 CReSIL 源码版本号全部更新到 1.2.1，重新构建 Docker 镜像和 Conda 包，并发布到 quay.io 和 Anaconda。源代码版本同步更新，便于备份到 GitHub。
- **Purpose**: 将 identify_wgls 模块的 bug 修复发布到生产环境，使用户可以直接通过 Docker/Conda 使用修复后的版本；同时统一所有版本号，便于 GitHub 版本管理。
- **Target Users**: CReSIL Nextflow 模块用户、直接使用 CReSIL 的研究人员

## Goals
- 所有版本号统一更新到 1.2.1（源代码 + Docker + Conda）
- 重新构建 Docker 镜像 `quay.io/bioinfortools/cresil:1.2.1` 并推送
- 重新构建 Conda 包 `cresil=1.2.1=py38_0` 并发布到 Anaconda
- 验证构建的包包含修复后的代码
- 更新 Nextflow 模块中的容器/conda 版本引用
- 源代码版本同步，便于 GitHub 备份

## Non-Goals (Out of Scope)
- 不修改 CReSIL 的功能逻辑（仅重新构建已修复的代码）
- 不添加新功能
- 不修改其他模块（trim, identify, annotate, visualize）的业务逻辑

## Background & Context
- 已修复 CReSIL identify_wgls.py 中的 `KeyError: 'name'` bug
- 修复内容：为 pybedtools `to_dataframe()` 调用添加明确的列名参数
- 当前发布版本：1.2.0 (build 0)
- Docker 基础镜像：`conda-builder-linux64:latest`
- Conda 构建环境：`conda_build` (位于 `/opt/homebrew/Caskroom/miniforge/base/envs/conda_build`)
- Apple Silicon Mac 需要使用 `--platform linux/amd64` 构建 Docker 镜像

## Functional Requirements
- **FR-1**: 所有版本号统一从 1.2.0 更新到 1.2.1
- **FR-2**: 构建 Docker 镜像 `quay.io/bioinfortools/cresil:1.2.1`（linux/amd64 平台）
- **FR-3**: 推送 Docker 镜像到 quay.io
- **FR-4**: 构建 Conda 包 1.2.1（linux-64 平台）
- **FR-5**: 发布 Conda 包到 Anaconda (anaconda.org/yangmingsi/cresil)
- **FR-6**: 更新 Nextflow 模块中的容器和 conda 版本引用

## Non-Functional Requirements
- **NFR-1**: Docker 镜像必须兼容 linux/amd64 平台
- **NFR-2**: Conda 包必须兼容 linux-64 平台
- **NFR-3**: 构建的包必须包含修复后的 identify_wgls.py 代码
- **NFR-4**: 版本号递增遵循语义化版本规范（Semantic Versioning）- patch 版本递增

## Constraints
- **Technical**: 
  - Apple Silicon Mac 构建 linux/amd64 镜像需要 Docker Buildx 或 `--platform` 参数
  - Conda 构建必须在 linux-64 环境中进行（使用 Docker 或现有 conda_build 环境）
- **Business**: 
  - 使用现有的 quay.io 账号 (bioinfortools)
  - 使用现有的 Anaconda 账号 (yangmingsi)
- **Dependencies**: 
  - Docker daemon 运行中
  - quay.io 登录凭证
  - Anaconda 登录凭证 (anaconda login)

## Assumptions
- Docker 和 Conda 构建环境已配置完成
- quay.io 和 Anaconda 的登录凭证有效
- 修复后的代码位于 `/Users/siyangming/nextflow_nf_core/bio.nf/cresil/`
- Nextflow 模块位于 `/Users/siyangming/nextflow_nf_core/bio.nf/modules/cresil/`

## Acceptance Criteria

### AC-1: 所有版本号统一更新到 1.2.1
- **Given**: CReSIL 源码目录
- **When**: 版本号从 1.2.0 更新到 1.2.1
- **Then**: 
  - `cresil/__init__.py` 中 `__version__` 为 `'1.2.1'`
  - `Dockerfile` 中 `LABEL version` 为 `"1.2.1"`
  - `conda-recipe/meta.yaml` 中 `version` 为 `"1.2.1"`
  - `README.md` 中的版本引用更新为 1.2.1
- **Verification**: `programmatic`

### AC-2: Docker 镜像构建成功
- **Given**: 修复后的 CReSIL 源码和 Dockerfile
- **When**: 执行 Docker build 命令（linux/amd64 平台）
- **Then**: 
  - 镜像构建成功，无错误
  - 镜像标签为 `quay.io/bioinfortools/cresil:1.2.1`
  - 镜像中 `cresil --version` 输出 `cresil 1.2.1`
- **Verification**: `programmatic`

### AC-3: Docker 镜像推送成功
- **Given**: 已构建的 Docker 镜像
- **When**: 执行 docker push 到 quay.io
- **Then**: 镜像成功推送到 `quay.io/bioinfortools/cresil:1.2.1`
- **Verification**: `programmatic`

### AC-4: Conda 包构建成功
- **Given**: 修复后的 CReSIL 源码和 conda-recipe
- **When**: 执行 conda build（linux-64 平台）
- **Then**: 
  - Conda 包构建成功，无错误
  - 包文件为 `cresil-1.2.1-py38_0.conda`
  - 包包含修复后的 identify_wgls.py
- **Verification**: `programmatic`

### AC-5: Conda 包发布成功
- **Given**: 已构建的 Conda 包
- **When**: 执行 anaconda upload
- **Then**: 包成功发布到 `anaconda.org/yangmingsi/cresil`，版本为 1.2.1
- **Verification**: `programmatic`

### AC-6: Nextflow 模块版本更新
- **Given**: CReSIL Nextflow 模块
- **When**: 更新模块中的容器和 conda 版本引用
- **Then**: 
  - 所有模块的 container 引用更新为 `quay.io/bioinfortools/cresil:1.2.1`
  - Singularity 镜像路径同步更新
  - 模块 stub 测试通过
- **Verification**: `programmatic`

### AC-7: 修复验证
- **Given**: 新构建的 Docker 镜像
- **When**: 使用测试数据运行 `cresil identify_wgls`
- **Then**: 不再出现 `KeyError: 'name'` 错误
- **Verification**: `programmatic`

## Open Questions
- 无（已确认全部更新到 1.2.1）