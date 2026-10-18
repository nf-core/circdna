# 创建 Nextflow 模块 GitHub 分支 - 产品需求文档

## Overview
- **Summary**: 为 cresil 和 ecc_finder 的 Nextflow 模块分别创建独立的 GitHub 分支，只上传模块代码和 .trae 规划文档，排除源代码目录。
- **Purpose**: 分离模块代码与源代码，便于审核和后续合并。
- **Target Users**: 项目开发者和审核者

## Goals
- 为 cresil 模块创建独立分支并推送
- 为 ecc_finder 模块创建独立分支并推送
- 分支中只包含 modules/ 和 .trae/ 目录，排除 cresil/ 和 ecc_finder/ 源代码目录
- 保留项目结构完整性，便于后续 pull request

## Non-Goals (Out of Scope)
- 不合并分支到 main
- 不修改任何代码文件
- 不创建 pull request（用户审核后自行创建）
- 不处理其他目录或文件

## Background & Context
- 当前仓库状态：所有文件都是 untracked
- 远程仓库：https://github.com/SiYangming/bio.nf.git
- 需要分离上传：modules/ 和 .trae/ 目录上传，cresil/ 和 ecc_finder/ 源代码不上传

## Functional Requirements
- **FR-1**: 创建分支 feature/cresil-nf-modules
- **FR-2**: 创建分支 feature/ecc-finder-nf-modules
- **FR-3**: 在两个分支中，只添加 modules/cresil、modules/ecc_finder 和 .trae 目录
- **FR-4**: 推送两个分支到远程仓库

## Non-Functional Requirements
- **NFR-1**: 分支命名符合 Git Flow 规范（feature/*）
- **NFR-2**: 提交信息清晰描述内容
- **NFR-3**: 保留完整的目录结构

## Constraints
- **Technical**: Git 版本控制，远程仓库已配置
- **Dependencies**: 需要网络连接推送分支

## Assumptions
- 用户已配置 GitHub 认证
- 用户有远程仓库的推送权限

## Acceptance Criteria

### AC-1: cresil 分支创建并推送
- **Given**: 当前在 main 分支，远程仓库可用
- **When**: 创建 feature/cresil-nf-modules 分支并推送
- **Then**: 远程分支 feature/cresil-nf-modules 存在，包含 modules/cresil 和 .trae 目录
- **Verification**: `programmatic`

### AC-2: ecc_finder 分支创建并推送
- **Given**: 当前在 main 分支，远程仓库可用
- **When**: 创建 feature/ecc-finder-nf-modules 分支并推送
- **Then**: 远程分支 feature/ecc-finder-nf-modules 存在，包含 modules/ecc_finder 和 .trae 目录
- **Verification**: `programmatic`

### AC-3: 源代码目录不上传
- **Given**: 分支已创建
- **When**: 检查远程分支内容
- **Then**: cresil/ 和 ecc_finder/ 源代码目录不在分支中
- **Verification**: `programmatic`

## Open Questions
- [ ] 分支命名是否符合用户偏好？当前方案使用 feature/* 前缀
