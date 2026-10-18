# 为所有流程添加 Docker 用户映射权限配置 - Product Requirement Document

## Overview
- **Summary**: 为当前目录下 8 个 Nextflow 流程中缺少 Docker 用户映射参数的配置补充 `-u $(id -u):$(id -g)`，确保容器内进程始终以宿主机用户身份运行，从源头避免 work 目录文件归属 root 导致无法删除的问题。
- **Purpose**: 解决 SPAdes 等工具在 Docker 容器内以 root 身份运行生成的文件无法用普通用户删除的问题，确保所有运行场景（本地测试、服务器运行、模块测试）都生效。
- **Target Users**: 在本地或服务器上运行 Nextflow 流程的用户。

## Goals
1. 所有 8 个流程在任何 Docker 运行场景下都以宿主机用户身份执行
2. bio.nf 的 11 个模块测试配置补充 `-u` 参数
3. nanoseq.nf 的主 docker profile 添加显式 `docker.runOptions`
4. 以最小改动实现最大覆盖，采用最简单、最直接的方案

## Non-Goals (Out of Scope)
- 不为缺失 server.config 的流程（nanoseq.nf、riboseq.nf、fetchngs.nf、rnaseq、bio.nf）新建 server.config（超出"最简单方案"范围，docker profile 已覆盖所有运行场景）
- 不修改 circdna.nf、circrna.nf、isoseq.nf、fetchngs.nf、riboseq.nf、rnaseq 的 nextflow.config docker profile（已有正确配置）
- 不修改 Singularity/Apptainer 配置（仅处理 Docker）

## Background & Context
### 问题回顾
服务器上运行 circdna.nf 时，`rm -rf work/` 失败，提示"权限不够"。根本原因是 Docker 容器默认以 root 身份运行，SPAdes 等工具直接创建的中间文件（如 `.bin_reads/`）归属于 root，普通用户无法删除。

### 解决方案对比
| 方案 | 覆盖范围 | 改动量 | 生效时机 | 推荐 |
|------|---------|--------|---------|:---:|
| **A: 在 nextflow.config 的 docker profile 中添加 `-u`** | 所有 `docker,test_local` 等 profile 运行 | 小 | 运行时注入 | ✅ 最简单 |
| B: 为每个流程创建 conf/server.config | 仅 `-profile server` | 大 | 服务器专用 | ❌ 不覆盖本地测试 |
| C: 在每个 modules.config 中逐 process 配置 | 所有 process | 极大 | 逐进程 | ❌ 维护成本高 |
| D: 依赖 Nextflow `fixOwnership=true` | work 目录根级 | 中 | 任务结束后 | ❌ 不覆盖容器内直接写入的中间文件 |

**结论**：方案 A（在 docker profile 中统一设置 `docker.runOptions`）是最简单、覆盖最全的方案，在本次对话之前已在 6 个流程中正确使用，只需补充 bio.nf 模块测试和 nanoseq.nf。

### 当前配置状态（审计于 2026-07-31）
| # | 流程 | nextflow.config docker profile | bio.nf 模块测试 | 备注 |
|---|------|:---:|:---:|---|
| 1 | circdna.nf | ✅ 已有 | N/A | |
| 2 | circrna.nf | ✅ 已有 | N/A | |
| 3 | bio.nf | ✅ 已有 | ❌ 11个缺失 | 模块测试只设 `--platform linux/amd64` |
| 4 | nanoseq.nf | ⚠️ 仅 userEmulation | N/A | 无显式 runOptions |
| 5 | isoseq.nf | ✅ 已有 | N/A | |
| 6 | fetchngs.nf | ✅ 已有 | N/A | |
| 7 | riboseq.nf | ✅ 已有 | N/A | |
| 8 | rnaseq | ✅ 已有 | N/A | |

## Functional Requirements
- **FR-1**: bio.nf 的所有模块测试配置文件（11 个）中 `docker.runOptions` 必须包含 `-u $(id -u):$(id -g)`，原有 `--platform linux/amd64` 保留。
- **FR-2**: nanoseq.nf 的 `nextflow.config` 的 docker profile 必须包含显式 `docker.runOptions = '-u $(id -u):$(id -g)'`，原有 `docker.userEmulation = true` 保留（双重保险）。

## Non-Functional Requirements
- **NFR-1 (向后兼容性)**: 任何修改不得改变现有流程的功能、输出格式或参数接口。
- **NFR-2 (最小改动)**: 仅修改确实缺失配置的文件，不修改已正确配置的 6 个流程。
- **NFR-3 (一致性)**: bio.nf 所有模块测试的 runOptions 格式必须一致（`-u $(id -u):$(id -g)` + 原有参数）。

## Constraints
- **Technical**: 必须使用 Nextflow 原生的 `docker.runOptions` 机制，不引入外部脚本或包装层。
- **Dependencies**: 无新依赖，仅修改现有配置文件。
- **Business**: 所有改动必须保持向后兼容，不得破坏现有运行方式。

## Assumptions
- 用户在 Mac/Linux 宿主机上运行，`$(id -u)` 和 `$(id -g)` 能正确解析为当前用户的 UID/GID。
- Docker daemon 已正确配置，允许使用 `-u` 参数（标准 Docker 默认允许）。
- circdna.nf、circrna.nf、isoseq.nf、fetchngs.nf、riboseq.nf、rnaseq 的 docker profile 已正确配置，无需验证（基于上一轮审计结果）。

## Acceptance Criteria

### AC-1: bio.nf 11 个模块测试配置均包含 `-u`
- **Given**: bio.nf/modules/ 下存在 11 个 `*/tests/nextflow.config` 文件且均有 `runOptions` 行
- **When**: 逐一读取这些文件的 runOptions 值
- **Then**: 每个 runOptions 都包含 `-u $(id -u):$(id -g) --platform linux/amd64`
- **Verification**: `programmatic`
- **Notes**: 原有 `--platform linux/amd64` 必须保留，不能丢失

### AC-2: nanoseq.nf docker profile 含显式 runOptions
- **Given**: nanoseq.nf/nextflow.config 存在 docker profile 块（L146-153）
- **When**: 读取 docker profile 的内容
- **Then**: 包含 `docker.runOptions = '-u $(id -u):$(id -g)'` 且 `docker.userEmulation = true` 仍保留
- **Verification**: `programmatic`

### AC-3: 不修改已正确配置的流程
- **Given**: circdna.nf, circrna.nf, isoseq.nf, fetchngs.nf, riboseq.nf, rnaseq 的 nextflow.config
- **When**: 检查 git status 中这 6 个流程的配置文件
- **Then**: 这些文件未发生任何改动
- **Verification**: `programmatic`（通过 git diff 验证）

### AC-4: 改动最小且一致
- **Given**: 本 PRD 涉及的所有改动文件
- **When**: 审阅改动内容
- **Then**: bio.nf 模块测试的改动完全一致（都是在 `--platform` 前插入 `-u`），nanoseq.nf 只增加一行 runOptions
- **Verification**: `human-judgment`
- **Notes**: 提供改动摘要列表和批量一致性验证脚本

## Open Questions
- 无（方案已明确，改动范围可控）
