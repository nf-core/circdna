# 移除已废弃的 docker.userEmulation 配置 - Spec

## Why

Nextflow 26.04.6 已移除 `docker.userEmulation` 配置项支持。当前所有流程运行时均产生两条警告：

```
WARN: Unrecognized config option 'docker.userEmulation'
WARN: Config setting `docker.userEmulation` is not supported anymore
```

该配置项（原 A+B+C 三层方案中的 B 层）已无实际作用，仅产生噪音警告。A 层（`-u $(id -u):$(id -g)`）和 C 层（`fixOwnership = true`）已足够保证容器内文件归属正确——用户在 circdna.nf 的 61 任务实测中已验证此点（全部成功，无权限问题）。

## What Changes

- **移除** 所有 20 个活跃配置文件中的 `userEmulation = true` / `docker.userEmulation = true` 行
- **保留** A 层 `docker.runOptions = '-u $(id -u):$(id -g)'`（不动）
- **保留** C 层 `docker.fixOwnership = true`（不动）
- **不修改** `nf-core-nanoseq_3.1.0/3_1_0/` 参考归档目录（已下载的原始参考代码，不用于实际运行）
- **不修改** `modules/modules/nf-core/` 下的第三方模块测试配置（这些是 nf-core 官方模块的原始文件，且均未使用 `userEmulation`）
- **更新** `/Users/siyangming/nextflow_nf_core/AGENTS.md` 第 11 节：将 A+B+C 三层方案修订为 A+C 两层方案，移除 B 层描述

## Impact

- **Affected specs**: `bio.nf/.trae/specs/fix-all-docker-user-mapping`（原 spec 的 Task 2 和 AC-2 明确保留 `userEmulation`，现已过时）
- **Affected code**: 20 个配置文件 + AGENTS.md
  - 8 个主流程 `nextflow.config`：circdna.nf、circrna.nf、nanoseq.nf、isoseq.nf、fetchngs.nf、riboseq.nf、rnaseq、bio.nf
  - 12 个 bio.nf 模块测试配置：flye、cresil(5个)、ecc_finder(4个)、eccsplorer、fastqdl
  - AGENTS.md Section 11（Docker 用户映射标准规范）

## ADDED Requirements

无新增需求。

## MODIFIED Requirements

### Requirement: Docker 用户映射配置

所有 Nextflow 流程的 Docker 配置 SHALL 采用 **A+C 两层方案**：

| 层级 | 参数 | 作用 | 状态 |
|------|------|------|------|
| **A（预防层）** | `docker.runOptions = '-u $(id -u):$(id -g)'` | 显式将宿主机 UID/GID 传入容器 | ✅ 保留 |
| ~~B（预防层补充）~~ | ~~`docker.userEmulation = true`~~ | ~~Nextflow 动态生成 /etc/passwd 挂载~~ | ❌ **已废弃，移除** |
| **C（兜底层）** | `docker.fixOwnership = true` | 任务结束后 chown 修复文件所有权 | ✅ 保留 |

#### Scenario: 运行流程时无 userEmulation 警告
- **WHEN** 用户执行 `nextflow run main.nf -profile test_local,docker`
- **THEN** 控制台输出中不出现 `WARN: Unrecognized config option 'docker.userEmulation'` 或 `WARN: Config setting 'docker.userEmulation' is not supported anymore`

#### Scenario: 容器内文件归属仍正确
- **WHEN** Docker 任务执行完毕
- **THEN** work 目录下的文件归属为当前用户（非 root），`rm -rf work/` 可成功执行
- **BECAUSE** A 层 `-u $(id -u):$(id -g)` 在容器启动时即以宿主机用户身份运行，C 层 `fixOwnership` 在任务结束后兜底修复

#### Scenario: 需要 /etc/passwd 的工具仍可工作
- **WHEN** 某些工具（如 GATK4）在容器内需要查询 /etc/passwd
- **THEN** 该工具的测试配置已通过 `runOptions` 手动挂载 `/etc/passwd`、`/etc/group` 等文件（参见 `modules/modules/nf-core/gatk4/*/tests/nextflow.config` 的现有做法）
- **AND** 全局 docker profile 不需要为此额外配置

## REMOVED Requirements

### Requirement: docker.userEmulation 双重保险

**Reason**: Nextflow 26.04.6 已移除该配置项支持，设置后仅产生警告无实际效果。A 层 + C 层已充分覆盖用户映射需求。
**Migration**: 直接删除所有配置文件中的 `userEmulation = true` 和 `docker.userEmulation = true` 行，无需替代配置。
