# 所有流程统一 Docker 用户映射 A+B+C 方案 - Product Requirement Document

## Overview
- **Summary**: 为所有 8 个 Nextflow 流程（含 bio.nf 11 个模块测试）统一采用 **A+B+C 三层 Docker 权限方案**：A=`docker.runOptions = '-u $(id -u):$(id -g)'`（预防层）、B=`docker.userEmulation = true`（预防层·补充用户名）、C=`docker.fixOwnership = true`（兜底层）。三者在同一位置共存不冲突。
- **Purpose**: 从源头彻底避免 work 目录文件归属 root 的问题，同时提供事后补救兜底，确保任何场景（含 Nextflow 崩溃、SPAdes 隐藏中间文件竞态）都不遗留 root 文件。
- **Target Users**: 开发/维护本仓库 Nextflow 流程的所有用户。

## Goals
1. 所有流程在**同一标准位置**（`nextflow.config → profiles → docker { }`，bio.nf 为顶层 `docker { }`）同时设置 A+B+C 三个参数
2. 删除 server.config 中与标准位置重复的 A/B/C 设置（保留 `enabled = true` 以确保 `-profile server` 仍能启用 Docker）
3. 将 A+B+C 标准写入 AGENTS.md，指导后续所有流程维护
4. 完成后只同步 circdna.nf 的修改到 GitHub（其他流程的修改保留在本地，不同步）

## Non-Goals (Out of Scope)
- 不在流程外创建一个全局 Nextflow 配置（用户明确要求"在每个流程的同一位置设置"）
- 不修改 Singularity/Apptainer/Podman 等其他容器引擎配置
- 不为缺失 server.config 的流程新建 conf/server.config
- 不采用方案 D（userns-remap，需改 Docker daemon，HPC 不允许）
- 不修改 emulate_amd64 / arm / gpu profile 中的 runOptions（它们已经以 `-u` 开头，只需在 docker profile 主块加 B+C）
- 不修改 nf-core vendored 模块测试配置（`modules/modules/nf-core/*/tests/`，非我们的代码）
- 不修改 nf-core-nanoseq_3.1.0 参考副本
- 不修改 rnaseq 模块/子工作流独立测试配置（独立 nf-test 配置，非主流程重复）
- 不修改 nanoseq.nf/conf/base.conf 中 `if (params.deepvariant_gpu)` 条件块（A+GPU 组合设置，非简单重复）

## Background & Context

### Nextflow Docker 权限方案全对比

Nextflow 原生提供了**三种**解决 Docker 容器文件权限问题的机制，外加一种通过 Docker daemon 的通用机制。

---

#### 方案 A：`docker.runOptions = '-u $(id -u):$(id -g)'`（显式 runOptions）

**工作原理**：在启动 `docker run` 时把 `-u UID:GID` 作为命令行参数传入 Docker daemon，容器内主进程从一开始就以普通用户身份运行。

**作用时机**：容器启动前（预防）

**覆盖范围**：容器内所有进程（Nextflow 脚本 + 工具子进程 + SPAdes 等在子目录直接写的中间文件）

**优点**：从源头解决，不产生 root 文件；nf-core 通用规范；纯 Docker CLI 机制，跨版本稳定

**缺点**：容器内用户只是一个数字 ID，没有用户名（方案 B 补充这个）

---

#### 方案 B：`docker.userEmulation = true`

**工作原理**：Nextflow 动态生成 `/etc/passwd` 和 `/etc/group` 临时挂载文件，把当前用户的 UID/GID/用户名/组名映射进容器，同时切换容器用户到该 UID。

**作用时机**：容器启动前（预防，同方案 A 但附带 passwd/group 挂载）

**覆盖范围**：与方案 A 相同（所有进程）

**优点**：解决方案 A 中"容器内用户没有用户名"的问题（某些软件校验 `$HOME`、`$USER`）；不依赖 shell 展开 `$(id -u)`

**缺点**：依赖 Nextflow 内部实现；查问题不如 `-u` 直观

---

#### 方案 C：`docker.fixOwnership = true`

**工作原理**：事后补救机制。容器内进程仍以 root 运行，任务结束后 Nextflow 对 work 目录执行 `chown -R UID:GID`。

**作用时机**：每个任务结束之后（补救）

**覆盖范围**：仅限 work 目录的顶层和 Nextflow 知道的输出文件

**关键漏洞**：SPAdes 等工具在隐藏目录（`.bin_reads/`）中延迟写入中间文件，存在时序竞态；Nextflow 崩溃时 chown 不执行

**优点**：兼容任何 Docker 镜像；配置简单

**缺点**：非原子、有竞态、大规模 chown 很慢、不能单独依赖

**适用场景**：**只能作为辅助/兜底，不能作为唯一方案**

---

#### 方案 D（不采用）：`userns-remap`（Docker daemon 级别）

需要改 Docker daemon 配置，HPC 服务器不允许，对本项目不适用。

---

### 方案对比总表

| 维度 | A: runOptions `-u` | B: userEmulation | C: fixOwnership | D: userns-remap |
|------|:---:|:---:|:---:|:---:|
| **作用时机** | 启动前（预防） | 启动前（预防） | 结束后（补救） | daemon 级别 |
| **覆盖范围** | 所有进程写的所有文件 | 所有进程写的所有文件 | 仅 work 目录已知文件 | 宿主机全局 |
| **覆盖 SPAdes 隐藏中间文件** | ✅ | ✅ | ❌ 竞态漏网 | ✅ |
| **Nextflow 崩溃后仍安全** | ✅ | ✅ | ❌ | ✅ |
| **依赖 Nextflow 内部机制** | ❌ 纯 Docker CLI | ✅ 是 | ✅ 是 | ❌ daemon |
| **跨版本稳定性** | ✅ 最高 | ⚠️ 中 | ⚠️ 中 | ✅ 最高 |
| **本项目推荐度** | ⭐⭐⭐⭐⭐ 必选 | ⭐⭐⭐⭐ 必选 | ⭐⭐⭐⭐ 必选 | ⭐ 不适用 |

### 四方案同时使用的兼容性分析

**结论：A + B + C 可同时使用，不冲突。**

| 组合 | 冲突？ | 原因 | 本项目实例 |
|------|:---:|------|:---:|
| **A + B** | ❌ | 两者设同一 UID，Docker 取最后一个 `-u` 生效，值相同等价；B 额外挂载 passwd/group | nanoseq 已在用 |
| **A + C** | ❌ | A 让文件归属正确，C 的 chown 变成 no-op | circdna server.config 已在用 |
| **B + C** | ❌ | 同 A+C | — |
| **A + B + C** | ❌ | 预防（A+B）+ 冗余兜底（C），文件归属始终正确 | — |
| **D + A/B** | ⚠️ | UID 映射语义冲突 | 不适用 |

A+B 同时用的 Docker 命令行示意：
```
docker run -u 1000:1000 \                          ← userEmulation 注入
           -u $(id -u):$(id -g) \                   ← runOptions 注入（同值）
           -v /tmp/.passwd_xxx:/etc/passwd:ro \     ← userEmulation 注入
           <image> ...
```
Docker 对重复 `-u` 取后者覆盖，值相同所以等价，不会报错。

---

### 用户需求澄清
用户要求："不是统一设置一个（全局），而是在每个流程的同一位置设置，使其生效" → 即**逐流程、同位置、同配置**。
用户决定："所有流程应用 A+B+C 方案"。

### 标准位置

对有 `profiles {}` 块的标准 nf-core 流程（circdna.nf、circrna.nf、isoseq.nf、fetchngs.nf、riboseq.nf、nanoseq.nf、rnaseq）：
- 在 `nextflow.config` → `profiles { docker { <HERE> } }` 块内设置 A+B+C

对无 `profiles` 块的（bio.nf 顶层、各模块 tests/nextflow.config）：
- 在配置文件顶层 `docker { <HERE> }` 块内设置 A+B+C

### 各流程当前状态与所需改动

| 流程 | 标准位置 | A 当前 | B 当前 | C 当前 | 需要改动 |
|------|---------|:---:|:---:|:---:|---------|
| circdna.nf | nextflow.config L154-163 | ✅ | ❌ | ❌ | 加 B+C |
| circrna.nf | nextflow.config L164-173 | ✅ | ❌ | ❌ | 加 B+C |
| isoseq.nf | nextflow.config L128-137 | ✅ | ❌ | ❌ | 加 B+C |
| fetchngs.nf | nextflow.config L94-103 | ✅ | ❌ | ❌ | 加 B+C |
| riboseq.nf | nextflow.config L171-180 | ✅ | ❌ | ❌ | 加 B+C |
| nanoseq.nf | nextflow.config L146-153 | ❌ | ✅ | ❌ | 加 A+C |
| rnaseq | nextflow.config L184-193 | ✅ | ❌ | ❌ | 加 B+C |
| bio.nf 顶层 | nextflow.config L5-8 | ✅ | ❌ | ❌ | 加 B+C |
| bio.nf 11测试 | tests/nextflow.config L1-4 | ❌ | ❌ | ❌ | 加 A+B+C |

server.config 状态（circdna/circrna/isoseq 有）—— A/B/C 为重复设置，需删除：

| 流程 | server.config docker 块 | A | B | C | 需要改动 |
|------|------------------------|:---:|:---:|:---:|---------|
| circdna.nf | L178-182 | ✅ | ❌ | ✅ | 删除 A+C 行，保留 enabled=true |
| circrna.nf | L244-248 | ✅ | ❌ | ✅ | 删除 A+C 行，保留 enabled=true |
| isoseq.nf | L214-218 | ✅ | ❌ | ✅ | 删除 A+C 行，保留 enabled=true |

**注意**：删除 server.config 中的 A/B/C 后，`-profile server` 单独使用时仍能启用 Docker（`enabled = true` 保留），但 A+B+C 权限设置需通过 `-profile docker,server` 组合使用才能生效。这是预期行为：标准位置是 docker profile，server profile 只负责服务器资源/执行器配置。

### 非标准位置审计结果

除 server.config 外，还发现以下位置存在 A/B/C 相关设置，经审计后**不删除**：

| 位置 | 内容 | 不删除原因 |
|------|------|-----------|
| nanoseq.nf/conf/base.config L11-14 | `if (params.deepvariant_gpu) { docker.runOptions = '-u ... --gpus all' }` | A+GPU 组合设置，非简单重复；删除会丢失 `--gpus all` |
| rnaseq 5 个模块/子工作流测试配置 | `docker.runOptions = '-u $(id -u):$(id -g)'` | 独立 nf-test 配置，非主流程重复 |
| isoseq 9 个测试配置 | `runOptions = '--platform linux/amd64'` | 仅平台模拟参数，不含 A/B/C |
| fetchngs 1 个模块测试 | `runOptions = '--platform linux/amd64'` | 同上 |
| modules/modules/nf-core/*/tests/ | 各种 runOptions/fixOwnership | nf-core vendored 代码，非我们的文件 |
| nf-core-nanoseq_3.1.0/ | 各种 | 参考副本，不修改 |
| emulate_amd64/arm/gpu profile | A + `--platform`/`--gpus` | 方案内的扩展写法，非重复 |

## Functional Requirements
- **FR-1**: 7 个标准流程（circdna/circrna/isoseq/fetchngs/riboseq/nanoseq/rnaseq）的 `nextflow.config → profiles → docker { }` 块必须同时包含 A（`docker.runOptions = '-u $(id -u):$(id -g)'`）、B（`docker.userEmulation = true`）、C（`docker.fixOwnership = true`）。
- **FR-2**: bio.nf 顶层 `nextflow.config` 的 `docker { }` 块必须包含 B 和 C（A 已有，含 `--platform`）。
- **FR-3**: bio.nf 11 个模块测试配置文件的 `docker { }` 块必须包含 A+B+C（runOptions 改为 `'-u $(id -u):$(id -g) --platform linux/amd64'`，加 B+C）。
- **FR-4**: 3 个 server.config（circdna/circrna/isoseq）的 `docker { }` 块必须删除与标准位置重复的 A（`runOptions`）和 C（`fixOwnership`）行，仅保留 `enabled = true`。
- **FR-5**: 根目录 AGENTS.md 新增"Docker 用户映射标准设置规范"章节，明确 A+B+C 标准位置和配置行。

## Non-Functional Requirements
- **NFR-1 (一致性)**: 所有流程在标准位置的 docker 配置块必须按统一顺序包含 A+B+C 三行。
- **NFR-2 (最小改动)**: 只添加缺失的参数行，不修改已有正确配置行，不改变其他配置项。
- **NFR-3 (向后兼容)**: 不删除任何已有配置，不改变 emulate_amd64/arm/gpu profile。

## Constraints
- **Technical**: 必须使用 Nextflow 原生 docker 配置项，不引入外部脚本。
- **Dependencies**: 无新依赖。
- **Business**: 改动不得引入破坏性变更。

## Assumptions
- 用户在 Mac/Linux 上运行 Docker，`$(id -u)` / `$(id -g)` 能被 shell 正确解析。
- A+B+C 三者共存不冲突（已在前文兼容性分析中验证）。

## Acceptance Criteria

### AC-1: 7 个标准流程 nextflow.config docker profile 含 A+B+C
- **Given**: circdna/circrna/isoseq/fetchngs/riboseq/nanoseq/rnaseq 的 `nextflow.config → profiles → docker { }` 块
- **When**: 读取该块内容
- **Then**: 同时存在 `docker.runOptions = '-u $(id -u):$(id -g)'`、`docker.userEmulation = true`、`docker.fixOwnership = true`
- **Verification**: `programmatic`（Grep 三个参数全部匹配）

### AC-2: bio.nf 顶层 + 11 模块测试含 A+B+C
- **Given**: bio.nf/nextflow.config 顶层 docker 块 + 11 个 tests/nextflow.config 的 docker 块
- **When**: 读取 docker 块内容
- **Then**: 顶层含 runOptions（`-u` + `--platform`）+ userEmulation + fixOwnership；11 个测试含相同三行
- **Verification**: `programmatic`（Grep 12 个文件全部匹配）

### AC-3: 3 个 server.config 删除重复 A/B/C
- **Given**: circdna/circrna/isoseq 的 conf/server.config 的 docker 块
- **When**: 读取该块内容
- **Then**: 仅包含 `enabled = true`，不再包含 `runOptions`、`fixOwnership`、`userEmulation`
- **Verification**: `programmatic`（Grep 确认三个关键词不存在于 server.config）

### AC-4: AGENTS.md 新增 A+B+C 标准规范
- **Given**: 根目录 AGENTS.md
- **When**: 搜索"Docker 用户映射标准设置规范"
- **Then**: 该节存在且明确说明 A+B+C 三层方案、标准位置、标准配置行、兼容性说明
- **Verification**: `human-judgment`（审阅该节内容）

### AC-5: emulate_amd64/arm/gpu profile 未被修改
- **Given**: 各流程 nextflow.config 中的 emulate_amd64 / arm / gpu profile
- **When**: 执行 git diff
- **Then**: 这些 profile 的 runOptions 行无改动
- **Verification**: `programmatic`（git diff 验证）

## Open Questions
- 无
