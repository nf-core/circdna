# ECCsplorer 集成策略对比分析：Docker 合并 vs 模块拆分

> 分析日期：2026-08-06
> 分析范围：`circdna.nf` 中 ECCsplorer 运行环境构建策略
> 更新：根据 `kavonrtep/repeatexplorer:2.3.8` 调研结果修正分析

---

## 一、背景与现状

### 1.1 ECCsplorer 功能模式

| 模式  | 功能                            | 是否需要 RepeatExplorer2 | 当前状态         |
| ----- | ------------------------------- | ------------------------ | ---------------- |
| `map` | segemehl 比对 + eccDNA 检测     | **不需要**               | 已实现，生产可用 |
| `clu` | RepeatExplorer2 聚类 + 候选筛选 | **需要** seqclust        | 骨架占位，待实现 |
| `all` | map + clu 全流程                | **需要** seqclust        | 骨架占位，待实现 |

### 1.2 运行环境现状

**A. 现有 ECCsplorer 镜像**：`quay.io/bioinfortools/eccsplorer:2022.01.1.1`

- 基于 `condaforge/mambaforge`，当前仅包含 ECCsplorer map 模式依赖
- RepeatExplorer2/seqclust：Dockerfile 尝试从 `repeatexplorer.org` 下载但 URL 已 404，**镜像中实际无 seqclust**
- 有优雅降级：下载失败不报错，map 模式正常；clu 模式运行时报错

**B. 现有 Conda 包**：`yangmingsi::eccsplorer=2022.01.1.1`

- `noarch: python`，run 依赖 ~40 个包
- **不含 RepeatExplorer2**（meta.yaml 注释明确注明）

**C. RepeatExplorer2 官方 Docker 镜像**：`kavonrtep/repeatexplorer:2.3.8`

> **重要更正**：此前误判为"第三方不可控镜像"。经调研确认：
>
> - `kavonrtep` = **Petr Novak**，RepeatExplorer2 的**主要作者**（Biology Centre CAS, Czech Republic）
> - [GitHub repo](https://github.com/kavonrtep/repex_tarean) 是官方镜像（2023-11-30 从 Bitbucket 迁移），有 4 stars、维护活跃
> - [Docker Hub](https://hub.docker.com/r/kavonrtep/repeatexplorer) 是开发者官方维护的镜像
> - 文档明确："Singularity is preferred option but RepeatExplorer2 can be also run using Docker container"
> - `seqclust` 位于容器内 `/repex_tarean/seqclust`（非 PATH），需完整路径调用

### 1.3 当前架构拆分程度（已较高）

```
modules/local/
├── eccsplorer/                          # map 模式（已实现）
├── eccsplorer_clu_prepare/              # clu 准备（骨架）
├── eccsplorer_clu_core/                 # clu 核心（骨架，将调用 seqclust）
├── eccsplorer_clu_candidates_plot/      # clu 可视化（骨架）
├── eccsplorer_prepare_read_length/      # eccPrepare（已实现）
├── eccsplorer_prepare_read_count/       # eccPrepare（已实现）
└── eccsplorer_prepare_prexing/          # eccPrepare（已实现）

subworkflows/local/
├── eccdna_mode/                         # map 模式主工作流
├── eccsplorer_cluster/                  # clu 串联骨架
├── eccsplorer_all/                      # all 模式骨架
└── ecc_preprocessing/                   # 预处理链（已实现）
```

---

## 二、关键问题回答

### Q1: `kavonrtep/repeatexplorer:2.3.8` 能否与 ECCsplorer 合并实现完整功能？

**技术上是可行的**，有两种合并策略：

#### 策略 A：以 ECCsplorer 镜像为底，加入 RepeatExplorer2

在当前 `eccsplorer:2022.01.1.1` 的 Dockerfile 中替换 RepeatExplorer2 安装步骤，从 Bitbucket/GitHub 源码编译安装：

```dockerfile
# 替换原来的 wget repeatexplorer.org (404) 为：
RUN git clone https://github.com/kavonrtep/repex_tarean.git /opt/repex_tarean
RUN cd /opt/repex_tarean && make
RUN ln -sf /opt/repex_tarean/seqclust /opt/conda/envs/eccsplorer/bin/seqclust
```

**可行性**：高。ECCsplorer 镜像已包含 seqclust 所需的 R 生态依赖。`make` 仅编译 3 个 C++ Louvain 二进制（gcc 已有），32-bit 依赖（`formatdb`/`mgblast`）可通过 conda-forge 的 `blast-legacy` 包提供 64-bit 版本替代，无需安装 `libc6:i386`。

#### 策略 B：以 repeatexplorer 镜像为底，加入 ECCsplorer

基于 `kavonrtep/repeatexplorer:2.3.8`，叠加安装 ECCsplorer 的 Python 依赖和源码。

**可行性**：中等。需要了解 repeatexplorer 镜像的 base image 和内部结构（可能是 Ubuntu + conda）。若两边的 conda 环境冲突，需要处理。

#### 合并后的完整功能验证

合并后理论上能覆盖 ECCsplorer 的三个模式：

- **map**：不需要 seqclust，纯 ECCsplorer.py 运行，现有镜像已支持
- **clu**：需要 `ECCSPLORER_CLU_PREPARE`（纯 Python） + `ECCSPLORER_CLU_CORE`（调用 seqclust） + `ECCSPLORER_CLU_CANDIDATES_PLOT`（R 可视化）
- **all**：map + clu 串联

关键验证点：ECCsplorer 的 `lib/eccClusterer.py` 中通过 `subprocess.Popen` 调用 `seqclust`，需要确保：

1. `seqclust` 在 PATH 中或路径正确
2. ECCsplorer 的 `lib/config.py` 中 `CONDA_ENV` 设置正确（已通过 sed 修复）
3. 所需的分类数据库文件（REXdb）存在

**结论：合并可行，但需要解决数据库文件、32-bit 依赖等细节问题。**

---

### Q2: RepeatExplorer2 能否自己构建 conda 包并上传？

**答案：理论上可行，但实际工程代价过高，且存在体积障碍。**

#### 本地源码调研关键发现（`repex_tarean-master/`）

通过 `file` 命令分析 `bin/` 目录的预编译二进制，精确识别了依赖状况：

| 文件                         | 类型                         | 32/64-bit  |
| ---------------------------- | ---------------------------- | ---------- |
| `seqclust`                   | **Python 3** 脚本（非 Perl） | N/A        |
| `cap3`                       | ELF 可执行文件               | **64-bit** |
| `runOGDFlayout`              | ELF 可执行文件               | **64-bit** |
| `formatdb`                   | ELF 可执行文件               | **32-bit** |
| `mgblast`                    | ELF 可执行文件               | **32-bit** |
| `blastplus_wrapper.py`       | Python 脚本                  | N/A        |
| `last_wrapper.py`            | Python 脚本                  | N/A        |
| `align_parsing.pl`           | Perl 脚本                    | N/A        |
| `select_and_sort_contigs.pl` | Perl 脚本                    | N/A        |

> **重要更正**：此前误判 `seqclust` 为 Perl 脚本。实际是 `#!/usr/bin/env python3`，核心聚类逻辑为 Python + 编译后的 C++ Louvain 算法。

此外：

- `environment.yml` 中已声明 `blast-legacy` 作为 conda 依赖——该包在 conda-forge 中提供 **64-bit** 版的 `formatdb` 和 `mgblast`，可替代 `bin/` 中的 32-bit 预编译版本
- `databases/` 目录大小 **253MB**，含预构建的 BLAST 索引和 R 分类树文件

#### 有利因素

1. `environment.yml` 声明了完整 conda 依赖链（python 3.7, R 生态, blast, diamond, last 等）
2. `compilers` 包（conda-forge）提供 gcc/g++，可编译 `louvain/` C++ 代码
3. `blast-legacy` 包提供 64-bit `formatdb`/`mgblast`，**可消除 32-bit 依赖**
4. 源码托管在 GitHub，维护活跃
5. 上游开发者 Petr Novak 支持 conda 安装方式

#### 实际障碍

| 障碍              | 详情                                                                                      |
| ----------------- | ----------------------------------------------------------------------------------------- |
| **253MB 数据库**  | REXdb 分类数据库（含预构建 BLAST 索引），远超 conda 包合理体积（通常 <100MB）             |
| **编译步骤**      | `make` 需编译 3 个 C++ Louvain 二进制，不能声明 `noarch`，需为 linux-64/osx-64 分别构建   |
| **Perl 辅助脚本** | `bin/` 中 2 个 Perl 脚本依赖 BioPerl 等系统 Perl 模块                                     |
| **bioconda 审核** | RepeatExplorer2 涉及 Python + C++ + R + Perl 多语言混装，审核复杂度远高于普通 Python/R 包 |
| **维护成本**      | 上游更新频繁（CHANGELOG 显示活跃开发），需要持续同步 conda recipe                         |

#### 务实替代方案

既然 `kavonrtep/repeatexplorer:2.3.8` 已是开发者官方 Docker 镜像，**最务实的做法不是构建 conda 包，而是：**

1. **Docker 方式**：直接使用官方 `kavonrtep/repeatexplorer:2.3.8` 作为 `clu_core` 模块的容器
2. **自建 Docker**：基于现有 `eccsplorer` Dockerfile，`git clone` + `make` 构建，推送至 `quay.io/bioinfortools/`（推荐，保持命名空间一致）
3. **conda 混合**：map 模式用 conda 包，clu 模式用 Docker 兜底

---

## 三、两种思路的重新对比

### 思路 1：合并为单一镜像

> 将 `kavonrtep/repeatexplorer:2.3.8`（官方）与现有 `eccsplorer:2022.01.1.1` 合并为一个镜像。

**更新后的分析**：

| 维度                 | 评估（修正后）                                                                       |
| -------------------- | ------------------------------------------------------------------------------------ |
| 技术可行性           | **高**。kavonrtep 是官方维护，GitHub 有完整源码，`make` 编译即可，比之前认为的更容易 |
| 构建复杂度           | **中等**（原"高"→降低为"中"）。github clone + make，而非之前猜测的"两个陌生环境融合" |
| 镜像体积             | 仍大（~5GB+），但不比两个独立镜像总和更差                                            |
| 维护成本             | 仍需关注两方的独立更新节奏                                                           |
| kavonrtep 镜像可控性 | **可控**（原"不可控"→修正）。官方镜像 + GitHub 源码，有明确更新历史和许可证          |

**核心矛盾不变**：

- **仍然违背最小粒度原则**：map 模式不需要 RepeatExplorer2，但每次加载完整大镜像
- **conda 包仍不可行**：合并镜像的思路无法解决 RepeatExplorer2 无 conda 包的根本问题

### 思路 2：彻底拆分（多镜像）

> ECCSPLORER (map) 用现有轻量镜像，ECCSPLORER_CLU_CORE 用 repeatexplorer 官方镜像或自建 eccsplorer-repex 镜像。

**更新后的分析**：

| 维度                      | 评估（修正后）                                                                                                 |
| ------------------------- | -------------------------------------------------------------------------------------------------------------- |
| clu_core 镜像来源         | 可直接使用 `kavonrtep/repeatexplorer:2.3.8` 官方镜像（`seqclust` 在 `/repex_tarean/seqclust`），无需等自建镜像 |
| 或自建 `eccsplorer-repex` | 基于 GitHub 源码 + 现有 Dockerfile 构建，保证 `bioinfortools` 命名空间一致                                     |
| map 模式零影响            | 继续用 `quay.io/bioinfortools/eccsplorer:2022.01.1.1`                                                          |

---

## 四、更新后的核心对比

| 对比维度             | 思路 1：Docker 合并                | 思路 2：模块拆分（多镜像）                         | 优胜方     |
| -------------------- | ---------------------------------- | -------------------------------------------------- | ---------- |
| **技术可行性**       | 高（github clone + make）          | 高（clu_core 可直接用官方镜像）                    | **平手**   |
| **map 模式镜像体积** | 大（~5GB）                         | 小（~2GB）                                         | **思路 2** |
| **clu 模式镜像来源** | 需自建合并镜像                     | 可直接用 `kavonrtep/repeatexplorer:2.3.8` 官方镜像 | **思路 2** |
| **conda 支持**       | 不可行（RepeatExplorer2 无法打包） | map 支持 conda，clu-core 用 Docker                 | **思路 2** |
| **符合项目规范**     | 违背（最小粒度原则）               | 完全符合                                           | **思路 2** |
| **与现有架构兼容**   | 破坏（需回退拆分）                 | 复用 80%                                           | **思路 2** |
| **维护成本**         | 中（合并 = 耦合）                  | 低（官方镜像独立更新）                             | **思路 2** |
| **更新灵活性**       | 低                                 | 高                                                 | **思路 2** |
| **部署复杂度**       | 单镜像                             | 两镜像                                             | **思路 1** |

---

## 五、结论

**思路 2（模块拆分 + 多镜像）仍然是更优策略**，且研究修正后方案更加可行：

1. **clu_core 模块可以直接使用 `kavonrtep/repeatexplorer:2.3.8` 官方镜像**——这是开发者 Petr Novak 本人维护的，不需要等待自建镜像
2. **Conda 包构建**：RepeatExplorer2 因 32-bit 依赖和 Perl 脚本特性，构建 conda 包的难度极高。最务实的策略是 **map 模式用 conda，clu 模式用 Docker 兜底**（Nextflow 原生支持 process 级混合）
3. **不需合并镜像**：思路 1 的单镜像优势不足以抵消其架构退化成本

### 推荐方案

```
ECCSPLORER (map)          → conda: yangmingsi::eccsplorer=2022.01.1.1
                             container: quay.io/bioinfortools/eccsplorer:2022.01.1.1

ECCSPLORER_CLU_PREPARE    → conda: yangmingsi::eccsplorer=2022.01.1.1 (复用)
                             container: quay.io/bioinfortools/eccsplorer:2022.01.1.1 (复用)

ECCSPLORER_CLU_CORE       → container: kavonrtep/repeatexplorer:2.3.8 (官方镜像)
                             或自建 quay.io/bioinfortools/eccsplorer-repex:2022.01.1.1
                             (seqclust 在 /repex_tarean/seqclust)

ECCSPLORER_CLU_CANDIDATES_PLOT → conda/container: 复用 eccsplorer 或 repeatexplorer
```

---

## 六、RepeatExplorer2 自建 Docker 镜像方案（可选）

如果需要保持 `bioinfortools` 命名空间一致，可自建 `eccsplorer-repex` 镜像：

```dockerfile
FROM quay.io/bioinfortools/eccsplorer:2022.01.1.1

# 从 GitHub 克隆官方 RepeatExplorer2 源码并编译 Louvain C++ 二进制
# 注：formatdb/mgblast 的 32-bit 依赖已由 conda-forge 的 blast-legacy 包（64-bit）替代
RUN git clone --depth 1 https://github.com/kavonrtep/repex_tarean.git /opt/repex_tarean \
    && cd /opt/repex_tarean \
    && make \
    && ln -sf /opt/repex_tarean/seqclust /opt/conda/envs/eccsplorer/bin/seqclust \
    && cp /opt/repex_tarean/bin/louvain_* /opt/conda/envs/eccsplorer/bin/
```

然后推送至 quay.io：`docker push quay.io/bioinfortools/eccsplorer-repex:2022.01.1.1`
