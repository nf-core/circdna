# ECCsplorer 四模式模块化拆分 Spec

## Why

当前 `ECCsplorer` 在 `circdna.nf` 中只深度接入了 `map` 主链，其他模式仍停留在文档设计层，导致 `ECCsplorer` 仍然以“大镜像 + 大黑盒 process”的方式存在，不利于缓存复用、资源分配、断点续跑和容器瘦身。

用户希望将 `PRExer`、`map`、`clu`、`all/comparative` 四种模式按标准 Nextflow 模块和子工作流彻底拆开，并尽可能以“小 process / 单软件”为最小执行单位，以提升执行效率并降低镜像体积。

## What Changes

- 将 `ECCsplorer` 四种模式拆分为标准的 Nextflow 本地模块与子工作流
- 将原本单体 `ECCSPLORER` 大 process 拆解为多个最小职责 process，优先做到“一步一工具”或“一类明确变换一步”
- 为 `PRExer`、`map`、`clu`、`all/comparative` 分别定义标准输入输出契约
- 增加一层“流程归属分析”，明确哪些能力适合迁移到 Nextflow，哪些能力更适合保留或迁移到 Snakemake
- 对仓库中已经存在的模块，优先迁移 `ECCsplorer` 的参数语义和参数集，而不是默认重写实现
- 为 `clu` 与 `all/comparative` 引入独立 subworkflow，而不是继续依赖黑盒 CLI 串联
- 将可复用的公共步骤抽成共享模块，例如 reads 预处理、格式转换、BLAST 注释、结果提取
- 重新规划容器/环境边界：
  - 轻量模块优先使用更小的独立环境
  - 重型步骤单独绑定专用容器或环境
  - 避免所有步骤共用一个超大 `eccsplorer` 镜像
- 为每个模块定义资源标签、缓存边界和 `-resume` 复用策略
- 将 `all` 模式重定义为由 Nextflow DAG 组合 `PRExer + map + clu + comparative` 的编排层，而不再直接调用上游黑盒总入口

## Confirmed Decisions

### 1. 四模式的 Nextflow / Snakemake 归属结论

| 模式          | Nextflow 归属                                                                  | Snakemake 归属                                                  | 结论                                                                    |
| ------------- | ------------------------------------------------------------------------------ | --------------------------------------------------------------- | ----------------------------------------------------------------------- |
| `PRExer`      | reads 质控、裁剪、过滤、格式标准化、供 `map/clu` 复用的统一预处理链            | 不承担 `PRExer` 主链，仅消费标准化产物做探索性分析              | `PRExer` 不再保留黑盒入口，整体以 Nextflow 标准预处理模块链替代         |
| `map`         | 输入标准化、核心比对、候选提取、junction/序列/统计提取、可选 BLAST 注释        | 消费 `map` 标准化产物做轻量筛选、评分、统计与图形分析           | `map` 主链明确归属 Nextflow                                             |
| `clu`         | cluster 输入准备、RepeatExplorer2/seqclust 等重型聚类执行、聚类结果标准化导出  | 消费 `clu` 标准候选与汇总字段做轻量探索性分析                   | `clu` 的重计算部分归属 Nextflow，下游解释层暴露给 Snakemake             |
| `all`         | 作为 DAG 编排层组合 `PRExer + map + clu + comparative`，负责并行调度与缓存复用 | 不承担总入口编排                                                | `all` 只保留为 Nextflow 编排层，不再调用 ECCsplorer 黑盒总入口          |
| `comparative` | 消费 `map/clu` 标准化输出，完成规则固定、可缓存、可独立重跑的标准比较产物生成  | 对 comparative 产物继续做探索性评分、统计汇总、阈值试探和可视化 | `comparative` 采取“标准化比较在 Nextflow，探索性解释在 Snakemake”的分层 |

### 2. 可复用现有模块与“参数优先迁移”策略

- 首选复用现有 `BAM_PREPROCESSING`、`FASTQC`、`TRIMGALORE/fastp`、`MOSDEPTH`、`ECCSPLORER`、`ECCSPLORER_WITH_CONTROL` 以及现有结果提取/导出逻辑，先迁移 `ECCsplorer` 的参数语义、默认值和开关，再判断是否必须深拆。
- `PRExer` 替代链优先将最小长度、质量截断、过滤和格式标准化语义映射到现有 reads 预处理模块或其 `ext.args`，避免第一步就重写实现。
- `map` 优先基于现有 `ECCSPLORER` / `ECCSPLORER_WITH_CONTROL` 补齐参数集，再将其中已确认需要独立缓存和独立资源配置的步骤拆成更小模块。
- 现有结果提取、注释和导出模块优先承接 `map/comparative` 的字段扩展，包括候选 BED、junction、序列、alignment stats、comparative 标准化产物和 handoff 导出。
- 参数优先迁移后的验证维度固定为：结果文件数量与关键字段、候选区间数量/摘要统计、运行时间、资源使用、`-resume` 缓存命中表现；只有差异分析表明现有边界无法承载原始语义时，才进入更深层重构。

### 3. 四模式的最小职责模块拆分

**`PRExer` 最小职责模块**

- `eccsplorer_prx_qc_trim`：reads 质控与裁剪，承接长度/质量类参数语义。
- `eccsplorer_prx_filter_normalize`：低复杂度/低质量过滤与配对完整性整理。
- `eccsplorer_prx_format_prepare`：统一输出供 `map/clu` 消费的标准 FASTQ 对与元数据。

**`map` 最小职责模块**

- `eccsplorer_input_normalize`：统一 BAM/FASTQ 输入，产出标准化 FASTQ 对和输入元数据。
- `eccsplorer_map_core`：执行 `map` 核心比对/检测。
- `eccsplorer_map_candidates`：提取高置信/低置信候选区域 BED。
- `eccsplorer_map_extract`：提取 junction reads、ecc 序列和 alignment 统计。
- `eccsplorer_map_blast`：对候选序列执行可选 BLAST 注释。

**`clu` 最小职责模块**

- `eccsplorer_clu_prepare`：准备 cluster 所需 reads、标签和标准输入目录。
- `eccsplorer_clu_core`：执行 RepeatExplorer2/seqclust 等重型聚类。
- `eccsplorer_clu_annotate`：整理 cluster 注释和分类汇总。
- `eccsplorer_clu_candidates_plot`：提取 `cluster_candidates.bed` 或等价候选，并生成 cluster 相关绘图产物。

**`all/comparative` 最小职责模块**

- `eccsplorer_all_orchestrator`：只负责编排 `PRExer + map + clu + comparative` 的 DAG。
- `eccsplorer_comparative_prepare`：消费 `map/clu` 标准化输出并准备 comparative 输入契约。
- `eccsplorer_comparative_core`：生成标准化 comparative 表和候选对比结果。
- `eccsplorer_comparative_export`：输出供 `handoff.tsv` / `samples.auto.yaml` 和 Snakemake 下游消费的比较产物。

### 4. BAM / FASTQ 统一接口

- 统一接口只在输入标准化层区分 BAM 与 FASTQ；后续 `map`、`clu`、`comparative` 统一消费标准化后的 FASTQ 对与元数据。
- 若输入为 BAM，`eccsplorer_input_normalize` 负责完成 `samtools fastq` 或等价转换，并记录 `input_kind=bam`、`source_bam` 等追踪字段。
- 若输入为 FASTQ，`eccsplorer_input_normalize` 仅做配对检查、命名规范化和元数据补齐，并记录 `input_kind=fastq`。
- 统一后的模块契约为“标准化 reads + meta + fasta”，确保 `map`/`clu` 的核心执行层不再感知原始输入差异。
- `ECCSPLORER` / `ECCSPLORER_WITH_CONTROL` 现有接口继续作为过渡承载层，先补齐参数并复用 BAM/FASTQ 标准化入口，再决定是否继续细拆。

### 5. `handoff.tsv` 与 `samples.auto.yaml` 扩展字段

- `handoff.tsv` 必须在现有 `map` 主链字段之外增加可选字段：`map_candidates_bed`、`map_lowconf_bed`、`junction_reads_txt`、`ecc_sequences_fasta`、`alignment_stats_txt`、`cluster_bed`、`cluster_summary_tsv`、`comparative_tsv`、`comparative_bed`、`comparative_summary_tsv`。
- `samples.auto.yaml` 对应增加同名可选键，并允许按样本声明字段缺失状态，避免未启用 `clu`/`comparative` 时伪造空路径。
- `eccdna.smk` 的消费边界明确为：识别字段是否存在，存在则消费 `cluster_*` / `comparative_*` 标准化产物做轻量统计与可视化，不存在则只运行 `map` 相关下游。
- `clu` 的最小交接要求是暴露 `cluster_bed` 或等价候选字段；`comparative` 的最小交接要求是暴露可直接用于下游统计与绘图的标准化比较表。

### 6. 资源标签、镜像瘦身与 `-resume` 验证

- 资源标签按步骤职责拆分：预处理/格式转换优先 `process_single` 或 `process_low`，`map` 核心比对优先 `process_medium`，`clu` 核心聚类优先 `process_high` 或 `process_long`，轻量提取/注释/导出回落到 `process_single` 或 `process_low`。
- 镜像瘦身原则是“轻量步骤小环境、重型步骤独立大镜像”：samtools/bedtools/python/BLAST 类模块不再绑定包含 RepeatExplorer2 的大镜像；只有 `clu` 核心重型模块保留重依赖容器。
- `all/comparative` 不得重新把轻量步骤并回超大环境，避免因为单个提取或导出步骤改动而重拉大镜像、重跑重计算。
- `-resume` 验证的最小方案固定为三类场景：仅改提取/注释模块时 `map` 核心应缓存命中；仅改 comparative 逻辑时 `map/clu` 应缓存命中；仅改 handoff 导出时上游检测与比较步骤应缓存命中。
- 首期验收至少要求在 `test_local` 或等价本地数据集上完成一次全跑和一次 `-resume` 重跑，确认拆分后 `CACHED` 边界符合预期。

### 7. 首期范围与后置范围

**首期必须落地**

- `PRExer` 的标准预处理替代链与参数映射规则。
- `map` 模式的输入标准化、核心执行、候选/提取/BLAST 拆分方案。
- BAM/FASTQ 统一接口与基于 `ECCSPLORER` / `ECCSPLORER_WITH_CONTROL` 的过渡复用策略。
- `handoff.tsv` / `samples.auto.yaml` 的扩展字段设计。
- 资源标签拆分、镜像瘦身边界与 `test_local + -resume` 最小可行验收方案。

**后置实现范围**

- `clu` 的完整 RepeatExplorer2/seqclust 重型执行落地与其专用大镜像细化。
- `all` 的完整编排实现和 `comparative` 的全部独立模块化执行。
- comparative 之后的探索性评分、统计汇总、可视化和阈值试探，这部分继续优先放在 Snakemake。
- 针对重型依赖的进一步细粒度环境拆分、性能对比和大样本基准测试。

## Impact

- Affected specs:
  - `/Users/siyangming/nextflow_nf_core/.trae/specs/design-ngs-analysis-and-resume-real-test`
  - `/Users/siyangming/nextflow_nf_core/circdna.nf/.trae/specs/add-eccsplorer-database-and-control`
  - `/Users/siyangming/nextflow_nf_core/ECCsplorer/.trae/specs/build-eccsplorer-from-source-and-bam-support`
- Affected code:
  - `circdna.nf/modules/local/eccsplorer/`
  - `circdna.nf/subworkflows/local/eccdna_mode/`
  - `circdna.nf/subworkflows/local/eccsplorer_cluster/`
  - `circdna.nf/workflows/circdna.nf`
  - `circdna.nf/conf/modules.config`
  - `circdna.nf/nextflow.config`
  - `circdna.nf/nextflow_schema.json`
  - `circdna.nf/modules.json`
  - `circdna.nf/bin/`
  - `ECCsplorer/` 中与构建、运行环境、入口脚本相关的文件

## ADDED Requirements

### Requirement: 必须显式给出 Nextflow / Snakemake 归属分析

系统 SHALL 对 `ECCsplorer` 四种模式及其拆分步骤给出明确的流程归属建议，而不是默认所有能力都迁入同一引擎。

#### Scenario: 适合迁移到 Nextflow 的能力

- **WHEN** 某步骤具有以下特征：
  - 二进制或重型计算工具驱动
  - 输入输出边界稳定
  - 资源差异明显，适合按 process 分配 CPU/内存
  - 需要强依赖 `-resume`、缓存与并行 DAG
- **THEN** 该步骤 SHALL 优先归属 `circdna.nf` / Nextflow

#### Scenario: 适合迁移到 Snakemake 的能力

- **WHEN** 某步骤具有以下特征：
  - 参数敏感、探索性强
  - 以 Python / R 脚本编排、统计、筛选、可视化为主
  - 经常需要重跑轻量步骤而不希望触发上游重计算
  - 更适合以规则链方式快速迭代
- **THEN** 该步骤 SHALL 优先归属 `eccdna.smk` / Snakemake

#### Scenario: PRExer / map / clu / all 的归属结论

- **WHEN** 用户审查四种模式的归属分析
- **THEN** 系统 SHALL 至少给出以下结论：
  - `PRExer` 的标准预处理替代链优先归属 Nextflow
  - `map` 的核心检测链优先归属 Nextflow
  - `clu` 的重型聚类执行优先归属 Nextflow，但其下游统计汇总字段需面向 Snakemake 暴露
  - `all` 作为编排层优先归属 Nextflow
  - `comparative` 中偏规则化、可缓存的比较步骤可放 Nextflow
  - `comparative` 之后的探索性评分、统计汇总和图形分析可继续归属 Snakemake

### Requirement: 已有模块优先迁移参数而非立即重写

系统 SHALL 对仓库中已经存在、且功能边界基本匹配的模块采用“先迁参数、后看差异、再决定是否深拆”的策略，以降低改动风险。

#### Scenario: 已有模块功能边界可复用

- **WHEN** 当前仓库里已经存在与 ECCsplorer 某步骤功能近似的模块
- **THEN** 系统 SHALL 优先复用该模块
- **AND** 先补齐 `ECCsplorer` 原始参数语义、默认值和可调选项
- **AND** 不应在第一步就强行重写实现

#### Scenario: 参数迁移后的对比验证

- **WHEN** 已有模块完成 ECCsplorer 参数迁移
- **THEN** 系统 SHALL 比较迁移前后在以下维度的差异：
  - 结果文件数量与关键字段
  - 候选区域数量或统计摘要
  - 运行时间、资源使用与缓存表现
- **AND** 只有当差异分析表明现有模块边界无法承载 ECCsplorer 语义时，才进入更深层的模块重构

#### Scenario: 参数优先迁移的适用对象

- **WHEN** 审查当前已存在模块
- **THEN** 系统 SHALL 优先考虑以下对象是否采用“参数优先迁移”策略：
  - `fastp` 或等价预处理模块承接 `PRExer` 语义
  - 现有 `ECCSPLORER` / `ECCSPLORER_WITH_CONTROL` 模块承接 `map` 参数细化
  - 现有结果提取、BLAST 注释或交接导出模块承接 `map / comparative` 参数扩展

### Requirement: PRExer 模式必须模块化替代

系统 SHALL 不再把 `PRExer` 视为必须整体调用的黑盒步骤，而是将其拆解为 Nextflow 可编排的标准预处理模块。

#### Scenario: PRExer 被替代为标准模块链

- **WHEN** 用户启用与 `PRExer` 等价的数据预处理路径
- **THEN** 系统 SHALL 使用标准 Nextflow 模块完成 reads 质控、裁剪、过滤和格式准备
- **AND** 每个步骤 SHALL 可独立缓存和 `-resume`
- **AND** 不要求继续调用原始 `PRExer` 总入口

#### Scenario: PRExer 参数语义保留

- **WHEN** 原始 `PRExer` 存在最小长度、质量截断、格式标准化等核心参数语义
- **THEN** 新模块链 SHALL 保留这些语义
- **AND** 参数应映射到明确的 Nextflow params 或模块 `ext.args`
- **AND** 若现有 `fastp` 模块已可承载这些语义，则 SHALL 优先扩展参数而不是重写模块

### Requirement: map 模式必须拆成最小职责模块

系统 SHALL 将 `map` 模式拆解为多个标准模块，而不是继续通过单个 `ECCSPLORER` process 一次性完成全部工作。

#### Scenario: map 模式分层

- **WHEN** 用户运行参考基因组驱动的 ECCsplorer 检测
- **THEN** 系统 SHALL 至少拆分以下职责边界：
  - 输入标准化 / FASTQ 准备
  - 核心比对与 mapping 结果生成
  - 高置信/低置信候选区域提取
  - junction reads / 序列 / 统计结果提取
  - 可选 BLAST 注释
- **AND** 每个模块 SHALL 只承担一种清晰职责

#### Scenario: BAM 与 FASTQ 双输入兼容

- **WHEN** 输入为 BAM 或 FASTQ
- **THEN** 系统 SHALL 在输入标准化层统一处理差异
- **AND** 下游 mapping 与结果提取模块 SHALL 尽量复用同一接口
- **AND** 已有 `ECCSPLORER` / `ECCSPLORER_WITH_CONTROL` 模块 SHALL 先扩展参数集，再评估是否必须继续拆小

### Requirement: clu 模式必须拆为标准子工作流

系统 SHALL 将 `clu` 模式实现为标准 Nextflow subworkflow，而不是保留在后续文档中。

#### Scenario: clu 子工作流最小拆分

- **WHEN** 用户启用 `clu` 模式
- **THEN** 系统 SHALL 至少拆分为以下原子过程：
  - cluster 输入准备
  - 核心聚类执行
  - 聚类结果注释
  - cluster 候选提取与绘图
- **AND** 重型聚类步骤 SHALL 与轻量注释/绘图步骤分离

#### Scenario: clu 输出标准化

- **WHEN** `clu` 子工作流完成
- **THEN** 系统 SHALL 产出标准的 `cluster_candidates.bed` 或等价候选输出
- **AND** 输出 SHALL 能接入 `handoff.tsv` 与 `samples.auto.yaml`

### Requirement: all/comparative 必须改为 DAG 编排层

系统 SHALL 将 `all` / `comparative` 重新定义为 Nextflow 编排层，而不是继续通过单一 ECCsplorer 总入口串联执行。

#### Scenario: all 模式由子流程组合

- **WHEN** 用户请求 `all` 模式
- **THEN** 系统 SHALL 组合执行 `PRExer`、`map`、`clu` 及 comparative 相关模块
- **AND** 各分支 SHALL 尽可能并行
- **AND** comparative 逻辑 SHALL 作为独立步骤消费 `map` 与 `clu` 的标准化输出

#### Scenario: comparative 独立缓存

- **WHEN** `map` 或 `clu` 的输入未变化，仅调整 comparative 逻辑或比较参数
- **THEN** 上游 `map` 与 `clu` 结果 SHALL 可复用缓存
- **AND** comparative 步骤 SHALL 能够单独重跑

#### Scenario: comparative 与下游分析分界

- **WHEN** comparative 需要继续产出评分、汇总统计、下游富集或图形展示
- **THEN** 标准化比较产物的生成 SHALL 可放在 Nextflow
- **AND** 其后的探索性分析、分布统计、差异分析与可视化 SHALL 优先交给 Snakemake

### Requirement: 模块边界必须支持镜像瘦身

系统 SHALL 通过模块边界重构运行环境，使轻量步骤不再依赖超大 ECCsplorer 镜像。

#### Scenario: 轻量模块使用小环境

- **WHEN** 某模块仅依赖 Python、samtools、bedtools、blast 或单一工具
- **THEN** 该模块 SHALL 使用更小的独立环境或容器
- **AND** 不应被迫绑定包含 RepeatExplorer2 等重依赖的大镜像

#### Scenario: 重型步骤隔离

- **WHEN** 某步骤必须依赖 RepeatExplorer2、seqclust 或其他重型依赖
- **THEN** 系统 SHALL 将其隔离为独立模块
- **AND** 仅该模块使用重型镜像或环境

### Requirement: 模块拆分必须提高运行效率

系统 SHALL 将缓存粒度、资源标签和执行边界设计为优先优化 `-resume`、增量执行和资源利用率。

#### Scenario: 断点续跑

- **WHEN** 用户执行 `-resume`
- **THEN** 未变更输入的轻量步骤和重型步骤 SHALL 尽可能显示为 `CACHED`
- **AND** 不应因单个下游提取或注释步骤变更而强制重跑整个 ECCsplorer 黑盒流程

#### Scenario: 资源按步骤分配

- **WHEN** 模块分别属于轻量预处理、重型聚类、轻量结果提取、轻量绘图
- **THEN** 每类模块 SHALL 使用独立资源标签
- **AND** 不应继续沿用所有 ECCsplorer 步骤统一 `process_high` 的粗粒度分配

## MODIFIED Requirements

### Requirement: ECCsplorer 在 circdna.nf 中的集成方式

现有 ECCsplorer 集成方式 SHALL 从“以单体 process 为主，其他模式待后续实现”修改为“以标准模块 + 子工作流为主，四种模式均有明确模块边界”。

#### Scenario: map 不再是唯一深度接入模式

- **WHEN** 用户查看 ECCsplorer 的主线集成状态
- **THEN** 系统 SHALL 不再只有 `map` 模式深度接入
- **AND** `PRExer`、`clu`、`all/comparative` 也 SHALL 具备明确的模块化实现路径

### Requirement: 交接清单字段

现有 `samples.auto.yaml` / `handoff.tsv` 要求 SHALL 从仅覆盖 `map` 主链输出，扩展为支持 `clu` 与 comparative 相关输出的可选字段。

#### Scenario: 启用 clu 时的交接

- **WHEN** 用户启用 `clu` 子工作流
- **THEN** 交接清单 SHALL 能包含 `cluster_bed` 或等价字段
- **AND** 下游 SHALL 能识别这些字段是否存在

#### Scenario: 启用 comparative 时的交接

- **WHEN** 用户启用 comparative 相关步骤
- **THEN** 交接清单 SHALL 能包含 comparative 的标准化比较产物字段
- **AND** 这些字段 SHALL 作为 Snakemake 下游统计与可视化的输入候选

## REMOVED Requirements

### Requirement: ECCsplorer 黑盒式大 process 为默认集成方式

**Reason**: 这种方式导致环境过重、缓存粒度过粗、资源浪费严重，无法体现 Nextflow 的模块化与增量执行优势。
**Migration**: 使用标准模块和子工作流替代单体 `ECCSPLORER` 执行边界。

### Requirement: all 模式通过原始总入口直接运行

**Reason**: 直接调用黑盒总入口会掩盖中间状态、阻断缓存复用，并且不利于 comparative 逻辑独立重跑。
**Migration**: 使用 Nextflow DAG 组合 `PRExer`、`map`、`clu` 和 comparative 步骤。
