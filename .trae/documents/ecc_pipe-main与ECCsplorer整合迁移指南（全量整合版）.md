# **ecc_pipe-main 与 ECCsplorer 整合到 eccdna.smk / circdna.nf 的实施与迁移指南（全量整合版）**

**日期：**2026-08-04
**规划范围：**circdna.nf（Nextflow 上游）、eccdna.smk（Snakemake 下游）、ecc_pipe-main、ECCsplorer
**核心架构原则：**Nextflow 负责上游计算密集型与二进制数据处理，Snakemake 负责下游统计、筛选、注释与可视化编排，采用分阶段平滑迁移。

## ---

**一、 架构总览与分工设计**

### **1.1 三层整合架构**

| 架构层级                   | 核心引擎 / 载体            | 核心职责与输入输出                                                                                                                                                      | 包含的具体模块与工具                                                                           |
| :------------------------- | :------------------------- | :---------------------------------------------------------------------------------------------------------------------------------------------------------------------- | :--------------------------------------------------------------------------------------------- |
| **1\. 上游数据处理层**     | circdna.nf (Nextflow)      | **输入：** FASTQ / BAM **职责：** 数据 QC、Trimming、Alignment、BAM 预处理、Mosdepth 深度计算、候选 eccDNA 检测 **输出：** 候选 BED \+ 标准交接清单 (samples.auto.yaml) | fastp (替代 PRExer) Circle-Map ECCsplorer (map) ECCsplorer (clu Subworkflow) EXPORT_HANDOFF    |
| **2\. 下游统计与可视化层** | eccdna.smk (Snakemake)     | **输入：** samples.auto.yaml **职责：** 候选合并与评分、格式标准化、分布分析、注释、差异 Burden/DEG 分析、CIRCOS 圈图 **输出：** 样本级与队列级分析报告、统计图表       | rules/score.smk rules/standardize.smk rules/distribution.smk rules/deg.smk rules/visualize.smk |
| **3\. 能力来源与迁移参考** | ecc_pipe-main & ECCsplorer | **职责：** 作为代码逻辑与算法借鉴来源，丢弃原 ecc_pipe_master.py 总控入口与单体黑盒运行方式。                                                                           | 迁移 Python/R 分析逻辑 借鉴 PRExer 质控参数 借鉴 Comparative 交叉验证思想                      |

### **1.2 完整数据流向图**

\[ 原始 FASTQ / BAM \]
│
▼
┌─────────────────────────────────────────────────────────────────────────┐
│ NEXTFLOW: circdna.nf (上游计算密集型) │
│ │
│ \[ fastp Module \] ◄── 借鉴 PRExer 质控参数 (length\>=50, Q20/30截断) │
│ │ │
│ ├───► \[ CIRCLE_MAP Module \] ─────────┐ │
│ ├───► \[ ECCSPLORER_MAP Module \] ─────┼──► \[ EXPORT_HANDOFF Module \]│
│ │ │ │ │
│ └───► \[ ECCSPLORER_CLUSTER Subwf \] ──┘ │ │
│ (PREP \-\> RUN \-\> PLOT) ▼ │
│ handoff.tsv / samples.auto.yaml│
└────────────────────────────────────────────────────────┬────────────────┘
│
▼
┌─────────────────────────────────────────────────────────────────────────┐
│ SNAKEMAKE: eccdna.smk (下游统计与可视化编排) │
│ │
│ 读取 samples.auto.yaml │
│ │ │
│ ▼ │
│ \[ candidate_merge & ecc_score \] ◄── 借鉴 Comparative 交叉验证思想 │
│ │ │
│ ▼ │
│ \[ standardize \] ──► \[ distribution \] ──► \[ deg \] ──► \[ visualize \] │
│ (吸收 ecc_pipe 的 Standardize/Distribution/DEG/Circlize 逻辑) │
└─────────────────────────────────────────────────────────────────────────┘

## ---

**二、 核心组件拆解与整合策略**

### **2.1 ECCsplorer 4种运行模式的全面整合**

| ECCsplorer 模式 | 核心功能与算法特性                                                                  | 整合策略与实现路径                                                                                                                                           | 运行环境 / 容器配置                                                                   |
| :-------------- | :---------------------------------------------------------------------------------- | :----------------------------------------------------------------------------------------------------------------------------------------------------------- | :------------------------------------------------------------------------------------ |
| **PRExer**      | 基于 Trimmomatic 和 Python wrapper 的数据预处理与修剪。                             | **丢弃并替代。** 借鉴其严格质控参数（最小 Read 长度 $\\ge 50\\,\\text{bp}$、Phred Q20/Q30 截断、Poly-N 截断），直接在 circdna.nf 最上游统一使用 fastp 替代。 | Nextflow 原生 fastp 模块 (C++ 极速引擎)                                               |
| **map**         | 基于参考基因组比对检测 Split-reads 和 Discordant reads。                            | **主链优先集成。** 作为 circdna.nf 的默认候选检测模块 ECCSPLORER_MAP，产出 eccsplorer_candidates.bed。                                                       | Docker / Apptainer / Conda (eccsplorer_map.yml)                                       |
| **clu**         | 基于 RepeatExplorer2 / CD-HIT 的无参考基因组/重复序列 _de novo_ 聚类。              | **拆分为 Subworkflow 条件加载。** 在 circdna.nf 中作为可选分支（--run_eccsplorer_cluster），解耦为 4 个独立原子 Process。                                    | 优先 Apptainer (Singularity) 镜像；或独立 Conda 环境 (environment_eccsplorer_clu.yml) |
| **all**         | 串联 PRExer $\\rightarrow$ map $\\rightarrow$ clu $\\rightarrow$ Comparative 分析。 | **拆解由 Nextflow DAG 重构。** 废弃黑盒运行，利用 Nextflow 实现双路并行处理；将末端 Comparer 交叉验证逻辑移至 eccdna.smk 的评分规则中。                      | Nextflow 引擎自动化 DAG 调度                                                          |

#### **ECCsplorer clu 模式的 Subworkflow 拆分架构：**

为了实现细粒度资源调度（CPU/内存按需分配）与断点续传（-resume），clu 模式拆解为 subworkflows/local/eccsplorer_cluster/main.nf，包含 4 个 Process：

> 1. **ECCSPLORER_PREP_CLUSTER** (轻量)：Reads 下采样与 FASTA 格式规范化（cpus \= 2, memory \= '8 GB'）。
> 2. **ECCSPLORER_RUN_CLUSTER** (计算密集)：运行 RepeatExplorer2 / CD-HIT 聚类（cpus \= 16, memory \= '32 GB' 到 '64 GB'）。
> 3. **ECCSPLORER_ANNOTATE_CLUSTER** (中度)：比对 RepeatMasker / BLAST 数据库（cpus \= 8, memory \= '16 GB'）。
> 4. **ECCSPLORER_PLOT_CLUSTER** (单核)：调用 eccDNA_Rcodes.py 导出拓扑图并提取 cluster_candidates.bed（cpus \= 1, memory \= '4 GB'）。

### **2.2 ecc_pipe-main 功能完全拆解与 Snakemake 吸收**

彻底弃用 ecc_pipe-main/ecc_pipe_master.py 入口，保留其算法与代码逻辑，按功能拆解注入 eccdna.smk 的规则与脚本库中：

| 原 ecc_pipe-main 模块                                     | 目标 eccdna.smk Rule / 脚本       | 功能重构与改造说明                                                                                       |
| :-------------------------------------------------------- | :-------------------------------- | :------------------------------------------------------------------------------------------------------- |
| ecc_pipe_master.py / QC.py                                | **废弃 / 被替代**                 | QC 与 Trimming 上移至 circdna.nf 的 fastp 模块；总控调度交由 Snakemake 引擎。                            |
| 02.Detect (Circle-Map/AA/CReSIL/FLED)                     | circdna.nf 对应检测模块           | Circle-Map 保持在 circdna.nf；AA/CReSIL/FLED 放入 Nextflow 后续扩展路线图。                              |
| Distribution.py (前半部分)                                | rules/standardize.smk             | 吸收格式转换逻辑，统一 Circle-Map、ECCsplorer 及合并结果的 BED 列契约。                                  |
| Distribution.py (后半部分)                                | rules/distribution.smk            | 封装单样本分析：染色体分布、长度分布、HOMER 基因组注释及数据库注释。                                     |
| DEG.py \+ deseq2.R / edger.R / limma.R / clusterprofile.R | rules/deg.smk \+ scripts/r/       | 构建队列 Burden 矩阵，进行组间差异 eccDNA 筛选、DESeq2/edgeR/limma 差异分析与 ClusterProfiler 富集分析。 |
| Circlize.py \+ circlize.R                                 | rules/visualize.smk \+ scripts/r/ | 根据输入的特定 ecc_id 或基因组区域，一键绘制并导出高画质 CIRCOS 圈图。                                   |

### **2.3 上下游自动化交接契约（Hand-off Contract）**

为消除手动维护 config/samples.yaml 的易错点，在 Nextflow 尾端新增 EXPORT_HANDOFF 模块，调用 Python 脚本导出下游专属配置：

> - **handoff.tsv**：人类可读的汇总表格，用于排查与日志审计。
> - **samples.auto.yaml**：自动生成的结构化配置文件，直接供 eccdna.smk 消费。

#### **交接清单 (samples.auto.yaml) 核心字段契约：**

samples:
SAMPLE_01:
group: "Tumor"
genome: "hg38"
eccsplorer_bed: "/path/to/eccdna_out/eccsplorer/SAMPLE_01_eccsplorer.bed"
circle_map_bed: "/path/to/eccdna_out/circlemap/SAMPLE_01_circlemap.bed"
cluster_bed: "/path/to/eccdna_out/cluster/SAMPLE_01_cluster.bed" \# 若开启 clu 模式
eccdna_mosdepth_bed: "/path/to/eccdna_out/mosdepth/SAMPLE_01.per-base.bed.gz"
gdna_mosdepth_bed: "/path/to/eccdna_out/mosdepth/CONTROL_01.per-base.bed.gz" \# 可选

## ---

**三、 目录结构设计**

### **3.1 上游 circdna.nf 改造后目录结构**

circdna.nf/
├── bin/
│ └── export_handoff_manifest.py \# \[新增\] 生成 handoff.tsv 与 samples.auto.yaml
├── conf/
│ ├── modules.config \# \[修改\] 增加导出与各个 Module 资源配置
│ └── profiles.config
├── modules/
│ └── local/
│ ├── fastp/main.nf \# \[新增\] 替换 PRExer 的质控模块
│ ├── circle_map/main.nf
│ ├── export_handoff/main.nf \# \[新增\] 交接清单导出模块
│ └── eccsplorer/
│ ├── map.nf \# \[修改\] ECCSPLORER_MAP 模块
│ ├── prep_cluster.nf \# \[新增\] clu 步骤 1: 数据准备
│ ├── run_cluster.nf \# \[新增\] clu 步骤 2: 聚类计算
│ ├── annotate_cluster.nf \# \[新增\] clu 步骤 3: 数据库注释
│ └── plot_cluster.nf \# \[新增\] clu 步骤 4: 绘图与 BED 提取
├── subworkflows/
│ └── local/
│ ├── eccdna_mode/main.nf \# \[修改\] 主检测子工作流
│ └── eccsplorer_cluster/main.nf \# \[新增\] clu 专用子工作流
├── environment_eccsplorer_clu.yml \# \[新增\] clu 模式专属 Conda 环境文件
├── nextflow.config \# \[修改\] 新增 handoff 与 clu 开关参数
└── workflows/
└── circdna.nf \# \[修改\] 工作流入口，挂载交接导出

### **3.2 下游 eccdna.smk 改造后目录结构**

eccdna.smk/
├── Snakefile \# \[修改\] 模块化 include 顶层入口
├── config/
│ ├── config.yaml \# \[修改\] 全局分析参数 (阈值/注释/DEG设置)
│ └── samples.template.yaml \# \[修改\] 仅作为 fallback 和示例
├── rules/
│ ├── score.smk \# \[新增\] 包含 candidate_merge 和 ecc_score
│ ├── standardize.smk \# \[新增\] 标准化不同软件结果
│ ├── distribution.smk \# \[新增\] 单样本分布与 HOMER 注释
│ ├── deg.smk \# \[新增\] 队列 Burden 矩阵与差异/富集分析
│ └── visualize.smk \# \[新增\] CIRCOS 圈图可视化
└── scripts/
├── python/
│ ├── merge_candidates.py \# \[修改\] 增强列定义与校验
│ ├── calculate_ecc_score.py \# \[修改\] 拆分阈值至 config.yaml
│ ├── standardize_bed.py \# \[新增\] 从 Distribution.py 迁移
│ ├── distribution_analysis.py \# \[新增\] 从 Distribution.py 迁移
│ └── generate_burden_matrix.py \# \[新增\] 从 DEG.py 迁移
└── r/
├── deseq2.R \# \[新增\] 从 ecc_pipe 迁移并参数化
├── edger.R \# \[新增\] 从 ecc_pipe 迁移并参数化
├── limma.R \# \[新增\] 从 ecc_pipe 迁移并参数化
├── clusterprofile.R \# \[新增\] 从 ecc_pipe 迁移并参数化
└── circlize.R \# \[新增\] 从 ecc_pipe 迁移并参数化

## ---

**四、 分阶段实施路线图 (Phased Implementation)**

### **Phase 1：打通自动化交接与主链基础 (第 1-2 周)**

> - **目标：** 在 Nextflow 尾端产出 samples.auto.yaml，消除 Snakemake 手工配置痛点。
> - **具体任务：**

1. 在 circdna.nf 中配置 fastp 模块，设置 min_length \= 50 及 Q20/Q30 裁剪，替代 PRExer。
2. 编写 circdna.nf/bin/export_handoff_manifest.py 与 modules/local/export_handoff/main.nf。
3. 修改 workflows/circdna.nf，在 mode \== 'eccdna' 结束时触发导出。
4. 改造 eccdna.smk/Snakefile，使其默认加载 samples.auto.yaml 并通过 snakemake \-n 测试。
   > - **完成标志：** 运行 circdna.nf 自动产出 samples.auto.yaml，eccdna.smk \-n 可无缝读取并正确解析路径。

### **Phase 2：升级 eccdna.smk 评分与 Distribution 层 (第 3-4 周)**

> - **目标：** 在 Snakemake 中承接标准化与单样本分布分析。
> - **具体任务：**

1. 将 eccdna.smk/Snakefile 拆解为 include: "rules/score.smk"、standardize.smk、distribution.smk。
2. 迁移 ecc_pipe-main/analysis_code/Distribution.py 中的标准化与统计逻辑至 scripts/python/。
3. 增加 HOMER 和数据库注释 rule。
   > - **完成标志：** 对单样本自动生成合并评分 BED、标准化分析 BED、染色体/长度分布图及 HOMER 注释文件。

### **Phase 3：接入 DEG 与 CIRCOS 可视化模块 (第 5-6 周)**

> - **目标：** 实现队列级别差异表达/Burden 分析与圈图绘制。
> - **具体任务：**

1. 新增 rules/deg.smk，将 DEG.py 逻辑改编为矩阵构建脚本 generate_burden_matrix.py。
2. 迁移并参数化 deseq2.R、edger.R、limma.R 和 clusterprofile.R。
3. 新增 rules/visualize.smk，迁移 circlize.R。
   > - **完成标志：** 能够基于多个样本自动输出矩阵、DESeq2 差异表格、火山图、GO/KEGG 富集结果及特定区域的 CIRCOS 圈图。

### **Phase 4：解耦与接入 ECCsplorer clu 模式 (第 7 周)**

> - **目标：** 在 circdna.nf 中以可选子工作流形式运行无参考基因组聚类。
> - **具体任务：**

1. 创建 environment_eccsplorer_clu.yml 或 Apptainer 容器配置。
2. 编写 modules/local/eccsplorer/ 下的 4 个子模块及 subworkflows/local/eccsplorer_cluster/main.nf。
3. 在 nextflow.config 中配置 params.run_eccsplorer_cluster。
4. 将 clu 导出的 cluster_candidates.bed 整合进 export_handoff_manifest.py。
   > - **完成标志：** 开启 \--run_eccsplorer_cluster 后，能独立并行跑完聚类分支，并将其 BED 文件加入交接清单供下游评分加权。

### **Phase 5：扩展检测器路线图 (后续迭代)**

> - **目标：** 吸收 AmpliconArchitect (AA)、CReSIL 及 FLED。
> - **优先级顺序：**

1. AmpliconArchitect (AA)：作为短读长扩增子扩展支路收回 circdna.nf。
2. CReSIL / FLED：作为长读长 (ONT/PacBio) 处理分支引入。

## ---

**五、 验证与端到端验收方案**

### **5.1 模块级验证（单元测试）**

> 1. **Phase 1 验证：**
>    nextflow run circdna.nf \-profile test,docker \--mode eccdna \--outdir ./results
>    snakemake \-s eccdna.smk/Snakefile \--configfile results/handoff/samples.auto.yaml \-n*检查点：* 确认 samples.auto.yaml 中的路径均真实存在，Snakemake 干跑（Dry-run）无语法或依赖错误。
> 2. **Phase 2 验证：**
>    snakemake \-s eccdna.smk/Snakefile \--configfile results/handoff/samples.auto.yaml \-R standardize_bed*检查点：* 检查产出的分析 BED 列定义是否符合标准，HOMER 注释日志无报错。
> 3. **Phase 3 验证：**
>    snakemake \-s eccdna.smk/Snakefile \--configfile results/handoff/samples.auto.yaml \-R run_deg*检查点：* 检查 burden_matrix.txt 维度，验证 DESeq2/edgeR 输出的 PDF 火山图与富集柱状图是否成功渲染。
> 4. **Phase 4 验证：**
>    nextflow run circdna.nf \-profile test,docker \--mode eccdna \--run_eccsplorer_cluster true \--outdir ./results_cluster*检查点：* 检查 4 个子模块节点是否均触发，且集群断点续传（-resume）能够跳过已完成的聚类节点。

### **5.2 端到端集成验收标准**

采用一组包含 Control 和 Tumor 的最小真实测试数据集，执行完整的端到端运行：

> 1. **上游运行：** 执行 circdna.nf，一键完成质控、比对、Circle-Map、ECCsplorer (map \+ clu) 及 Mosdepth 覆盖度计算。
> 2. **清单检查：** 确认自动生成了标准的 samples.auto.yaml，且正确配对样本的 group 信息。
> 3. **下游运行：** 运行 eccdna.smk，一键消费 samples.auto.yaml。
> 4. **最终交付物验收：**

- results/scored_beds/：高质量的候选 eccDNA 评分与合并结果。
- results/distribution/：单样本染色体/长度分布图及基因组注释表。
- results/deg/：队列差异 Burden 表格、DESeq2/edgeR/limma 差异基因集与 ClusterProfiler 富集图。
- results/visualize/：重点候选 eccDNA 区域的高分辨率 CIRCOS 圈图。

## ---

**六、 历史文档合并记录与冲突决议**

本节用于吸收 `.trae/documents` 中其余历史方案文档的有效信息，并明确与当前“全量整合版”发生冲突时的判定规则。

### **6.1 已合并的历史文档**

本指南已吸收以下历史文档中的有效内容：

- `ecc_pipe_eccsplorer_nf_upstream_smk_downstream_plan.md`
- `ecc_pipe_reevaluation_plan.md`
- `ecc_pipe_revisit_plan.md`
- `ecc_pipe_revisit_plan_v2.md`

吸收方式不是机械拼接原文，而是将其核心判断、边界分析、优先级结论与争议点整理为本节附录内容；避免旧方案和现行方案并列造成歧义。

### **6.2 冲突处理原则**

如历史文档与本指南存在冲突，统一以当前文件为准。具体优先级如下：

1. **架构归属冲突**
   以“`circdna.nf` 负责上游、`eccdna.smk` 负责下游”这一主原则为准。

2. **模块归属冲突**
   凡历史文档主张“继续把下游分析收回 `circdna.nf analysis mode`”者，均视为已被本指南否决的旧路线。

3. **实施顺序冲突**
   以本指南的 `Phase 1 -> Phase 5` 顺序为准；旧文档中的阶段 A/B/C、P0/P1/P2、v4.5.0 等分期建议仅保留为历史讨论背景。

4. **细节参数冲突**
   若旧文档中出现 `count_type`、`peak_path`、GO/KEGG barplot、GSEA 注释、多样本比较等实现建议，与本指南不一致时，应将其视为“待迁移能力池”，而非当前主架构的最终归属。

### **6.3 来自历史文档的有效补充结论**

以下判断来自旧文档，虽原归属方案已调整，但这些内容本身仍有保留价值，因此并入本指南：

#### **A. 下游能力池清单（已保留）**

旧文档集中识别出的 7 类高价值下游能力，仍然具有实施价值：

1. 多样本比较分析
2. GO/KEGG 富集可视化
3. GSEA 结果基因符号注释
4. 基本统计摘要与汇总信息
5. `count_type` 参数与计数口径明确化
6. `peak_path` / 扩展 BED 输入兼容
7. 数据库注释（SNP / Enhancer / SuperEnhancer / eQTL）

**当前决议：**
这些能力保留为待迁移能力池，但优先在 `eccdna.smk` 下游规则中吸收，而不是默认放回 `circdna.nf analysis mode`。

#### **B. 条件性保留的物种与数据库判断**

旧文档对数据库注释和部分分析功能做过适用性判断，这些结论继续有效：

- SNP / Enhancer / SuperEnhancer / eQTL 数据库注释优先面向 `hg38`
- 植物或其他非人类物种不应强行套用人类数据库注释
- 相关注释能力在下游规则中应当做成**条件执行**，由配置显式开启

#### **C. 长读长检测路线仍属后续扩展**

历史文档多次提到以下能力，但并未改变其阶段定位：

- `AmpliconArchitect`
- `CReSIL`
- `FLED`
- `CReSIL WGS` 模式

**当前决议：**
它们仍属于后续阶段扩展路线图，不进入首期主链验收门槛。

### **6.4 已否决但需留档的旧路线**

为了避免后续重复讨论，下列路线已被正式否决，但因其曾出现在历史文档中，仍保留在本节中作为“已决策不采用”的归档内容：

#### **路线 1：把下游分析全部并回 `circdna.nf`**

旧文档曾提出：

- 让 `circdna.nf analysis mode` 成为 Distribution / DEG / Visualize 的唯一承载者
- `eccdna.smk` 维持为仅包含 `candidate_merge + ecc_score` 的薄评分层

**否决原因：**

- 与当前已经确认的“Nextflow 上游、Snakemake 下游”职责划分冲突
- `ecc_pipe-main` 的 Python/R 探索分析更适合交由 Snakemake 编排
- 会导致 `circdna.nf` 主工作流再次膨胀

#### **路线 2：维持 `eccdna.smk` 为长期超薄评分层**

旧文档曾主张：

- `eccdna.smk` 不吸收 Distribution / DEG / Visualize
- 后续分析继续依赖 `circdna.nf analysis mode` 或手工 notebook

**否决原因：**

- 无法形成真正稳定的“Nextflow -> Snakemake”端到端闭环
- 会导致迁移后仍保留大量手工探索步骤
- 与本次整合目标“让 Snakemake 接手下游分析编排”不一致

#### **路线 3：继续保留 `ecc_pipe_master.py` 作为运行入口**

**否决原因：**

- 单体总控结构不利于模块化迁移
- 与 Nextflow / Snakemake 的分层编排设计冲突
- 会重新引入目录耦合和参数散落问题

### **6.5 历史文档的最终归档原则**

自本次合并后：

- “全量整合版”作为唯一对外有效的总规划文档
- 其余历史文档删除，不再并列保留
- 若后续需要追溯旧判断，应以本节中的“历史合并记录与冲突决议”为准，而不是恢复旧文档
