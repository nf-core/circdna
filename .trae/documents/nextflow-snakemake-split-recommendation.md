# Nextflow ↔ Snakemake 工作流划分建议

> **目的**：基于 circdna.nf 三种模式（reference / eccdna / integrated）的当前架构与计算特性，给出"哪些保留在 Nextflow、哪些迁移到 Snakemake"的建议。
>
> **状态**：建议稿（待用户确认后再细化迁移计划）

---

## 一、划分原则

| 维度         | Nextflow 适合                    | Snakemake 适合                  |
| ------------ | -------------------------------- | ------------------------------- |
| 计算量       | 重计算（高 CPU/内存/长时）       | 轻计算（< 4 CPU、< 8 GB、秒级） |
| 频次         | 每样本一次的标准化批处理         | 一次性/跨样本汇总               |
| 参数敏感度   | 固定参数、确定性输出             | 参数可调、需迭代探索            |
| 依赖关系     | 强依赖、必须串行、重算成本高     | 弱依赖、可独立重算、成本低      |
| 容器化需求   | 需要 Docker/Singularity 统一环境 | 可用轻量 Python/R 环境          |
| 输入数据规模 | 大文件（FASTQ/BAM）              | 小文件（BED/TXT/TSV）           |

---

## 二、三种模式步骤特性分析

### 2.1 reference 模式

```
FASTQ → INPUT_CHECK → CAT_FASTQ → FASTQC → TRIMGALORE
                                          ↓
                BWA_INDEX → BWA_MEM (process_high, 12C/72G)
                                   ↓
                SAMTOOLS_INDEX_BAM → SAMTOOLS_FAIDX (一次性)
                                   ↓
                BAM_STATS_SAMTOOLS → BAM_MARKDUPLICATES_PICARD (内存密集)
                                   ↓
                SAMTOOLS_VIEW_FILTER → SAMTOOLS_SORT_FILTERED → SAMTOOLS_INDEX_FILTERED
                                   ↓
                MOSDEPTH (process_medium, 6C/36G) → mosdepth_bed / global.dist.txt
```

**特性**：纯标准化"上游加工"流水线，所有步骤确定性输出，无探索性环节。

### 2.2 eccdna 模式

```
BAM_PREPROCESSING (同上)
        │
        ├── MOSDEPTH → mosdepth_bed
        ├── ECCSPLORER → eccsplorer_bed
        └── CIRCLE_MAP_PIPELINE:
                SAMTOOLS_SORT_QNAME_CM → CIRCLEMAP_READEXTRACTOR
                  → SAMTOOLS_SORT_RE + SAMTOOLS_INDEX_RE
                  → CIRCLEMAP_REPEATS + CIRCLEMAP_REALIGN (process_high + 96h)
                  → circle_map_bed
        │
        ▼
CANDIDATE_MERGE (process_single, Python, --max-distance 可调) → merged_bed
```

**特性**：在 reference 基础上增加两条并行检测支路。`CIRCLEMAP_REALIGN` 是全流程最重计算步骤（process_high + 96h）。`CANDIDATE_MERGE` 是参数敏感的轻量 Python 脚本。

### 2.3 integrated 模式

```
ch_trimmed_reads
   ├── filter(datatype=='gdna')   → REFERENCE_MODE → reference_mosdepth_bed
   └── filter(datatype=='eccdna') → ECCDNA_MODE    → eccdna_mosdepth_bed + merged_bed
                                                                   ↓
                                              ECC_SCORE (process_single, Python)
                                                  score = w1*JR + w2*log2(D_ecc/D_gdna) - w3*TE_penalty
                                                                   ↓
                                              scored_bed (含 grade: High/Medium/Low)
```

**特性**：先并行运行 reference + eccdna 两套子流程，再以两者 mosdepth + merged_bed 作为输入，调用 `ECC_SCORE` 计算综合评分。`ECC_SCORE` 的 `w1/w2/w3` 权重 + High/Medium/Low 分级阈值（10/5）均为可调参数，**天然适合迭代调优**。

---

## 三、各模块特性标签

### 3.1 重计算模块（建议保留 Nextflow）

| 模块                        | Label              | 资源    | 计算类型      | 频次         | 重算成本 |
| --------------------------- | ------------------ | ------- | ------------- | ------------ | -------- |
| `BWA_MEM`                   | process_high       | 12C/72G | CPU+IO 密集   | 每样本       | **极高** |
| `BAM_MARKDUPLICATES_PICARD` | process_medium     | 6C/36G  | 内存密集      | 每样本       | 高       |
| `MOSDEPTH`                  | process_medium     | 6C/36G  | CPU+IO 中等   | 每样本       | 中       |
| `CIRCLEMAP_REALIGN`         | process_high + 96h | 12C/72G | CPU+IO 重计算 | 每样本       | **极高** |
| `UNICYCLER`                 | 默认 + 96h         | -       | CPU 长时      | 每样本       | **极高** |
| `AMPLICONSUITE`             | process_low + 96h  | -       | CPU+内存 长时 | 每样本       | **极高** |
| `BWA_INDEX`                 | process_single     | 1C/6G   | CPU 密集      | 每参考组一次 | 中       |
| `ECCSPLORER`                | process_medium     | 6C/36G  | (stub 占位)   | 每样本       | 中       |

### 3.2 中间加工模块（建议保留 Nextflow）

| 模块                             | Label          | 资源   | 计算类型        | 频次     |
| -------------------------------- | -------------- | ------ | --------------- | -------- |
| `SAMTOOLS_FAIDX/INDEX/VIEW/SORT` | process_single | 1C/6G  | IO 轻量         | 每样本   |
| `BAM_STATS_SAMTOOLS`             | process_single | 1C/6G  | IO 轻量         | 每样本   |
| `CIRCLEMAP_READEXTRACTOR`        | process_low    | 2C/12G | CPU 单线程      | 每样本   |
| `CIRCLEMAP_REPEATS`              | process_low    | 2C/12G | CPU 轻量        | 每样本   |
| `CIRCLEFINDER`                   | process_low    | 2C/12G | awk 串行 IO+CPU | 每样本   |
| `MULTIQC`                        | process_single | 1C/6G  | IO 轻量         | 全局一次 |

### 3.3 探索性轻量模块（建议迁移 Snakemake）

| 模块              | Label          | 资源  | 计算类型       | 参数敏感                      | 重算成本 |
| ----------------- | -------------- | ----- | -------------- | ----------------------------- | -------- |
| `CANDIDATE_MERGE` | process_single | 1C/6G | IO 轻量 Python | `--max-distance` (默认 100bp) | **秒级** |
| `ECC_SCORE`       | process_single | 1C/6G | IO 轻量 Python | `w1/w2/w3` + 分级阈值         | **秒级** |

---

## 四、划分建议（核心结论）

### 4.1 保留在 Nextflow 的部分

**所有重计算 + 中间加工步骤** —— 这些步骤具有"输入大文件、输出确定性、参数固定、重算成本高"的特点：

```
Nextflow 保留范围:
├── BAM 预处理全链路（BWA_INDEX/MEM → SAMTOOLS_* → MARKDUPLICATES → MOSDEPTH）
├── ECCSPLORER（eccDNA 候选检测）
├── CIRCLE_MAP_PIPELINE 全流程（含 REALIGN 96h 重计算）
├── CIRCLEFINDER（legacy 模式）
├── UNICYCLER（长读组装）
├── AMPLICONSUITE（96h 重计算）
└── MULTIQC（全局汇总报告）
```

### 4.2 迁移到 Snakemake 的部分

**两个参数敏感的轻量 Python 评分步骤 + 探索性可视化**：

```
Snakemake 迁移范围:
├── CANDIDATE_MERGE   （bin/merge_candidates.py，--max-distance 调参）
├── ECC_SCORE         （bin/calculate_ecc_score.py，w1/w2/w3 + 阈值网格搜索）
└── 探索性可视化（新增）:
    ├── mosdepth 深度分布对比（gDNA vs eccDNA）
    ├── scored_bed 分级统计 / 阈值敏感性曲线
    ├── eccDNA 候选长度 / 染色体分布
    ├── 工具一致性 Venn 图（ECCsplorer vs Circle-Map）
    └── 跨样本候选汇总（当前 Nextflow 无此步骤）
```

### 4.3 契约接口（Nextflow 产物 → Snakemake 输入）

| Nextflow 产物          | 路径模式                                            | Snakemake 用途       |
| ---------------------- | --------------------------------------------------- | -------------------- |
| mosdepth 深度分布      | `{outdir}/mosdepth/{sample}/*.global.dist.txt`      | 深度可视化、阈值探索 |
| mosdepth 区域深度      | `{outdir}/mosdepth/{sample}/*.regions.bed.gz`       | ECC_SCORE 输入       |
| Circle-Map realign bed | `{outdir}/circlemap/realign/*.bed`                  | CANDIDATE_MERGE 输入 |
| ECCsplorer bed         | `{outdir}/eccdna_mode/eccsplorer/*_candidates.bed`  | CANDIDATE_MERGE 输入 |
| (可选) merged_bed      | `{outdir}/eccdna_mode/candidate_merge/*_merged.bed` | ECC_SCORE 备用输入   |
| BAM 标准化产物         | `{outdir}/bwa/*.sorted.bam` + `.bai`                | 备用 IGV 检查        |

---

## 五、迁移理由详述

### 5.1 为什么 CANDIDATE_MERGE 适合迁移

1. **计算成本极低**：process_single（1 CPU/6GB），仅读取两个小 BED 文件，重算秒级
2. **参数明确可调**：`--max-distance`（默认 100bp）直接影响两个工具候选的合并结果，适合做参数扫描
3. **依赖独立**：仅依赖 ECCsplorer + Circle-Map 的 BED 输出，不依赖 BAM/FASTQ
4. **重算无损**：调整 max_distance 不需要重新触发上游 96h 的 CIRCLEMAP_REALIGN

### 5.2 为什么 ECC_SCORE 适合迁移

1. **计算成本极低**：process_single（1 CPU/6GB），纯 Python 脚本，输入均为小文本 BED
2. **多参数可调**：`w1/w2/w3` 三个权重 + High/Medium/Low 两个分级阈值（10/5）
3. **输入独立于参考基因组**：仅需 merged_bed + 两个 mosdepth bed + repeat_bed（可选）
4. **输出可解释**：scored_bed 含 grade 列，可直接评估不同权重对候选分级的影响
5. **天然适合网格搜索**：在 Snakemake 中可用 `expand` 或 configfile 做参数网格搜索

### 5.3 为什么 MOSDEPTH 保留在 Nextflow

虽然 mosdepth 输出的 `global.dist.txt` 适合在 Snakemake 中做探索性可视化，但 mosdepth 本身是 process_medium（6C/36G）的中等计算步骤，且其输出是 ECC_SCORE 的必要输入。将其保留在 Nextflow 中可确保：

- BAM → mosdepth 的强依赖链在 Nextflow 内部完整
- mosdepth 产物作为契约接口供 Snakemake 探索使用

### 5.4 为什么可视化/统计适合放 Snakemake

当前 Nextflow 流程**缺少**以下探索性环节，它们均基于已有 BED/BAM 产物，不触发新比对：

- eccDNA 候选长度分布、染色体分布统计
- ECC_SCORE 分级分布、阈值敏感性曲线
- 跨样本候选重叠热图
- mosdepth 深度分布对比（gDNA vs eccDNA）
- 工具一致性 Venn 图（ECCsplorer vs Circle-Map）

这些是典型的"低密集计算 + 探索性"任务，与 Snakemake 的定位完全契合。

---

## 六、划分边界图

```
┌─────────────────────── Nextflow (重计算 + 标准化) ──────────────────────┐
│                                                                        │
│  BWA_INDEX → BWA_MEM → SAMTOOLS_* → MARKDUPLICATES → MOSDEPTH         │
│       │                                                                │
│       ├── ECCSPLORER ──┐                                              │
│       │                 └── (原始 bed 输出)                            │
│       └── CIRCLE_MAP ──── (原始 bed 输出)                              │
│                                                                        │
│  [契约接口]：mosdepth_bed + eccsplorer_bed + circle_map_bed            │
└────────────────────────────────────────────────────────────────────────┘
                                │
                                ▼
┌─────────────────────── Snakemake (轻分析 + 探索) ───────────────────────┐
│                                                                        │
│  CANDIDATE_MERGE (max_distance 调参)                                   │
│       │                                                                │
│       ▼                                                                │
│  ECC_SCORE (w1/w2/w3 网格搜索)                                         │
│       │                                                                │
│       ├── 深度分布可视化 (gDNA vs eccDNA)                              │
│       ├── 候选分级统计 / 阈值敏感性                                     │
│       ├── 工具一致性 Venn 图                                            │
│      └── 跨样本候选汇总 (新增)                                         │
│                                                                        │
└────────────────────────────────────────────────────────────────────────┘
```

---

## 七、关键文件清单

### Nextflow 侧（保留）

- [modules/local/eccsplorer/main.nf](file:///Users/siyangming/nextflow_nf_core/circdna.nf/modules/local/eccsplorer/main.nf)
- [modules/local/circle_map_pipeline/main.nf](file:///Users/siyangming/nextflow_nf_core/circdna.nf/modules/local/circle_map_pipeline/main.nf)
- [modules/local/mosdepth/main.nf](file:///Users/siyangming/nextflow_nf_core/circdna.nf/modules/local/mosdepth/main.nf)
- [subworkflows/local/bam_preprocessing/main.nf](file:///Users/siyangming/nextflow_nf_core/circdna.nf/subworkflows/local/bam_preprocessing/main.nf)

### Snakemake 侧（迁移目标）

- [modules/local/candidate_merge/main.nf](file:///Users/siyangming/nextflow_nf_core/circdna.nf/modules/local/candidate_merge/main.nf)
- [modules/local/ecc_score/main.nf](file:///Users/siyangming/nextflow_nf_core/circdna.nf/modules/local/ecc_score/main.nf)
- [bin/merge_candidates.py](file:///Users/siyangming/nextflow_nf_core/circdna.nf/bin/merge_candidates.py)
- [bin/calculate_ecc_score.py](file:///Users/siyangming/nextflow_nf_core/circdna.nf/bin/calculate_ecc_score.py)

### 三种模式入口

- [workflows/circdna.nf](file:///Users/siyangming/nextflow_nf_core/circdna.nf/workflows/circdna.nf) 第 271-329 行（模式路由）
- [subworkflows/local/reference_mode/main.nf](file:///Users/siyangming/nextflow_nf_core/circdna.nf/subworkflows/local/reference_mode/main.nf)
- [subworkflows/local/eccdna_mode/main.nf](file:///Users/siyangming/nextflow_nf_core/circdna.nf/subworkflows/local/eccdna_mode/main.nf)
- [subworkflows/local/integrated_mode/main.nf](file:///Users/siyangming/nextflow_nf_core/circdna.nf/subworkflows/local/integrated_mode/main.nf)

---

## 八、待用户决策的关键问题

在制定具体迁移实施计划前，需用户确认以下决策点：

1. **迁移粒度**：

   - 选项 A：仅迁移 `ECC_SCORE`（最小改动，integrated 模式收益最大）
   - 选项 B：迁移 `CANDIDATE_MERGE` + `ECC_SCORE`（推荐，覆盖两个参数敏感步骤）
   - 选项 C：迁移上述两者 + 新增探索性可视化（最完整，需新建 Snakemake 项目）

2. **Nextflow 侧处理方式**：

   - 选项 A：保留 Nextflow 中的 CANDIDATE_MERGE/ECC_SCORE 模块作为"默认实现"，Snakemake 作为"探索副本"
   - 选项 B：从 Nextflow 中移除这两个模块，eccdna_mode/integrated_mode 在产出原始 bed 后即终止，下游完全交给 Snakemake

3. **Snakemake 项目位置**：

   - 选项 A：在 circdna.nf 仓库内新建 `snakemake/` 子目录
   - 选项 B：新建独立仓库（如 `circdna-explore.smk`）

4. **是否需要参数网格搜索**：

   - 是否需要 Snakemake 支持 `w1/w2/w3` 网格搜索，还是仅支持单组参数调优

5. **是否新增可视化模块**：
   - 是否需要在 Snakemake 中新增深度分布、候选统计、工具一致性等可视化规则

---

## 九、验证步骤（迁移实施后）

1. **契约接口验证**：Nextflow 产出的 mosdepth_bed/eccsplorer_bed/circle_map_bed 路径与 Snakemake 期望一致
2. **结果一致性验证**：相同输入下，Snakemake 实现的 CANDIDATE_MERGE/ECC_SCORE 输出与原 Nextflow 模块完全一致
3. **参数扫描验证**：在 Snakemake 中调整 `--max-distance` 和 `w1/w2/w3`，确认输出随参数变化符合预期
4. **三种模式覆盖验证**：
   - reference 模式：Nextflow 独立产出 mosdepth_bed（无 Snakemake 下游）
   - eccdna 模式：Nextflow 产出原始 bed → Snakemake 接管 CANDIDATE_MERGE
   - integrated 模式：Nextflow 产出原始 bed + mosdepth → Snakemake 接管 CANDIDATE_MERGE + ECC_SCORE

---

_本建议基于当前 v3.2.1 版本的代码分析，若用户确认迁移意图，将基于上述划分进一步细化实施计划。_
