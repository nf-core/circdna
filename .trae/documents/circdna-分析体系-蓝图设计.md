# 植物 circRNA/eccDNA 多平台分析体系完整构建计划

## 1. 总体架构

建议最终形成一个：

> **三层架构 + 五个分析模块**

```
                    ┌──────────────────────────┐
                    │     Metadata Layer       │
                    │  circRNA / circDNA 分离  │
                    └────────────┬─────────────┘
                                 │
              ┌──────────────────┴──────────────────┐
              │                                     │
              ▼                                     ▼
     ┌─────────────────┐                   ┌─────────────────┐
     │   circRNA 层     │                   │   eccDNA 层     │
     │ circRNA.nf       │                   │ circDNA.nf      │
     └────────┬────────┘                   └────────┬────────┘
              │                                     │
              │                         ┌───────────┴───────────┐
              │                         │                       │
              │                         ▼                       ▼
              │                  ┌──────────────┐      ┌──────────────┐
              │                  │  Illumina    │      │ Long-read    │
              │                  │ eccDNA       │      │ PacBio/ONT   │
              │                  └──────────────┘      └──────────────┘
              │                         │                       │
              └─────────────────┬───────┴───────────────────────┘
                                │
                                ▼
                  ┌──────────────────────────┐
                  │  Unified Result Layer    │
                  │  circRNA / eccDNA        │
                  │  annotation / abundance  │
                  │  structure / sequence    │
                  └────────────┬─────────────┘
                               │
                               ▼
                  ┌──────────────────────────┐
                  │     Integration Layer    │
                  │     Snakemake建议         │
                  │                          │
                  │  DNA ↔ RNA ↔ Genome     │
                  │  eccDNA ↔ circRNA        │
                  │  TE ↔ SV ↔ Expression   │
                  └──────────────────────────┘
```

最终建议拆成：

| 模块     | 基础流程          | 主要任务            |
| -------- | ----------------- | ------------------- |
| Module 1 | `circrna.nf`      | circRNA             |
| Module 2 | `circdna.nf`      | Illumina eccDNA     |
| Module 3 | `circdna.nf` 扩展 | PacBio/ONT eccDNA   |
| Module 4 | 新增 Integration  | circRNA + eccDNA    |
| Module 5 | Benchmark         | 33 个代表数据集测试 |

---

# 2. Metadata 总体设计原则

最重要的设计原则是：

> **cirRNA metadata 与 eccDNA metadata 不共用同一个 schema。**

因为两者实验设计完全不同。

---

# 3. circRNA Metadata

建议建立：

```
metadata/
└── circRNA/
    ├── Oryza_sativa.tsv
    ├── Triticum_aestivum.tsv
    ├── Arabidopsis_thaliana.tsv
    ├── Artemisia_annua.tsv
    └── ...
```

每个物种一个 metadata。

例如：

```
metadata/circRNA/Oryza_sativa.tsv
```

建议字段：

| 字段           | 说明                             |
| -------------- | -------------------------------- |
| sample_id      | 样本唯一 ID                      |
| species        | 物种                             |
| strain         | 品种/品系                        |
| tissue         | 组织                             |
| treatment      | 处理                             |
| condition      | 实验条件                         |
| replicate      | 生物学重复                       |
| sample_type    | circRNA / total RNA              |
| run_id         | ERR/SRR/CRR                      |
| study_id       | SRA Study                        |
| project_id     | BioProject                       |
| platform       | Illumina                         |
| strategy       | RNA-Seq                          |
| layout         | SINGLE/PAIRED                    |
| read_length    | 150                              |
| library_type   | stranded/unstranded              |
| strandedness   | RF/RR/FR/unknown                 |
| enrichment     | rRNA depletion / polyA / RNase R |
| source         | TRANSCRIPTOMIC                   |
| reference      | 参考基因组                       |
| annotation     | GTF/GFF                          |
| data_source    | SRA/ENA/GSA/CNGB                 |
| fq1            | R1                               |
| fq2            | R2                               |
| analysis_group | control/treatment                |
| batch          | 批次                             |
| notes          | 备注                             |

---

## 3.1 当前 circRNA 数据分类

### A 类：真正适合 circRNA.nf

例如：

```
Oryza RNA-seq
ERR10889820
Illumina
RNA-Seq
TRANSCRIPTOMIC
cDNA
PAIRED
```

以及：

```
Wheat RNA-seq
ERR12724677
Illumina
RNA-Seq
TRANSCRIPTOMIC
RANDOM PCR
PAIRED
```

这两类数据：

```
RNA-seq
      │
      ▼
QC
      │
      ▼
Adapter trimming
      │
      ▼
Genome alignment
      │
      ▼
Back-splice junction detection
      │
      ├── CIRI2/CIRIquant
      ├── find_circ
      ├── DCC
      └── 其他
      │
      ▼
circRNA consensus
      │
      ▼
Annotation
      │
      ▼
Quantification
      │
      ▼
Differential circRNA
```

建议作为 `circrna.nf` 的标准输入。

---

### B 类：普通 WGS

例如：

```
ERR11535563
Artemisia gDNA
Illumina WGS
genome-skimming
```

不应该直接进入 circRNA.nf。

因为：

```
gDNA ≠ RNA
```

所以：

```
circRNA.nf
    │
    ├── RNA-seq       → YES
    ├── total RNA     → YES
    ├── rRNA depletion → YES
    ├── polyA RNA     → 视实验而定
    │
    └── gDNA/WGS      → NO
```

---

# 4. eccDNA Metadata

eccDNA 单独建立：

```
metadata/
└── eccDNA/
    ├── Oryza_sativa.tsv
    ├── Triticum_aestivum.tsv
    ├── Arabidopsis_thaliana.tsv
    ├── Amaranthus_palmeri.tsv
    ├── Alopecurus_myosuroides.tsv
    └── ...
```

推荐字段：

| 字段                  | 说明                    |
| --------------------- | ----------------------- |
| sample_id             | 样本 ID                 |
| species               | 物种                    |
| strain                | 品系                    |
| tissue                | 组织                    |
| treatment             | 处理                    |
| condition             | 条件                    |
| replicate             | 重复                    |
| run_id                | Run                     |
| study_id              | SRA Study               |
| project_id            | BioProject              |
| archive               | SRA/ENA/GSA/CNGB        |
| platform              | Illumina/PacBio/ONT     |
| sequencing_generation | short_read/long_read    |
| strategy              | WGS/AMPLICON/OTHER      |
| source                | GENOMIC                 |
| selection             | RANDOM/RANDOM PCR/other |
| layout                | SINGLE/PAIRED           |
| read_length           | 150/long                |
| library_name          | 文库                    |
| eccdna_enrichment     | yes/no                  |
| exonuclease_treatment | yes/no                  |
| plasmid_safe          | yes/no                  |
| RCA                   | yes/no                  |
| phi29                 | yes/no                  |
| debranching           | yes/no                  |
| T7_treatment          | yes/no                  |
| linearization         | yes/no                  |
| CIDER_seq             | yes/no                  |
| mobilome_seq          | yes/no                  |
| circle_seq            | yes/no                  |
| eccDNA_confidence     | high/medium/low         |
| reference             | genome                  |
| annotation            | GTF/GFF                 |
| fq1                   | R1                      |
| fq2                   | R2                      |
| bam                   | BAM                     |
| fast5/pod5            | ONT raw                 |
| ccs                   | PacBio HiFi             |
| notes                 | 备注                    |

---

# 5. eccDNA 数据分层

建议不要简单按照：

```
二代
三代
```

划分，而应该同时考虑：

```
测序平台
+
实验富集方法
+
RCA
+
是否真正 eccDNA enriched
```

因此推荐：

```
eccDNA
│
├── Level A
│   └── Enriched eccDNA + Illumina
│
├── Level B
│   └── Enriched eccDNA + PacBio
│
├── Level C
│   └── Enriched eccDNA + ONT
│
├── Level D
│   └── WGS / genome-skimming
│
└── Level E
    └── RNA-seq
```

其中：

```
A+B+C
```

是核心 eccDNA 数据。

```
D
```

是辅助验证。

```
E
```

不进入 eccDNA 主流程。

---

# 6. 二代 eccDNA 流程

## 6.1 基础流程

以现有 `circdna.nf` 为基础：

```
FASTQ
  │
  ▼
FastQC
  │
  ▼
fastp
  │
  ▼
Bowtie2/BWA
  │
  ▼
Alignment
  │
  ▼
BAM QC
  │
  ▼
ECCsplorer / circle detection
  │
  ▼
Candidate eccDNA
  │
  ▼
Breakpoint detection
  │
  ▼
Reference annotation
  │
  ├── Gene
  ├── TE
  ├── Repeat
  ├── CDS
  └── Intergenic
  │
  ▼
Quantification
  │
  ▼
Differential eccDNA
```

---

## 6.2 建议增加的模块

现有：

```
circdna.nf
```

建议新增：

```
modules/
├── fastp
├── bwa
├── bowtie2
├── samtools
├── eccsplorer
├── circle_map
├── circle_finder
├── annotation
└── quantification
```

如果当前流程已有部分模块，则只补缺失模块。

---

## 6.3 二代 eccDNA 推荐工具组合

建议：

```
Primary:
ECCsplorer

Secondary:
Circle-Map
Circle_finder

Validation:
BLAST
minimap2
```

如果目标是：

> **发现 + 定量**

建议：

```
ECCsplorer
      +
Circle-Map
      +
统一 consensus
```

如果目标是：

> **结构解析**

则：

```
Illumina
   │
   └── 主要负责
       ├── abundance
       ├── breakpoint
       └── candidate detection
```

不建议让 Illumina 独立承担复杂 eccDNA 全长结构解析。

---

# 7. 三代 eccDNA 流程

这是整个项目最需要扩展的部分。

建议从：

```
circdna.nf
```

派生出：

```
circdna_longread.nf
```

或者：

```
circdna.nf
│
├── shortread mode
│
├── pacbio mode
│
└── ont mode
```

我更推荐后者：

```
params.mode = 'shortread'
params.mode = 'pacbio'
params.mode = 'ont'
```

---

# 8. PacBio eccDNA Pipeline

主要包括：

```
PacBio BAM
   │
   ▼
BAM → FASTQ
   │
   ▼
CCS/HiFi
   │
   ▼
Quality filtering
   │
   ▼
DeConcat
   │
   ▼
TideHunter
   │
   ▼
CircleSeeker
   │
   ▼
Candidate circular molecules
   │
   ▼
minimap2
   │
   ▼
Reference alignment
   │
   ▼
Circular structure reconstruction
   │
   ▼
eccDNA annotation
```

---

## 8.1 CIDER-seq

对于：

```
SRR16958693
```

属于：

```
PacBio
AMPLICON
RANDOM PCR
CIDER-seq
```

建议：

```
CIDER-seq
   │
   ▼
DeConcat
   │
   ▼
TideHunter
   │
   ▼
CircleSeeker
   │
   ▼
eccDNA candidates
```

---

## 8.2 RCA-PacBio

例如：

```
SRR26069818
```

设计：

```
circular DNA
→ RCA
→ linearization
→ PacBio
```

建议：

```
PacBio
   │
   ▼
DeConcat
   │
   ▼
TideHunter
   │
   ▼
CircleSeeker
```

与 CIDER-seq 使用同一个主干流程。

因此：

> **SRR16958693 与 SRR26069818 可以归入同一个 PacBio long-read eccDNA 分析大类。**

但 metadata 中必须保留：

```
CIDER_seq = yes/no
RCA = yes/no
linearization = yes/no
```

不能把两者实验方法直接标成完全相同。

---

# 9. PacBio HiFi WGS

例如：

```
ERR11838731
SRR30359583
```

属于：

```
PacBio HiFi WGS
```

不是标准 CIDER-seq。

建议：

```
PacBio HiFi WGS
      │
      ▼
CCS/HiFi
      │
      ▼
QC
      │
      ▼
minimap2
      │
      ▼
Genome alignment
      │
      ▼
SV / circular structure
      │
      ├── Sniffles2
      ├── cuteSV
      └── SVIM
      │
      ▼
Candidate eccDNA
      │
      ▼
Circle validation
```

这里建议作为：

> **eccDNA 结构验证 / genomic origin analysis**

而不是主 eccDNA discovery pipeline。

---

# 10. ONT eccDNA Pipeline

建议：

```
ONT FASTQ
    │
    ▼
NanoPlot
    │
    ▼
Filtlong
    │
    ▼
minimap2
    │
    ▼
Genome alignment
    │
    ▼
eccDNA candidate detection
    │
    ├── ecc_finder
    ├── SVIM
    ├── Sniffles2
    └── cuteSV
    │
    ▼
Circular molecule reconstruction
    │
    ▼
TE / gene annotation
```

---

# 11. ONT 数据进一步分三类

### 第一类：真正 eccDNA enriched

例如：

```
ERR12724336
ERR6326020
SRR24335762
SRR24693334
SRR36603439
```

进入：

```
ONT eccDNA mode
```

---

### 第二类：RCA / T7 / 特殊处理

例如：

```
SRR28004411
SRR31773424
```

需要 metadata 标记：

```
RCA = yes/no
T7 = yes/no
```

然后：

```
ONT special mode
```

---

### 第三类：WGS

例如：

```
ERR12723706
```

属于：

```
ONT WGS
```

主要用于：

```
genome
+
SV
+
eccDNA origin
```

而不是直接等价于 eccDNA-enriched sequencing。

---

# 12. 推荐的三代统一架构

```
                  Long-read eccDNA
                         │
              ┌──────────┴──────────┐
              │                     │
           PacBio                  ONT
              │                     │
       ┌──────┴──────┐       ┌──────┴──────┐
       │             │       │             │
     CIDER/RCA     WGS    eccDNA-enriched  WGS
       │             │       │             │
       ▼             ▼       ▼             ▼
   DeConcat       minimap2  ecc_finder    minimap2
   TideHunter     SV tools  SV tools      SV tools
   CircleSeeker
       │             │       │             │
       └─────────────┴───────┴─────────────┘
                         │
                         ▼
                 Unified eccDNA
                 candidate table
```

---

# 13. 最终统一 eccDNA 数据结构

无论二代还是三代，最终都输出：

```
eccDNA_ID
species
sample_id
chromosome
start
end
length
strand
support_reads
support_method
platform
technology
circularity_score
breakpoint_confidence
gene
TE
repeat
GC
copy_number
abundance
```

这样可以实现：

```
Illumina
PacBio
ONT
```

三者进入同一个 downstream analysis。

---

# 14. circRNA + eccDNA 整合流程

这里我建议**不要继续用 Nextflow 作为最高层 workflow orchestration**。

你的想法：

> 用 Snakemake 构建第三个整合流程

我认为是合理的。

最终：

```
                 Snakemake
                     │
        ┌────────────┴────────────┐
        │                         │
        ▼                         ▼
   Nextflow                    Nextflow
   circRNA.nf                  circDNA.nf
        │                         │
        ▼                         ▼
    circRNA                    eccDNA
        │                         │
        └────────────┬────────────┘
                     ▼
               Integration
                     │
        ┌────────────┼────────────┐
        │            │            │
        ▼            ▼            ▼
      DNA           RNA          Genome
        │            │            │
        └────────────┼────────────┘
                     ▼
              Multi-omics result
```

---

# 15. 为什么推荐 Snakemake 做第三层

你的项目目前已经有：

```
circrna.nf
circdna.nf
```

如果继续用 Nextflow：

```
Nextflow
    ├── circRNA
    ├── circDNA
    └── integration
```

会导致：

```
pipeline
    └── pipeline
        └── pipeline
```

管理复杂度增加。

而 Snakemake 更适合作为：

> **Analysis Orchestrator**

例如：

```
Snakefile

rule circRNA:
    input:
        circRNA_metadata
    output:
        circRNA_results

rule circDNA:
    input:
        circDNA_metadata
    output:
        eccDNA_results

rule integration:
    input:
        circRNA_results,
        eccDNA_results,
        genome_annotation
    output:
        integrated_results
```

---

# 16. Snakemake 推荐结构

```
integration/
│
├── Snakefile
├── config/
│   ├── config.yaml
│   ├── circRNA.yaml
│   └── eccDNA.yaml
│
├── metadata/
│   ├── circRNA/
│   │   ├── Oryza_sativa.tsv
│   │   └── ...
│   │
│   └── eccDNA/
│       ├── Oryza_sativa.tsv
│       └── ...
│
├── rules/
│   ├── circRNA.smk
│   ├── eccDNA.smk
│   ├── annotation.smk
│   ├── overlap.smk
│   ├── expression.smk
│   └── integration.smk
│
├── scripts/
│   ├── map_circRNA_eccDNA.py
│   ├── annotate_TE.py
│   └── integrate_results.py
│
└── results/
```

---

# 17. 三个流程的职责

## Pipeline 1：circRNA

```
circrna.nf
```

负责：

```
RNA-seq
↓
circRNA detection
↓
annotation
↓
quantification
↓
differential analysis
```

---

## Pipeline 2：eccDNA

```
circdna.nf
```

扩展：

```
circdna.nf

├── shortread
│
├── pacbio
│
└── ont
```

负责：

```
eccDNA detection
↓
structure
↓
breakpoint
↓
copy number
↓
annotation
```

---

## Pipeline 3：Integration

```
Snakemake
```

负责：

```
circRNA
    │
    ├── genomic origin
    │
    ├── gene
    │
    ├── TE
    │
    └── expression
    │
    +
    │
eccDNA
    │
    ├── genomic origin
    │
    ├── gene
    │
    ├── TE
    │
    └── abundance
    │
    ▼
Integration
```

---

# 18. 最终建议的分析关系

建议重点研究：

```
Genome
   │
   ├──────────────┐
   │              │
   ▼              ▼
eccDNA          circRNA
   │              │
   │              │
   ▼              ▼
TE activity    Host gene
   │              │
   └──────┬───────┘
          ▼
      Expression
          │
          ▼
      Phenotype
```

尤其建议增加：

```
eccDNA ↔ circRNA
```

的三类关联：

### ① Genomic overlap

```
eccDNA genomic region
        ∩
circRNA host gene
```

---

### ② TE-mediated

```
TE-derived eccDNA
        │
        ▼
TE-associated circRNA
```

---

### ③ Stress response

```
Stress
 │
 ├── eccDNA formation
 │
 ├── TE activation
 │
 └── circRNA regulation
```

这个方向对于植物材料尤其有价值。

---

# 19. 建议的最终目录

```
plant_circle_pipeline/
│
├── circRNA/
│   └── circrna.nf
│
├── circDNA/
│   └── circdna.nf
│
├── integration/
│   └── Snakefile
│
├── metadata/
│   │
│   ├── circRNA/
│   │   ├── Oryza_sativa.tsv
│   │   ├── Triticum_aestivum.tsv
│   │   └── ...
│   │
│   └── eccDNA/
│       ├── Oryza_sativa.tsv
│       ├── Triticum_aestivum.tsv
│       └── ...
│
├── benchmark/
│   ├── shortread/
│   ├── pacbio/
│   ├── ont/
│   └── rna/
│
├── configs/
│
├── containers/
│
├── docs/
│
└── README.md
```

---

# 20. 我建议的开发优先级

### Phase 1：先稳定 circRNA

基于：

```
circrna.nf
```

完成：

```
RNA-seq metadata
↓
QC
↓
alignment
↓
circRNA detection
↓
annotation
↓
quantification
```

优先测试：

```
ERR10889820
ERR12724677
```

---

### Phase 2：完善 Illumina eccDNA

基于：

```
circdna.nf
```

完善：

```
Illumina
↓
ECCsplorer
↓
Circle-Map
↓
Consensus
```

测试：

```
ERR11535564
ERR6004146
CRR3168890
CRR3168885
CRR1082975
```

---

### Phase 3：加入 PacBio

增加：

```
DeConcat
TideHunter
CircleSeeker
```

重点测试：

```
SRR16958693
SRR26069818
```

然后加入：

```
ERR11838731
SRR30359583
```

作为 HiFi WGS 模式。

---

### Phase 4：加入 ONT

重点：

```
eccDNA enriched
↓
ecc_finder
↓
SV tools
```

测试：

```
ERR12724336
ERR6326020
SRR24335762
SRR24693334
SRR36603439
```

---

### Phase 5：Snakemake Integration

最后才做：

```
circRNA
    +
eccDNA
    +
Genome
    +
TE
```

不要一开始就做整合层。

---

# 21. 最终推荐的技术路线

我最终建议你把整个项目定义成：

> **Plant Circular Nucleic Acid Analysis Framework**

其中：

```
                    Plant Circular Nucleic Acid
                              │
               ┌──────────────┴──────────────┐
               │                             │
             circRNA                       eccDNA
               │                             │
         circrna.nf                     circdna.nf
               │                             │
         Illumina RNA              ┌─────────┴─────────┐
                                   │                   │
                               Illumina             Long-read
                                   │                   │
                               ECCsplorer       ┌───────┴───────┐
                                                │               │
                                             PacBio            ONT
                                                │               │
                                        DeConcat/TideHunter   ecc_finder
                                        CircleSeeker          SV tools
                                                │               │
                                                └───────┬───────┘
                                                        │
                                                        ▼
                                             Unified eccDNA table
                                                        │
                          ┌─────────────────────────────┴───────────────────────────┐
                          │                                                         │
                          ▼                                                         ▼
                     circRNA results                                         eccDNA results
                          │                                                         │
                          └──────────────────────┬──────────────────────────────────┘
                                                 ▼
                                           Snakemake
                                           Integration
                                                 │
                                  ┌──────────────┼──────────────┐
                                  ▼              ▼              ▼
                              Genome          TE          Expression
                                  │              │              │
                                  └──────────────┼──────────────┘
                                                 ▼
                                      Circular DNA-RNA atlas
```

**核心建议是：不要把 33 个数据集当成 33 个完全独立的流程开发对象，而是把它们作为一个跨平台 benchmark matrix。** 其中最重要的 6 个 benchmark 组合应当是：

| Benchmark       | 数据                      | 目的                      |
| --------------- | ------------------------- | ------------------------- |
| RNA-Seq         | ERR10889820               | circRNA.nf                |
| RNA-Seq         | ERR12724677               | circRNA.nf                |
| Illumina eccDNA | ERR11535564               | circDNA short-read        |
| PacBio CIDER    | SRR16958693               | PacBio circular consensus |
| PacBio RCA      | SRR26069818               | PacBio RCA eccDNA         |
| ONT eccDNA      | ERR12724336 / SRR24335762 | ONT eccDNA                |

这样做的好处是，最终可以明确回答三个问题：

1. **circRNA.nf 能否稳定处理不同植物 RNA-seq 数据？**
2. **circdna.nf 能否统一处理 Illumina、PacBio、ONT 三种 eccDNA 数据？**
3. **同一物种中，eccDNA 和 circRNA 是否存在共同的基因组来源、TE 来源或应激响应机制？**

其中，我认为目前最值得优先推进的是 **`circdna.nf` 的 PacBio/ONT 模式改造**。现有 `circdna.nf` 如果本身已经具备 Nextflow 的 `main.nf → subworkflows → modules` 架构，就没有必要另起炉灶；可以直接在现有架构上增加 `SHORTREAD / PACBIO / ONT` 三个 execution mode，并通过统一的 `eccDNA consensus TSV/BED/FASTA` 作为输出接口。这样后续的 Snakemake 整合层就只需要读取标准化结果，而不需要理解每种测序技术的底层差异。

如果按照这个方案实施，**下一步最合适的工作就是直接针对你现有的 `circdna.nf` 和 `circrna.nf` 仓库目录结构，逐个文件制定修改清单**，包括 `main.nf`、`nextflow.config`、`conf/`、`subworkflows/`、`modules/`、Docker/Conda 环境、metadata schema，以及 33 个 benchmark 数据如何映射到每一个 process。
