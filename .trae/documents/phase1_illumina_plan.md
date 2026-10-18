# 阶段一：二代测序（Illumina）流程打通与基础架构搭建 实施计划

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 完成 circdna.nf 主框架开发，打通 Illumina 数据的 Reference Mode、eccDNA Mode 与 Integrated Mode，并构建基础的 ECC_SCORE 评价评分系统（v1.0）。

**Architecture:** 基于现有 nf-core/circdna 代码库，采用 Nextflow DSL2 模块化架构。通过 `params.mode` 控制三种模式（reference / eccdna / integrated）的调度。新增 mosdepth 深度分析、ECCsplorer 模块、Candidate Merge 脚本和 ECC_SCORE 计算模块。将现有 trimgalore 替换为 fastp，BWA 替换为 BWA-MEM2 以提升性能。

**Tech Stack:** Nextflow DSL2, Python, BEDTools, nf-core modules (fastp, bwamem2, mosdepth, samtools, etc.)

---

## 仓库调研结论

### 现有架构

- **主入口**: [main.nf](file:///Users/siyangming/nextflow_nf_core/circdna.nf/main.nf) → 调用 `CIRCDNA` workflow
- **主工作流**: [workflows/circdna.nf](file:///Users/siyangming/nextflow_nf_core/circdna.nf/workflows/circdna.nf) - 包含输入检查、QC、比对、多个circDNA识别工具
- **子工作流**: `subworkflows/local/` 下已有 input_check, bam_preprocessing, circle_map_pipeline, circle_finder_pipeline, ampliconarchitect_pipeline, unicycler_pipeline
- **现有工具**: trimgalore (质控), BWA (比对), circexplorer2, circle_map, circle_finder, ampliconarchitect, unicycler
- **Samplesheet**: 支持 FASTQ 和 BAM 两种格式，已支持 lane 列实现增量缓存

### 可用 nf-core 模块（在 modules 仓库中）

- `fastp` - 二代测序质量控制
- `bwamem2/index` 和 `bwamem2/mem` - BWA-MEM2 快速比对
- `mosdepth` - 快速测序深度计算

### 需要新增

- `params.mode` 顶层路由控制
- Reference Mode (gDNA背景分析): mosdepth深度分析、重复序列注释加载
- ECCsplorer 模块
- Candidate Merge 脚本（合并 ECCsplorer + Circle-Map 结果）
- Integrated Mode: ECC_SCORE v1.0 计算、深度比计算、候选分级

---

## 文件结构总览

### 新增文件

| 文件路径                                  | 职责                                              |
| ----------------------------------------- | ------------------------------------------------- |
| `modules/nf-core/fastp/`                  | fastp 质控模块（从 nf-core/modules 复制）         |
| `modules/nf-core/bwamem2/`                | BWA-MEM2 索引+比对模块（从 nf-core/modules 复制） |
| `modules/nf-core/mosdepth/`               | mosdepth 深度计算模块（从 nf-core/modules 复制）  |
| `modules/local/eccsplorer/`               | ECCsplorer 本地模块                               |
| `modules/local/candidate_merge/`          | 候选区间合并模块（Python脚本）                    |
| `modules/local/ecc_score/`                | ECC_SCORE 计算模块（Python脚本）                  |
| `subworkflows/local/reference_mode/`      | Reference Mode 子工作流                           |
| `subworkflows/local/eccdna_mode/`         | eccDNA Mode 子工作流                              |
| `subworkflows/local/integrated_mode/`     | Integrated Mode 子工作流                          |
| `bin/merge_candidates.py`                 | 候选区间合并 Python 脚本                          |
| `bin/calculate_ecc_score.py`              | ECC_SCORE 计算 Python 脚本                        |
| `bin/annotate_repeats.py`                 | 重复序列注释脚本                                  |
| `.trae/documents/phase1_illumina_plan.md` | 本计划文档                                        |

### 修改文件

| 文件路径                                                                                                                              | 修改内容                                             |
| ------------------------------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------- |
| [workflows/circdna.nf](file:///Users/siyangming/nextflow_nf_core/circdna.nf/workflows/circdna.nf)                                     | 新增 mode 路由控制，调用三个 mode 子工作流           |
| [nextflow.config](file:///Users/siyangming/nextflow_nf_core/circdna.nf/nextflow.config)                                               | 新增 mode, repeat_gff, ecc_score_weights 等参数      |
| [conf/modules.config](file:///Users/siyangming/nextflow_nf_core/circdna.nf/conf/modules.config)                                       | 新增 fastp, bwamem2, mosdepth, eccsplorer 等模块配置 |
| [conf/base.config](file:///Users/siyangming/nextflow_nf_core/circdna.nf/conf/base.config)                                             | 新增进程资源配置（如需）                             |
| [conf/test_local.config](file:///Users/siyangming/nextflow_nf_core/circdna.nf/conf/test_local.config)                                 | 新增测试参数                                         |
| [bin/check_samplesheet.py](file:///Users/siyangming/nextflow_nf_core/circdna.nf/bin/check_samplesheet.py)                             | 新增 datatype, platform, protocol 字段支持           |
| [subworkflows/local/input_check/main.nf](file:///Users/siyangming/nextflow_nf_core/circdna.nf/subworkflows/local/input_check/main.nf) | 支持新字段解析                                       |
| [CHANGELOG.md](file:///Users/siyangming/nextflow_nf_core/circdna.nf/CHANGELOG.md)                                                     | 版本更新记录                                         |

---

## Sprint 1: 基础框架与 Reference Mode

### Task 1.1: Samplesheet 扩展与多模式解析器

**Files:**

- Modify: `bin/check_samplesheet.py`
- Modify: `subworkflows/local/input_check/main.nf`
- Modify: `assets/samplesheet.csv`

**目标:** 支持多模式 samplesheet，新增 `datatype`、`platform`、`protocol` 等字段

**步骤:**

- [ ] **Step 1: 修改 check_samplesheet.py 支持新字段**

在 [check_samplesheet.py](file:///Users/siyangming/nextflow_nf_core/circdna.nf/bin/check_samplesheet.py) 中扩展：

1. 新增可选字段：`datatype`（gdna / eccdna）、`platform`（illumina / pacbio / ont）、`protocol`（short_read / long_read）
2. 保留向后兼容：缺失字段时使用默认值（datatype=eccdna, platform=illumina, protocol=short_read）
3. 输出文件增加对应列

```python
def check_samplesheet(file_in, file_out, input_format):
    # 扩展 HEADER 定义，支持可选列
    HEADER_FASTQ = ["sample", "fastq_1", "fastq_2"]
    OPTIONAL_HEADERS = ["lane", "datatype", "platform", "protocol"]

    # 检测可选列是否存在
    # 验证 datatype 取值为 gdna / eccdna
    # 验证 platform 取值为 illumina / pacbio / ont
    # 验证 protocol 取值为 short_read / long_read
    # 输出时保留所有输入列
```

- [ ] **Step 2: 修改 input_check/main.nf 解析新字段**

在 [input_check/main.nf](file:///Users/siyangming/nextflow_nf_core/circdna.nf/subworkflows/local/input_check/main.nf) 的 `create_fastq_channels` 函数中：

```groovy
def create_fastq_channels(LinkedHashMap row) {
    def meta = [:]
    meta.id           = row.sample
    meta.single_end   = row.containsKey('single_end') ? (row.single_end ? row.single_end.toBoolean() : false) : (!row.fastq_2 || row.fastq_2.isEmpty())
    meta.datatype     = row.containsKey('datatype') ? row.datatype : 'eccdna'
    meta.platform     = row.containsKey('platform') ? row.platform : 'illumina'
    meta.protocol     = row.containsKey('protocol') ? row.protocol : 'short_read'
    if (row.containsKey('lane') && row.lane) {
        meta.lane = row.lane
    }
    // ... rest of function
}
```

- [ ] **Step 3: 更新 assets/samplesheet.csv 示例**

更新示例 samplesheet，展示新字段用法：

```csv
sample,fastq_1,fastq_2,datatype,platform,protocol
SAMPLE_GDNA,gdna_R1.fastq.gz,gdna_R2.fastq.gz,gdna,illumina,short_read
SAMPLE_ECCDNA,eccdna_R1.fastq.gz,eccdna_R2.fastq.gz,eccdna,illumina,short_read
```

- [ ] **Step 4: 验证 samplesheet 解析**

运行测试验证 samplesheet 解析正确：

```bash
cd /Users/siyangming/nextflow_nf_core/circdna.nf
python bin/check_samplesheet.py testdatasets/samplesheet/samplesheet_local.csv /tmp/test_out.csv FASTQ
```

Expected: 输出文件包含正确的列，缺失字段有默认值

---

### Task 1.2: 顶层路由控制结构（params.mode）

**Files:**

- Modify: `workflows/circdna.nf`
- Modify: `nextflow.config`

**目标:** 建立 `params.mode` 路由机制，支持 reference / eccdna / integrated 三种模式

**步骤:**

- [ ] **Step 1: 在 nextflow.config 中新增 mode 参数**

在 [nextflow.config](file:///Users/siyangming/nextflow_nf_core/circdna.nf/nextflow.config) 的 params 块中添加：

```groovy
// Mode options
mode                       = 'eccdna'     // 'reference', 'eccdna', 'integrated'
repeat_gff                 = null          // EDTA/RepeatMasker GFF3/BED 文件路径
ecc_score_w1               = 1.0           // Junction reads 权重
ecc_score_w2               = 1.0           // Depth ratio 权重
ecc_score_w3               = 0.5           // TE repeat penalty 权重
```

- [ ] **Step 2: 在 circdna.nf 中新增 mode 校验和路由逻辑**

在 [workflows/circdna.nf](file:///Users/siyangming/nextflow_nf_core/circdna.nf/workflows/circdna.nf) 开头添加 mode 校验：

```groovy
// Mode validation
def valid_modes = ['reference', 'eccdna', 'integrated']
if (!valid_modes.contains(params.mode)) {
    exit 1, "Invalid mode: ${params.mode}. Valid modes: ${valid_modes.join(', ')}"
}

// Integrated mode requires both gdna and eccdna samples
if (params.mode == 'integrated') {
    // 将在后续 Task 中实现样本分离逻辑
}
```

- [ ] **Step 3: 搭建 mode 子工作流调用框架**

在 `CIRCDNA` workflow 中，根据 mode 条件调用不同子工作流：

```groovy
// 按 datatype 分离样本
if (params.mode == 'reference') {
    // 仅处理 gdna 样本
    REFERENCE_MODE ( ... )
} else if (params.mode == 'eccdna') {
    // 仅处理 eccdna 样本
    ECCDNA_MODE ( ... )
} else if (params.mode == 'integrated') {
    // 同时处理 gdna 和 eccdna 样本
    REFERENCE_MODE ( ... )
    ECCDNA_MODE ( ... )
    INTEGRATED_MODE ( ... )
}
```

注意：此步骤仅搭建框架，子工作流实现在后续 Sprint 中完成。

---

### Task 1.3: fastp 替换 trimgalore

**Files:**

- Create: `modules/nf-core/fastp/` (从 nf-core/modules 复制)
- Modify: `workflows/circdna.nf`
- Modify: `conf/modules.config`
- Modify: `modules.json`

**目标:** 用 fastp 替换 trimgalore 作为二代测序质控工具

**步骤:**

- [ ] **Step 1: 从 nf-core/modules 复制 fastp 模块**

```bash
cd /Users/siyangming/nextflow_nf_core/circdna.nf
# 从 modules 仓库复制 fastp 模块
cp -r ../modules/modules/nf-core/fastp modules/nf-core/
```

- [ ] **Step 2: 修改 workflows/circdna.nf 引入 fastp**

替换 TRIMGALORE 相关的 include 和调用：

```groovy
// 替换 include
// include { TRIMGALORE } from '../modules/nf-core/trimgalore/main'
include { FASTP } from '../modules/nf-core/fastp/main'

// 替换调用逻辑
if (!params.skip_trimming) {
    FASTP (
        ch_cat_fastq
    )
    ch_trimmed_reads = FASTP.out.reads
    ch_fastp_html = FASTP.out.html
    ch_fastp_json = FASTP.out.json
    ch_versions = ch_versions.mix(FASTP.out.versions)
}
```

- [ ] **Step 3: 更新 modules.config**

添加 fastp 配置，移除 trimgalore 配置：

```groovy
process {
    withName: 'FASTP' {
        ext.args = '--qualified_quality_phred 15 --length_required 50'
        publishDir = [
            [
                path: { "${params.outdir}/fastp" },
                mode: params.publish_dir_mode,
                pattern: '*.{html,json}',
                saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
            ],
            [
                path: { "${params.outdir}/fastp" },
                mode: params.publish_dir_mode,
                pattern: '*.fq.gz',
                saveAs: { filename -> filename.equals('versions.yml') ? null : filename },
                enabled: params.save_trimmed
            ]
        ]
    }
}
```

- [ ] **Step 4: 更新 modules.json**

运行 `nf-core modules list` 或手动更新 modules.json

---

### Task 1.4: BWA-MEM2 替换 BWA

**Files:**

- Create: `modules/nf-core/bwamem2/` (从 nf-core/modules 复制)
- Modify: `subworkflows/local/bam_preprocessing/main.nf`
- Modify: `workflows/circdna.nf`
- Modify: `conf/modules.config`

**目标:** 用 BWA-MEM2 替换 BWA 提升比对速度

**步骤:**

- [ ] **Step 1: 从 nf-core/modules 复制 bwamem2 模块**

```bash
cd /Users/siyangming/nextflow_nf_core/circdna.nf
cp -r ../modules/modules/nf-core/bwamem2 modules/nf-core/
```

- [ ] **Step 2: 修改 bam_preprocessing subworkflow**

在 [bam_preprocessing/main.nf](file:///Users/siyangming/nextflow_nf_core/circdna.nf/subworkflows/local/bam_preprocessing/main.nf) 中：

```groovy
// 替换 include
// include { BWA_MEM } from '../../../modules/nf-core/bwa/mem/main'
include { BWAMEM2_MEM } from '../../../modules/nf-core/bwamem2/mem/main'

// 替换调用
if (run_bwa) {
    BWAMEM2_MEM (
        trimmed_reads,
        bwa_index,  // 注意：bwamem2 的 index 格式不同
        fasta_meta,
        channel.value(true)
    )
    ch_bam_sorted = BWAMEM2_MEM.out.bam
    // ...
}
```

- [ ] **Step 3: 更新 BWA INDEX 为 BWAMEM2 INDEX**

在主工作流中替换 BWA_INDEX 为 BWAMEM2_INDEX

- [ ] **Step 4: 更新 modules.config**

添加 BWAMEM2 相关配置

---

### Task 1.5: mosdepth 深度分析模块

**Files:**

- Create: `modules/nf-core/mosdepth/` (从 nf-core/modules 复制)
- Modify: `conf/modules.config`

**目标:** 集成 mosdepth 用于快速全基因组测序深度计算

**步骤:**

- [ ] **Step 1: 从 nf-core/modules 复制 mosdepth 模块**

```bash
cd /Users/siyangming/nextflow_nf_core/circdna.nf
cp -r ../modules/modules/nf-core/mosdepth modules/nf-core/
```

- [ ] **Step 2: 在 modules.config 中添加 mosdepth 配置**

```groovy
process {
    withName: 'MOSDEPTH' {
        ext.args = '--fast-mode --by 1000'  // 1kb sliding window
        publishDir = [
            path: { "${params.outdir}/mosdepth/${meta.id}" },
            mode: params.publish_dir_mode,
            saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
        ]
    }
}
```

---

### Task 1.6: Reference Mode 子工作流

**Files:**

- Create: `subworkflows/local/reference_mode/main.nf`
- Create: `subworkflows/local/reference_mode/meta.yml`
- Modify: `workflows/circdna.nf`

**目标:** 实现 Reference Mode (gDNA 背景分析): QC + Alignment + Coverage Matrix + Genomic Context

**步骤:**

- [ ] **Step 1: 创建 reference_mode subworkflow 骨架**

```groovy
// subworkflows/local/reference_mode/main.nf
include { MOSDEPTH } from '../../../modules/nf-core/mosdepth/main'

workflow REFERENCE_MODE {
    take:
    reads          // channel: [ val(meta), [ reads ] ]
    bwa_index      // channel: [ "bwamem2_index", index_dir ]
    fasta_meta     // channel: [ val(meta), path(fasta) ]
    repeat_gff     // channel: [ path(repeat_gff) ] 或 channel.empty()

    main:
    ch_versions = channel.empty()

    // 调用 BAM 预处理（比对、排序、去重）
    BAM_PREPROCESSING (
        reads,
        bwa_index,
        fasta_meta,
        true
    )
    ch_versions = ch_versions.mix(BAM_PREPROCESSING.out.versions)

    // mosdepth 深度分析
    MOSDEPTH (
        BAM_PREPROCESSING.out.bam_sorted.join(BAM_PREPROCESSING.out.bam_sorted_bai),
        BAM_PREPROCESSING.out.fasta_fai
    )
    ch_versions = ch_versions.mix(MOSDEPTH.out.versions)

    emit:
    bam_sorted        = BAM_PREPROCESSING.out.bam_sorted
    bam_sorted_bai    = BAM_PREPROCESSING.out.bam_sorted_bai
    mosdepth_bed      = MOSDEPTH.out.regions_bed  // per-base depth bed
    mosdepth_summary  = MOSDEPTH.out.summary
    versions          = ch_versions
}
```

- [ ] **Step 2: 在主工作流中集成 REFERENCE_MODE 调用**

在 [workflows/circdna.nf](file:///Users/siyangming/nextflow_nf_core/circdna.nf/workflows/circdna.nf) 中：

```groovy
include { REFERENCE_MODE } from '../subworkflows/local/reference_mode/main'
```

并在 mode 路由逻辑中调用

- [ ] **Step 3: 添加重复序列注释加载功能**

创建 Python 脚本加载 EDTA/RepeatMasker GFF3/BED：

```python
# bin/annotate_repeats.py
#!/usr/bin/env python
"""加载重复序列注释并生成 TE 区域 BED 文件"""
import sys
import argparse

def parse_gff(gff_file, output_bed):
    """解析 GFF3 文件，输出 TE 区域 BED"""
    with open(gff_file) as fin, open(output_bed, 'w') as fout:
        for line in fin:
            if line.startswith('#'):
                continue
            cols = line.strip().split('\t')
            if len(cols) < 9:
                continue
            chrom = cols[0]
            start = cols[3]
            end = cols[4]
            strand = cols[6]
            feature_type = cols[2]
            # 提取重复序列类型
            attrs = cols[8]
            # 输出 BED6 格式
            fout.write(f"{chrom}\t{start}\t{end}\t{feature_type}\t0\t{strand}\n")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('input_gff', help='Input GFF3 file')
    parser.add_argument('output_bed', help='Output BED file')
    args = parser.parse_args()
    parse_gff(args.input_gff, args.output_bed)

if __name__ == '__main__':
    main()
```

- [ ] **Step 4: 验证 Reference Mode 运行**

使用测试数据运行 Reference Mode：

```bash
nextflow run main.nf -profile test_local,docker --mode reference --input test_reference.csv
```

Expected: 成功运行至 mosdepth 步骤，输出深度文件

---

## Sprint 2: Illumina eccDNA 工具链

### Task 2.1: ECCsplorer 模块封装

**Files:**

- Create: `modules/local/eccsplorer/main.nf`
- Create: `modules/local/eccsplorer/meta.yml`
- Create: `modules/local/eccsplorer/environment.yml`
- Modify: `conf/modules.config`

**目标:** 将 ECCsplorer 封装为独立的 Nextflow Process，处理 Junction 读长与 Reference/Non-reference 环状候选

**步骤:**

- [ ] **Step 1: 创建 ECCsplorer 模块目录和 environment.yml**

```yaml
# modules/local/eccsplorer/environment.yml
name: eccsplorer
channels:
  - bioconda
  - conda-forge
dependencies:
  - python=3.9
  - pysam
  - bedtools
  - numpy
```

- [ ] **Step 2: 编写 ECCsplorer main.nf**

```groovy
process ECCSPLORER {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "quay.io/biocontainers/eccsplorer:latest"  // 如无则自建

    input:
    tuple val(meta), path(bam), path(bai)
    path fasta

    output:
    tuple val(meta), path("*_candidates.bed"), emit: candidates_bed
    tuple val(meta), path("*_junction_reads.txt"), emit: junction_reads
    path "versions.yml", emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # ECCsplorer 分析脚本
    # 识别 junction reads 和 circular DNA candidates
    # 输出 BED 格式的候选区间

    echo 'eccsplorer 0.1.0' > versions.yml

    # 实际调用 eccsplorer 或自定义脚本
    python \$NXF_SCRIPTS/eccsplorer.py \\
        --bam ${bam} \\
        --fasta ${fasta} \\
        --output ${prefix}
    """
}
```

注意：ECCsplorer 具体实现需要根据工具文档调整。如果 eccsplorer 没有现成容器，需要先确认是否在 bioconda 上有包。

- [ ] **Step 3: 创建 meta.yml**

```yaml
name: eccsplorer
description: Detects extrachromosomal circular DNA candidates from BAM files
keywords:
  - circular dna
  - eccdna
  - junction reads
tools:
  - eccsplorer:
      description: ECCsplorer detects eccDNA candidates using junction reads
      homepage: https://github.com/MWSchmid/ECCsplorer
      documentation: https://github.com/MWSchmid/ECCsplorer
      licence: ["MIT"]
input:
  - meta:
      type: map
      description: Sample metadata
  - bam:
      type: file
      description: Sorted BAM file
  - fasta:
      type: file
      description: Reference genome FASTA
output:
  - candidates_bed:
      type: file
      description: BED file with eccDNA candidates
  - versions:
      type: file
      description: File containing software versions
authors:
  - "@siyangming"
```

- [ ] **Step 4: 在 modules.config 中添加配置**

```groovy
process {
    withName: 'ECCSPLORER' {
        ext.args = ''
        ext.prefix = { "${meta.id}.eccsplorer" }
        publishDir = [
            path: { "${params.outdir}/eccsplorer/${meta.id}" },
            mode: params.publish_dir_mode,
            saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
        ]
    }
}
```

---

### Task 2.2: Circle-Map 结果标准化

**Files:**

- Modify: `subworkflows/local/circle_map_pipeline/main.nf`
- Create: `bin/standardize_circle_map.py` (如需)

**目标:** 确保 Circle-Map 输出标准化的 BED 格式，便于后续合并

**步骤:**

- [ ] **Step 1: 检查 Circle-Map 现有输出格式**

查看现有 circle_map 模块的输出，确认是否为 BED 格式

- [ ] **Step 2: 如需要，添加标准化脚本**

如果 Circle-Map 输出格式不是标准 BED，编写转换脚本

- [ ] **Step 3: 更新 circle_map_pipeline emit**

确保 emit 标准化的 BED 文件（含 junction reads 计数）

---

### Task 2.3: eccDNA Mode 子工作流

**Files:**

- Create: `subworkflows/local/eccdna_mode/main.nf`
- Create: `subworkflows/local/eccdna_mode/meta.yml`
- Modify: `workflows/circdna.nf`

**目标:** 实现 eccDNA Mode: QC + Alignment + ECCsplorer + Circle-Map

**步骤:**

- [ ] **Step 1: 创建 eccdna_mode subworkflow**

```groovy
// subworkflows/local/eccdna_mode/main.nf
include { ECCSPLORER } from '../../../modules/local/eccsplorer/main'
include { CIRCLE_MAP_PIPELINE } from '../circle_map_pipeline/main'
include { MOSDEPTH } from '../../../modules/nf-core/mosdepth/main'

workflow ECCDNA_MODE {
    take:
    reads
    bwa_index
    fasta_meta
    run_eccsplorer
    run_circle_map

    main:
    ch_versions = channel.empty()

    // BAM 预处理
    BAM_PREPROCESSING (
        reads,
        bwa_index,
        fasta_meta,
        true
    )
    ch_versions = ch_versions.mix(BAM_PREPROCESSING.out.versions)

    // mosdepth 深度分析（eccDNA 样本）
    MOSDEPTH (
        BAM_PREPROCESSING.out.bam_sorted.join(BAM_PREPROCESSING.out.bam_sorted_bai),
        BAM_PREPROCESSING.out.fasta_fai
    )
    ch_versions = ch_versions.mix(MOSDEPTH.out.versions)

    // ECCsplorer
    if (run_eccsplorer) {
        ECCSPLORER (
            BAM_PREPROCESSING.out.bam_sorted.join(BAM_PREPROCESSING.out.bam_sorted_bai).map { meta, bam, bai -> [meta, bam, bai] },
            fasta_meta.map { _meta, fasta -> fasta }
        )
        ch_versions = ch_versions.mix(ECCSPLORER.out.versions)
    }

    // Circle-Map
    if (run_circle_map) {
        CIRCLE_MAP_PIPELINE (
            BAM_PREPROCESSING.out.bam_sorted,
            BAM_PREPROCESSING.out.bam_sorted_bai,
            BAM_PREPROCESSING.out.fasta_fai,
            true,  // run_realign
            false  // run_repeats
        )
        ch_versions = ch_versions.mix(CIRCLE_MAP_PIPELINE.out)
    }

    emit:
    bam_sorted        = BAM_PREPROCESSING.out.bam_sorted
    bam_sorted_bai    = BAM_PREPROCESSING.out.bam_sorted_bai
    mosdepth_bed      = MOSDEPTH.out.regions_bed
    mosdepth_summary  = MOSDEPTH.out.summary
    eccsplorer_bed    = run_eccsplorer ? ECCSPLORER.out.candidates_bed : channel.empty()
    eccsplorer_junc   = run_eccsplorer ? ECCSPLORER.out.junction_reads : channel.empty()
    circle_map_bed    = run_circle_map ? CIRCLE_MAP_PIPELINE.out.bed : channel.empty()
    versions          = ch_versions
}
```

- [ ] **Step 2: 在主工作流中集成 ECCDNA_MODE 调用**

- [ ] **Step 3: 验证 eccDNA Mode 运行**

---

## Sprint 3: 二代候选区 Merge

### Task 3.1: Candidate Merge Python 脚本

**Files:**

- Create: `bin/merge_candidates.py`
- Create: `modules/local/candidate_merge/main.nf`
- Create: `modules/local/candidate_merge/meta.yml`

**目标:** 合并 ECCsplorer 与 Circle-Map 产出，进行区间 Merge 与重叠度整理

**步骤:**

- [ ] **Step 1: 编写 merge_candidates.py**

```python
#!/usr/bin/env python
"""
合并多个 eccDNA 候选区间 BED 文件
支持 ECCsplorer, Circle-Map 等多个工具的结果合并
"""
import argparse
import sys
import os

def parse_bed(bed_file, source_name):
    """解析 BED 文件，返回区间列表
    BED6+ format: chrom, start, end, name, score, strand, [extra...]
    """
    candidates = []
    with open(bed_file) as f:
        for line in f:
            if line.startswith('#') or line.strip() == '':
                continue
            cols = line.strip().split('\t')
            if len(cols) < 3:
                continue
            chrom = cols[0]
            start = int(cols[1])
            end = int(cols[2])
            name = cols[3] if len(cols) > 3 else f"{source_name}_{len(candidates)+1}"
            score = cols[4] if len(cols) > 4 else '0'
            strand = cols[5] if len(cols) > 5 else '.'
            candidates.append({
                'chrom': chrom,
                'start': start,
                'end': end,
                'name': name,
                'score': score,
                'strand': strand,
                'source': source_name,
                'junction_reads': int(score) if score.isdigit() else 0
            })
    return candidates

def merge_candidates(candidates_list, max_distance=100):
    """
    合并重叠或邻近的候选区间
    使用 BEDTools merge 风格的合并逻辑
    """
    if not candidates_list:
        return []

    # 按染色体和起始位置排序
    all_candidates = []
    for source, cands in candidates_list:
        for c in cands:
            c['source'] = source
            all_candidates.append(c)

    all_candidates.sort(key=lambda x: (x['chrom'], x['start']))

    merged = []
    current = None
    sources = set()
    total_junction = 0

    for cand in all_candidates:
        if current is None:
            current = cand.copy()
            sources = {cand['source']}
            total_junction = cand.get('junction_reads', 0)
        elif (cand['chrom'] == current['chrom'] and
              cand['start'] <= current['end'] + max_distance):
            # 重叠或邻近，合并
            current['end'] = max(current['end'], cand['end'])
            sources.add(cand['source'])
            total_junction += cand.get('junction_reads', 0)
        else:
            # 新区间
            current['num_sources'] = len(sources)
            current['sources'] = ','.join(sorted(sources))
            current['junction_reads'] = total_junction
            merged.append(current)
            current = cand.copy()
            sources = {cand['source']}
            total_junction = cand.get('junction_reads', 0)

    if current is not None:
        current['num_sources'] = len(sources)
        current['sources'] = ','.join(sorted(sources))
        current['junction_reads'] = total_junction
        merged.append(current)

    return merged

def write_bed(merged, output_file):
    """输出合并后的 BED 文件"""
    with open(output_file, 'w') as f:
        f.write('#chrom\tstart\tend\tname\tscore\tstrand\tsources\tnum_tools\tjunction_reads\n')
        for i, m in enumerate(merged):
            name = f"ecc_candidate_{i+1}"
            score = m['junction_reads']
            f.write(f"{m['chrom']}\t{m['start']}\t{m['end']}\t{name}\t{score}\t{m['strand']}\t{m['sources']}\t{m['num_sources']}\t{m['junction_reads']}\n")

def main():
    parser = argparse.ArgumentParser(description='Merge eccDNA candidates from multiple tools')
    parser.add_argument('--eccsplorer', help='ECCsplorer BED output')
    parser.add_argument('--circle_map', help='Circle-Map BED output')
    parser.add_argument('--output', required=True, help='Output merged BED file')
    parser.add_argument('--max-distance', type=int, default=100,
                       help='Maximum distance for merging adjacent candidates (default: 100)')
    args = parser.parse_args()

    candidates_list = []

    if args.eccsplorer and os.path.exists(args.eccsplorer):
        cands = parse_bed(args.eccsplorer, 'ECCsplorer')
        candidates_list.append(('ECCsplorer', cands))

    if args.circle_map and os.path.exists(args.circle_map):
        cands = parse_bed(args.circle_map, 'Circle-Map')
        candidates_list.append(('Circle-Map', cands))

    if not candidates_list:
        print("Warning: No input files found", file=sys.stderr)
        with open(args.output, 'w') as f:
            f.write('# No candidates\n')
        return

    merged = merge_candidates(candidates_list, args.max_distance)
    write_bed(merged, args.output)
    print(f"Merged {sum(len(c) for _, c in candidates_list)} candidates into {len(merged)} consensus candidates")

if __name__ == '__main__':
    main()
```

- [ ] **Step 2: 创建 candidate_merge Nextflow 模块**

```groovy
process CANDIDATE_MERGE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/../../environment.yml"  // 共享 bedtools 环境
    container "quay.io/biocontainers/bedtools:2.31.0--hf5e1c6e_1"

    input:
    tuple val(meta), path(eccsplorer_bed), path(circle_map_bed)

    output:
    tuple val(meta), path("*_merged_candidates.bed"), emit: merged_bed
    path "versions.yml", emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo 'candidate_merge 1.0' > versions.yml

    python ${projectDir}/bin/merge_candidates.py \\
        --eccsplorer ${eccsplorer_bed} \\
        --circle_map ${circle_map_bed} \\
        --output ${prefix}_merged_candidates.bed
    """
}
```

- [ ] **Step 3: 单元测试 merge_candidates.py**

```bash
cd /Users/siyangming/nextflow_nf_core/circdna.nf
python bin/merge_candidates.py --help
```

---

### Task 3.2: 在 eccDNA Mode 中集成分并

**Files:**

- Modify: `subworkflows/local/eccdna_mode/main.nf`

**步骤:**

- [ ] **Step 1: 在 eccdna_mode 中调用 CANDIDATE_MERGE**

```groovy
include { CANDIDATE_MERGE } from '../../../modules/local/candidate_merge/main'

// 在 ECCsplorer 和 Circle-Map 运行后
ch_merge_input = ECCSPLORER.out.candidates_bed
    .join(CIRCLE_MAP_PIPELINE.out.bed)
    .map { meta, ecc_bed, cm_bed -> [meta, ecc_bed, cm_bed] }

CANDIDATE_MERGE (ch_merge_input)
ch_versions = ch_versions.mix(CANDIDATE_MERGE.out.versions)
```

- [ ] **Step 2: 更新 emit 输出 merged_bed**

---

## Sprint 4: Integrated v1.0 与测试验证

### Task 4.1: 深度比计算 (Depth Ratio)

**Files:**

- Create: `bin/calculate_depth_ratio.py`
- Create: `modules/local/depth_ratio/main.nf`

**目标:** 计算 eccDNA 样本与 gDNA 样本的深度比 log2(eccDNA_depth + 1 / gDNA_depth + 1)

**步骤:**

- [ ] **Step 1: 编写 calculate_depth_ratio.py**

```python
#!/usr/bin/env python
"""
计算 eccDNA 与 gDNA 样本的测序深度比
输入: mosdepth 输出的 per-base bed 或 window bed
输出: 深度比 BED 文件
"""
import argparse
import sys

def load_depth(bed_file):
    """加载 mosdepth bed 文件，返回 {chrom: {pos: depth}}"""
    depths = {}
    with open(bed_file) as f:
        for line in f:
            if line.startswith('#') or line.strip() == '':
                continue
            cols = line.strip().split('\t')
            chrom = cols[0]
            start = int(cols[1])
            end = int(cols[2])
            depth = float(cols[3]) if len(cols) > 3 else 0.0
            if chrom not in depths:
                depths[chrom] = []
            depths[chrom].append((start, end, depth))
    return depths

def calculate_ratio(eccdna_depths, gdna_depths, output_file):
    """计算 log2 深度比"""
    import math

    with open(output_file, 'w') as fout:
        fout.write('#chrom\tstart\tend\teccdna_depth\tgdna_depth\tlog2_ratio\n')

        # 遍历所有染色体
        all_chroms = set(eccdna_depths.keys()) | set(gdna_depths.keys())

        for chrom in sorted(all_chroms):
            ecc_regions = eccdna_depths.get(chrom, [])
            gdna_regions = gdna_depths.get(chrom, [])

            # 简化：假设窗口对齐，直接逐行计算
            # 实际生产环境需要更复杂的区间匹配逻辑
            for (e_start, e_end, e_depth), (g_start, g_end, g_depth) in zip(ecc_regions, gdna_regions):
                if e_start != g_start or e_end != g_end:
                    continue  # 跳过不对齐的窗口
                ratio = math.log2((e_depth + 1) / (g_depth + 1))
                fout.write(f"{chrom}\t{e_start}\t{e_end}\t{e_depth:.2f}\t{g_depth:.2f}\t{ratio:.4f}\n")

def main():
    parser = argparse.ArgumentParser(description='Calculate depth ratio between eccDNA and gDNA')
    parser.add_argument('--eccdna', required=True, help='eccDNA mosdepth BED file')
    parser.add_argument('--gdna', required=True, help='gDNA mosdepth BED file')
    parser.add_argument('--output', required=True, help='Output depth ratio BED file')
    args = parser.parse_args()

    eccdna_depths = load_depth(args.eccdna)
    gdna_depths = load_depth(args.gdna)

    calculate_ratio(eccdna_depths, gdna_depths, args.output)

if __name__ == '__main__':
    main()
```

- [ ] **Step 2: 创建 depth_ratio Nextflow 模块**

---

### Task 4.2: ECC_SCORE v1.0 计算模块

**Files:**

- Create: `bin/calculate_ecc_score.py`
- Create: `modules/local/ecc_score/main.nf`

**目标:** 实现 ECC*SCORE v1.0:
$$\text{ECC_SCORE}*{\text{v1}} = w*1 \cdot \text{Junction_Reads} + w_2 \cdot \log_2\left(\frac{\text{Depth}*{\text{eccDNA}} + 1}{\text{Depth}\_{\text{gDNA}} + 1}\right) - w_3 \cdot \text{TE_Repeat_Penalty}$$

**步骤:**

- [ ] **Step 1: 编写 calculate_ecc_score.py**

```python
#!/usr/bin/env python
"""
ECC_SCORE v1.0 计算
对每个候选 eccDNA 区间计算综合评分
"""
import argparse
import sys
import os

def load_repeat_regions(repeat_bed):
    """加载 TE/重复序列区域"""
    repeats = {}
    if not repeat_bed or not os.path.exists(repeat_bed):
        return repeats

    with open(repeat_bed) as f:
        for line in f:
            if line.startswith('#') or line.strip() == '':
                continue
            cols = line.strip().split('\t')
            chrom = cols[0]
            start = int(cols[1])
            end = int(cols[2])
            te_type = cols[3] if len(cols) > 3 else 'Unknown'
            if chrom not in repeats:
                repeats[chrom] = []
            repeats[chrom].append((start, end, te_type))
    return repeats

def calculate_te_overlap(cand_chrom, cand_start, cand_end, repeats):
    """计算候选区间与 TE 区域的重叠比例"""
    if cand_chrom not in repeats:
        return 0.0, 'None'

    cand_length = cand_end - cand_start
    overlap_total = 0
    te_types = set()

    for te_start, te_end, te_type in repeats[cand_chrom]:
        overlap_start = max(cand_start, te_start)
        overlap_end = min(cand_end, te_end)
        if overlap_start < overlap_end:
            overlap_total += overlap_end - overlap_start
            te_types.add(te_type)

    overlap_ratio = overlap_total / cand_length if cand_length > 0 else 0
    return overlap_ratio, ','.join(sorted(te_types)) if te_types else 'None'

def get_depth_for_region(depth_file, chrom, start, end):
    """从深度文件中获取指定区间的平均深度"""
    # 简化实现：实际需要根据 mosdepth 输出格式调整
    total_depth = 0
    total_bases = 0

    with open(depth_file) as f:
        for line in f:
            if line.startswith('#') or line.strip() == '':
                continue
            cols = line.strip().split('\t')
            d_chrom = cols[0]
            d_start = int(cols[1])
            d_end = int(cols[2])
            depth = float(cols[3]) if len(cols) > 3 else 0.0

            if d_chrom != chrom:
                continue
            if d_end <= start:
                continue
            if d_start >= end:
                break

            overlap_start = max(start, d_start)
            overlap_end = min(end, d_end)
            overlap_bases = overlap_end - overlap_start
            total_depth += depth * overlap_bases
            total_bases += overlap_bases

    return total_depth / total_bases if total_bases > 0 else 0.0

def main():
    parser = argparse.ArgumentParser(description='Calculate ECC_SCORE v1.0 for eccDNA candidates')
    parser.add_argument('--candidates', required=True, help='Merged candidates BED file')
    parser.add_argument('--eccdna-depth', required=True, help='eccDNA sample mosdepth BED')
    parser.add_argument('--gdna-depth', required=True, help='gDNA sample mosdepth BED')
    parser.add_argument('--repeat-bed', help='TE/Repeat BED file for penalty calculation')
    parser.add_argument('--output', required=True, help='Output scored BED file')
    parser.add_argument('--w1', type=float, default=1.0, help='Weight for junction reads (default: 1.0)')
    parser.add_argument('--w2', type=float, default=1.0, help='Weight for depth ratio (default: 1.0)')
    parser.add_argument('--w3', type=float, default=0.5, help='Weight for TE repeat penalty (default: 0.5)')
    args = parser.parse_args()

    import math

    repeats = load_repeat_regions(args.repeat_bed)

    with open(args.candidates) as fin, open(args.output, 'w') as fout:
        header = fin.readline()
        fout.write('#chrom\tstart\tend\tname\tscore\tstrand\tsources\tnum_tools\tjunction_reads\teccdna_depth\tgdna_depth\tlog2_depth_ratio\tte_overlap_ratio\tte_types\tecc_score\tgrade\n')

        for line in fin:
            if line.startswith('#') or line.strip() == '':
                continue
            cols = line.strip().split('\t')
            chrom = cols[0]
            start = int(cols[1])
            end = int(cols[2])
            name = cols[3]
            strand = cols[5] if len(cols) > 5 else '.'
            sources = cols[6] if len(cols) > 6 else 'Unknown'
            num_tools = cols[7] if len(cols) > 7 else '1'
            junction_reads = int(cols[8]) if len(cols) > 8 and cols[8].isdigit() else 0

            # 计算深度比
            ecc_depth = get_depth_for_region(args.eccdna_depth, chrom, start, end)
            gdna_depth = get_depth_for_region(args.gdna_depth, chrom, start, end)
            log2_ratio = math.log2((ecc_depth + 1) / (gdna_depth + 1))

            # 计算 TE 重叠惩罚
            te_overlap, te_types = calculate_te_overlap(chrom, start, end, repeats)
            te_penalty = te_overlap * 10  # 归一化惩罚项

            # 计算 ECC_SCORE
            ecc_score = (args.w1 * junction_reads +
                        args.w2 * log2_ratio -
                        args.w3 * te_penalty)

            # 分级
            if ecc_score >= 10:
                grade = 'High'
            elif ecc_score >= 5:
                grade = 'Medium'
            else:
                grade = 'Low'

            fout.write(f"{chrom}\t{start}\t{end}\t{name}\t{ecc_score:.2f}\t{strand}\t{sources}\t{num_tools}\t{junction_reads}\t{ecc_depth:.2f}\t{gdna_depth:.2f}\t{log2_ratio:.4f}\t{te_overlap:.4f}\t{te_types}\t{ecc_score:.2f}\t{grade}\n")

if __name__ == '__main__':
    main()
```

- [ ] **Step 2: 创建 ecc_score Nextflow 模块**

---

### Task 4.3: Integrated Mode 子工作流

**Files:**

- Create: `subworkflows/local/integrated_mode/main.nf`
- Create: `subworkflows/local/integrated_mode/meta.yml`
- Modify: `workflows/circdna.nf`

**目标:** 实现 Integrated Mode: gDNA + eccDNA 联合分析，输出 ECC_SCORE 分级 catalog

**步骤:**

- [ ] **Step 1: 创建 integrated_mode subworkflow**

```groovy
// subworkflows/local/integrated_mode/main.nf
include { DEPTH_RATIO } from '../../../modules/local/depth_ratio/main'
include { ECC_SCORE } from '../../../modules/local/ecc_score/main'

workflow INTEGRATED_MODE {
    take:
    reference_mosdepth   // gDNA 深度结果
    eccdna_mosdepth      // eccDNA 深度结果
    eccdna_merged_bed    // eccDNA 合并候选
    repeat_bed           // TE 重复序列 BED
    w1, w2, w3           // ECC_SCORE 权重

    main:
    ch_versions = channel.empty()

    // 样本配对（按样本ID匹配 gDNA 和 eccDNA）
    // 这里需要根据实际 samplesheet 设计配对逻辑

    // 计算 ECC_SCORE
    ECC_SCORE (
        eccdna_merged_bed,
        eccdna_mosdepth,
        reference_mosdepth,
        repeat_bed
    )
    ch_versions = ch_versions.mix(ECC_SCORE.out.versions)

    emit:
    scored_bed    = ECC_SCORE.out.scored_bed
    versions      = ch_versions
}
```

- [ ] **Step 2: 样本配对逻辑**

实现样本配对：根据样本命名规则或 samplesheet 中的配对信息，将 gDNA 和 eccDNA 样本配对

---

### Task 4.4: 完整集成与测试验证

**Files:**

- Modify: `workflows/circdna.nf`
- Modify: `conf/test_local.config`
- Create: `testdatasets/samplesheet/test_integrated.csv`

**步骤:**

- [ ] **Step 1: 完整集成所有模式到主工作流**

在 [workflows/circdna.nf](file:///Users/siyangming/nextflow_nf_core/circdna.nf/workflows/circdna.nf) 中完成三种模式的完整调用

- [ ] **Step 2: 创建测试 samplesheet**

创建包含 gDNA 和 eccDNA 样本的测试用 samplesheet

- [ ] **Step 3: Reference Mode 测试**

```bash
nextflow run main.nf -profile test_local,docker --mode reference
```

- [ ] **Step 4: eccDNA Mode 测试**

```bash
nextflow run main.nf -profile test_local,docker --mode eccdna
```

- [ ] **Step 5: Integrated Mode 测试**

```bash
nextflow run main.nf -profile test_local,docker --mode integrated
```

- [ ] **Step 6: 更新 CHANGELOG.md**

版本号从 v3.2.0 升级到 v4.0.0（因为这是重大架构变更）

---

## 潜在依赖与注意事项

### 工具可用性

- **ECCsplorer**: 需要确认是否有 bioconda 包或 Docker 镜像。如没有，需要自定义容器
- **BWA-MEM2**: nf-core 已有模块，索引格式与 BWA 不同，需要注意迁移
- **mosdepth**: nf-core 已有模块，直接可用

### 样本配对策略

- Integrated Mode 需要明确 gDNA 和 eccDNA 样本的配对关系
- 策略选项：
  1. 样本命名约定（如 `sampleA_gdna`, `sampleA_eccdna`）
  2. samplesheet 中新增 `pair_id` 列
  3. 单样本 + `--gdna_input` 和 `--eccdna_input` 两个参数

### 向后兼容

- 保留 `--circle_identifier` 参数的兼容性
- 无 `--mode` 参数时默认行为不变（eccdna 模式）
- 旧版 samplesheet（无新字段）应能正常运行

### 性能考虑

- mosdepth 比 samtools depth 快很多，但大基因组仍需注意内存
- ECC_SCORE 计算是单进程，候选多时可能较慢

### 风险处理

- **ECCsplorer 无现成容器**: 先实现 Python 版本的简化版 junction read 检测，后续再集成完整 ECCsplorer
- **BWA-MEM2 索引格式问题**: 保留 BWA 作为备选，通过参数切换
- **样本配对复杂**: 第一版先支持单样本配对，多样本配对后续迭代

---

## 风险与缓解策略

| 风险                          | 影响 | 概率 | 缓解策略                                            |
| ----------------------------- | ---- | ---- | --------------------------------------------------- |
| ECCsplorer 无可用 bioconda 包 | 高   | 中   | 先实现简化版 junction-read 检测脚本，后续集成完整版 |
| BWA-MEM2 索引不兼容           | 中   | 低   | 保留 BWA 模块作为 fallback，通过参数切换            |
| 样本配对逻辑复杂              | 中   | 高   | v1.0 仅支持单样本配对，多样本后续迭代               |
| mosdepth 内存不足             | 中   | 低   | 调整 `--fast-mode` 和窗口大小参数                   |
| 向后兼容性破坏                | 高   | 低   | 默认 mode='eccdna'，旧参数全部保留                  |

---

## 验证标准

### Sprint 1 完成标准

- [ ] Samplesheet 支持 datatype/platform/protocol 字段
- [ ] `--mode reference` 能运行完 fastp + bwamem2 + mosdepth
- [ ] 旧版 samplesheet 仍能正常运行（向后兼容）

### Sprint 2 完成标准

- [ ] ECCsplorer 模块可运行并输出 BED
- [ ] Circle-Map 输出标准化 BED
- [ ] `--mode eccdna` 能运行完所有 eccDNA 工具

### Sprint 3 完成标准

- [ ] merge_candidates.py 能正确合并多个工具的结果
- [ ] 输出包含来源工具数、junction reads 数的 BED 文件

### Sprint 4 完成标准

- [ ] ECC_SCORE v1.0 计算正确
- [ ] `--mode integrated` 能输出带评分和分级的候选 catalog
- [ ] 三种模式均通过测试验证
- [ ] CHANGELOG 更新到 v4.0.0
