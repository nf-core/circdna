# circdna.nf 三代长读长 eccDNA 分析集成计划 (circdnalr) - 方案 B

## 1. 执行摘要

在 `circdna.nf` v3.1.0 中集成 PacBio / ONT 三代长读长 eccDNA 分析流程。**方案 B：内部构建预处理子流程**，直接引用 nf-core 独立模块（PBCCS, LIMA, CHOPPER, PYCHOPPER）按需组装，而非调用完整的 `isoseq.nf` 或 `nanoseq.nf`。核心分析采用 **CReSIL（mapping）+ FLED（mapping）+ Flye（assembly）** 三引擎架构，并增加结果收敛与过滤子流程。

## 2. 方案决策理由

### 为什么不调用完整的外部流程？

**isoseq.nf（PacBio）**：

- 核心设计目标是**全长转录本分析（FLNC）**，包含 `ISOSEQ_REFINE`（丢弃无 polyA 的 reads）和 `GSTAMA_POLYACLEANUP`（清理 polyA 尾巴）。
- eccDNA 是基因组环状 DNA，**没有 polyA 尾巴**，直接套用会导致大量真实 eccDNA reads 被误删。
- 最终输出是 `FASTA`，而 FLED/CReSIL 期望 `FASTQ` 输入。

**nanoseq.nf（ONT）**：

- DNA protocol 分支仅包含 `QCAT`（demultiplexing）→ `NANOLYSE`（去污染）→ `QCFASTQ`（质控）。
- **关键缺失**：`CHOPPER`（质量过滤）和 `PYCHOPPER`（接头切除）仅在 cDNA/directRNA 分支中使用，DNA 分支完全没有 read cleaning/adapter removal。
- eccDNA ONT RCA 扩增数据**必须切除接头和引物**，nanoseq DNA 分支无法满足。

### 方案 B 的优势

| 维度       | 方案 A（外部流程）       | 方案 B（内部组装）              |
| ---------- | ------------------------ | ------------------------------- |
| 技术可行性 | ❌ polyA 过滤误伤 eccDNA | ✅ 无转录本特异性过滤           |
| 操作复杂度 | ❌ 需管理两套样本表/流程 | ✅ 单一样本表、单一命令         |
| 数据兼容性 | ❌ FASTA vs FASTQ 不匹配 | ✅ 输出直接是三引擎需要的 FASTQ |
| 维护成本   | ⚠️ 需维护格式转换脚本    | ✅ 模块来自 nf-core 标准库      |
| 灵活性     | ❌ 无法跳过不需要的步骤  | ✅ entrypoint 参数灵活控制      |

## 3. 当前状态分析

- **Branch**: `circdnalr` 已创建（当前分支）。
- **版本**: `manifest.version` 已改为 `3.1.0`，但 `params` 块未加入新参数。
- **CHANGELOG**: 已更新 v3.1.0 条目，但代码层面尚未实现。
- **主流程** (`workflows/circdna.nf`):
  - 仅有 `--input_format` (FASTQ/BAM) 和 `--circle_identifier`。
  - 无 `protocol`/`entrypoint` 分支，无长读预处理逻辑。
- **样本表** (`assets/schema_input.json` + `bin/check_samplesheet.py`):
  - 仅支持 `sample,fastq_1,fastq_2` (FASTQ) 或 `sample,bam` (BAM)。
- **nf-core 模块状态**:
  - ✅ 已安装：`flye` (`modules/modules/nf-core/flye/main.nf`)
  - ✅ 可用但未安装到本流程：`pbccs`, `lima`, `chopper`, `pychopper`, `minimap2/align`, `samtools/sort`, `samtools/index`, `bedtools/intersect`
- **bio.nf 模块**:
  - ✅ 已存在：`CRESIL_IDENTIFY` (`modules/cresil/identify/main.nf`)
  - ✅ 已存在：`FLYE` (`modules/flye/main.nf`)
  - ❌ 缺失：`FLED`（需在 bio.nf 新建后复制）
- **容器**: Docker registry 指向 `quay.io`，已配置 wave/apptainer 等 profiles。FLED 容器需自构建。

## 4. 待修改文件清单与详细步骤

### 4.1 参数与 Schema

**文件**: `nextflow.config`

- **What**: 在 `params` 块新增：
  ```groovy
  params {
      protocol              = 'short_read'    // short_read | pacbio | ont
      entrypoint            = 'cleaned_fastq' // cleaned_fastq | raw_fastq | subreads | hifi_bam (pacbio)
      primers               = null             // Primers FASTA for Lima/Pychopper
      long_read_identifier  = 'cresil,fled,flye'
      min_read_support      = 2
      blacklist_bed         = null
      repeats_bed           = null
      save_long_read_intermediate = false
      skip_long_read_qc     = false
  }
  ```
- **Why**: `protocol` 区分测序平台（参考 nanoseq.nf），`entrypoint` 控制预处理深度（参考 isoseq.nf）。
- **How**: 追加到现有 `params` 块，默认 `short_read` 保证二代兼容。

**文件**: `nextflow_schema.json`

- **What**: 为新增参数添加 JSON Schema 定义：
  - `protocol`: enum `["short_read", "pacbio", "ont"]`
  - `entrypoint`: enum `["cleaned_fastq", "raw_fastq", "subreads", "hifi_bam"]`
  - `long_read_identifier`: string, pattern `^[a-z,]+$`
  - `min_read_support`: integer, minimum 1
  - `blacklist_bed`, `repeats_bed`: string, format file-path
- **Why**: nf-schema 插件验证用户输入合法性。

**文件**: `assets/schema_input.json`

- **What**: 扩展样本表 schema，增加可选列：
  - `input_bam` (type string, optional): 原始 BAM 入口（subreads/hifi_bam）。
  - `entrypoint` (type string, optional): 覆盖全局 entrypoint。
- **Why**: 长读数据允许以原始 BAM 或清洗后 FASTQ 为入口。

### 4.2 样本表解析扩展

**文件**: `bin/check_samplesheet.py`

- **What**: 扩展脚本，支持长读样本表：
  - `protocol == 'short_read'`: 保持现有逻辑（FASTQ/BAM）。
  - `protocol == 'pacbio' 或 'ont'`:
    - 接受 `sample,fastq_1`（清洗后 FASTQ）或 `sample,input_bam`（原始 BAM）。
    - 允许 `fastq_2` 为空（长读天然单端）。
- **Why**: 当前脚本不支持长读列。
- **How**: 新增第4个命令行参数 `PROTOCOL`，根据 protocol 切换校验逻辑。

**文件**: `modules/local/samplesheet_check/main.nf`

- **What**: 修改 script 块，传入 `params.protocol`：
  ```groovy
  check_samplesheet.py \
      $samplesheet \
      samplesheet.valid.csv \
      $params.input_format \
      $params.protocol
  ```

**文件**: `subworkflows/local/input_check/main.nf`

- **What**: 增加长读解析逻辑：
  - `protocol == 'short_read'`: 保持现有逻辑。
  - `protocol == 'pacbio' 或 'ont'`: 解析 `fastq_1`（可选）和 `input_bam`（可选），构建 channel `[meta, fastq, input_bam, entrypoint]`。

### 4.3 安装/复制模块

**通过 `nf-core modules install` 安装的模块**（安装后自动更新 `modules.json`）：

- `pbccs` — 从 isoseq.nf 底层复用
- `lima` — 从 isoseq.nf 底层复用
- `chopper` — 从 nanoseq.nf 底层复用（cDNA/directRNA 分支使用）
- `pychopper` — 从 nanoseq.nf 底层复用（cDNA/directRNA 分支使用）
- `minimap2/align` — 长读比对（需保留 SA 标签）
- `samtools/sort` — BAM 排序
- `samtools/index` — BAM 索引
- `bedtools/intersect` — 假阳性过滤

**从 bio.nf 复制的模块**：

- `modules/local/cresil_identify/main.nf` — 复制自 `bio.nf/modules/cresil/identify/main.nf`
- `modules/local/fled/main.nf` — 在 bio.nf 新建后复制（见 4.3.3）

**复用 nf-core 已安装模块**：

- `FLYE` — 已安装在 `modules/modules/nf-core/flye/main.nf`，直接 include 使用

#### 4.3.1 CReSIL 模块

**文件**: `modules/local/cresil_identify/main.nf`

- **What**: 复制 `bio.nf/modules/cresil/identify/main.nf`，调整输入输出格式适配 circdna.nf。
- **容器**: `quay.io/bioinfortools/cresil:1.2.1`（与 bio.nf 一致）。

#### 4.3.2 Flye 模块

**文件**: `modules/modules/nf-core/flye/main.nf`

- **What**: 直接在主流程中 include nf-core 的 `FLYE` 模块，无需复制到 local。
- **Why**: nf-core 版本更完善（wave 容器支持），减少维护成本。

#### 4.3.3 FLED 模块（需新建）

**文件**: `bio.nf/modules/fled/main.nf`（先在 bio.nf 新建）

```groovy
process FLED {
    tag "$meta.id"
    label 'process_high'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fled:1.0.0--pyhdfd78af_0' :
        'quay.io/bioinfortools/fled:1.0.0' }"
    input:
    tuple val(meta), path(fastq)
    tuple val(meta2), path(fasta)
    output:
    tuple val(meta), path("*.full_eccDNA.txt"), emit: full_eccdna
    tuple val(meta), path("*.break_eccDNA.txt"), emit: break_eccdna
    tuple val(meta), path("*.consensus.fasta"), emit: consensus
    tuple val("${task.process}"), val('fled'), eval("fled --version"), topic: versions, emit: versions_fled
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    fled -i ${fastq} -r ${fasta} -o ${prefix} -t ${task.cpus} ${args}
    """
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.full_eccDNA.txt ${prefix}.break_eccDNA.txt ${prefix}.consensus.fasta
    """
}
```

**文件**: `bio.nf/modules/fled/environment.yml`

```yaml
name: fled
channels:
  - conda-forge
  - bioconda
dependencies:
  - python>=3.8
  - pyspoa
  - numpy
  - pandas
  - pysam
```

**复制到 circdna.nf**: `modules/local/fled/main.nf` 和 `modules/local/fled/environment.yml`

### 4.4 新建长读子流程

在 `subworkflows/local/` 下新建以下子流程：

#### A. `long_read_preprocessing/main.nf`

**目标**: 根据 protocol + entrypoint 组装预处理链。

```groovy
workflow LONG_READ_PREPROCESSING {
    take:
    ch_samplesheet   // [meta, fastq, input_bam, entrypoint]
    ch_primers       // [meta, primers_fasta]
    ch_fasta_meta    // [meta, fasta]

    main:
    ch_versions = Channel.empty()
    ch_cleaned_fastq = Channel.empty()

    // PacBio 分支
    if (params.protocol == 'pacbio') {
        if (params.entrypoint == 'subreads') {
            PBCCS(ch_samplesheet.map { meta, fastq, bam, ep -> [meta, bam] })
            ch_versions = ch_versions.mix(PBCCS.out.versions)
            LIMA(PBCCS.out.bam, ch_primers)
            ch_versions = ch_versions.mix(LIMA.out.versions)
            ch_cleaned_fastq = LIMA.out.bam
        } else if (params.entrypoint == 'hifi_bam') {
            LIMA(ch_samplesheet.map { meta, fastq, bam, ep -> [meta, bam] }, ch_primers)
            ch_versions = ch_versions.mix(LIMA.out.versions)
            ch_cleaned_fastq = LIMA.out.bam
        } else {
            ch_cleaned_fastq = ch_samplesheet.map { meta, fastq, bam, ep -> [meta, fastq] }
        }
    }

    // ONT 分支
    if (params.protocol == 'ont') {
        if (params.entrypoint == 'raw_fastq') {
            CHOPPER(ch_samplesheet.map { meta, fastq, bam, ep -> [meta, fastq] })
            ch_versions = ch_versions.mix(CHOPPER.out.versions)
            PYCHOPPER(CHOPPER.out.fastq)
            ch_versions = ch_versions.mix(PYCHOPPER.out.versions)
            ch_cleaned_fastq = PYCHOPPER.out.fastq
        } else {
            ch_cleaned_fastq = ch_samplesheet.map { meta, fastq, bam, ep -> [meta, fastq] }
        }
    }

    emit:
    cleaned_fastq = ch_cleaned_fastq
    versions = ch_versions
}
```

**Why**: 只保留 eccDNA 需要的步骤：PBCCS（生成 HiFi）→ LIMA（切除接头）或 CHOPPER → PYCHOPPER。**跳过 isoseq.nf 的 polyA 过滤**。

#### B. `long_read_mapping/main.nf`

**目标**: 长读比对，输出带 SA 标签的 BAM。

```groovy
workflow LONG_READ_MAPPING {
    take:
    ch_cleaned_fastq  // [meta, fastq]
    ch_fasta_meta     // [meta, fasta]

    main:
    def map_preset = params.protocol == 'pacbio' ? '-x map-pb' : '-x map-ont'
    MINIMAP2_ALIGN(
        ch_cleaned_fastq,
        ch_fasta_meta,
        Channel.value(true),   // bam_format
        Channel.value("bai"),  // index_format
        Channel.value(false),  // cigar_paf
        Channel.value(false),  // cigar_bam
        Channel.value(map_preset)
    )
    SAMTOOLS_SORT_BAM(MINIMAP2_ALIGN.out.bam, Channel.value([]), Channel.value('bai'))
    SAMTOOLS_INDEX_BAM(SAMTOOLS_SORT_BAM.out.bam)

    emit:
    bam_sorted = SAMTOOLS_SORT_BAM.out.bam
    bam_bai    = SAMTOOLS_INDEX_BAM.out.bai
}
```

**Why**: CReSIL 需要带 SA 标签的比对 BAM。

#### C. `long_read_cresil/main.nf`

**目标**: 封装 CReSIL 分析。

```groovy
workflow LONG_READ_CRESIL {
    take:
    ch_bam_sorted  // [meta, bam]
    ch_bam_bai     // [meta, bai]
    ch_fasta_meta  // [meta, fasta]
    ch_fasta_fai   // [meta, fai]

    main:
    CRESIL_IDENTIFY(ch_fasta_meta, ch_fasta_fai, ch_bam_sorted, ch_bam_bai)

    emit:
    eccdna_results = CRESIL_IDENTIFY.out.identify
    versions       = CRESIL_IDENTIFY.out.versions_cresil
}
```

#### D. `long_read_fled/main.nf`

**目标**: 封装 FLED 分析。

```groovy
workflow LONG_READ_FLED {
    take:
    ch_cleaned_fastq  // [meta, fastq]
    ch_fasta_meta     // [meta, fasta]

    main:
    FLED(ch_cleaned_fastq, ch_fasta_meta)

    emit:
    full_eccdna  = FLED.out.full_eccdna
    break_eccdna = FLED.out.break_eccdna
    consensus    = FLED.out.consensus
    versions     = FLED.out.versions_fled
}
```

**Why**: FLED 直接从 FASTQ 进行单分子解折叠，不依赖 BAM。

#### E. `long_read_flye/main.nf`

**目标**: 封装 Flye 从头组装 + 环状 Contig 提取 + 回比对。

```groovy
workflow LONG_READ_FLYE {
    take:
    ch_cleaned_fastq  // [meta, fastq]
    ch_fasta_meta     // [meta, fasta]

    main:
    def flye_mode = params.protocol == 'pacbio' ? '--pacbio-hifi' : '--nano-hq'
    FLYE(ch_cleaned_fastq, Channel.value(flye_mode))

    EXTRACT_CIRCULAR_CONTIG(FLYE.out.fasta, FLYE.out.txt)

    MINIMAP2_ALIGN(
        EXTRACT_CIRCULAR_CONTIG.out.circular_fasta,
        ch_fasta_meta,
        Channel.value(true),
        Channel.value("bai"),
        Channel.value(false),
        Channel.value(false),
        Channel.value('-x map-pb')
    )
    SAMTOOLS_SORT_BAM(MINIMAP2_ALIGN.out.bam, Channel.value([]), Channel.value('bai'))
    SAMTOOLS_INDEX_BAM(SAMTOOLS_SORT_BAM.out.bam)

    emit:
    circular_fasta = EXTRACT_CIRCULAR_CONTIG.out.circular_fasta
    mapped_bam     = SAMTOOLS_SORT_BAM.out.bam
    mapped_bai     = SAMTOOLS_INDEX_BAM.out.bai
}
```

**注意**: `EXTRACT_CIRCULAR_CONTIG` 是新的本地模块，解析 `assembly_info.txt` 提取环状 contig。

#### F. `long_read_filtering/main.nf`

**目标**: 统一三引擎结果格式并过滤假阳性。

```groovy
workflow LONG_READ_FILTERING {
    take:
    ch_cresil_results  // [meta, eccDNA_final.txt]
    ch_fled_break      // [meta, break_eccDNA.txt]
    ch_flye_mapped_bam // [meta, bam]
    ch_blacklist_bed   // path(bed)
    ch_repeats_bed     // path(bed)

    main:
    CONVERT_CRESIL_TO_BED(ch_cresil_results)
    CONVERT_FLED_TO_BED(ch_fled_break)
    BAM_TO_BED(ch_flye_mapped_bam)

    ch_all_bed = CONVERT_CRESIL_TO_BED.out.bed
        .mix(CONVERT_FLED_TO_BED.out.bed)
        .mix(BAM_TO_BED.out.bed)
        .groupTuple()

    if (params.blacklist_bed) {
        BEDTOOLS_FILTER_BLACKLIST(ch_all_bed, ch_blacklist_bed)
        ch_filtered = BEDTOOLS_FILTER_BLACKLIST.out.bed
    } else {
        ch_filtered = ch_all_bed
    }

    if (params.repeats_bed) {
        BEDTOOLS_FILTER_REPEATS(ch_filtered, ch_repeats_bed)
        ch_filtered = BEDTOOLS_FILTER_REPEATS.out.bed
    }

    FILTER_MIN_READ_SUPPORT(ch_filtered, params.min_read_support)

    emit:
    final_bed = FILTER_MIN_READ_SUPPORT.out.bed
}
```

### 4.5 集成到主流程

**文件**: `workflows/circdna.nf`

- **参数校验**（在现有校验之后）：
  ```groovy
  if (params.protocol != 'short_read' && params.protocol != 'pacbio' && params.protocol != 'ont') {
      exit 1, "Invalid protocol: ${params.protocol}"
  }
  if (params.protocol != 'short_read') {
      def lri = params.long_read_identifier.split(',')
      if (!lri.any { it in ['cresil', 'fled', 'flye'] }) {
          exit 1, "Invalid long_read_identifier"
      }
  }
  ```
- **Include 新增子流程**（文件顶部）：
  ```groovy
  include { LONG_READ_PREPROCESSING } from '../subworkflows/local/long_read_preprocessing/main'
  include { LONG_READ_MAPPING       } from '../subworkflows/local/long_read_mapping/main'
  include { LONG_READ_CRESIL        } from '../subworkflows/local/long_read_cresil/main'
  include { LONG_READ_FLED          } from '../subworkflows/local/long_read_fled/main'
  include { LONG_READ_FLYE          } from '../subworkflows/local/long_read_flye/main'
  include { LONG_READ_FILTERING     } from '../subworkflows/local/long_read_filtering/main'
  ```
- **主流程分支**（在现有二代逻辑之后，MultiQC 之前）：

  ```groovy
  if (params.protocol == 'pacbio' || params.protocol == 'ont') {
      LONG_READ_PREPROCESSING(INPUT_CHECK.out.reads,
          params.primers ? channel.fromPath(params.primers).map { [[id: it.baseName], it] } : channel.empty(),
          ch_fasta_meta)
      ch_cleaned_fastq = LONG_READ_PREPROCESSING.out.cleaned_fastq

      LONG_READ_MAPPING(ch_cleaned_fastq, ch_fasta_meta)
      ch_lr_bam = LONG_READ_MAPPING.out.bam_sorted

      def lri = params.long_read_identifier.split(',')
      if ('cresil' in lri) {
          LONG_READ_CRESIL(ch_lr_bam, LONG_READ_MAPPING.out.bam_bai, ch_fasta_meta, ch_fasta_fai)
          ch_cresil_results = LONG_READ_CRESIL.out.eccdna_results
      }
      if ('fled' in lri) {
          LONG_READ_FLED(ch_cleaned_fastq, ch_fasta_meta)
          ch_fled_break = LONG_READ_FLED.out.break_eccdna
      }
      if ('flye' in lri) {
          LONG_READ_FLYE(ch_cleaned_fastq, ch_fasta_meta)
          ch_flye_bam = LONG_READ_FLYE.out.mapped_bam
      }

      if (lri.any()) {
          LONG_READ_FILTERING(
              lri.contains('cresil') ? ch_cresil_results : channel.empty(),
              lri.contains('fled') ? ch_fled_break : channel.empty(),
              lri.contains('flye') ? ch_flye_bam : channel.empty(),
              params.blacklist_bed ? channel.fromPath(params.blacklist_bed) : channel.empty(),
              params.repeats_bed ? channel.fromPath(params.repeats_bed) : channel.empty()
          )
      }
  }
  ```

### 4.6 配置更新

**文件**: `conf/modules.config`

- 添加 publishDir 配置：
  - `PBCCS`: `reports/pbccs/`, enabled: `params.save_long_read_intermediate`
  - `LIMA`: `reports/lima/`, enabled: `params.save_long_read_intermediate`
  - `CHOPPER`: `reports/chopper/`, enabled: `params.save_long_read_intermediate`
  - `PYCHOPPER`: `reports/pychopper/`, enabled: `params.save_long_read_intermediate`
  - `MINIMAP2_ALIGN`: ext.args 添加 `-x map-pb` / `-x map-ont` 和 `--secondary=yes`
  - `CRESIL_IDENTIFY`: `cresil/`, enabled: true
  - `FLED`: `fled/`, enabled: true
  - `FLYE`: `flye/`, enabled: true

**文件**: `conf/base.config`

- 添加资源分配：
  - `FLYE`: cpus 16, memory 64.GB, time 48.h
  - `MINIMAP2_ALIGN`: cpus 8, memory 32.GB, time 24.h
  - `PBCCS`: cpus 16, memory 64.GB, time 24.h
  - `FLED`: cpus 8, memory 32.GB, time 24.h

### 4.7 新增本地脚本

**文件**: `bin/convert_cresil_to_bed.py` — CReSIL txt → BED
**文件**: `bin/convert_fled_to_bed.py` — FLED txt → BED
**文件**: `bin/extract_circular_contig.py` — 提取 Flye 环状 contig
**文件**: `bin/filter_min_read_support.py` — read support 过滤

### 4.8 容器与 Conda 环境

**FLED 容器构建**（无现成 biocontainer）：

1. 在 `bio.nf` 新建 branch `fled`，构建 `modules/fled/environment.yml`。
2. 构建 Docker: `docker build -t quay.io/bioinfortools/fled:1.0.0 .`
3. 推送: `docker push quay.io/bioinfortools/fled:1.0.0`
4. conda 包上传至 `anaconda.org/yangmingsi/fled`。

**其他工具**: 使用 nf-core / biocontainers 镜像（PBCCS, LIMA, CHOPPER, PYCHOPPER, flye, minimap2, samtools, bedtools）。

### 4.9 文档与变更记录

**文件**: `CHANGELOG.md` — 补充 v3.1.0 条目。
**文件**: `CHANGES&FIX/20260720.md` — 详细记录变更。
**文件**: `README.md` / `docs/usage.md` — 更新使用文档。

## 5. 假设与决策

| 决策项     | 选择                                                    | 理由                                                       |
| ---------- | ------------------------------------------------------- | ---------------------------------------------------------- |
| 预处理方案 | 方案 B：内部组装                                        | isoseq polyA 过滤误伤 eccDNA；nanoseq DNA 分支缺少接头切除 |
| protocol   | `short_read` / `pacbio` / `ont`                         | 与 nanoseq.nf 模式一致，默认 short_read 保证兼容           |
| entrypoint | `cleaned_fastq` / `raw_fastq` / `subreads` / `hifi_bam` | 参考 isoseq.nf，控制预处理深度                             |
| 模块来源   | nf-core 标准模块 > bio.nf > 新建 local                  | 减少维护成本，符合用户偏好                                 |
| FLED 容器  | 自构建上传                                              | 无现成 biocontainer                                        |
| 结果过滤   | 统一 BED 后过滤                                         | 三引擎输出格式差异大，BED 是最佳通用格式                   |

## 6. 验证步骤

1. **语法检查**: `nextflow run . -stub` 确认模块加载正确。
2. **Schema 验证**: `nextflow run . --help` 和 `--validate_params`。
3. **样本表校验**: 短读表和长读表分别验证。
4. **nf-core lint**: `nf-core lint` 检查规范。
5. **Stub 测试**:
   - `--protocol pacbio --entrypoint cleaned_fastq --long_read_identifier cresil,fled,flye`
   - `--protocol ont --entrypoint raw_fastq --long_read_identifier flye`
   - `--protocol short_read --circle_identifier circle_map_realign`（回归测试）
6. **容器测试**: 确认所有容器可获取。
7. **实际数据测试**: 用小样本运行完整流程。

## 7. 修改文件清单

### 新增文件

- `subworkflows/local/long_read_preprocessing/main.nf`
- `subworkflows/local/long_read_mapping/main.nf`
- `subworkflows/local/long_read_cresil/main.nf`
- `subworkflows/local/long_read_fled/main.nf`
- `subworkflows/local/long_read_flye/main.nf`
- `subworkflows/local/long_read_filtering/main.nf`
- `modules/local/cresil_identify/main.nf`
- `modules/local/fled/main.nf`
- `modules/local/fled/environment.yml`
- `modules/local/extract_circular_contig/main.nf`
- `modules/local/convert_cresil_to_bed/main.nf`
- `modules/local/convert_fled_to_bed/main.nf`
- `modules/local/bam_to_bed/main.nf`
- `modules/local/filter_min_read_support/main.nf`
- `bin/convert_cresil_to_bed.py`
- `bin/convert_fled_to_bed.py`
- `bin/extract_circular_contig.py`
- `bin/filter_min_read_support.py`

### 修改文件

- `nextflow.config`
- `nextflow_schema.json`
- `assets/schema_input.json`
- `bin/check_samplesheet.py`
- `modules/local/samplesheet_check/main.nf`
- `subworkflows/local/input_check/main.nf`
- `workflows/circdna.nf`
- `conf/modules.config`
- `conf/base.config`
- `modules.json`
- `CHANGELOG.md`
- `CHANGES&FIX/20260720.md`
- `README.md` / `docs/usage.md`

### 外部依赖（bio.nf）

- `bio.nf` 新建 branch `fled`，构建 `modules/fled/` 模块和容器。
