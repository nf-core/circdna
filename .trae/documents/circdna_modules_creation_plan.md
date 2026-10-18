# circDNA 模块与子工作流创建计划

## 概述

在 `/Users/siyangming/nextflow_nf_core/circdna.nf` 项目中创建 5 个 Nextflow DSL2 文件，遵循现有代码风格。

## 创建文件列表

### 1. modules/local/eccsplorer/main.nf

- **进程名**: ECCSPLORER
- **输入**: `tuple val(meta), path(bam), path(bai), path(fasta)`
- **输出**:
  - `tuple val(meta), path("*_candidates.bed")` (emit: candidates_bed)
  - `tuple val(meta), path("*_junction_reads.txt")` (emit: junction_reads)
  - `tuple val("${task.process}"), val('eccsplorer'), eval(...), emit: versions, topic: versions`
- **标签**: `process_medium`
- **容器**: `python:3.9-slim`
- **Stub 块**: 生成假的 BED 文件和 junction reads 文件
- **脚本**: stub 版本（先输出假结果）

### 2. modules/local/candidate_merge/main.nf

- **进程名**: CANDIDATE_MERGE
- **输入**: `tuple val(meta), path(eccsplorer_bed), path(circle_map_bed)`
- **输出**:
  - `tuple val(meta), path("*_merged_candidates.bed")` (emit: merged_bed)
  - `tuple val("${task.process}"), val('candidate_merge'), eval(...), emit: versions, topic: versions`
- **脚本调用**: `python ${projectDir}/bin/merge_candidates.py --eccsplorer ${eccsplorer_bed} --circle_map ${circle_map_bed} --output ${prefix}_merged_candidates.bed`
- **标签**: `process_single`
- **容器**: `python:3.9-slim`
- **Stub 块**: 生成假的 merged BED 文件
- **注意**: 支持输入文件为空的情况

### 3. modules/local/ecc_score/main.nf

- **进程名**: ECC_SCORE
- **输入**: `tuple val(meta), path(candidates_bed), path(eccdna_depth_bed), path(gdna_depth_bed), path(repeat_bed)`
- **输出**:
  - `tuple val(meta), path("*_scored.bed")` (emit: scored_bed)
  - `tuple val("${task.process}"), val('ecc_score'), eval(...), emit: versions, topic: versions`
- **脚本调用**: `python ${projectDir}/bin/calculate_ecc_score.py --candidates ${candidates_bed} --eccdna-depth ${eccdna_depth_bed} --gdna-depth ${gdna_depth_bed} --repeat-bed ${repeat_bed} --output ${prefix}_scored.bed --w1 ${w1} --w2 ${w2} --w3 ${w3}`
- **权重配置**: 通过 `task.ext.w1`, `task.ext.w2`, `task.ext.w3` 配置，默认值 1.0, 1.0, 0.5
- **标签**: `process_single`
- **容器**: `python:3.9-slim`
- **Stub 块**: 生成假的 scored BED 文件

### 4. subworkflows/local/eccdna_mode/main.nf

- **工作流名**: ECCDNA_MODE
- **Include**: BAM_PREPROCESSING, MOSDEPTH, CIRCLE_MAP_PIPELINE, CANDIDATE_MERGE
- **输入**:
  - reads (channel: [ val(meta), [ reads ] ])
  - bwa_index (channel: [ "bwa_index", index_dir ])
  - fasta_meta (channel: [ val(meta), path(fasta) ])
  - run_circle_map (boolean)
- **流程**: BAM预处理 -> mosdepth深度分析 -> (可选)Circle-Map -> Candidate Merge
- **输出**:
  - bam_sorted
  - bam_sorted_bai
  - mosdepth_bed
  - mosdepth_summary
  - circle_map_bed
  - merged_bed
  - versions
- **注意**: 没有eccsplorer时，candidate_merge也要能工作（输入空文件）

### 5. subworkflows/local/integrated_mode/main.nf

- **工作流名**: INTEGRATED_MODE
- **Include**: ECC_SCORE
- **输入**:
  - reference_mosdepth_bed
  - eccdna_mosdepth_bed
  - eccdna_merged_bed
  - repeat_bed (可以为空)
  - w1, w2, w3
- **流程**: 计算 ECC_SCORE
- **输出**:
  - scored_bed
  - versions

## 代码风格遵循

- 参考 `modules/nf-core/mosdepth/main.nf` 的 module 写法
- 参考 `subworkflows/local/reference_mode/main.nf` 的 subworkflow 写法
- 使用 `def prefix = task.ext.prefix ?: "${meta.id}"` 模式
- versions 输出使用 tuple 格式（参考 mosdepth）
- 包含 `when: task.ext.when == null || task.ext.when` 块
- stub 块使用 `touch` 或简单 echo 生成假文件
- 相对路径引用正确
