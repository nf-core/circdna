# 修复参考基因组通道类型错误 Spec

## Why

流程第一次运行只处理一个样本就显示完成，使用 `-resume` 后才能逐个处理后续样本。根本原因是参考基因组的 `fasta`/`fai` 通道为 queue channel（队列通道），被首个样本消费后即空，导致后续样本的 SAMTOOLS_SORT 等进程"无输入可处理"而被跳过。Queue channel 只能被消费一次，而所有样本共享的单一参考基因组文件必须使用 value channel（值通道），可被无限次消费。

## What Changes

- `subworkflows/local/bam_preprocessing/main.nf`：将 `ch_fasta` 与 `ch_fasta_fai` 通过 `.first()` 转为 value channel；新增 `fasta_fai` 与 `fai` emit，供下游子流程复用
- `subworkflows/local/circle_finder_pipeline/main.nf`：新增 `fasta_fai` 入参（value channel），替换原先对 SAMTOOLS_SORT 传入的 `channel.empty()`
- `subworkflows/local/circle_map_pipeline/main.nf`：将入参 `fasta` 改为 `fasta_fai`（value channel）；SAMTOOLS_SORT 使用 `ch_fasta_fai`；CIRCLEMAP_REALIGN 的 fasta 路径从 `fasta_fai.map { meta, fasta, fai -> fasta }` 提取
- `workflows/circdna.nf`：调用 CIRCLE_FINDER_PIPELINE、CIRCLE_MAP_PIPELINE、AMPLICONARCHITECT_PIPELINE 时传入 `BAM_PREPROCESSING.out.fasta_fai`；移除不再需要的 `SAMTOOLS_FAIDX` import；移除对不存在的 `SAMTOOLS_FAIDX.out.fai` 的引用

## Impact

- Affected specs: 无（独立 bug 修复）
- Affected code:
  - [subworkflows/local/bam_preprocessing/main.nf](file:///Users/siyangming/nextflow_nf_core/circdna.nf/subworkflows/local/bam_preprocessing/main.nf)
  - [subworkflows/local/circle_finder_pipeline/main.nf](file:///Users/siyangming/nextflow_nf_core/circdna.nf/subworkflows/local/circle_finder_pipeline/main.nf)
  - [subworkflows/local/circle_map_pipeline/main.nf](file:///Users/siyangming/nextflow_nf_core/circdna.nf/subworkflows/local/circle_map_pipeline/main.nf)
  - [workflows/circdna.nf](file:///Users/siyangming/nextflow_nf_core/circdna.nf/workflows/circdna.nf)

## ADDED Requirements

### Requirement: 参考基因组文件作为 value channel 传递

当参考基因组文件（单一文件，所有样本共享）需要被多个样本的 process 使用时，系统 SHALL 通过对单元素 queue channel 调用 `.first()` 将其转为 value channel，确保每个需要它的 process 都能读取。

#### Scenario: 多样本一次性全部处理

- **WHEN** 用户使用 `conf/test_local.config` 运行流程（3 个样本 circdna_1/2/3，不启用 unicycler）
- **THEN** 一次性生成 `results_testdata` 中所有样本的结果，无需 `-resume`

#### Scenario: value channel 无限消费

- **WHEN** BAM_PREPROCESSING 内 SAMTOOLS_FAIDX 产出 fasta/fai 并转为 value channel
- **THEN** BAM_STATS_SAMTOOLS、BAM_MARKDUPLICATES_PICARD、SAMTOOLS_VIEW_FILTER、SAMTOOLS_SORT_FILTERED 均能读取同一份 fasta/fai，且下游 CIRCLE_FINDER_PIPELINE、CIRCLE_MAP_PIPELINE、AMPLICONARCHITECT_PIPELINE 通过 `BAM_PREPROCESSING.out.fasta_fai` 复用

## MODIFIED Requirements

### Requirement: BAM_PREPROCESSING 子流程

BAM_PREPROCESSING SHALL 将 SAMTOOLS_FAIDX 产出的 `fa` 与 `fai` 通过 `.first()` 转为 value channel，并通过新增的 `fasta_fai`（`[meta, fasta, fai]`）与 `fai` emit 暴露给下游子流程。BAM_STATS_SAMTOOLS SHALL 接收 value channel 形式的 `ch_fasta` 而非 queue channel 形式的 `ch_fasta_fai`。

### Requirement: 主流程 CIRCDNA 的 ch_fasta_meta 与 ch_bwa_index 通道（追加）

主流程 CIRCDNA SHALL 将 `ch_fasta_meta` 与 `ch_bwa_index` 通过 `.first()` 转为 value channel，确保它们被 BWA_MEM 等多个样本共享的 process 广播消费。

**问题根因**：BWA_MEM 的输入签名为三个 `tuple val(meta), path(...)` 输入（reads/index/fasta）。当 `index`/`fasta` 为 queue channel（`.collect()` 产生）时，首个样本消费后 channel 即空，后续样本因无 index/fasta 被静默跳过。46 个样本的 Arabidopsis 运行只有 ERR9968690 产出完整结果即此问题。

**修复**：

- `ch_fasta_meta = ch_fasta.map{ fasta -> [[id: fasta.baseName], fasta] }.first()`（原为 `.collect()`）
- `ch_bwa_index = ch_bwa_index.map{ index -> ["bwa_index", index] }.first()`（原为 `.collect()`）
- `ch_fasta_fai = fasta_meta.join(SAMTOOLS_FAIDX.out.fai).map{...}.first()`（原为 `.collect()`，spec 已要求但代码未生效）

### Requirement: CIRCLE_FINDER_PIPELINE 子流程

CIRCLE_FINDER_PIPELINE SHALL 新增 `fasta_fai` 入参（value channel），并在调用 SAMTOOLS_SORT_QNAME_CF 时传入 `ch_fasta_fai`，替代原先的 `channel.empty()`。

### Requirement: CIRCLE_MAP_PIPELINE 子流程

CIRCLE_MAP_PIPELINE SHALL 将入参 `fasta` 更名为 `fasta_fai`（value channel `[meta, fasta, fai]`）；SAMTOOLS_SORT_QNAME_CM 与 SAMTOOLS_SORT_RE SHALL 使用 `ch_fasta_fai`；CIRCLEMAP_REALIGN 的 fasta 路径 SHALL 通过 `fasta_fai.map { meta, fasta, fai -> fasta }` 提取。

### Requirement: 主流程 CIRCDNA 调用

主流程 SHALL 移除 `SAMTOOLS_FAIDX` 的 import（其调用已封装在 BAM_PREPROCESSING 内）；SHALL 在调用 AMPLICONARCHITECT_PIPELINE、CIRCLE_FINDER_PIPELINE、CIRCLE_MAP_PIPELINE 时传入 `BAM_PREPROCESSING.out.fasta_fai`；SHALL 移除对不存在的 `SAMTOOLS_FAIDX.out.fai` 的引用（原先 AMPLICONARCHITECT_PIPELINE 处的 `ch_fasta_meta.join(SAMTOOLS_FAIDX.out.fai)`）。
