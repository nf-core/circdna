# Fix Circle-Map Realign Missing Index Files Spec

## Why

CIRCLEMAP_REALIGN 进程因缺少两个索引文件而失败：(1) qname-sorted BAM 的 `.bai` 索引未传入进程输入，(2) bgzip 压缩的参考基因组 FASTA（`.fa.gz`）缺少 `.gzi` 索引，导致 pysam 无法随机访问。

## What Changes

- 在 CIRCLE_MAP_PIPELINE 子工作流中将 `ch_qname_sorted_bai` 并入 CIRCLEMAP_REALIGN 的输入通道，并更新输入 tuple 结构
- 在 CIRCLEMAP_REALIGN 模块的 script 块中，调用 Circle-Map 前先运行 `samtools faidx` 以在任务工作目录中生成 `.gzi` 索引

## Impact

- Affected specs: 无
- Affected code:
  - `subworkflows/local/circle_map_pipeline/main.nf`（通道 join + 入参传递）
  - `modules/local/circlemap/realign/main.nf`（script 块 + input tuple）

## MODIFIED Requirements

### Requirement: CIRCLEMAP_REALIGN 输入包含完整的 BAM 索引

CIRCLEMAP_REALIGN SHALL 接收包含全部四组 BAM/BAI 的输入 tuple：`[meta, re_bam, re_bai, qname_bam, qname_bai, sbam, sbai]`（原为 `[meta, re_bam, re_bai, qname_bam, sbam, sbai]`）。

#### Scenario: Realign 正常执行

- **WHEN** CIRCLE_MAP_PIPELINE 将 qname BAI 并入 CIRCLEMAP_REALIGN 输入
- **THEN** Circle-Map Realign 命令可正常通过 pysam 随机访问 qname-sorted BAM

### Requirement: CIRCLEMAP_REALIGN 自行创建 FASTA 的 .gzi 索引

CIRCLEMAP_REALIGN script 块 SHALL 在调用 `circle_map.py Realign` 之前执行 `samtools faidx $fasta`，确保 bgzip 压缩的 FASTA 文件在任务工作目录中有 `.gzi` 索引，供 pysam.FastaFile 使用。

#### Scenario: bgzip FASTA 被正确索引

- **WHEN** 参考基因组为 `.fa.gz` 格式
- **THEN** samtools faidx 在任务目录创建 `.fai` + `.gzi`，pysam.FastaFile 可正常打开并随机访问
