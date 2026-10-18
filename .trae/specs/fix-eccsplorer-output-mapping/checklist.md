# Checklist

## 阶段一：circdna.nf ECCSPLORER 模块修复

### main.nf output 块
- [x] output 块新增 `tuple val(meta), path("*_lowconf_ecc_regions.bed"), emit: lowconf_ecc_regions`
- [x] output 块新增 `path "*_alignment_stats.txt", emit: alignment_stats`
- [x] output 块新增 `path "*_ecc_sequences.fasta", emit: ecc_sequences`
- [x] output 块新增 `path "*_eccpipe_results", emit: eccpipe_results`
- [x] 保留原有 `candidates_bed`、`junction_reads`、`versions` 输出声明

### main.nf script 块（BAM 模式）
- [x] `*_candidates.bed` 的源文件改为 `${prefix}_output/eccpipe_results/mapping_results/*_hiconf-ECC-REGIONS.bed`，使用 `cp` 而非 `mv`
- [x] `*_junction_reads.txt` 的源文件改为 `${prefix}_output/eccpipe_results/mapping_results/*.trns.txt`，使用 `cp`
- [x] 新增 `cp ${prefix}_output/eccpipe_results/mapping_results/*_lowconf-ECC-regions.bed ${prefix}_lowconf_ecc_regions.bed 2>/dev/null || touch ${prefix}_lowconf_ecc_regions.bed`
- [x] 新增 `cp ${prefix}_output/eccpipe_results/mapping_results/*_alignment-stats.txt ${prefix}_alignment_stats.txt 2>/dev/null || touch ${prefix}_alignment_stats.txt`
- [x] 新增 `cp ${prefix}_output/eccpipe_results/mapping_results/*_ECC-SEQUENCES.fasta ${prefix}_ecc_sequences.fasta 2>/dev/null || touch ${prefix}_ecc_sequences.fasta`
- [x] 新增 `cp -r ${prefix}_output/eccpipe_results ${prefix}_eccpipe_results`
- [x] 保留 `|| touch` 兜底逻辑（高置信度区域可能为空）

### main.nf script 块（FASTQ 模式）
- [x] 同 BAM 模式的所有修改

### main.nf stub 块
- [x] 新增 `mkdir -p ${prefix}_eccpipe_results`
- [x] 新增 `touch ${prefix}_lowconf_ecc_regions.bed`
- [x] 新增 `touch ${prefix}_alignment_stats.txt`
- [x] 新增 `touch ${prefix}_ecc_sequences.fasta`
- [x] 保留原有 `touch ${prefix}_candidates.bed` 和 `touch ${prefix}_junction_reads.txt`

### meta.yml
- [x] output 部分新增 `lowconf_ecc_regions` 描述（低置信度 eccDNA 区域）
- [x] output 部分新增 `alignment_stats` 描述（segemehl 比对统计）
- [x] output 部分新增 `ecc_sequences` 描述（提取的 eccDNA 序列）
- [x] output 部分新增 `eccpipe_results` 描述（完整 ECCsplorer 结果树）

### modules.config
- [x] ECCSPLORER 的 publishDir 配置无需特别修改（默认无 pattern 会发布所有输出声明的文件/目录），已验证新增输出会被正确发布

## 阶段二：bio.nf 同步

### bio.nf/modules/eccsplorer/main.nf
- [x] output 块与 circdna.nf 一致（新增 4 个输出声明）
- [x] script 块 BAM 模式与 circdna.nf 一致
- [x] script 块 FASTQ 模式与 circdna.nf 一致
- [x] stub 块与 circdna.nf 一致

### bio.nf/modules/eccsplorer/meta.yml
- [x] output 部分与 circdna.nf 一致（新增 4 个输出描述）

## 阶段三：版本与文档

### nextflow.config
- [x] `manifest.version` 从 `'4.2.0'` 改为 `'4.2.1'`

### CHANGELOG.md
- [x] 顶部新增 `## v4.2.1 - [2026-08-03]` 版本段
- [x] Enhancements & fixes 部分新增"修复 ECCSPLORER 模块输出文件映射错误"条目，说明根因与修复方案
- [x] 新增条目说明新增的 4 个输出通道

## 阶段四：验证

### stub 模式运行
- [x] `nextflow run main.nf -profile test_local,docker -stub --outdir /tmp/circdna_stub_eccsplorer` 执行成功（exit code 0，52 任务全部 succeeded）
- [x] ECCSPLORER 任务无 failed（3 of 3 ✔）
- [x] `/tmp/circdna_stub_eccsplorer/eccsplorer/` 目录存在 `*_candidates.bed`（0 字节占位）
- [x] `/tmp/circdna_stub_eccsplorer/eccsplorer/` 目录存在 `*_junction_reads.txt`（0 字节占位）
- [x] `/tmp/circdna_stub_eccsplorer/eccsplorer/` 目录存在 `*_lowconf_ecc_regions.bed`（0 字节占位）
- [x] `/tmp/circdna_stub_eccsplorer/eccsplorer/` 目录存在 `*_alignment_stats.txt`（0 字节占位）
- [x] `/tmp/circdna_stub_eccsplorer/eccsplorer/` 目录存在 `*_ecc_sequences.fasta`（0 字节占位）
- [x] `/tmp/circdna_stub_eccsplorer/eccsplorer/` 目录存在 `*_eccpipe_results/`（空目录占位）
