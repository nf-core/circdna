# Tasks

## 阶段一：修复 circdna.nf ECCSPLORER 模块

- [x] Task 1: 修复 `circdna.nf/modules/local/eccsplorer/main.nf` 的输出映射
  - [x] SubTask 1.1: 修改 output 块，新增 `lowconf_ecc_regions`、`alignment_stats`、`ecc_sequences`、`eccpipe_results` 输出声明
  - [x] SubTask 1.2: 修改 script 块（BAM 模式），将 `mv ${prefix}_output/*_candidates.bed` 改为 `cp ${prefix}_output/eccpipe_results/mapping_results/*_hiconf-ECC-REGIONS.bed ${prefix}_candidates.bed`，junction_reads 改为 `*.trns.txt`，新增 lowconf/stats/sequences 的 cp 命令，新增 `cp -r ${prefix}_output/eccpipe_results ${prefix}_eccpipe_results`
  - [x] SubTask 1.3: 修改 script 块（FASTQ 模式），同 SubTask 1.2
  - [x] SubTask 1.4: 修改 stub 块，新增 `mkdir -p ${prefix}_eccpipe_results`、`touch ${prefix}_lowconf_ecc_regions.bed`、`touch ${prefix}_alignment_stats.txt`、`touch ${prefix}_ecc_sequences.fasta`

- [x] Task 2: 更新 `circdna.nf/modules/local/eccsplorer/meta.yml`
  - [x] SubTask 2.1: 在 output 部分新增 `lowconf_ecc_regions`、`alignment_stats`、`ecc_sequences`、`eccpipe_results` 的描述

- [x] Task 3: 更新 `circdna.nf/conf/modules.config` 的 ECCSPLORER publishDir
  - [x] SubTask 3.1: 保持现有 publishDir 配置（已使用无 pattern 的默认发布，会自动发布所有输出声明的文件/目录），验证新增输出会被正确发布

## 阶段二：同步修复 bio.nf ECCSPLORER 模块

- [x] Task 4: 同步修复 `bio.nf/modules/eccsplorer/main.nf`
  - [x] SubTask 4.1: 应用与 Task 1 相同的修改（output 块、script 块 BAM 模式、script 块 FASTQ 模式、stub 块）

- [x] Task 5: 同步更新 `bio.nf/modules/eccsplorer/meta.yml`
  - [x] SubTask 5.1: 应用与 Task 2 相同的修改

## 阶段三：版本管理与文档

- [x] Task 6: 更新 `circdna.nf/nextflow.config` 版本号
  - [x] SubTask 6.1: 将 `version = '4.2.0'` 改为 `version = '4.2.1'`

- [x] Task 7: 更新 `circdna.nf/CHANGELOG.md`
  - [x] SubTask 7.1: 在文件顶部新增 `## v4.2.1 - [2026-08-03]` 版本段
  - [x] SubTask 7.2: 在 Enhancements & fixes 部分新增条目："修复 ECCSPLORER 模块输出文件映射错误"，详细说明根因与修复方案
  - [x] SubTask 7.3: 新增条目说明新增的输出通道（lowconf_ecc_regions、alignment_stats、ecc_sequences、eccpipe_results）

## 阶段四：验证

- [x] Task 8: stub 模式验证
  - [x] SubTask 8.1: 执行 `nextflow run main.nf -profile test_local,docker -stub --outdir /tmp/circdna_stub_eccsplorer` 验证 stub 模式通过（52 任务全部成功，exit code 0）
  - [x] SubTask 8.2: 验证 stub 产出文件存在（3 样本 × 6 文件/目录 = 18 项全部存在：`*_candidates.bed`、`*_junction_reads.txt`、`*_eccpipe_results/`、`*_lowconf_ecc_regions.bed`、`*_alignment_stats.txt`、`*_ecc_sequences.fasta`）

# Task Dependencies

- Task 1、Task 2、Task 3 可顺序执行（同一文件组的修改）
- Task 4、Task 5 依赖 Task 1、Task 2 完成（同步上游模块）
- Task 6、Task 7 依赖 Task 1-5 完成（所有代码修改就绪后更新版本与文档）
- Task 8 依赖 Task 1-7 全部完成

# 并行化说明

- Task 4+5（bio.nf 同步）可与 Task 3（modules.config 更新）并行修改（不同文件）
- Task 6（版本号）与 Task 7（CHANGELOG）建议串行（CHANGELOG 引用版本号）
