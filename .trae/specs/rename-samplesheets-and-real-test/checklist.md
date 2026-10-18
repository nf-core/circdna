# Checklist

## 阶段一：Samplesheets 重命名

### 文件重命名

- [x] `samplesheet.csv` → `test_online.csv`
- [x] `samplesheet_bam_test.csv` → `test_local_bam.csv`
- [x] `samplesheet_local.csv` → `test_local_eccdna.csv`
- [x] `samplesheet_local_with_gdna.csv` → `test_local_gdna.csv`
- [x] `samplesheet_integrated.csv` → `test_local_integrated.csv`
- [x] 原 `samplesheet_*` 文件已不存在（git mv 或 mv 后删除原文件）

## 阶段二：更新引用文件

### conf/test_local.config

- [x] `input` 路径已更新为 `test_local_eccdna.csv`

### conf/test_local_gdna.config

- [x] 注释中文件名已更新为 `test_local_gdna.csv`
- [x] `input` 路径已更新为 `test_local_gdna.csv`

### scripts/test_incremental_cache.py

- [x] `ORIGINAL_SAMPLESHEET` 路径已更新
- [x] backup 文件名已更新

### AGENTS.md

- [x] 所有 `samplesheet_local.csv` 引用已更新为 `test_local_eccdna.csv`
- [x] 所有 `samplesheet_integrated.csv` 引用已更新为 `test_local_integrated.csv`

### testdatasets/README.md

- [x] `samplesheet_local.csv` 引用已更新为 `test_local_eccdna.csv`

## 阶段三：stub 模式快速验证

- [x] `nextflow run main.nf -profile test_local,docker -stub` 执行成功（52 tasks, 1m 42s）
- [x] `nextflow run main.nf -profile test_local_gdna,docker -stub` 执行成功（61 tasks, 1m 44s）

## 阶段四：真实模式测试

### test_local（无 control）

- [x] `nextflow run main.nf -profile test_local,docker` 执行成功（61 tasks, 1h26m58s）
- [x] 无 failed 任务
- [x] ECCSPLORER 任务全部完成（3/3）

### test_local_gdna（有 control）

- [x] `nextflow run main.nf -profile test_local_gdna,docker` 执行成功（70 tasks, 全部 COMPLETED）
- [x] 无 failed 任务
- [x] ECCSPLORER_WITH_CONTROL 任务全部完成（3/3）

## 阶段五：验证真实产出

### test_local 产出

- [x] `results/test_local/eccsplorer/*_candidates.bed` 非空（78/102/119 行）
- [x] `results/test_local/eccsplorer/*_blast.m6` 非空（85-445 bytes，BLAST 注释生效）
- [x] `results/test_local/eccsplorer/*_junction_reads.txt` 非空（144KB-409KB）
- [x] `results/test_local/eccsplorer/*_ecc_sequences.fasta` 非空（256KB-436KB）
- [x] `results/test_local/multiqc/multiqc_report.html` 存在（2.4MB）

### test_local_gdna 产出

- [x] `results/test_local_gdna/eccsplorer/*_candidates.bed` 非空（78/102/119 行）
- [x] `results/test_local_gdna/eccsplorer/*_blast.m6` 非空（55/78、68/102、66/119 个非空文件）
- [x] `results/test_local_gdna/eccsplorer/*_junction_reads.txt` 非空（144KB-408KB）
- [x] `results/test_local_gdna/eccsplorer/*_ecc_sequences.fasta` 非空（256KB-436KB）
- [x] `results/test_local_gdna/multiqc/multiqc_report.html` 存在（2.4MB）

## 阶段六：对比分析

- [x] 已统计 `results/test_local/eccsplorer/` 各样本 candidates.bed 行数（78/102/119）
- [x] 已统计 `results/test_local_gdna/eccsplorer/` 各样本 candidates.bed 行数（78/102/119）
- [x] 已对比两组数据，记录 gDNA control 对候选数量的影响
  - 两组候选数量完全相同
  - 原因：ECCsplorer comparison/clustering 阶段失败（RepeatExplorer 聚类数据库未配置）
  - 模块脚本 `|| true` 使流程继续，复制了未过滤的 mapping 结果
  - 结论：gDNA control 比对功能未实际生效，需后续配置 RepeatExplorer 数据库
