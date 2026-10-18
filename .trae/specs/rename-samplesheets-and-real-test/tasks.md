# Tasks

## 阶段一：Samplesheets 重命名

- [x] Task 1: 重命名 5 个测试 samplesheets 文件
  - [x] SubTask 1.1: `samplesheet.csv` → `test_online.csv`（在线 URL 数据，无 data_type）
  - [x] SubTask 1.2: `samplesheet_bam_test.csv` → `test_local_bam.csv`（本地 BAM 输入）
  - [x] SubTask 1.3: `samplesheet_local.csv` → `test_local_eccdna.csv`（本地 FASTQ，eccDNA only）
  - [x] SubTask 1.4: `samplesheet_local_with_gdna.csv` → `test_local_gdna.csv`（本地 FASTQ，含 gDNA control）
  - [x] SubTask 1.5: `samplesheet_integrated.csv` → `test_local_integrated.csv`（本地 FASTQ，integrated 模式）

## 阶段二：更新引用文件

- [x] Task 2: 更新配置文件中的 samplesheet 引用

  - [x] SubTask 2.1: `conf/test_local.config` line 18: `samplesheet_local.csv` → `test_local_eccdna.csv`
  - [x] SubTask 2.2: `conf/test_local_gdna.config` line 6: 注释中 `samplesheet_local_with_gdna.csv` → `test_local_gdna.csv`
  - [x] SubTask 2.3: `conf/test_local_gdna.config` line 20: `input` 路径更新为 `test_local_gdna.csv`

- [x] Task 3: 更新脚本中的 samplesheet 引用

  - [x] SubTask 3.1: `scripts/test_incremental_cache.py` line 18: `ORIGINAL_SAMPLESHEET` 路径更新
  - [x] SubTask 3.2: `scripts/test_incremental_cache.py` line 131: backup 文件名更新

- [x] Task 4: 更新文档中的 samplesheet 引用
  - [x] SubTask 4.1: `AGENTS.md` 中所有 `samplesheet_local.csv` → `test_local_eccdna.csv`
  - [x] SubTask 4.2: `AGENTS.md` 中所有 `samplesheet_integrated.csv` → `test_local_integrated.csv`
  - [x] SubTask 4.3: `testdatasets/README.md` line 17: `samplesheet_local.csv` → `test_local_eccdna.csv`

## 阶段三：stub 模式快速验证（重命名后）

- [x] Task 5: stub 验证重命名后流程仍正常
  - [x] SubTask 5.1: 执行 `nextflow run main.nf -profile test_local,docker -stub` 验证通过（52 tasks, 1m 42s）
  - [x] SubTask 5.2: 执行 `nextflow run main.nf -profile test_local_gdna,docker -stub` 验证通过（61 tasks, 1m 44s）

## 阶段四：真实模式测试

- [x] Task 6: 执行 `test_local` 真实模式（无 control）

  - [x] SubTask 6.1: 清理旧产物 `results/test_local/`（如已存在）
  - [x] SubTask 6.2: 执行 `nextflow run main.nf -profile test_local,docker`（真实模式，含注释库）
  - [x] SubTask 6.3: 记录执行结果（成功，61 tasks，1h26m58s，CPU hours 7.9）

- [x] Task 7: 执行 `test_local_gdna` 真实模式（有 control）
  - [x] SubTask 7.1: 清理旧产物 `results/test_local_gdna/`（如已存在）
  - [x] SubTask 7.2: 执行 `nextflow run main.nf -profile test_local_gdna,docker`（真实模式，含注释库 + gDNA control）
  - [x] SubTask 7.3: 记录执行结果（成功，70 tasks，全部 COMPLETED，无 FAILED）

## 阶段五：验证真实产出

- [x] Task 8: 验证 `test_local` 产出

  - [x] SubTask 8.1: 验证 `results/test_local/eccsplorer/*_candidates.bed` 非空（78/102/119 行）
  - [x] SubTask 8.2: 验证 `results/test_local/eccsplorer/*_blast.m6` 非空（85-445 bytes，BLAST 注释生效）
  - [x] SubTask 8.3: 验证 `results/test_local/eccsplorer/*_junction_reads.txt` 非空（144KB-409KB）
  - [x] SubTask 8.4: 验证 `results/test_local/eccsplorer/*_ecc_sequences.fasta` 非空（256KB-436KB）
  - [x] SubTask 8.5: 验证 `results/test_local/multiqc/multiqc_report.html` 存在（2.4MB）

- [x] Task 9: 验证 `test_local_gdna` 产出
  - [x] SubTask 9.1: 验证 `results/test_local_gdna/eccsplorer/*_candidates.bed` 非空（78/102/119 行）
  - [x] SubTask 9.2: 验证 `results/test_local_gdna/eccsplorer/*_blast.m6` 非空（55/78、68/102、66/119 个非空）
  - [x] SubTask 9.3: 验证 `results/test_local_gdna/eccsplorer/*_junction_reads.txt` 非空（144KB-408KB）
  - [x] SubTask 9.4: 验证 `results/test_local_gdna/eccsplorer/*_ecc_sequences.fasta` 非空（256KB-436KB）
  - [x] SubTask 9.5: 验证 `results/test_local_gdna/multiqc/multiqc_report.html` 存在（2.4MB）

## 阶段六：对比分析

- [x] Task 10: 对比有无 gDNA control 的结果差异
  - [x] SubTask 10.1: 统计 `results/test_local/eccsplorer/` 各样本 candidates.bed 行数（无 control：78/102/119）
  - [x] SubTask 10.2: 统计 `results/test_local_gdna/eccsplorer/` 各样本 candidates.bed 行数（有 control：78/102/119）
  - [x] SubTask 10.3: 对比两组数据，记录 gDNA control 对候选数量的影响
    - **发现**：两组候选数量完全相同，原因如下：
      - ECCSPLORER_WITH_CONTROL 命令正确传入 gDNA fastq（4 个 fastq 输入）
      - ECCsplorer mapping 阶段成功完成，生成 TR_hiconf-ECC-REGIONS.bed
      - ECCsplorer comparison/clustering 阶段失败：`comp_cl_tab_eccCandidates_list.csv not found`（RepeatExplorer 聚类数据库未配置）
      - 模块脚本中的 `|| true` 使流程继续执行，`cp` 复制了未过滤的 mapping 结果
      - 因此 `test_local_gdna` 的 candidates.bed 实际是 eccDNA-only mapping 结果，非 control-filtered 结果
      - **结论**：gDNA control 比对功能目前因 RepeatExplorer 聚类缺失而未实际生效，需后续配置 RepeatExplorer 数据库才能完整运行 comparison 模式

# Task Dependencies

- Task 2-4（更新引用）依赖 Task 1（重命名完成）
- Task 5（stub 验证）依赖 Task 2-3（配置引用更新完成）
- Task 6-7（真实测试）依赖 Task 5（stub 验证通过）
- Task 8-9（产出验证）依赖 Task 6-7（真实运行完成）
- Task 10（对比分析）依赖 Task 8-9（两组产出都验证完成）

# 并行化说明

- Task 2、3、4 可并行（更新不同文件）
- Task 6 和 Task 7 建议串行（本地资源有限，避免 Docker 容器竞争）
- Task 8 和 Task 9 可并行（验证不同目录的产出）
