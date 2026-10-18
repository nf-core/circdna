# Tasks

- [x] Task 1: 修改 `bam_preprocessing/main.nf` 中 `ch_fasta_fai` 的构造方式

  - [x] SubTask 1.1: 将 `.join().map()` 改为 `.join().map().collect().map { it[0] }`
  - [x] SubTask 1.2: 添加注释说明为什么需要 `.collect().map { it[0] }` 而不是 `.first()`

- [x] Task 2: 使用 test_local 验证无警告且流程正常

  - [x] SubTask 2.1: 用 test_local,docker 运行（全新 work-dir，无 -resume）
  - [x] SubTask 2.2: 确认 `.nextflow.log` 中无 "first is useless" 警告（WARN 计数 = 0）
  - [x] SubTask 2.3: 确认所有 3 个样本都完成处理（CIRCLEMAP_REALIGN 3 of 3 ✔，CANDIDATE_MERGE 3 of 3 ✔，64 任务全部成功）

- [x] Task 3: 选择性提交到 GitHub
  - [x] SubTask 3.1: stash 本地其他修改
  - [x] SubTask 3.2: 只提交 `bam_preprocessing/main.nf` 的修改
  - [x] SubTask 3.3: 推送到 GitHub master 分支（commit 32cdcbc）
  - [x] SubTask 3.4: 恢复本地完整修改

# Task Dependencies

- Task 2 depends on Task 1
- Task 3 depends on Task 2
