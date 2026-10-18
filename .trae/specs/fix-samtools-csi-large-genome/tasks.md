# Tasks

- [x] Task 1: 修复 `conf/large_genome.config` 的 SAMTOOLS_INDEX 匹配范围

  - [x] 将 `withName: 'SAMTOOLS_INDEX'` 改为 `withName: '.*SAMTOOLS_INDEX.*'`（正则匹配，覆盖全部 SAMTOOLS_INDEX 实例）
  - [x] 验证正则不会误匹配其他进程（如 SAMTOOLS_SORT 等）

- [x] Task 2: 更新 `SERVER_RUN_GUIDE.md` 标注 Tragopogon_porrifolius hap1 为大基因组

  - [x] hap1 命令附加 `-c circdna.nf/conf/large_genome.config`
  - [x] 标题/说明中标注 hap1 大基因组（hap2 保持不变）
  - [x] 注意事项"大基因组"列表补充 Tragopogon_porrifolius hap1

- [x] Task 3: 更新 `README.md` 说明大基因组 CSI 用法

  - [x] Usage 部分添加大基因组物种附加 `-c conf/large_genome.config` 的提示（CSI 索引，支持 >512Mb 染色体）

- [x] Task 4: 更新 `CHANGELOG.md`（**不修改 nextflow.config，不 bump 版本**）

  - [x] 在现有 v4.4.1 条目下追加记录：large_genome.config 匹配修复（覆盖全部 SAMTOOLS_INDEX 实例）+ hap1 大基因组标注 + README 更新

- [x] Task 5: Git 提交范围限定（同步 GitHub）
  - [x] 提交前先 `git status` 确认工作区变更清单
  - [x] 仅 `git add` 本次大基因组相关文件：`conf/large_genome.config`、`SERVER_RUN_GUIDE.md`、`README.md`、`CHANGELOG.md`
  - [x] 提交并推送，确认其他无关更改未被包含

# Task Dependencies

- [Task 4] depends on [Task 1]、[Task 2]、[Task 3]（变更记录需包含全部改动）
- [Task 5] depends on [Task 1]、[Task 2]、[Task 3]、[Task 4]（提交需包含全部改动）
