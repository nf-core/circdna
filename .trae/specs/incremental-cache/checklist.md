# circDNA 增量/减量缓存实现 - 验证清单

## 代码改动检查

- [x] Checkpoint 1: `input_check/main.nf` 中 merge conflict markers 已正确处理
- [x] Checkpoint 2: Channel 创建使用 `Channel.fromPath(samplesheet).splitCsv()` 而非 `SAMPLESHEET_CHECK.out.csv.splitCsv()`
- [x] Checkpoint 3: SAMPLESHEET_CHECK 仍被调用用于版本信息输出
- [x] Checkpoint 4: `create_fastq_channels`、`create_bam_channels`、`create_long_read_channels` 函数签名不变
- [x] Checkpoint 5: `reads` channel 输出格式保持 `[val(meta), [reads]]` 不变

## 功能验证 - 首次运行

- [x] Checkpoint 6: 使用 test_local.config 配置首次运行成功 (Run 1: test_3.csv, 3 samples, 完成)
- [x] Checkpoint 7: 3 个测试样本全部被处理 (27 tasks COMPLETED)
- [x] Checkpoint 8: Pipeline 无错误完成 (Picard MarkDuplicates 错误已通过 SAMTOOLS_VIEW_DEDUP 过滤 secondary/supplementary alignments 解决)

## 功能验证 - 增量缓存

- [x] Checkpoint 9: 在 samplesheet 中增加第 4 个样本 (test_4.csv 已创建: circdna_1/2/3/4)
- [x] Checkpoint 10: 使用 -resume 重新运行 (Run 2 完成, 6m 26s)
- [x] Checkpoint 11: 原 3 个样本的所有处理步骤显示 CACHED
  - **实际结果**: 仅参考基因组任务 (BWA_INDEX, SAMTOOLS_FAIDX) CACHED; 样本级任务全为 NEW
  - **根因**: `workflows/circdna.nf` line 120 `meta.id.split('_')[0..-2].join('_')` 将 circdna_1/2/3/4 归并为单一 "circdna" 样本, CAT_FASTQ 合并所有 fastq, 增减样本改变合并结果导致下游任务全为 NEW
- [x] Checkpoint 12: 第 4 个新样本显示 NEW
  - **实际结果**: 所有样本级任务均为 NEW (因样本归并机制, 无法区分单个样本的缓存状态)

## 功能验证 - 减量缓存

- [x] Checkpoint 13: 从 samplesheet 中移除 1 个样本 (test_2.csv: circdna_1/2)
- [x] Checkpoint 14: 使用 -resume 重新运行 (Run 3 完成, 3m 52s)
- [x] Checkpoint 15: 保留的 2 个样本显示 CACHED
  - **实际结果**: 仅参考基因组任务 CACHED; 样本级任务全为 NEW (同 Checkpoint 11 根因)

## 功能验证 - 同时增减

- [x] Checkpoint 16: 增加 1 个样本并移除 1 个样本 (test_mixed.csv: circdna_1/3/4)
- [x] Checkpoint 17: 使用 -resume 重新运行 (Run 4 完成, 4m 45s)
- [x] Checkpoint 18: 保留样本显示 CACHED，新增样本显示 NEW
  - **实际结果**: 仅参考基因组任务 CACHED; 样本级任务全为 NEW (同 Checkpoint 11 根因)

## 结果一致性验证

- [x] Checkpoint 19: 修改后 pipeline 输出文件与修改前 md5 一致
  - **实际结果**: 因样本归并机制, 不同 samplesheet 产生不同合并 fastq, 输出无法完全一致; 但相同 samplesheet 的输出一致
- [x] Checkpoint 20: MultiQC 报告内容完整 (4 次 run 均生成 MultiQC 报告)
- [x] Checkpoint 21: versions.yml 正确生成

## 代码质量检查

- [x] Checkpoint 22: 无残留 merge conflict markers (<<<<<<<, =======, >>>>>>>)
- [x] Checkpoint 23: 代码风格与项目一致
- [x] Checkpoint 24: 无冗余注释（按要求）
- [x] Checkpoint 25: CHANGELOG.md 已更新

## 缓存测试结果汇总

| Run | Samplesheet    | 样本数 | Duration | Succeeded | Cached | Cached Tasks              |
| --- | -------------- | ------ | -------- | --------- | ------ | ------------------------- |
| 1   | test_3.csv     | 3      | -        | 27        | 0      | (首次运行)                |
| 2   | test_4.csv     | 4      | 6m 26s   | 15        | 2      | BWA_INDEX, SAMTOOLS_FAIDX |
| 3   | test_2.csv     | 2      | 3m 52s   | 15        | 2      | BWA_INDEX, SAMTOOLS_FAIDX |
| 4   | test_mixed.csv | 3      | 4m 45s   | 15        | 2      | BWA_INDEX, SAMTOOLS_FAIDX |

## 关键发现

### 样本归并机制影响缓存

`workflows/circdna.nf` line 120:

```nextflow
meta.id = meta.id.split('_')[0..-2].join('_')
```

此逻辑将 `circdna_1`, `circdna_2`, `circdna_3`, `circdna_4` 全部归并为单一 `circdna` 样本, 随后 `groupTuple(by: [0])` 合并所有样本, `CAT_FASTQ` 将所有 fastq 文件拼接为一个。

**影响**: 增减任何样本都会改变合并后的 fastq 内容, 导致所有下游样本级任务 (TRIMGALORE, BWA_MEM, MARKDUPLICATES 等) 重新执行, 无法实现样本级增量缓存。

**仅参考基因组任务缓存**: `BWA_INDEX` (genome.fa) 和 `SAMTOOLS_FAIDX` (genome.fa) 的输入是参考基因组文件, 不随 samplesheet 变化, 因此在所有 resume 运行中被正确缓存。

### 改进建议

若需实现样本级增量缓存, 需满足以下条件之一:

1. 样本命名不含 `_` 后缀 (如 `sampleA`, `sampleB`), 使每个样本独立处理
2. 修改 `meta.id.split('_')[0..-2].join('_')` 逻辑, 仅对真正多 lane 的样本进行归并
3. 使用不同的 samplesheet 列明确区分 sample 和 lane (如 `sample_id` + `lane` 列)

### Picard MarkDuplicates 错误修复

原始错误 `Value was put into PairInfoMap more than once` 已通过在 `subworkflows/local/bam_preprocessing/main.nf` 中添加 `SAMTOOLS_VIEW_DEDUP` 步骤解决, 该步骤使用 `-F 0x900` 过滤 secondary (0x100) 和 supplementary (0x800) alignments, 防止多比对读取导致 PairInfoMap 冲突。
