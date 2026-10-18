# circDNA 样本级增量缓存 - 验证清单

## 代码改动检查

- [x] Checkpoint 1: `bin/check_samplesheet.py` 支持可选 lane 列
  - [x] 三列 header（无 lane）仍正常工作
  - [x] 四列 header（有 lane）正常解析
  - [x] lane 列为空时报错
- [x] Checkpoint 2: `subworkflows/local/input_check/main.nf` 的 create_fastq_channels 支持 lane 列
  - [x] 有 lane 列时 meta.lane 被正确设置
  - [x] 无 lane 列时 meta.lane 为 null
  - [x] meta.id 直接使用 row.sample，不做 split 处理
- [x] Checkpoint 3: `workflows/circdna.nf` 样本归并逻辑更新
  - [x] 有 lane 列时直接按 sample 列值分组（无 split）
  - [x] 无 lane 列时保留 \_T\d+ 后缀归并（向后兼容）
  - [x] groupTuple 按 meta.id 分组（通过 id 字段提取）
  - [x] CAT_FASTQ 逻辑不变

## 功能验证 - 无 lane 列（向后兼容）

- [x] Checkpoint 4: 3 个独立样本（无 lane 列）全部独立处理
- [x] Checkpoint 5: 每个样本有独立的 TRIMGALORE, BWA_MEM, MARKDUPLICATES 等任务
- [x] Checkpoint 6: 输出文件以 sample 列值命名
- [x] Checkpoint 7: Pipeline 无错误完成

## 功能验证 - 有 lane 列（多 lane 合并）

- [x] Checkpoint 8: 2 个样本各 2 个 lane → CAT_FASTQ 任务数 == 2
- [x] Checkpoint 9: 每个 CAT_FASTQ 输入包含对应样本的所有 lane fastq
- [x] Checkpoint 10: 下游任务（TRIMGALORE, BWA_MEM 等）任务数 == 2（样本数）
- [x] Checkpoint 11: Pipeline 无错误完成

## 缓存验证 - 增量（新增样本）

- [x] Checkpoint 12: Run 1（3 样本）成功完成
- [x] Checkpoint 13: Run 2（4 样本，-resume）中原 3 样本的样本级任务全部 CACHED
- [x] Checkpoint 14: Run 2 中新样本的任务全部 NEW
- [x] Checkpoint 15: Run 2 中参考基因组任务（BWA_INDEX, SAMTOOLS_FAIDX）CACHED

## 缓存验证 - 减量（删除样本）

- [x] Checkpoint 16: Run 3（2 样本，-resume from Run 2）中保留的 2 样本 CACHED
- [x] Checkpoint 17: Run 3 中删除的样本不执行

## 缓存验证 - 混合增减

- [x] Checkpoint 18: Run 4（混合增减，-resume from Run 3）中保留样本 CACHED
- [x] Checkpoint 19: Run 4 中新增样本 NEW
- [x] Checkpoint 20: Run 4 中删除样本不执行

## 结果一致性验证

- [x] Checkpoint 21: 同一样本相同输入 → 输出文件一致（缓存命中证明）
- [x] Checkpoint 22: MultiQC 报告包含所有样本
- [x] Checkpoint 23: versions.yml 正确生成

## 代码质量检查

- [x] Checkpoint 24: 代码风格与项目一致
- [x] Checkpoint 25: 无冗余注释（按要求）
