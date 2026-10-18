# circDNA 样本级增量缓存 - The Implementation Plan (Decomposed and Prioritized Task List)

## [x] Task 1: 修改 check_samplesheet.py 支持可选 lane 列

- **Priority**: high
- **Depends On**: None
- **Description**:
  - 修改 `bin/check_samplesheet.py` 支持可选的 `lane` 列
  - FASTQ 模式下，header 允许 `["sample", "fastq_1", "fastq_2"]` 或 `["sample", "fastq_1", "fastq_2", "lane"]`
  - 有 `lane` 列时验证其非空
  - 输出时保留 `lane` 列信息
  - 保持 BAM 模式不变
- **Acceptance Criteria Addressed**: AC-6
- **Test Requirements**:
  - `programmatic` TR-1.1: 三列 samplesheet（无 lane）验证通过，行为不变
  - `programmatic` TR-1.2: 四列 samplesheet（有 lane）验证通过
  - `programmatic` TR-1.3: lane 列为空时报错
- **Notes**: SAMPLESHEET_CHECK 的输出当前未被使用（INPUT_CHECK 直接从原始 CSV 读取），但仍需保持校验脚本正确性

## [x] Task 2: 修改 INPUT_CHECK 支持 lane 列

- **Priority**: high
- **Depends On**: Task 1
- **Description**:
  - 修改 `subworkflows/local/input_check/main.nf` 的 `create_fastq_channels` 函数
  - 检测 CSV 是否包含 `lane` 列
  - 有 `lane` 列时，将 lane 值存入 `meta.lane`
  - 无 `lane` 列时，`meta.lane` 为 null
  - `meta.id` 始终使用 `row.sample`（不再做 split 处理）
- **Acceptance Criteria Addressed**: AC-1, AC-2, AC-6
- **Test Requirements**:
  - `programmatic` TR-2.1: 三列 samplesheet 时 meta.id == row.sample，meta.lane == null
  - `programmatic` TR-2.2: 四列 samplesheet 时 meta.id == row.sample，meta.lane == row.lane
  - `programmatic` TR-2.3: BAM 模式行为不变

## [x] Task 3: 修改 circdna.nf 样本归并逻辑

- **Priority**: high
- **Depends On**: Task 2
- **Description**:
  - 移除 `workflows/circdna.nf` 中 line 120 的 `meta.id.split('_')[0..-2].join('_')`
  - 样本归并改为直接使用 `meta.id`（即 sample 列值）
  - `groupTuple(by: [0])` 仍按 meta.id 分组
  - 同一样本多行（多个 lane）时被归并，单行（单 lane 或无 lane 列）时保持独立
  - CAT_FASTQ 逻辑保持不变
- **Acceptance Criteria Addressed**: AC-1, AC-2, AC-7
- **Test Requirements**:
  - `programmatic` TR-3.1: 无 lane 列 + 3 个独立样本 → 3 个独立样本流程
  - `programmatic` TR-3.2: 有 lane 列 + 同一样本 2 个 lane → 1 个样本 + CAT_FASTQ 合并
  - `programmatic` TR-3.3: 输出文件命名以 sample 列值为前缀

## [x] Task 4: 创建带 lane 列的测试 samplesheet

- **Priority**: high
- **Depends On**: Task 1
- **Description**:
  - 创建 `testdatasets/samplesheet/test_3_lane.csv`（3 个样本，无 lane 列对照）
  - 创建 `testdatasets/samplesheet/test_multilane.csv`（2 个样本，每个 2 个 lane）
  - 创建 `testdatasets/samplesheet/test_incr_lane_3.csv`（3 个独立样本，用于 Run 1）
  - 创建 `testdatasets/samplesheet/test_incr_lane_4.csv`（4 个独立样本，用于 Run 2 增量）
  - 创建 `testdatasets/samplesheet/test_incr_lane_2.csv`（2 个独立样本，用于 Run 3 减量）
  - 创建 `testdatasets/samplesheet/test_incr_lane_mixed.csv`（3 个样本混合，用于 Run 4）
- **Acceptance Criteria Addressed**: AC-3, AC-4, AC-5
- **Test Requirements**:
  - `programmatic` TR-4.1: 所有 samplesheet 格式正确，可被 INPUT_CHECK 解析

## [x] Task 5: 端到端缓存验证测试

- **Priority**: high
- **Depends On**: Task 3, Task 4
- **Description**:
  - Run 1: 3 个样本首次运行（test_incr_lane_3.csv）
  - Run 2: 4 个样本增量缓存测试（test_incr_lane_4.csv，-resume from Run 1）
  - Run 3: 2 个样本减量缓存测试（test_incr_lane_2.csv，-resume from Run 2）
  - Run 4: 混合增减缓存测试（test_incr_lane_mixed.csv，-resume from Run 3）
  - 每次运行使用独立的 work-dir 和 outdir
  - 通过 trace 文件验证 CACHED/NEW 状态
- **Acceptance Criteria Addressed**: AC-3, AC-4, AC-5
- **Test Requirements**:
  - `programmatic` TR-5.1: Run 1 成功完成，3 个样本全部处理
  - `programmatic` TR-5.2: Run 2 中原 3 个样本的所有样本级任务 CACHED，新样本 NEW
  - `programmatic` TR-5.3: Run 3 中保留的 2 个样本 CACHED
  - `programmatic` TR-5.4: Run 4 中保留样本 CACHED，新增样本 NEW
  - `programmatic` TR-5.5: 参考基因组任务（BWA_INDEX, SAMTOOLS_FAIDX）始终 CACHED

## [x] Task 6: 多 lane 合并验证测试

- **Priority**: medium
- **Depends On**: Task 3, Task 4
- **Description**:
  - 使用 test_multilane.csv 运行 pipeline（2 个样本，每个 2 个 lane）
  - 验证 CAT_FASTQ 正确合并同一样本的不同 lane
  - 验证下游任务按样本数（2 个）执行，而非按 lane 数（4 个）
- **Acceptance Criteria Addressed**: AC-1
- **Test Requirements**:
  - `programmatic` TR-6.1: CAT_FASTQ 任务数 == 样本数（2）
  - `programmatic` TR-6.2: 每个 CAT_FASTQ 输入包含对应样本的所有 lane fastq
  - `programmatic` TR-6.3: 下游任务（TRIMGALORE, BWA_MEM 等）任务数 == 样本数

## [x] Task 7: 向后兼容验证

- **Priority**: medium
- **Depends On**: Task 3, Task 4
- **Description**:
  - 使用旧格式 samplesheet（无 lane 列）运行 pipeline
  - 验证每个样本独立处理
  - 验证输出与修改前一致
- **Acceptance Criteria Addressed**: AC-2, AC-7
- **Test Requirements**:
  - `programmatic` TR-7.1: 无 lane 列时每行作为独立样本
  - `programmatic` TR-7.2: 输出文件命名正确
  - `human-judgement` TR-7.3: 行为符合 nf-core 管道标准

# Task Dependencies

- Task 2 depends on Task 1
- Task 3 depends on Task 2
- Task 4 depends on Task 1
- Task 5 depends on Task 3, Task 4
- Task 6 depends on Task 3, Task 4
- Task 7 depends on Task 3, Task 4
