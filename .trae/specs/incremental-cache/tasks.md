# circDNA 增量/减量缓存实现 - 实施计划

## [x] Task 1: 修改 input_check/main.nf 实现分离验证和 Channel 创建

- **Priority**: high
- **Depends On**: None
- **Description**:
  - 修改 `subworkflows/local/input_check/main.nf` 文件
  - 将 SAMPLESHEET_CHECK 的调用与 channel 创建解耦
  - Channel 创建改为直接使用 `file(samplesheet).splitCsv()` 从原始 CSV 解析
  - SAMPLESHEET_CHECK 仍然执行，但其输出不用于 channel 创建
  - 同时解决文件中存在的 merge conflict markers (<<<<<<<, =======, >>>>>>>)
  - 确保 FASTQ、BAM、long-read 三种模式均正确处理
- **Implementation Details**:

  ```nextflow
  workflow INPUT_CHECK {
      take:
      samplesheet

      main:
      // 1. 独立验证 (仅用于版本信息)
      SAMPLESHEET_CHECK ( samplesheet )
      ch_versions = SAMPLESHEET_CHECK.out.versions

      // 2. 直接从原始 CSV 创建 channel (关键改动)
      if ( params.input_format == "FASTQ" ) {
          file(samplesheet).splitCsv(header:true, sep:',')
              .map { row -> create_fastq_channels(row) }
              .set { reads }
      } else if ( params.input_format == "BAM" ) {
          file(samplesheet).splitCsv(header:true, sep:',')
              .map { row -> create_bam_channels(row) }
              .set { reads }
      } else if ( params.protocol in ["pacbio", "ont"] ) {
          file(samplesheet).splitCsv(header:true, sep:',')
              .map { row -> create_long_read_channels(row) }
              .set { reads }
      }

      emit:
      reads
      versions = ch_versions
  }
  ```

- **Acceptance Criteria Addressed**: AC-1, AC-2, AC-3, AC-4, AC-5
- **Test Requirements**:
  - `programmatic` TR-1.1: 使用 test_local.config 首次运行 3 个样本成功完成
  - `programmatic` TR-1.2: 修改 samplesheet 增加 1 个样本后 -resume 运行，原 3 个样本 CACHED
  - `programmatic` TR-1.3: 修改 samplesheet 减少 1 个样本后 -resume 运行，保留样本 CACHED
  - `programmatic` TR-1.4: 同时增减样本后 -resume 运行，保留 CACHED + 新增 NEW
  - `programmatic` TR-1.5: 输出结果与修改前完全一致
- **Notes**:
  - 需要先处理 merge conflict markers
  - 保持 `create_fastq_channels`、`create_bam_channels`、`create_long_read_channels` 函数不变
  - `file(samplesheet).splitCsv()` 中 `file()` 会确保路径存在性检查

## [x] Task 2: 创建缓存验证测试脚本

- **Priority**: high
- **Depends On**: Task 1
- **Description**:
  - 创建测试脚本用于验证增量/减量缓存行为
  - 基于 test_local.config 的 3 个测试样本
  - 生成不同版本的 samplesheet（增加/减少样本）
  - 自动运行 pipeline 并检查任务缓存状态
- **Acceptance Criteria Addressed**: AC-1, AC-2, AC-3
- **Test Requirements**:
  - `programmatic` TR-2.1: 脚本能自动创建测试 samplesheet 变体
  - `programmatic` TR-2.2: 脚本能正确解析 nextflow log 中的 CACHED/NEW 状态
  - `programmatic` TR-2.3: 脚本能给出明确的测试通过/失败结论
- **Notes**:
  - 测试脚本应放在 `scripts/` 目录
  - 可使用 `nextflow log -f 'CACHED|NEW'` 过滤日志

## [x] Task 3: 端到端验证测试

- **Priority**: high
- **Depends On**: Task 1, Task 2
- **Description**:
  - 执行完整的增量缓存验证流程 (Run 1-4 全部完成)
  - 执行完整的减量缓存验证流程
  - 记录每个步骤的 nextflow log 输出
  - 对比修改前后的 pipeline 输出一致性
  - **关键发现**: 样本归并机制 (circdna.nf line 120) 导致样本级缓存不可用, 仅参考基因组任务缓存
- **Acceptance Criteria Addressed**: AC-1, AC-2, AC-3, AC-4, AC-5 (部分, 受样本归并机制限制)
- **Test Requirements**:
  - `programmatic` TR-3.1: 首次运行完成，3 个样本全部处理 ✅ (Run 1: 27 tasks)
  - `programmatic` TR-3.2: 增加样本运行，原样本 CACHED ⚠️ (仅 BWA_INDEX/SAMTOOLS_FAIDX CACHED, 样本级任务全 NEW)
  - `programmatic` TR-3.3: 减少样本运行，保留样本 CACHED ⚠️ (同上)
  - `programmatic` TR-3.4: 同时增减运行，缓存正确 ⚠️ (同上)
  - `programmatic` TR-3.5: 输出文件 md5 与修改前一致 ⚠️ (因样本归并, 不同 samplesheet 输出不同)
- **Notes**:
  - 使用 `conda activate nextflow` 环境
  - 使用 test_local.config 配置
  - 测试脚本: `scripts/test_incremental_cache.py`
  - 测试数据: `testdatasets/samplesheet/test_*.csv`
  - Run 1: test_3.csv (3 samples), Run 2: test_4.csv (4 samples), Run 3: test_2.csv (2 samples), Run 4: test_mixed.csv (3 samples)
  - 详见 checklist.md 中的缓存测试结果汇总表

## [x] Task 4: 代码审查和文档更新

- **Priority**: medium
- **Depends On**: Task 1, Task 3
- **Description**:
  - 审查修改后的代码，确保符合 Nextflow 最佳实践
  - 检查是否需要更新 CHANGELOG.md
  - 确认所有 merge conflict markers 已正确处理 ✅
  - 验证代码风格与项目一致 ✅
  - 测试待 shell 执行后完成
- **Acceptance Criteria Addressed**: AC-6
- **Test Requirements**:
  - `human-judgement` TR-4.1: 代码改动最小且逻辑清晰
  - `human-judgement` TR-4.2: 无 merge conflict markers 残留
  - `human-judgement` TR-4.3: 代码风格与项目一致
- **Notes**:
  - 参考 `.trae` 中已有 spec 的变更记录格式
  - 可能需要更新 CHANGELOG.md
