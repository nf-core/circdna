# circDNA 样本级增量缓存 - Product Requirement Document

## Overview

- **Summary**: 在 samplesheet 中引入 `sample` + `lane` 两列，明确区分样本 ID 和测序 lane，替代当前基于命名约定 `meta.id.split('_')[0..-2]` 的样本归并机制，实现真正的样本级增量缓存。
- **Purpose**: 解决当前 samplesheet 中 `circdna_1/2/3` 被错误归并为单一 `circdna` 样本的问题，使每个样本独立处理，支持 `-resume` 时样本级任务缓存。
- **Target Users**: nf-core/circdna 管道用户，尤其是需要批量处理多样本并利用缓存加速的研究者。

## Goals

- 引入 `lane` 列，与 `sample` 列配合明确标识样本与 lane 的关系
- 移除基于命名约定的隐式样本归并 (`meta.id.split('_')[0..-2].join('_')`)
- 同一样本的不同 lane 仍可通过 CAT_FASTQ 合并
- 不同样本之间相互独立，支持 `-resume` 时的样本级缓存
- 保持向后兼容：无 `lane` 列时，每行视为独立样本

## Non-Goals (Out of Scope)

- 不修改 BAM 输入模式 (input_format == "BAM") 的样本归并逻辑
- 不修改 long-read 模式 (protocol in ["pacbio", "ont"]) 的样本处理
- 不改变 SAMPLESHEET_CHECK 的输出格式（当前不使用其输出）
- 不修改下游分析流程的逻辑
- 不添加新的 pipeline 参数

## Background & Context

### 当前机制

[workflows/circdna.nf:118-121](file:///Users/siyangming/nextflow_nf_core/circdna.nf/workflows/circdna.nf#L118-L121) 中存在隐式样本归并逻辑：

```nextflow
.map { meta, fastq ->
    meta.id = meta.id.split('_')[0..-2].join('_')
    [ meta, fastq ] }
```

此逻辑假设样本名中最后一个 `_` 后缀是 lane 编号（如 `sampleA_L1`, `sampleA_L2`），将其后缀去掉后归并。

### 问题

1. **独立样本被错误归并**：`circdna_1`, `circdna_2`, `circdna_3` 被归并为 `circdna`，CAT_FASTQ 合并所有 fastq
2. **缓存失效**：增减任何样本都会改变合并后的 fastq，所有下游任务重新执行
3. **语义不明确**：依赖命名约定容易出错，用户难以理解

### 验证结果

在 [incremental-cache spec](file:///Users/siyangming/nextflow_nf_core/circdna.nf/.trae/specs/incremental-cache/checklist.md) 的 4 次缓存测试中，仅参考基因组任务（BWA_INDEX, SAMTOOLS_FAIDX）被缓存，样本级任务全部为 NEW，证实了此问题。

## Functional Requirements

### FR-1: 支持 lane 列

- samplesheet 可选包含 `lane` 列
- 有 `lane` 列时：同 `sample` + 不同 `lane` 的多行归并为同一样本
- 无 `lane` 列时：每行视为独立样本（与当前 nf-core 标准行为一致）

### FR-2: 移除隐式归并

- 删除 `meta.id.split('_')[0..-2].join('_')` 隐式归并逻辑
- 样本 ID 直接使用 `sample` 列的值
- lane 信息保存在 `meta.lane` 中供参考

### FR-3: 多 lane 合并

- 同一样本有多个 lane 时，仍通过 CAT_FASTQ 合并 fastq
- 合并逻辑与当前一致（按 lane 顺序拼接）

### FR-4: 向后兼容

- 无 `lane` 列的 samplesheet 仍可正常工作
- `sample` 列值直接作为样本 ID
- 每个样本独立处理，不进行归并

### FR-5: 更新 samplesheet 校验

- `check_samplesheet.py` 支持可选 `lane` 列
- 有 `lane` 列时验证其非空
- 无 `lane` 列时保持当前行为

## Non-Functional Requirements

### NFR-1: 缓存性能

- `-resume` 时，未修改样本的所有样本级任务应显示 CACHED
- 新增样本的任务显示 NEW
- 删除样本的任务不执行

### NFR-2: 输出一致性

- 同一样本输入数据相同的情况下，输出文件 md5 一致

## Constraints

- **技术**: Nextflow DSL2, nf-core 管道规范
- **依赖**: 当前的 CAT_FASTQ 模块、INPUT_CHECK 子工作流
- **兼容性**: 必须兼容现有无 `lane` 列的 samplesheet

## Assumptions

- 用户理解 `lane` 列的含义（同一样本的多条测序 lane）
- lane 的顺序由 samplesheet 中行的顺序决定
- 测试数据可通过添加 `lane` 列重新生成 samplesheet

## Acceptance Criteria

### AC-1: 有 lane 列时样本归并正确

- **Given**: samplesheet 包含 `sample`, `fastq_1`, `fastq_2`, `lane` 四列，同一样本有多个 lane
- **When**: pipeline 运行
- **Then**: 同一样本的不同 lane 被 CAT_FASTQ 合并为一个 fastq 后进入下游分析
- **Verification**: `programmatic`
- **Notes**: 验证 CAT_FASTQ 输入包含对应 lane 的所有 fastq

### AC-2: 无 lane 列时样本独立

- **Given**: samplesheet 只有 `sample`, `fastq_1`, `fastq_2` 三列（旧格式）
- **When**: pipeline 运行
- **Then**: 每行作为独立样本处理，不进行归并
- **Verification**: `programmatic`

### AC-3: 增量缓存 - 新增样本

- **Given**: Run A 使用 3 个样本成功完成
- **When**: Run B 使用 `-resume` 从 Run A 恢复，samplesheet 新增 1 个样本
- **Then**: 原有 3 个样本的所有样本级任务显示 CACHED，新样本任务显示 NEW
- **Verification**: `programmatic`
- **Notes**: 参考基因组任务（BWA_INDEX, SAMTOOLS_FAIDX）应始终 CACHED

### AC-4: 增量缓存 - 减少样本

- **Given**: Run A 使用 4 个样本成功完成
- **When**: Run B 使用 `-resume` 从 Run A 恢复，samplesheet 减少为 2 个样本
- **Then**: 保留的 2 个样本的所有样本级任务显示 CACHED
- **Verification**: `programmatic`

### AC-5: 增量缓存 - 混合增减

- **Given**: Run A 使用 3 个样本成功完成 (sampleA, sampleB, sampleC)
- **When**: Run B 使用 `-resume` 从 Run A 恢复，samplesheet 为 (sampleA, sampleC, sampleD)
- **Then**: sampleA 和 sampleC 显示 CACHED，sampleD 显示 NEW，sampleB 不执行
- **Verification**: `programmatic`

### AC-6: check_samplesheet.py 兼容 lane 列

- **Given**: samplesheet 包含可选的 `lane` 列
- **When**: SAMPLESHEET_CHECK 运行
- **Then**: 验证通过，lane 列信息保留
- **Verification**: `programmatic`

### AC-7: 输出结果正确性

- **Given**: 相同输入数据
- **When**: 分别使用旧格式（无 lane 列）和新格式（有 lane 列，每样本单 lane）运行
- **Then**: 输出结果一致
- **Verification**: `programmatic`
- **Notes**: 验证关键输出文件的 md5 或内容一致

## Open Questions

- [ ] BAM 模式是否也需要 lane 列支持？（当前 BAM 模式没有 CAT_FASTQ 归并）
- [ ] long-read 模式是否需要 lane 列支持？
