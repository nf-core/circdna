# Fix CReSIL identify_wgls KeyError Bug Spec

## Why
CReSIL `identify_wgls` 模块在运行时出现 `KeyError: 'name'` 错误。**根本原因**是代码中比较 strand 时使用整数 `1` 和 `-1`，但 trim 输出的 strand 值是字符串 `'+'` 和 `'-'`。这导致 `df_breaks_bed` 为空 DataFrame，进而 `breaks_bdg.to_dataframe()` 返回没有列的空 DataFrame，访问 `'name'` 列时抛出 KeyError。

此外，还发现了 bedtools merge 排序问题：`filtered_region_depth_bed.merge()` 要求输入文件已排序，否则会报错。

## What Changes
在 `/Users/siyangming/nextflow_nf_core/bio.nf/cresil/cresil/cli/identify_wgls.py` 中修复以下问题：

1. **Strand 比较修复**（主要 bug）：
   - 将 `trim_sup['strand'] == 1` 改为 `trim_sup['strand'] == '+'`
   - 将 `trim_sup['strand'] == -1` 改为 `trim_sup['strand'] == '-'`
   - 修复位置：第 430-431 行（linkage mode）、第 517-518 行（depth mode）
   - 修复 `check_breakpoint_direction` 函数中的 strand 比较（第 198-221 行）
   - 修复 `majority_strand` 函数中的 strand 计数（第 58 行）

2. **Bedtools merge 排序修复**：
   - 在 `merge()` 前添加 `.sort()`
   - 修复位置：第 511-512 行

## Impact
- Affected specs: CReSIL identify_wgls 模块
- Affected code: `/Users/siyangming/nextflow_nf_core/bio.nf/cresil/cresil/cli/identify_wgls.py`
- 受影响的用户: 所有使用 CReSIL identify_wgls 功能的用户
- 影响: 修复后 identify_wgls 模块可以正常运行，不再出现 KeyError 和排序错误

## ADDED Requirements

### Requirement: Strand 值使用正确的类型
CReSIL identify_wgls 代码 SHALL 使用字符串 `'+'` 和 `'-'` 来比较 strand 值，与 trim 模块输出的格式保持一致。

#### Scenario: 正链 reads 过滤
- **WHEN** 过滤正链（plus strand）的 trimmed reads
- **THEN** 应使用 `trim_sup['strand'] == '+'` 进行比较
- **AND** 能正确匹配所有正链 reads

#### Scenario: 负链 reads 过滤
- **WHEN** 过滤负链（minus strand）的 trimmed reads
- **THEN** 应使用 `trim_sup['strand'] == '-'` 进行比较
- **AND** 能正确匹配所有负链 reads

#### Scenario: breakpoint 方向检查
- **WHEN** 调用 `check_breakpoint_direction` 函数
- **THEN** 函数中的 strand 比较应使用 `'+'` 和 `'-'`
- **AND** 能正确判断断点方向

#### Scenario: 多数链计算
- **WHEN** 调用 `majority_strand` 函数
- **THEN** 函数应统计 `'+'` 和 `'-'` 的数量
- **AND** 返回正确的多数链方向

### Requirement: Bedtools merge 输入排序
CReSIL identify_wgls 代码 SHALL 在调用 `merge()` 前确保 BED 文件已排序。

#### Scenario: depth mode 区域合并
- **WHEN** 调用 `filtered_region_depth_bed.merge()`
- **THEN** 应先调用 `.sort()` 确保输入已排序
- **AND** 不会出现 "Sorted input specified" 错误

## MODIFIED Requirements

### Requirement: genome_coverage 输出处理
**原有理解**: 认为是 `to_dataframe()` 列名不匹配导致 KeyError

**实际根本原因**: strand 类型不匹配导致 `df_breaks_bed` 为空，进而 `to_dataframe()` 返回无列的空 DataFrame

**修复后**: 修复 strand 比较，确保 `breaks_bdg` 有数据，`to_dataframe()` 返回的 DataFrame 包含默认列名 `['chrom', 'start', 'end', 'name']`

## REMOVED Requirements
无需移除任何要求。

## Root Cause Analysis
| 问题 | 原因 | 修复 |
|------|------|------|
| `KeyError: 'name'` | strand 比较使用 `1`/`-1`（int），但实际是 `'+'`/`'-'`（str），导致空 DataFrame | 改为 `'+'`/`'-'` 比较 |
| `BEDToolsError: Sorted input` | bedtools merge 要求输入排序 | 添加 `.sort()` 调用 |
