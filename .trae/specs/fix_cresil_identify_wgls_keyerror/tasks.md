# Fix CReSIL identify_wgls KeyError Bug - The Implementation Plan

## [x] Task 1: 修复 strand 比较错误（根本原因）
- **Priority**: high
- **Depends On**: None
- **Description**: 
  - 修复 `identify_wgls.py` 中 strand 比较类型不匹配的问题
  - 代码使用整数 `1`/`-1` 比较，但 trim 输出的 strand 是字符串 `'+'`/`'-'`
  - 这导致 `df_breaks_bed` 为空，进而 `to_dataframe()` 返回无列的空 DataFrame，引发 `KeyError: 'name'`
  - 修复位置：
    - `majority_strand` 函数：`'-1'`/`'1'` → `'-'`/`'+'`
    - `check_breakpoint_direction` 函数：`1`/`-1` → `'+'`/`'-'`
    - linkage mode（第 430-431 行）：`1`/`-1` → `'+'`/`'-'`
    - depth mode（第 517-518 行）：`1`/`-1` → `'+'`/`'-'`
- **Acceptance Criteria Addressed**: AC-1, AC-2
- **Test Requirements**:
  - `programmatic` TR-1.1: `breaks_bdg.to_dataframe()` 返回有数据的 DataFrame，列名为 `['chrom', 'start', 'end', 'name']`
  - `programmatic` TR-1.2: `df_breaks_bed` 不为空，包含正链和负链的 reads
  - `programmatic` TR-1.3: 运行 identify_wgls 不再出现 `KeyError: 'name'`
- **Notes**: 这是 KeyError 的根本原因，不是 to_dataframe() 列名问题

## [x] Task 2: 修复 bedtools merge 排序问题
- **Priority**: high
- **Depends On**: None
- **Description**: 
  - depth mode 中 `filtered_region_depth_bed.merge()` 要求输入文件已排序
  - 在 `merge()` 前添加 `.sort()` 调用
  - 修复位置：第 511-512 行
- **Acceptance Criteria Addressed**: AC-3
- **Test Requirements**:
  - `programmatic` TR-2.1: 运行 identify_wgls depth mode 不再出现 "Sorted input specified" 错误
- **Notes**: bedtools merge 默认假设输入已排序

## [x] Task 3: 验证修复
- **Priority**: high
- **Depends On**: Task 1, Task 2
- **Description**: 
  - 构建 Docker 镜像 quay.io/bioinfortools/cresil:1.2.1
  - 使用测试数据运行 identify_wgls
  - 确认 KeyError 和排序错误均已解决
- **Acceptance Criteria Addressed**: AC-1, AC-2, AC-3
- **Test Requirements**:
  - `programmatic` TR-3.1: Docker 镜像构建成功
  - `programmatic` TR-3.2: 运行 identify_wgls 不再出现 KeyError
  - `programmatic` TR-3.3: 运行 identify_wgls 不再出现排序错误
- **Notes**: 测试数据较小，可能不会走到完整流程，但 KeyError 应已解决

# Task Dependencies
- Task 1 和 Task 2 独立，可以并行执行
- Task 3 依赖于 Task 1 和 Task 2
