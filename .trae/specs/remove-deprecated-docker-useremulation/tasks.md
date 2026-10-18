# 移除已废弃的 docker.userEmulation 配置 - Tasks

## [ ] Task 1: 移除 8 个主流程 nextflow.config 中的 userEmulation 行

- **Priority**: high
- **Depends On**: None
- **Description**: 逐一删除以下 8 个文件中的 `docker.userEmulation = true`（或 `userEmulation = true`）行，不修改其他任何内容。

  | # | 文件 | 行号 | 当前内容 |
  |---|------|------|----------|
  | 1 | `circdna.nf/nextflow.config` | 165 | `        docker.userEmulation   = true` |
  | 2 | `circrna.nf/nextflow.config` | 173 | `        docker.userEmulation    = true` |
  | 3 | `nanoseq.nf/nextflow.config` | 149 | `        docker.userEmulation   = true` |
  | 4 | `isoseq.nf/nextflow.config` | 137 | `        docker.userEmulation    = true` |
  | 5 | `fetchngs.nf/nextflow.config` | 103 | `        docker.userEmulation   = true` |
  | 6 | `riboseq.nf/nextflow.config` | 180 | `        docker.userEmulation    = true` |
  | 7 | `rnaseq/nextflow.config` | 193 | `        docker.userEmulation    = true` |
  | 8 | `bio.nf/nextflow.config` | 8 | `    userEmulation = true` |

  - **操作方式**：对每个文件使用 Edit 工具删除该行（含行尾换行符）
  - **注意**：`bio.nf/nextflow.config` 的缩进为 4 空格（顶层 `docker {}` 块），其余 7 个为 8 空格（`profiles { docker {} }` 块内）
- **Acceptance Criteria**: `grep -r 'userEmulation' circdna.nf/nextflow.config circrna.nf/nextflow.config nanoseq.nf/nextflow.config isoseq.nf/nextflow.config fetchngs.nf/nextflow.config riboseq.nf/nextflow.config rnaseq/nextflow.config bio.nf/nextflow.config` 无输出
- **Test**: TR-1.1 `grep -c 'userEmulation'` 对上述 8 个文件均返回 0

## [ ] Task 2: 移除 12 个 bio.nf 模块测试配置中的 userEmulation 行

- **Priority**: high
- **Depends On**: None（与 Task 1 无依赖，可并行）
- **Description**: 逐一删除以下 12 个文件中的 `    userEmulation = true` 行（4 空格缩进），不修改其他任何内容。

  | # | 文件 |
  |---|------|
  | 1 | `bio.nf/modules/flye/tests/nextflow.config` |
  | 2 | `bio.nf/modules/cresil/identify_wgls/tests/nextflow.config` |
  | 3 | `bio.nf/modules/cresil/trim/tests/nextflow.config` |
  | 4 | `bio.nf/modules/cresil/visualize/tests/nextflow.config` |
  | 5 | `bio.nf/modules/cresil/annotate/tests/nextflow.config` |
  | 6 | `bio.nf/modules/cresil/identify/tests/nextflow.config` |
  | 7 | `bio.nf/modules/eccsplorer/tests/nextflow.config` |
  | 8 | `bio.nf/modules/ecc_finder/map_ont/tests/nextflow.config` |
  | 9 | `bio.nf/modules/ecc_finder/map_sr/tests/nextflow.config` |
  | 10 | `bio.nf/modules/ecc_finder/asm_sr/tests/nextflow.config` |
  | 11 | `bio.nf/modules/ecc_finder/asm_ont/tests/nextflow.config` |
  | 12 | `bio.nf/modules/fastqdl/download/tests/nextflow.config` |

  - **操作方式**：对每个文件使用 Edit 工具删除 `    userEmulation = true` 行
  - **注意**：保留同文件中的 `runOptions` 和 `fixOwnership` 行不变
- **Acceptance Criteria**: `grep -r 'userEmulation' bio.nf/modules/*/tests/nextflow.config` 无输出
- **Test**: TR-2.1 `grep -c 'userEmulation'` 对上述 12 个文件均返回 0；TR-2.2 `grep -c 'runOptions\|fixOwnership'` 确认这些行仍在

## [ ] Task 3: 更新 AGENTS.md 第 11 节为 A+C 两层方案

- **Priority**: medium
- **Depends On**: Task 1, Task 2
- **Description**: 修改 `/Users/siyangming/nextflow_nf_core/AGENTS.md` 的 Section 11（Docker 用户映射 A+B+C 标准设置规范）：
  1. 将标题 `## 11. Docker 用户映射 A+B+C 标准设置规范` 改为 `## 11. Docker 用户映射 A+C 标准设置规范`
  2. 在 11.1 三层方案定义表格中，移除 B 层行，仅保留 A 层和 C 层
  3. 移除 11.1 中关于 "三者共存不冲突" 的说明，改为 "A+C 两层共存：A 在容器启动时即以宿主机用户身份运行，C 在任务结束后兜底修复（文件已正确归属时为 no-op）"
  4. 在 11.3 标准配置行中，移除 `docker.userEmulation = true` / `userEmulation = true` 行
  5. 在 11.5 检查清单中，移除 "同时设置 A+B+C" 改为 "同时设置 A+C"
  6. 在 11.4 禁止重复设置部分，移除对 B 的引用
  7. 在 Section 11 开头或 11.1 中添加简要说明："`docker.userEmulation` 在 Nextflow 26.04+ 已被移除，原 B 层不再需要"
- **Acceptance Criteria**: AGENTS.md Section 11 中不再出现 `userEmulation` 字样（除废弃说明外）
- **Test**: TR-3.1 阅读 Section 11 确认方案已从 A+B+C 修订为 A+C；TR-3.2 `grep 'userEmulation' AGENTS.md` 仅在废弃说明上下文中出现

## [ ] Task 4: 全局验证无残留 userEmulation

- **Priority**: medium
- **Depends On**: Task 1, Task 2, Task 3
- **Description**: 运行全局 grep 确认所有活跃配置文件中不再有 `userEmulation`（参考归档 `nf-core-nanoseq_3.1.0/` 除外）。
  - 验证命令：`grep -r 'userEmulation' --include='*.config' . | grep -v 'nf-core-nanoseq_3.1.0'` 应无输出
  - 同时确认 AGENTS.md 中的 `userEmulation` 仅出现在废弃说明上下文
- **Acceptance Criteria**: 活跃配置中零 `userEmulation` 残留
- **Test**: TR-4.1 grep 结果为空（排除 nf-core-nanoseq_3.1.0 参考归档）

## Task Dependencies

- Task 1 与 Task 2 无依赖，可并行执行
- Task 3 依赖 Task 1 + Task 2 完成
- Task 4 依赖 Task 1 + Task 2 + Task 3 完成
