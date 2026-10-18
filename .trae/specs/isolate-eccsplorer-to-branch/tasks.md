# Tasks

- [x] Task 1: 盘点 `master` 中所有 ECCsplorer 归属内容

  - [x] SubTask 1.1: 列出代码层归属文件和目录
  - [x] SubTask 1.2: 列出配置、样本表、测试数据和脚本归属文件
  - [x] SubTask 1.3: 列出 `.trae/specs` 与文档中需跟随分支迁移的内容

- [x] Task 2: 定义 `master` 与 `ECCsplorer` 分支的职责边界

  - [x] SubTask 2.1: 定义 `master` 保留哪些能力
  - [x] SubTask 2.2: 定义 `ECCsplorer` 分支承载哪些能力
  - [x] SubTask 2.3: 明确今后 ECCsplorer 相关开发与测试都只在分支进行

- [x] Task 3: 设计 `master` 的回退方案

  - [x] SubTask 3.1: 定义如何将 `circdna.nf` 恢复到加入 ECCsplorer 之前的状态
  - [x] SubTask 3.2: 定义哪些文件应从主线移除
  - [x] SubTask 3.3: 定义如何避免误删非 ECCsplorer 功能

- [x] Task 4: 设计 `ECCsplorer` 分支的承接方案

  - [x] SubTask 4.1: 定义如何保留现有 ECCsplorer 模块、参数、测试与数据
  - [x] SubTask 4.2: 定义如何让 `modularize-eccsplorer-modes` 转移到该分支继续实施
  - [x] SubTask 4.3: 定义分支中的版本、文档和测试策略

- [x] Task 5: 设计分支切换与迁移顺序

  - [x] SubTask 5.1: 定义先保护现有工作区、再整理分支、最后回退主线的顺序
  - [x] SubTask 5.2: 定义如何处理当前脏工作区与未提交变更
  - [x] SubTask 5.3: 定义 branch/worktree 的推荐操作路径

- [x] Task 6: 设计测试与验收策略

  - [x] SubTask 6.1: 定义 `master` 的最小验证范围（本次至少执行 `git status --short --branch`、`git rev-parse HEAD` 与 `nextflow config -profile test_local`）
  - [x] SubTask 6.2: 定义 `ECCsplorer` 分支的最小验证范围（本次至少执行 `git status` 与 `nextflow config -profile test_local`）
  - [x] SubTask 6.3: 定义 `test_local`、`-resume` 和数据库准备在哪个分支执行

- [x] Task 7: 明确与现有 spec 的关系
  - [x] SubTask 7.1: 说明 `modularize-eccsplorer-modes` 在分支隔离后的实施位置
  - [x] SubTask 7.2: 说明 `design-ngs-analysis-and-resume-real-test` 中哪些内容转移到分支
  - [x] SubTask 7.3: 明确哪些旧任务在 `master` 上停止执行

# Task Dependencies

- Task 2 依赖 Task 1
- Task 3 依赖 Task 1、Task 2
- Task 4 依赖 Task 1、Task 2
- Task 5 依赖 Task 3、Task 4
- Task 6 依赖 Task 3、Task 4、Task 5
- Task 7 依赖 Task 2、Task 4

# 并行化说明

- Task 3 与 Task 4 可在职责边界明确后并行整理
- Task 6 与 Task 7 可在迁移顺序确定后并行补充
