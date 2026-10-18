# 迁移 bio.nf 模块至 modules 仓库计划

## 1. 摘要 (Summary)
本计划旨在将 `/Users/siyangming/nextflow_nf_core/bio.nf` 中的软件模块及子工作流（subworkflows），按照 nf-core 的标准目录结构迁移至 `/Users/siyangming/nextflow_nf_core/modules` 仓库中。迁移过程中将以“每个软件/子工作流”为单位，分别进行 `git commit` 提交，以保持代码提交历史的清晰和模块化。

## 2. 当前状态分析 (Current State Analysis)
- **源路径**：`/Users/siyangming/nextflow_nf_core/bio.nf`
  - 包含 `modules/`（如 `flair`, `kingfisher`, `longshot`, `minimap2`, `samtools`）
  - 包含 `subworkflows/`（如 `flair_pipeline`, `minimap2_longshot_flair`）
- **目标仓库**：`/Users/siyangming/nextflow_nf_core/modules`
  - 这是一个标准的 nf-core 模块仓库。
  - 模块的存放规范路径为：`modules/nf-core/<软件名>/`
  - 子工作流的存放规范路径为：`subworkflows/nf-core/<子工作流名>/`

## 3. 拟议更改 (Proposed Changes)

### 3.1 迁移 Modules（按软件分别提交）
将依次遍历源目录下的每个软件文件夹，复制到目标仓库并独立提交：
1. **flair**：复制 `bio.nf/modules/flair` -> `modules/modules/nf-core/flair`，执行 `git add` 并提交 `feat(modules): add flair modules`。
2. **kingfisher**：复制 `bio.nf/modules/kingfisher` -> `modules/modules/nf-core/kingfisher`，执行 `git add` 并提交 `feat(modules): add kingfisher modules`。
3. **longshot**：复制 `bio.nf/modules/longshot` -> `modules/modules/nf-core/longshot`，执行 `git add` 并提交 `feat(modules): add longshot module`。
4. **minimap2**：复制 `bio.nf/modules/minimap2` -> `modules/modules/nf-core/minimap2`，执行 `git add` 并提交 `feat(modules): add minimap2 modules`。
5. **samtools**：复制 `bio.nf/modules/samtools` -> `modules/modules/nf-core/samtools`，执行 `git add` 并提交 `feat(modules): add samtools modules`。

### 3.2 迁移 Subworkflows（按工作流分别提交）
将依次遍历源目录下的每个子工作流文件夹，复制到目标仓库并独立提交：
1. **flair_pipeline**：复制 `bio.nf/subworkflows/flair_pipeline` -> `modules/subworkflows/nf-core/flair_pipeline`，执行 `git add` 并提交 `feat(subworkflows): add flair_pipeline`。
2. **minimap2_longshot_flair**：复制 `bio.nf/subworkflows/minimap2_longshot_flair` -> `modules/subworkflows/nf-core/minimap2_longshot_flair`，执行 `git add` 并提交 `feat(subworkflows): add minimap2_longshot_flair`。

*(执行方案：可以编写一个简单的 bash 脚本自动化完成上述 cp 和 git 操作)*

## 4. 假设与决策 (Assumptions & Decisions)
- **范围假设**：除了软件 `modules` 外，同在 `bio.nf` 下的 `subworkflows` 也属于用户需要迁移和管理的模块资产，因此一并纳入迁移计划。
- **目标路径假设**：由于目标仓库是 nf-core 格式，模块统一存放在 `modules/nf-core/` 和 `subworkflows/nf-core/` 命名空间下。
- **提交策略决策**：所有提交都将直接在目标仓库的当前分支（通常为 `main` 或 `master`）上进行。如果后续需要推送到远程或提 PR，用户可自行基于当前分支进行推送或分支拆分。

## 5. 验证步骤 (Verification Steps)
- 迁移完成后，在目标仓库 `/Users/siyangming/nextflow_nf_core/modules` 中运行 `git log -n 7`，验证是否生成了按软件拆分的 commit 记录。
- 运行 `tree modules/nf-core` 和 `tree subworkflows/nf-core`（或通过 `ls`）确认文件已按预期路径正确存放。
