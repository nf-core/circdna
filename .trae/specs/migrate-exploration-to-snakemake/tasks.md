# Tasks

## 阶段一：创建 eccdna.smk Snakemake 仓库

- [x] Task 1: 创建 eccdna.smk 项目目录结构与基础文件

  - [x] SubTask 1.1: 在 `/Users/siyangming/nextflow_nf_core/eccdna.smk/` 创建项目结构（config/, scripts/, results/, logs/）
  - [x] SubTask 1.2: 复制 `circdna.nf/bin/merge_candidates.py` → `eccdna.smk/scripts/merge_candidates.py`
  - [x] SubTask 1.3: 复制 `circdna.nf/bin/calculate_ecc_score.py` → `eccdna.smk/scripts/calculate_ecc_score.py`
  - [x] SubTask 1.4: 创建 `config/config.yaml`（含 w1/w2/w3, max_distance, nextflow_outdir 等参数）
  - [x] SubTask 1.5: 创建 `config/samples.yaml`（样本配置模板，指向 Nextflow 产物路径）
  - [x] SubTask 1.6: 创建 `environment.yml`（Python 3.12 依赖）
  - [x] SubTask 1.7: 创建 `README.md`（使用说明、契约接口、参数说明）
  - [x] SubTask 1.8: 创建 `.gitignore`

- [x] Task 2: 编写 Snakemake 工作流 Snakefile
  - [x] SubTask 2.1: 创建 `eccdna.smk/Snakefile`，定义 `rule all`
  - [x] SubTask 2.2: 实现 `rule candidate_merge`，调用 scripts/merge_candidates.py
  - [x] SubTask 2.3: 实现 `rule ecc_score`，调用 scripts/calculate_ecc_score.py
  - [x] SubTask 2.4: 在 Snakefile 中加载 config/config.yaml 和 config/samples.yaml
  - [x] SubTask 2.5: 验证 Snakemake 工作流语法（--lint 或 --dry-run）

## 阶段二：迁移脚本至 eccdna.smk 并初始化 GitHub 仓库

- [~] Task 3: 初始化 eccdna.smk Git 仓库并推送至 GitHub
  - [x] SubTask 3.1: 在 `/Users/siyangming/nextflow_nf_core/eccdna.smk/` 执行 `git init`
  - [x] SubTask 3.2: 设置远程仓库 URL `origin → https://github.com/SiYangming/eccdna.smk.git`（仓库需用户手动创建或授权 MCP 后创建）
  - [x] SubTask 3.3: 执行 `git add . && git commit -m "feat: initial commit with candidate_merge and ecc_score rules"`
  - [ ] SubTask 3.4: 执行 `git push -u origin main`（阻塞：GitHub 仓库尚未创建，需用户手动创建或重新授权 GitHub MCP）

## 阶段三：从 circdna.nf 移除迁移模块

- [x] Task 4: 删除迁移模块文件

  - [x] SubTask 4.1: 删除 `circdna.nf/modules/local/candidate_merge/` 目录（含 main.nf, meta.yml）
  - [x] SubTask 4.2: 删除 `circdna.nf/modules/local/ecc_score/` 目录（含 main.nf, meta.yml）
  - [x] SubTask 4.3: 删除 `circdna.nf/subworkflows/local/integrated_mode/` 目录（含 main.nf, meta.yml 如有）
  - [x] SubTask 4.4: 删除 `circdna.nf/bin/merge_candidates.py`
  - [x] SubTask 4.5: 删除 `circdna.nf/bin/calculate_ecc_score.py`

- [x] Task 5: 修改 Nextflow 工作流文件

  - [x] SubTask 5.1: 修改 `subworkflows/local/eccdna_mode/main.nf`：移除 `include { CANDIDATE_MERGE }` 行；移除 `ch_merge_input` 和 `CANDIDATE_MERGE(ch_merge_input)` 调用；移除 `merged_bed` emit；保留 `eccsplorer_bed` 和 `circle_map_bed` emit
  - [x] SubTask 5.2: 修改 `workflows/circdna.nf`：移除 `INTEGRATED_MODE` include（如有）；移除 `else if (params.mode == 'integrated')` 整个分支（第 295-328 行）；移除 `INTEGRATED_MODE(...)` 调用；移除 `integrated` 从 `valid_modes` 列表（第 32-35 行）；移除 integrated 相关的 ch_mosdepth_multiqc 处理
  - [x] SubTask 5.3: 检查并修改 `workflows/circdna.nf` 顶部 include 语句，移除 INTEGRATED_MODE 引用

- [x] Task 6: 修改 Nextflow 配置文件
  - [x] SubTask 6.1: 修改 `nextflow.config`：移除 `ecc_score_w1/w2/w3` 三个参数定义；移除 `mode` 参数注释中的 'integrated' 选项；更新 manifest version 从 '3.2.1' 为 '4.0.0'
  - [x] SubTask 6.2: 修改 `conf/test_integrated.config`：删除该文件（integrated 模式已移除）
  - [x] SubTask 6.3: 检查 `conf/modules.config` 是否有 CANDIDATE_MERGE/ECC_SCORE 配置，如有则移除
  - [x] SubTask 6.4: 检查 `conf/test_local.config` 是否有 ecc_score_w1/w2/w3 配置，如有则移除
  - [x] SubTask 6.5: 检查 `conf/test_integrated.config` 是否在 `nextflow.config` profiles 中被引用，如有则移除引用

## 阶段四：更新版本号、CHANGELOG 与 Schema

- [x] Task 7: 更新 circdna.nf 版本号与 CHANGELOG
  - [x] SubTask 7.1: 在 `CHANGELOG.md` 顶部新增 `## v4.0.0 - [2026-08-02]` 条目，标注 BREAKING CHANGES，描述 integrated 模式移除、CANDIDATE_MERGE/ECC_SCORE 迁移至 eccdna.smk
  - [x] SubTask 7.2: 确认 `nextflow.config` 第 318 行 `version = '4.0.0'`（在 Task 6.1 中完成）
  - [x] SubTask 7.3: 运行 `conda activate nextflow && nf-core schema build` 重新生成 `nextflow_schema.json`（移除 ecc_score_w1/w2/w3 与 integrated mode 选项）
  - [x] SubTask 7.4: 检查 `modules.json` 是否引用 candidate_merge/ecc_score，如有则移除

## 阶段五：验证与提交 circdna.nf 至 GitHub

- [x] Task 8: 本地验证 circdna.nf 修改

  - [x] SubTask 8.3: 确认 `nextflow run main.nf -profile test_local --mode integrated` 会因 mode 验证失败而退出（预期行为）
  - [x] SubTask 8.1: 执行 `nextflow run main.nf -profile test_local,docker --mode reference -stub-run` 验证 reference 模式（成功启动，路由正确）
  - [x] SubTask 8.2: 执行 `nextflow run main.nf -profile test_local,docker --mode eccdna -stub-run` 验证 eccdna 模式（成功启动，无 CANDIDATE_MERGE/ECC_SCORE 步骤）
  - [x] SubTask 8.4: `nf-core lint` 跳过（nf-core 3.5.2 在 Python 3.14 下有 bug）

- [x] Task 9: 提交 circdna.nf v4.0.0 至 GitHub
  - [x] SubTask 9.1: 在 `circdna.nf/` 执行 `git status` 确认未提交修改范围
  - [x] SubTask 9.2: 执行 `git add -A` 暂存所有修改（排除 scripts/nextflow.log）
  - [x] SubTask 9.3: 执行 `git commit -m "release: v4.0.0 — migrate CANDIDATE_MERGE and ECC_SCORE to eccdna.smk"` 提交（commit 7c01758）
  - [x] SubTask 9.4: 执行 `git tag -a v4.0.0 -m "Release v4.0.0 — BREAKING: integrated mode removed, exploration migrated to eccdna.smk"` 打标签
  - [x] SubTask 9.5: 执行 `git push origin master --tags` 推送至 GitHub（成功）

# Task Dependencies

- Task 2 依赖 Task 1（需要项目结构才能编写 Snakefile）
- Task 3 依赖 Task 1 + Task 2（需要完整项目才能初始化 git）
- Task 5 + Task 6 可并行（修改不同文件）
- Task 7 依赖 Task 5 + Task 6（schema 重新生成需要配置已更新）
- Task 8 依赖 Task 5 + Task 6 + Task 7
- Task 9 依赖 Task 8

# 并行化说明

- 阶段一（Task 1-3）与阶段三（Task 4-6）原则上可并行，但为避免脚本文件在迁移过程中丢失，建议先完成阶段一（脚本已复制到新仓库），再执行阶段三（从 Nextflow 删除）
- 阶段三中的 Task 4（删除文件）与 Task 5/6（修改文件）可并行执行
