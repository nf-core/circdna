# 将探索性分析迁移至 Snakemake 工作流 Spec

## Why

circdna.nf 当前的 `CANDIDATE_MERGE` 和 `ECC_SCORE` 两个模块是参数敏感的轻量 Python 步骤（process_single，秒级重算），但被嵌入到包含 `CIRCLEMAP_REALIGN`（process_high, 96h）等重计算步骤的 Nextflow 流程中。这导致：

1. 调整 `--max-distance`、`w1/w2/w3` 等探索性参数时，无法避免重新触发上游重计算
2. 参数网格搜索难以在 Nextflow 中实现
3. 探索性可视化（深度分布、候选统计、工具一致性）在 Nextflow 中缺失

将这两个步骤迁移至独立的 Snakemake 工作流（eccdna.smk），可以解耦"重计算"与"轻探索"，让 Nextflow 专注标准化批处理，Snakemake 专注参数调优与可视化。

## What Changes

### Nextflow 侧（circdna.nf）— **BREAKING**

- **BREAKING**: 移除 `modules/local/candidate_merge/` 目录及 CANDIDATE_MERGE process
- **BREAKING**: 移除 `modules/local/ecc_score/` 目录及 ECC_SCORE process
- **BREAKING**: 移除 `bin/merge_candidates.py`（迁移至 eccdna.smk 仓库）
- **BREAKING**: 移除 `bin/calculate_ecc_score.py`（迁移至 eccdna.smk 仓库）
- **BREAKING**: 移除 `subworkflows/local/integrated_mode/` 目录（ECC_SCORE 是其唯一实质步骤）
- 修改 `subworkflows/local/eccdna_mode/main.nf`：移除 CANDIDATE_MERGE 调用，eccdna 模式产出原始 `eccsplorer_bed` + `circle_map_bed` 后即终止
- 修改 `workflows/circdna.nf`：移除 INTEGRATED_MODE 调用及 integrated 分支（第 295-328 行）；integrated 模式不再可用，仅保留 reference + eccdna 两种模式
- 修改 `nextflow.config`：移除 `ecc_score_w1/w2/w3` 参数；移除 `mode` 参数中的 `integrated` 选项；版本 3.2.1 → 4.0.0
- 修改 `nextflow_schema.json`：通过 `nf-core schema build` 重新生成
- 修改 `conf/test_integrated.config`：标记为废弃或删除（integrated 模式已移除）
- 修改 `CHANGELOG.md`：新增 v4.0.0 条目，标注 BREAKING CHANGES
- 修改 `conf/modules.config`：移除 CANDIDATE_MERGE/ECC_SCORE 相关配置（如有）

### Snakemake 侧（eccdna.smk 新仓库）

- 新建 GitHub 仓库 `SiYangming/eccdna.smk`
- 创建 `Snakefile` 包含以下 rules：
  - `rule all`: 聚合最终目标
  - `rule candidate_merge`: 调用 `merge_candidates.py` 合并 ECCsplorer + Circle-Map 候选
  - `rule ecc_score`: 调用 `calculate_ecc_score.py` 计算综合评分
- 复制 `bin/merge_candidates.py` 和 `bin/calculate_ecc_score.py` 至 `eccdna.smk/scripts/`
- 创建 `config/config.yaml` 配置参数（w1/w2/w3, max_distance, 路径）
- 创建 `config/samples.yaml` 样本配置（指向 Nextflow 产物路径）
- 创建 `README.md` 说明使用方法
- 创建 `environment.yml` 定义 Python 依赖
- 初始化 git 仓库并推送至 GitHub

### 契约接口（Nextflow 产物 → Snakemake 输入）

| Nextflow 产物          | 路径模式                                           | Snakemake 用途       |
| ---------------------- | -------------------------------------------------- | -------------------- |
| mosdepth 区域深度      | `{outdir}/mosdepth/{sample}/*.regions.bed.gz`      | ECC_SCORE 输入       |
| mosdepth 深度分布      | `{outdir}/mosdepth/{sample}/*.global.dist.txt`     | 深度可视化           |
| Circle-Map realign bed | `{outdir}/circlemap/realign/*.bed`                 | candidate_merge 输入 |
| ECCsplorer bed         | `{outdir}/eccdna_mode/eccsplorer/*_candidates.bed` | candidate_merge 输入 |

## Impact

- **Affected specs**: 无（首个针对此迁移的 spec）
- **Affected code**:
  - `workflows/circdna.nf`（移除 integrated 分支）
  - `subworkflows/local/eccdna_mode/main.nf`（移除 CANDIDATE_MERGE）
  - `subworkflows/local/integrated_mode/main.nf`（删除）
  - `modules/local/candidate_merge/`（删除）
  - `modules/local/ecc_score/`（删除）
  - `bin/merge_candidates.py`（迁移）
  - `bin/calculate_ecc_score.py`（迁移）
  - `nextflow.config`（版本 + 参数）
  - `nextflow_schema.json`（重新生成）
  - `CHANGELOG.md`（v4.0.0 条目）
  - `conf/test_integrated.config`（废弃）
  - `conf/modules.config`（清理）
- **Affected modes**:
  - `reference` 模式：无变化
  - `eccdna` 模式：产出原始 bed 后终止，不再产出 merged_bed
  - `integrated` 模式：**BREAKING** 完全移除，由 Snakemake 工作流替代

## ADDED Requirements

### Requirement: Snakemake eccdna.smk 工作流

系统 SHALL 提供独立的 Snakemake 工作流仓库 `eccdna.smk`，接管原 Nextflow 中的 CANDIDATE_MERGE 和 ECC_SCORE 两个探索性分析步骤。

#### Scenario: Snakemake 接管候选合并

- **WHEN** 用户提供 Nextflow 产出的 ECCsplorer bed 和 Circle-Map bed 路径
- **THEN** Snakemake 通过 `rule candidate_merge` 调用 `scripts/merge_candidates.py` 生成 merged_candidates.bed
- **AND** 用户可通过 `config/config.yaml` 调整 `max_distance` 参数

#### Scenario: Snakemake 接管 ECC_SCORE 评分

- **WHEN** 用户提供 merged_candidates.bed + gDNA mosdepth bed + eccDNA mosdepth bed 路径
- **THEN** Snakemake 通过 `rule ecc_score` 调用 `scripts/calculate_ecc_score.py` 生成 scored.bed
- **AND** 用户可通过 `config/config.yaml` 调整 `w1/w2/w3` 权重

#### Scenario: 参数网格搜索（可选，本次不实现）

- **WHEN** 用户在 `config/config.yaml` 中配置 `grid_search` 启用
- **THEN** Snakemake 通过 `expand` 生成多组参数组合的 scored.bed
- **AND** 输出对比表汇总不同参数组合下的分级结果

注：参数网格搜索为可选增强项，本次迁移仅实现单组参数支持，网格搜索留作后续增强。

## MODIFIED Requirements

### Requirement: circdna.nf eccdna 模式

eccdna 模式 SHALL 在产出 ECCsplorer bed 和 Circle-Map bed 后终止，不再执行候选合并。原 `merged_bed` 输出通道 SHALL 移除。下游合并分析由独立的 Snakemake 工作流（eccdna.smk）接管。

### Requirement: circdna.nf 模式选择

`params.mode` SHALL 仅接受 `reference` 和 `eccdna` 两个值。`integrated` 模式 SHALL 被移除。原 `integrated` 分支（workflows/circdna.nf 第 295-328 行）SHALL 删除。

## REMOVED Requirements

### Requirement: circdna.nf integrated 模式

**Reason**: ECC_SCORE 评分步骤迁移至 Snakemake 工作流，integrated 模式在 Nextflow 中失去实质内容
**Migration**: 用户应使用 eccdna.smk 仓库的 Snakemake 工作流执行综合评分分析。Nextflow 仅负责产出 gDNA mosdepth + eccDNA mosdepth + 原始候选 bed，作为 Snakemake 输入

### Requirement: circdna.nf CANDIDATE_MERGE 模块

**Reason**: 候选合并步骤迁移至 Snakemake 工作流，便于参数调优
**Migration**: 用户应使用 eccdna.smk 的 `rule candidate_merge`，通过 `config/config.yaml` 配置 `max_distance` 参数

### Requirement: circdna.nf ECC_SCORE 模块

**Reason**: ECC_SCORE 评分步骤迁移至 Snakemake 工作流，便于权重调优和参数扫描
**Migration**: 用户应使用 eccdna.smk 的 `rule ecc_score`，通过 `config/config.yaml` 配置 `w1/w2/w3` 权重

### Requirement: circdna.nf ecc_score_w1/w2/w3 参数

**Reason**: 这些参数仅被 ECC_SCORE 使用，ECC_SCORE 迁移后参数失去用途
**Migration**: 用户应在 eccdna.smk 的 `config/config.yaml` 中配置这些权重
