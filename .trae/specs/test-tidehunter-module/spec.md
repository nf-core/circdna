# slim 补全：REPEX_TAREAN(clu) 模块 + 未完成项完成 + 模块独立测试 Spec

## Why

1. eccsplorer_slim 尚无 clu（RepeatExplorer2 聚类）能力。之前 slim 流程设计（modularize-eccsplorer-modes）中 **Task 14（clu 标准骨架）、Task 15（测试）** 与 implement-eccsplorer-clu-all 的 **Task 5（clu 真实测试）** 均未完成
2. TIDEHUNTER 模块（ecc_finder 工具链 ONT 串联重复检测）CLI 实现错误（`-f ref -b bam -o bed` vs 官方 `TideHunter in.fa/fq > cons.fa`），从未测试
3. 其余 slim 设计项已实际完成：map 黑盒拆分（segemehl/haarz/peak/candidate 等模块）、eccsplorer_pipeline 重命名、circle_identifier 统一（含 eccsplorer_map/eccsplorer_clu）、冗余子工作流清理、slim Docker 镜像构建推送

## What Changes

- **新建 REPEATEXPLORER2 模块（bio.nf）**：`bio.nf/modules/repeatexplorer2/main.nf`，独立 RepeatExplorer2 镜像 `docker.1ms.run/kavonrtep/repeatexplorer:2.3.8`（从完整 eccsplorer 镜像拆分，仅含 RepeatExplorer2 + 内置 `/repex_tarean/databases/`）。命令：`/repex_tarean/seqclust --paired --prefix_length <N> --output_dir seqclust --taxon <T> --cpu <N> reads.fa --cleanup --keep_names --options ILLUMINA`
- **新建 clu_prepare**：提取 `eccPrepare.py prepare_for_clustering` 为 `clu_prepare.py`（读长统一 → 子采样 → 前缀拼接 → REPEATEXPLORER_READY.fa），加入 eccsplorer_slim 镜像 + `modules/local/eccsplorer_slim/clu_prepare/main.nf`
- **新建 clu_candidates**：提取 eccClusterer/eccComparer 候选逻辑为 `clu_candidates.py`（解析 seqclust 输出 → cluster_candidates.csv / comparative_cluster_table.csv），加入 slim 镜像 + `modules/local/eccsplorer_slim/clu_candidates/main.nf`
- **新建 `subworkflows/local/eccsplorer_clu_slim/main.nf`**：组合 clu_prepare → repex_tarean → clu_candidates；circle_identifier 新增 `eccsplorer_clu_slim`
- **修复 TIDEHUNTER 模块**：输入长读 reads（fasta/fastq），输出 consensus fasta，`TideHunter ${reads} > ${prefix}.fasta`
- **独立测试**：`tests/tidehunter/test.nf`（模拟串联重复长读）、`tests/eccsplorer_clu/test.nf`（circdna_1 treatment + gdna_1 control），输出到 `results/slim_run/tidehunter_test/` 与 `results/slim_run/eccsplorer_clu_test/`
- 回归验证现有 slim 流水线（-resume）

## Impact

- Affected specs: modularize-eccsplorer-modes、implement-eccsplorer-clu-all（Task 5）
- Affected code:
  - `bio.nf/modules/repeatexplorer2/main.nf`（新建，镜像 docker.1ms.run/kavonrtep/repeatexplorer:2.3.8）
  - `circdna.nf/modules/local/tidehunter/main.nf`（修复）
  - `circdna.nf/modules/local/eccsplorer_slim/clu_prepare/main.nf`（新建）
  - `circdna.nf/modules/local/eccsplorer_slim/clu_candidates/main.nf`（新建）
  - `circdna.nf/subworkflows/local/eccsplorer_clu_slim/main.nf`（新建）
  - `circdna.nf/workflows/circdna.nf`（circle_identifier 解析）
  - `circdna.nf/conf/modules.config`（新模块 publishDir）
  - `ECCsplorer/bin/clu_prepare.py`、`ECCsplorer/bin/clu_candidates.py`（新建）+ eccsplorer_slim 镜像重建
  - `circdna.nf/tests/tidehunter/test.nf`、`circdna.nf/tests/eccsplorer_clu/test.nf`（新建）

## ADDED Requirements

### Requirement: REPEATEXPLORER2 模块执行 RepeatExplorer2 聚类

系统 SHALL 提供独立 REPEATEXPLORER2 模块，直接调用 seqclust（RepeatExplorer2）完成串联/分散重复聚类，使用内置数据库。

#### Scenario: 成功案例

- **WHEN** 对 REPEATEXPLORER_READY.fa 运行 REPEATEXPLORER2
- **THEN** 生成 seqclust 聚类结果目录（clusters/superclusters/annotations）

### Requirement: slim clu 流程可独立运行

系统 SHALL 提供 eccsplorer_clu_slim 子工作流（prepare → 聚类 → 候选提取），可通过 circle_identifier 触发。

#### Scenario: 成功案例

- **WHEN** 运行 `tests/eccsplorer_clu/test.nf`（circdna_1 + gdna_1）
- **THEN** 生成 cluster_candidates.csv / comparative_cluster_table.csv 且非空

### Requirement: TIDEHUNTER 模块正确调用

系统 SHALL 以官方 CLI 运行 TideHunter，输入长读 reads，输出 consensus fasta。

#### Scenario: 成功案例

- **WHEN** 对模拟串联重复长读运行 TIDEHUNTER
- **THEN** 生成非空 `${prefix}.fasta` 含重复单元 consensus

## MODIFIED Requirements

### Requirement: 之前 slim 设计未完成项

modularize-eccsplorer-modes Task 14/15 与 implement-eccsplorer-clu-all Task 5 由本次实现并标记完成。

#### Scenario: 盘点确认

- **WHEN** 检查相关 spec tasks.md
- **THEN** clu 骨架、clu 测试、TIDEHUNTER 修复与测试均勾选完成
