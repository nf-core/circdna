# 修复 ECCSPLORER 模块输出文件映射 Spec

## Why

`circdna.nf/modules/local/eccsplorer/main.nf` 的 script 块使用 `mv ${prefix}_output/*_candidates.bed` 和 `mv ${prefix}_output/*_junction_reads.txt` 重命名 ECCsplorer 输出，但 ECCsplorer 从不生成这两个文件名。实际输出文件位于 `${prefix}_output/eccpipe_results/mapping_results/` 下，文件名为 `TR_hiconf-ECC-REGIONS.bed`、`TR_lowconf-ECC-regions.bed`、`TR.trns.txt` 等。

由于 `mv` 的 glob 匹配不到任何文件，`|| touch` 兜底创建了 0 字节空文件，导致 `circdna_test_local_eccdna/eccsplorer/` 下 6 个文件全部为空。实际工作目录中 ECCsplorer 正常生成了 78 个高置信度 eccDNA 区域和 2388 个低置信度区域，但未被发布到输出目录。

## What Changes

- 修复 `circdna.nf/modules/local/eccsplorer/main.nf` 的 script 块：将 `*_candidates.bed` 和 `*_junction_reads.txt` 的 glob 替换为 ECCsplorer 真实输出文件名
- 新增输出声明：`eccpipe_results` 目录（完整 ECCsplorer 结果树）、`lowconf_ecc_regions`（低置信度区域）、`alignment_stats`（比对统计）、`ecc_sequences`（eccDNA 序列）
- 更新 `circdna.nf/modules/local/eccsplorer/meta.yml`：补充新增输出的描述
- 更新 `circdna.nf/conf/modules.config` 的 ECCSPLORER publishDir：发布所有新增输出
- 同步修复 `bio.nf/modules/eccsplorer/main.nf` 和 `bio.nf/modules/eccsplorer/meta.yml`（保持上游模块集合一致）
- 更新 `circdna.nf/CHANGELOG.md`：v4.2.1 版本下新增条目
- 版本 bump：v4.2.0 → v4.2.1（PATCH，bug fix）

### 文件名映射关系

| 模块输出名（保持不变） | ECCsplorer 实际源文件 | 说明 |
|---|---|---|
| `${prefix}_candidates.bed` | `*_hiconf-ECC-REGIONS.bed` | 高置信度 eccDNA 候选区域 |
| `${prefix}_junction_reads.txt` | `*.trns.txt` | junction/translocation reads |

### 新增输出

| 输出名 | 源文件/目录 | 说明 |
|---|---|---|
| `lowconf_ecc_regions` | `*_lowconf-ECC-regions.bed` | 低置信度 eccDNA 区域 |
| `alignment_stats` | `*_alignment-stats.txt` | segemehl 比对统计 |
| `ecc_sequences` | `*_ECC-SEQUENCES.fasta` | 提取的 eccDNA 序列 |
| `eccpipe_results` | `${prefix}_output/eccpipe_results/` | 完整结果树（含可视化、per-candidate 分析、HTML 报告） |

## Impact

- Affected specs: 无（本次为 bug fix，不改变模块接口契约）
- Affected code:
  - `circdna.nf/modules/local/eccsplorer/main.nf`（核心修复）
  - `circdna.nf/modules/local/eccsplorer/meta.yml`（补充输出描述）
  - `circdna.nf/conf/modules.config`（publishDir 更新）
  - `bio.nf/modules/eccsplorer/main.nf`（同步修复）
  - `bio.nf/modules/eccsplorer/meta.yml`（同步更新）
  - `circdna.nf/CHANGELOG.md`（版本记录）
  - `circdna.nf/nextflow.config`（版本号 bump）

## MODIFIED Requirements

### Requirement: ECCSPLORER 模块输出

ECCSPLORER 模块 SHALL 正确映射 ECCsplorer 工具的实际输出文件到模块声明的输出通道，并发布完整的分析结果到输出目录。

#### Scenario: 高置信度 eccDNA 区域输出
- **WHEN** ECCsplorer 在 mapping 模式下运行完成
- **THEN** `${prefix}_candidates.bed` SHALL 包含 `*_hiconf-ECC-REGIONS.bed` 的完整内容（非 0 字节空文件）
- **AND** 文件 SHALL 被发布到 `${params.outdir}/eccsplorer/`

#### Scenario: Junction reads 输出
- **WHEN** ECCsplorer 在 mapping 模式下运行完成
- **THEN** `${prefix}_junction_reads.txt` SHALL 包含 `*.trns.txt` 的完整内容（非 0 字节空文件）
- **AND** 文件 SHALL 被发布到 `${params.outdir}/eccsplorer/`

#### Scenario: 完整 eccpipe_results 输出
- **WHEN** ECCsplorer 在 mapping 模式下运行完成
- **THEN** 完整的 `eccpipe_results/` 目录树 SHALL 被发布到 `${params.outdir}/eccsplorer/${prefix}_eccpipe_results/`
- **AND** 目录 SHALL 包含 `mapping_results/`（含 `TR_hiconf-ECC-REGIONS.bed`、`TR_lowconf-ECC-regions.bed`、`TR.bed`、`TR.trns.txt`、`TR_alignment-stats.txt`、`TR_ECC-SEQUENCES.fasta`、`TR_chr_manhattan-plot.png` 等）
- **AND** 目录 SHALL 包含 `eccMap_summary.html`

#### Scenario: 低置信度区域输出
- **WHEN** ECCsplorer 在 mapping 模式下运行完成
- **THEN** `${prefix}_lowconf_ecc_regions.bed` SHALL 包含 `*_lowconf-ECC-regions.bed` 的完整内容
- **AND** 文件 SHALL 被发布到 `${params.outdir}/eccsplorer/`

#### Scenario: 无候选 eccDNA 时的兜底
- **WHEN** ECCsplorer 运行完成但未检测到高置信度 eccDNA 区域（`*_hiconf-ECC-REGIONS.bed` 为空或不存在）
- **THEN** `${prefix}_candidates.bed` SHALL 为 0 字节空文件（通过 `touch` 兜底）
- **AND** 流程 SHALL 正常完成，不报错

#### Scenario: stub 模式
- **WHEN** 流程以 stub 模式运行
- **THEN** SHALL 创建 0 字节的 `${prefix}_candidates.bed` 和 `${prefix}_junction_reads.txt` 占位文件
- **AND** SHALL 创建空的 `${prefix}_eccpipe_results` 目录占位
- **AND** SHALL 创建 0 字节的 `${prefix}_lowconf_ecc_regions.bed`、`${prefix}_alignment_stats.txt`、`${prefix}_ecc_sequences.fasta` 占位文件

### Requirement: bio.nf 上游模块同步

`bio.nf/modules/eccsplorer/main.nf` SHALL 与 `circdna.nf/modules/local/eccsplorer/main.nf` 保持同步，确保上游模块集合与下游流程使用相同的输出映射逻辑。

#### Scenario: bio.nf 模块同步
- **WHEN** circdna.nf 的 ECCSPLORER 模块输出映射被修复
- **THEN** bio.nf 的 ECCSPLORER 模块 SHALL 应用相同的修复
- **AND** bio.nf 的 `meta.yml` SHALL 同步更新输出描述
