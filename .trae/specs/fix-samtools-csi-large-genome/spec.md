# 修复 SAMTOOLS_INDEX BAI 上限与大基因组运行 Spec

## Why

Tragopogon_porrifolius hap1 参考序列 OZ388095 长度 608,949,837 bp，超过 BAI 索引坐标上限 2^29-1 (536,870,911)。`SAMTOOLS_INDEX_BAM` 默认生成 BAI 索引时以 `Numerical result out of range` 失败（exit code 1），触发 Execution cancelled。`conf/large_genome.config` 已配置 `ext.args = '-c'`（CSI 索引），但其 `withName: 'SAMTOOLS_INDEX'` 精确匹配仅覆盖无 alias 的实例（BAM_MARKDUPLICATES_PICARD 内部），无法匹配 alias 后的 `SAMTOOLS_INDEX_BAM`/`SAMTOOLS_INDEX_FILTERED`/`SAMTOOLS_INDEX_RE`，导致加了 `-c` 配置后大基因组索引仍走 BAI 而失败。

## What Changes

- **核心修复（仅影响大基因组，全部收口至 `conf/large_genome.config`）**：将 `withName: 'SAMTOOLS_INDEX'` 改为正则 `'.*SAMTOOLS_INDEX.*'`，使 `-c`（CSI 索引）覆盖 circdna.nf 全部 SAMTOOLS_INDEX 实例（SAMTOOLS_INDEX_BAM、BAM_MARKDUPLICATES_PICARD:SAMTOOLS_INDEX、SAMTOOLS_INDEX_FILTERED、SAMTOOLS_INDEX_RE）
- **不修改 `nextflow.config`**：本次变更仅作用于显式附加 `-c large_genome.config` 的大基因组运行，对所有流程默认行为无影响，故不做 manifest 版本 bump
- `SERVER_RUN_GUIDE.md` 中标注 Tragopogon_porrifolius hap1 为大基因组，命令附加 `-c circdna.nf/conf/large_genome.config`（hap2 保持不变）
- 更新 `README.md`：在 Usage 部分说明大基因组物种需附加 `-c conf/large_genome.config`（启用 CSI 索引）
- `CHANGELOG.md` 在现有 v4.4.1 条目下追加本次配置修复记录（不 bump 版本）
- **Git 提交范围限定**：提交（push 至 GitHub）时仅暂存本次大基因组相关文件（`conf/large_genome.config`、`SERVER_RUN_GUIDE.md`、`README.md`、`CHANGELOG.md`），工作区中其他无关更改一律不提交

## Impact

- Affected specs: 无（新增）
- Affected code: `conf/large_genome.config`、`SERVER_RUN_GUIDE.md`、`README.md`、`CHANGELOG.md`
- **不影响**：`nextflow.config`（manifest version 不变，所有流程默认行为不变）

## ADDED Requirements

### Requirement: large_genome.config 覆盖全部 SAMTOOLS_INDEX 实例

`conf/large_genome.config` SHALL 使用正则匹配使 `ext.args = '-c'`（CSI 索引）对 circdna.nf 中所有 SAMTOOLS_INDEX 进程生效，包括 alias 名称实例（`SAMTOOLS_INDEX_BAM`、`SAMTOOLS_INDEX_FILTERED`、`SAMTOOLS_INDEX_RE`）与无 alias 实例（`BAM_MARKDUPLICATES_PICARD:SAMTOOLS_INDEX`）。

#### Scenario: 大基因组样本索引成功

- **WHEN** 用户以 `-profile server -c circdna.nf/conf/large_genome.config` 运行参考序列超过 536,870,911 bp 的物种（如 Tragopogon_porrifolius hap1）
- **THEN** 全部索引步骤 SHALL 生成 `.csi` 索引而非失败
- **AND** pipeline 不再以 `samtools index: failed to create index ... out of range` 退出

#### Scenario: 小基因组样本不受影响

- **WHEN** 用户运行参考序列较短的物种（不附加 large_genome.config）
- **THEN** 索引仍 SHALL 正常生成（`.bai`，与原行为一致，nextflow.config 未改动）

### Requirement: Tragopogon_porrifolius hap1 标注大基因组

`SERVER_RUN_GUIDE.md` SHALL 将 Tragopogon_porrifolius hap1 标注为大基因组，命令附加 `-c circdna.nf/conf/large_genome.config`；hap2 命令保持不变。

#### Scenario: 用户复制 hap1 命令

- **WHEN** 用户从指南复制 Tragopogon_porrifolius hap1 运行命令
- **THEN** 命令 SHALL 包含 `-c circdna.nf/conf/large_genome.config`
- **AND** hap2 命令不附加该配置

### Requirement: README 说明大基因组 CSI 用法

`README.md` SHALL 在运行说明中提示：参考序列超过 BAI 上限（约 512 Mb/染色体）的大基因组物种需附加 `-c conf/large_genome.config` 以启用 CSI 索引。

### Requirement: Git 提交范围限定

本次变更提交（同步至 GitHub）SHALL 仅包含大基因组相关文件，工作区中的其他无关更改 SHALL NOT 被暂存或提交。

#### Scenario: 提交本次修复

- **WHEN** 用户提交并推送本次大基因组修复
- **THEN** 暂存区 SHALL 仅包含 `conf/large_genome.config`、`SERVER_RUN_GUIDE.md`、`README.md`、`CHANGELOG.md`
- **AND** 其他未相关文件的更改 SHALL 保持未提交状态

## MODIFIED Requirements

无

## REMOVED Requirements

无
