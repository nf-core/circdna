# ecc_finder + ampliconsuite 接入 ECCDNA_MODE Spec

## Why

`ecc_finder`（MAP_SR / ASM_SR 模式）是 NGS 短读 eccDNA 检测工具，适合加入 `ECCDNA_MODE`。`ampliconsuite` 当前仅通过 legacy 模式运行，应迁移到新流程。三者均适合二代测序。

**circle_identifier 独立控制**：MAP_SR 和 ASM_SR 是两种不同的检测方式，各自独立开关。

**NGS 适用性**：

- `ECC_FINDER_MAP_SR`：是（split-read+discordant read，需 FASTQ+ref）
- `ECC_FINDER_ASM_SR`：是（de novo 组装，需 FASTQ）
- `AMPLICONSUITE`：是（WGS ecDNA 扩增子，需 sorted BAM）

## What Changes

- `ECC_FINDER_MAP_SR`（模块已存在）接入 ECCDNA_MODE，由 `ecc_finder_map_sr` 控制
- `ECC_FINDER_ASM_SR` 从 bio.nf 拷贝并接入，由 `ecc_finder_asm_sr` 控制
- `AMPLICONSUITE` 加入 ECCDNA_MODE，由 `ampliconsuite` 控制（保留 legacy 兼容）
- 新建分支 `integrate-eccfinder-ampliconsuite`
- test_local 验证

## Impact

- 新增: `modules/local/ecc_finder_asm_sr/`（从 bio.nf 拷贝）
- 已有: `modules/local/ecc_finder_map_sr/`（仅接入流程）
- 修改: `eccdna_mode/main.nf`、`workflows/circdna.nf`
- 修改: `conf/modules.config`、`conf/test_local.config`
- **NOT BREAKING**

## ADDED Requirements

### Requirement: ecc_finder_map_sr 通过 circle_identifier 独立控制

The system SHALL run ECC_FINDER_MAP_SR when `circle_identifier` 含 `ecc_finder_map_sr`.

**Scenario**: 输入 sorted BAM + FASTA → CSV + FASTA → `results/<outdir>/ecc_finder_map_sr/`

### Requirement: ecc_finder_asm_sr 通过 circle_identifier 独立控制

The system SHALL run ECC_FINDER_ASM_SR when `circle_identifier` 含 `ecc_finder_asm_sr`.

**Scenario**: 输入 paired FASTQ → FASTA → `results/<outdir>/ecc_finder_asm_sr/`

### Requirement: ampliconsuite 接入 ECCDNA_MODE

The system SHALL run AMPLICONSUITE when `circle_identifier` 含 `ampliconsuite`.

**Scenario**: 输入 sorted BAM + mosek + aa_data_repo → `results/<outdir>/ampliconsuite/`

### Requirement: test_local

**Scenario**: `-profile test_local,docker -resume --circle_identifier eccsplorer,ecc_finder_map_sr,ecc_finder_asm_sr,ampliconsuite` 通过
