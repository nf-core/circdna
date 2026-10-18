# Build Slim Images and Modules Spec

## Why

当前 ECCsplorer（~1.5GB）和 ecc_finder（~2GB）Docker 镜像体积大，因为它们打包了所有外部工具（bwa、segemehl、unicycler 等）。原子化拆分方案已将这些外部工具调用替换为 nf-core 标准模块。需要构建仅含 Python/R 脚本依赖的最小化 slim 镜像，同时保护原始源码的分支隔离。

## What Changes

- 在 ecc_finder 和 ECCsplorer 的 GitHub 仓库中创建 `slim` 分支，用于存放 slim 版本的源码修改
- master 分支保留 slime 前的原版 Docker/conda 配置
- 构建最小化 conda 包：`yangmingsi::eccsplorer_slim=1.0.0`、`yangmingsi::ecc_finder_slim=1.0.0`
- 构建最小化 Docker 镜像：`quay.io/bioinfortools/eccsplorer_slim:1.0.0`、`quay.io/bioinfortools/ecc_finder_slim:1.0.0`
- 从 ecc_finder 源码提取 map-sr.py 中的 merge_score 逻辑为独立脚本
- 从 ecc_finder 源码提取 asm-sr.py 中的 asm_filter 逻辑为独立脚本
- 从 ECCsplorer 源码提取 6 个专有分析步骤为独立 Python/R 脚本
- 在 bio.nf 中创建对应的 Nextflow 模块定义
- 完成 test_local.config 的 stub 测试验证

## Impact

- Affected specs: integrate-eccfinder-ampliconsuite, modularize-eccsplorer-modes
- Affected code:
  - `bio.nf/modules/eccsplorer_slim/*` (新增)
  - `bio.nf/modules/ecc_finder_slim/*` (新增)
  - `ECCsplorer/` (slim 分支)
  - `ecc_finder/` (slim 分支)
- Affected repos: `crimBubble/ECCsplorer` fork, `njaupanpan/ecc_finder` fork

## Branch Strategy

| 分支     | 内容                                                        | 用途            |
| -------- | ----------------------------------------------------------- | --------------- |
| `master` | 原版 ECCsplorer / ecc_finder 完整源码 + Docker + conda 配置 | 原版黑盒模式    |
| `slim`   | 精简后的 Python/R 脚本 + 最小化 Dockerfile + conda recipe   | slim 原子化模式 |

## ADDED Requirements

### Requirement: Slim Branch Isolation

The system SHALL maintain a separate `slim` branch in the ecc_finder and ECCsplorer GitHub repositories for slim-version source code modifications. The master branch SHALL retain unmodified pre-slim Docker/conda configurations.

#### Scenario: Create slim branch

- **WHEN** developer needs to modify ecc_finder or ECCsplorer source for slim version
- **THEN** changes go to a new `slim` branch created from current master
- **AND** master branch Dockerfile and conda-recipe remain unchanged

### Requirement: ECCsplorer Slim Conda Package

The system SHALL provide a minimal conda package `yangmingsi::eccsplorer_slim=1.0.0` containing only Python(numpy/scipy/biopython) + R(ggplot2) + 6 custom analysis scripts.

#### Scenario: Build conda package

- **WHEN** `conda build conda-recipe/` is run in ECCsplorer slim branch
- **THEN** a noarch conda package is produced
- **AND** the package can be tested with `conda create -n test --use-local eccsplorer_slim`

### Requirement: ECCsplorer Slim Docker Image

The system SHALL provide a minimal Docker image `quay.io/bioinfortools/eccsplorer_slim:1.0.0` (~400MB) based on condaforge/mambaforge with only Python + R dependencies.

#### Scenario: Build Docker image

- **WHEN** `docker build` is run with the slim Dockerfile
- **THEN** a minimal image is created containing only the slim conda environment + custom scripts
- **AND** the image excludes segemehl, haarz, samtools, bedtools, blast, trimmomatic

### Requirement: ecc_finder Slim Conda Package

The system SHALL provide a minimal conda package `yangmingsi::ecc_finder_slim=1.0.0` containing only Python(numpy/pandas/matplotlib/pybedtools) + bedtools + 2 custom analysis scripts.

### Requirement: ecc_finder Slim Docker Image

The system SHALL provide a minimal Docker image `quay.io/bioinfortools/ecc_finder_slim:1.0.0` (~250MB) based on condaforge/mambaforge.

### Requirement: Custom Script Extraction from ecc_finder

The system SHALL extract the merge_score logic from `map-sr.py` and the asm_filter logic from `asm-sr.py` into standalone Python scripts that accept explicit input/output file paths.

#### Scenario: Merge score script

- **WHEN** merge_score.py is called with enrichment BED, split-read BED, and reference FASTA
- **THEN** it outputs candidates.csv and candidates.fasta

### Requirement: Custom Script Extraction from ECCsplorer

The system SHALL extract 6 analysis steps from ECCsplorer source into standalone scripts:

- `peak_detect.py`: coverage.tsv → peak_regions.bed (scipy find_peaks)
- `candidate_extract.py`: SR.bed + peak_all.bed + peak_DR.bed → candidates.bed (bedtools intersect)
- `coverage_profile.py`: candidates.bed + all.bed + ref.fa → per-candidate coverage.tsv
- `normalize.R`: raw coverage + mapped bases → RPM + fold enrichment
- `visualize.R`: normalized data → PNG plots (manhattan + candidate)
- `html_report.py`: candidate data → HTML summary

### Requirement: bio.nf Module Definitions

The system SHALL create Nextflow module main.nf files under `bio.nf/modules/eccsplorer_slim/` and `bio.nf/modules/ecc_finder_slim/` for each custom script, referencing the respective slim Docker images.

### Requirement: Stub Test Verification

The system SHALL pass `nextflow run main.nf -profile test_local -stub` with circle_identifier set to slim mode.
