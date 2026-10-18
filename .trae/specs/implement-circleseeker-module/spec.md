# CircleSeeker 模块构建与流程接入 Spec

## Why

circdna.nf 长读长（PacBio HiFi）eccDNA 检测目前已有 CRESIL、FLED、Flye、ECC_FINDER 四个引擎，缺少 CircleSeeker（HiFi 优化的综合 eccDNA 检测工具，支持 CtcReads-Caller 与 SplitReads-Caller 双证据检测）。CircleSeeker 已发布 bioconda 包（1.1.2）与 biocontainers 镜像（`quay.io/biocontainers/circleseeker:1.1.2--pyhdfd78af_0`），nf-core/modules 中尚无该模块，按通用规范第 12 章需在目标流程 `modules/local/` 中自建 nf-core 标准模块并接入。

## What Changes

- 新增 nf-core 标准自定义模块 `modules/local/circleseeker/`（main.nf + meta.yml + environment.yml + nextflow.config + tests/），支持 docker/conda 双引擎，含 stub 模式
- 新增 `bin/circleseeker_to_bed.py`：将 CircleSeeker 的 `*_merged_output.csv` 转换为带 read support 列的 BED（供现有 `LONG_READ_FILTERING` 复用）
- 新增子流程 `subworkflows/local/circleseeker_pipeline/main.nf`（+ meta.yml）：运行 CircleSeeker → 产出 BED → 接入现有长读候选过滤
- 在 `workflows/circdna.nf` 长读分支（`params.protocol in ["pacbio","ont"]`）新增 `circleseeker` 作为 `long_read_identifier` 选项，并复用 `LONG_READ_FILTERING` 做支持度/黑名单/重复序列过滤
- `conf/modules.config` 新增 CIRCLESEEKER 的 publishDir（`${params.outdir}/long_read/circleseeker/${meta.id}`）
- `conf/test_pacbio_lr.config`、`conf/test_nanopore_lr.config` 的 `long_read_identifier` 加入 `circleseeker` 并补充资源覆盖
- 版本 bump：4.1.0 → 4.2.0（MINOR，新功能），CHANGELOG 更新
- 测试：nf-test 模块测试（真实 + stub，docker 拉取 biocontainer 镜像在本地 testdata 上运行）

## Impact

- Affected specs: 长读长（TGS）检测能力
- Affected code:
  - `modules/local/circleseeker/*`（新增）
  - `bin/circleseeker_to_bed.py`（新增）
  - `subworkflows/local/circleseeker_pipeline/*`（新增）
  - `workflows/circdna.nf`（长读分支接线）
  - `conf/modules.config`、`conf/test_pacbio_lr.config`、`conf/test_nanopore_lr.config`
  - `nextflow.config`（版本号）、`CHANGELOG.md`
- 参考外部仓库：`/Users/siyangming/nextflow_nf_core/CircleSeeker`（使用方式与测试数据）

## ADDED Requirements

### Requirement: CircleSeeker nf-core 标准模块

系统 SHALL 提供符合 nf-core 模块规范的 `modules/local/circleseeker/` 模块：

- main.nf 定义 process `CIRCLESEEKER`，包含 `conda`/`container` 指令、`task.ext.args`、`task.ext.prefix`、stub 模式与 `versions.yml` emit
- container 指向 `quay.io/biocontainers/circleseeker:1.1.2--pyhdfd78af_0`（singularity 回退到 galaxy depot）；environment.yml 引用 `bioconda::circleseeker=1.1.2`
- 输入：`tuple val(meta), path(reads)`（FASTA 或 FASTQ，gz/plain 均可）+ `tuple val(meta2), path(fasta)`（参考基因组 FASTA）
- 脚本内将 gz 解压、FASTQ 转 FASTA（awk），调用 `circleseeker -i reads.fa -r ref.fa -o . -p <prefix> -t <cpus>`，并调用 `bin/circleseeker_to_bed.py` 生成 BED
- 输出 emit：`merged`（`<prefix>_merged_output.csv`）、`bed`（`<prefix>.circleseeker.bed`）、`summary`（`<prefix>_summary.txt`）、`report`（`<prefix>_report.html`）、`versions`

#### Scenario: 模块真实运行成功

- **WHEN** 以 docker 运行 nf-test 模块测试（本地 testdata：`test_reads.fa` + `test_ref.fa`）
- **THEN** 进程成功，产出 merged CSV / BED / summary / report / versions.yml，snapshot 匹配

#### Scenario: stub 模式

- **WHEN** 以 `-stub` 运行 nf-test
- **THEN** 进程成功并产出预期的空/占位输出文件

### Requirement: CircleSeeker 流水线接入

系统 SHALL 将 CircleSeeker 接入长读长检测分支：

- 新增 `subworkflows/local/circleseeker_pipeline/main.nf`，take `reads` + `genome_fasta`，emit `merged`、`bed`、`versions`
- `workflows/circdna.nf` 长读分支解析 `long_read_identifier` 中新增 `circleseeker`；启用时将子流程产出的 `bed` 输入 `LONG_READ_FILTERING`（复用现有支持度/黑名单/重复序列过滤），versions 汇入 MultiQC
- `conf/modules.config` 配置 CIRCLESEEKER publishDir 至 `${params.outdir}/long_read/circleseeker/${meta.id}`
- 资源标签使用 `process_high`；测试 profile 中覆盖为 2 CPU / 6 GB

#### Scenario: 长读流程接线正确

- **WHEN** 在 `long_read_identifier` 中包含 `circleseeker` 运行长读分支
- **THEN** CircleSeeker 子流程运行、BED 进入过滤流程、versions 汇入 MultiQC 汇总

#### Scenario: 未启用 circleseeker

- **WHEN** `long_read_identifier` 不包含 `circleseeker`
- **THEN** 不执行 CircleSeeker，现有行为不受影响

### Requirement: 测试与文档

系统 SHALL 完成测试与版本管理：

- `modules/local/circleseeker/tests/main.nf.test`（真实 + stub 两个用例）、`tests/nextflow.config`（docker 启用 + `--platform linux/amd64` + 资源覆盖）、`main.nf.test.snap`
- 测试数据使用模块自带 `testdata/`（`test_reads.fa`、`test_ref.fa`），不引用其他仓库绝对路径
- 更新 `conf/test_pacbio_lr.config` 与 `conf/test_nanopore_lr.config`（`long_read_identifier` 增加 circleseeker + withName 资源覆盖）
- 版本 4.1.0 → 4.2.0，CHANGELOG 增加对应条目（新模块 + 新依赖 circleseeker 1.1.2）

#### Scenario: nf-test 全通过

- **WHEN** 在 `circleseeker` 分支运行 `nf-test test modules/local/circleseeker/tests/main.nf.test`
- **THEN** 真实与 stub 用例均通过

## MODIFIED Requirements

### Requirement: 长读检测引擎集合

`long_read_identifier` 支持的引擎由 `cresil,fled,flye,eccfinder` 扩展为包含 `circleseeker`。默认值 `'cresil,fled,flye,eccfinder'` 保持不变（circleseeker 为显式启用），测试 profile 中显式加入。
