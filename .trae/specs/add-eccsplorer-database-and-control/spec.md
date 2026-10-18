# 为 ECCsplorer 添加注释数据库与 gDNA control 支持 Spec

## Why

当前 ECCSPLORER 模块未启用 BLAST 注释（`-d` 参数未传入），导致所有候选 eccDNA 的 `*_blast.m6` 为 0 字节空文件，无法识别候选序列的重复元件来源。同时 ECCsplorer 支持传入 gDNA（genomic DNA）对照组进行背景去除，能显著降低假阳性，但当前流程未启用该能力。用户已在 `/Users/siyangming/Library/CloudStorage/OneDrive-个人/Project_backup/PublicDB/` 下载了 RepBase 31.07 和 Dfam（RepeatMasker 格式）注释库，需要将这些库接入流程并测试 gDNA control 对结果的影响。

现有 samplesheets 目录下的 CSV 已使用 `data_type` 列（列名带下划线，值为 camelCase 的 `eccDNA` 或 `gDNA`），见 `circdna_ngs_clean.csv`、`circdna_Arabidopsis_thaliana_eccDNA.csv`。本 spec 严格适配该现有格式：**如果 samplesheet 中包含 `data_type=gDNA` 的行，则自动将 gDNA 样本作为 ECCsplorer 的 control（dataset B）启用；如果只有 `data_type=eccDNA` 的行，则不启用 control**。流程不新增任何命令行 control 参数，control 完全由 samplesheet 驱动。

为建立 eccDNA 样本与 gDNA control 的对应关系，**新增可选的 `group` 分组列**：同一 `group` 值下的 eccDNA 样本与 gDNA 样本配对（如 `group=A` 下有 1 个 eccDNA 样本 `circdna_1` 和 1 个 gDNA 样本 `gdna_1`，则 `circdna_1` 使用 `gdna_1` 作为 control）。不再依赖 sample 名的数字后缀匹配。

## What Changes

### 1. 注释数据库参数支持

- 新增 `params.eccsplorer_database` 参数（默认 null），指向未压缩的 FASTA 注释库
- 在 `conf/modules.config` 的 ECCSPLORER `ext.args` 中条件性传递 `-d ${params.eccsplorer_database}`
- 在 `conf/test_local.config` 与 `conf/test_local_gdna.config` 中设置默认测试库路径（合并后的 Dfam + RepBase 物种库）

### 2. gDNA control 支持（完全由 samplesheet `data_type` + `group` 列驱动）

- **不新增** `eccsplorer_control_fastq_1` / `eccsplorer_control_fastq_2` 等命令行参数
- 适配现有 samplesheet 的 `data_type` 列（列名带下划线，值为 `eccDNA` 或 `gDNA`）
- **新增可选的 `group` 分组列**：用于建立 eccDNA 样本与 gDNA control 的配对关系
- 流程从 samplesheet 中分离 `eccDNA` 样本（dataset A）和 `gDNA` 样本（dataset B control）
- **匹配规则**：按 `group` 列值匹配，同一 group 下的 eccDNA 样本使用该 group 下的 gDNA 样本作为 control
  - 一个 group 下有 1 个 gDNA + N 个 eccDNA → 这 N 个 eccDNA 共用该 gDNA control
  - 一个 group 下无 gDNA → 该 group 的 eccDNA 不启用 control
  - 一个 group 下有多个 gDNA → 报错（一个 group 仅允许 1 个 gDNA control）
  - 无 `group` 列但 samplesheet 含 gDNA：若 gDNA 仅 1 个，所有 eccDNA 共用该 control（向后兼容）；若 gDNA 多于 1 个，报错要求用户提供 `group` 列
- ECCSPLORER 模块新增可选 control 输入通道，script 块条件性附加 control 位置参数（ECCsplorer 的位置参数 3、4）

### 3. `check_samplesheet.py` 适配现有 `data_type` 列格式 + 新增 `group` 列

- 现状：`check_samplesheet.py` 仅识别 `datatype`（无下划线）列名，且值仅接受小写 `gdna`/`eccdna`
- 修改：同时识别 `data_type` 和 `datatype` 两种列名；值归一化为小写后校验，接受 `eccDNA`/`gDNA`/`eccdna`/`gdna` 任一写法
- 新增：识别可选的 `group` 列（任意非空字符串值），输出 samplesheet 保留该列
- 输出 samplesheet 统一保留原始列名（若输入是 `data_type` 则输出也是 `data_type`）
- meta map 中以 `data_type` 键存储归一化后的小写值（`eccdna`/`gdna`），以 `group` 键存储分组值，供 subworkflow 判断

### 4. gDNA 测试数据

- 现状：`testdatasets/testdata/` 已有 `gdna_1_R1.fastq.gz`、`gdna_1_R2.fastq.gz`，但缺少 `gdna_2`、`gdna_3`
- 复制 `circdna_2` 的 R1/R2 → `gdna_2_R1.fastq.gz` / `gdna_2_R2.fastq.gz`
- 复制 `circdna_3` 的 R1/R2 → `gdna_3_R1.fastq.gz` / `gdna_3_R2.fastq.gz`
- 现有 `samplesheets/samplesheet_local.csv` 补充 `data_type` 列（3 个样本均为 `eccDNA`，无 gDNA，无 group 列，用于无 control 对比测试）
- 新增 `samplesheets/samplesheet_local_with_gdna.csv`（3 个 eccDNA 样本 + 3 个 gDNA 样本，`data_type` 列区分，`group` 列建立配对：circdna_1↔gdna_1 同属 group=A，circdna_2↔gdna_2 同属 group=B，circdna_3↔gdna_3 同属 group=C）

### 5. 注释数据库预处理脚本

- 新增 `scripts/prepare_eccsplorer_database.sh`：将 Dfam RepeatMasker.lib.gz 解压 + 与 RepBase 对应物种库合并，生成单一的 `eccsplorer_db.fa`

### 6. 测试配置

- 新增 `conf/test_local_gdna.config`：eccdna 模式 + gDNA control（由 samplesheet `data_type` + `group` 列驱动）+ 注释数据库，输出到 `circdna.nf/results/test_local_gdna/`
- 修改 `conf/test_local.config`：samplesheet 补充 `data_type=eccDNA` 列，新增注释数据库默认路径，输出到 `circdna.nf/results/test_local/`

### 7. 输出目录规范

- 所有 test profile 的 `--outdir` 统一设置为 `circdna.nf/results/<profile_name>/`（如 `results/test_local/`、`results/test_local_gdna/`）
- 不再使用 `/tmp/` 或 `circdna_test_local_*` 等旧式输出目录

## Impact

- Affected specs:
  - `fix-eccsplorer-output-mapping`（输出映射保持不变，新增 BLAST 注释输出会被自动发布）
- Affected code:
  - `circdna.nf/nextflow.config`（新增 `eccsplorer_database` 参数；版本号 4.2.3 → 4.3.0）
  - `circdna.nf/conf/modules.config`（ECCSPLORER 的 ext.args 动态构造 `-d`）
  - `circdna.nf/conf/test_local.config`（设置默认数据库路径、输出目录改为 `results/test_local/`、samplesheet 补充 data_type 列）
  - `circdna.nf/conf/test_local_gdna.config`（新增 profile）
  - `circdna.nf/modules/local/eccsplorer/main.nf`（input 块新增可选 control 通道，script 块条件性附加 control 参数）
  - `circdna.nf/subworkflows/local/eccdna_mode/main.nf`（根据 data_type + group 分离 eccDNA/gDNA，路由 control 通道）
  - `circdna.nf/workflows/circdna.nf`（将 control 通道传入 ECCDNA_MODE）
  - `circdna.nf/samplesheets/samplesheet_local.csv`（补充 data_type=eccDNA 列，无 group 列）
  - `circdna.nf/samplesheets/samplesheet_local_with_gdna.csv`（新增，含 data_type + group 列）
  - `circdna.nf/scripts/prepare_eccsplorer_database.sh`（新增）
  - `circdna.nf/bin/check_samplesheet.py`（兼容 `data_type` 列名与 camelCase 值，新增 `group` 列识别）
  - `circdna.nf/CHANGELOG.md`（v4.3.0 版本记录，MINOR bump：新增功能）
  - `bio.nf/modules/eccsplorer/main.nf`（同步 input/script 修改）
  - `bio.nf/modules/eccsplorer/meta.yml`（同步更新）

## ADDED Requirements

### Requirement: ECCsplorer BLAST 注释数据库支持

系统 SHALL 允许用户通过 `--eccsplorer_database` 参数传入未压缩的 FASTA 注释库，ECCSPLORER 模块 SHALL 通过 `-d` 参数传递给 ECCsplorer 工具，用于对候选 eccDNA 序列进行重复元件/功能注释。

#### Scenario: 启用注释数据库

- **WHEN** 用户运行时传入 `--eccsplorer_database /path/to/anno.fa`
- **THEN** ECCSPLORER 模块的 `ext.args` SHALL 包含 `-d /path/to/anno.fa`
- **AND** 所有候选 eccDNA 的 `*_blast.m6` 文件 SHALL 包含 BLAST 命中结果（非 0 字节，前提是数据库中有匹配序列）

#### Scenario: 未启用注释数据库（默认）

- **WHEN** 用户未传入 `--eccsplorer_database` 或值为 null
- **THEN** ECCSPLORER 模块 SHALL 不传递 `-d` 参数
- **AND** `*_blast.m6` 文件 SHALL 为 0 字节（保持当前行为）

### Requirement: 基于现有 samplesheet `data_type` 列 + 新增 `group` 列的 gDNA control 自动启用

系统 SHALL 通过现有 samplesheet 的 `data_type` 列（列名带下划线，值为 `eccDNA` 或 `gDNA`）自动决定是否启用 gDNA control，并通过新增的可选 `group` 列建立 eccDNA 样本与 gDNA control 的配对关系，**无需用户新增任何命令行参数**。

#### Scenario: samplesheet 包含 gDNA 样本且提供 `group` 列（启用 control）

- **WHEN** samplesheet 中同时存在 `data_type=eccDNA` 和 `data_type=gDNA` 的行，且包含 `group` 列
- **THEN** 流程 SHALL 将同一 `group` 值下的 gDNA 样本 reads 作为该 group 下 eccDNA 样本的 dataset B（control）
- **AND** 一个 group 下有 1 个 gDNA + N 个 eccDNA 时，这 N 个 eccDNA SHALL 共用该 gDNA control
- **AND** 一个 group 下无 gDNA 时，该 group 的 eccDNA SHALL 不启用 control
- **AND** 一个 group 下有多个 gDNA 时，流程 SHALL 报错（一个 group 仅允许 1 个 gDNA control）
- **AND** ECCsplorer 命令 SHALL 包含 gDNA control 的 R1/R2 作为位置参数 3、4

#### Scenario: samplesheet 包含 gDNA 样本但无 `group` 列（向后兼容，仅 1 个 gDNA）

- **WHEN** samplesheet 中同时存在 eccDNA 和 gDNA 行，无 `group` 列，且 gDNA 样本仅 1 个
- **THEN** 流程 SHALL 让所有 eccDNA 样本共用该唯一的 gDNA control
- **AND** ECCsplorer 命令 SHALL 包含 gDNA control 的 R1/R2 作为位置参数 3、4

#### Scenario: samplesheet 包含多个 gDNA 样本但无 `group` 列（报错）

- **WHEN** samplesheet 中存在多个 gDNA 样本，无 `group` 列
- **THEN** 流程 SHALL 报错，提示用户添加 `group` 列以明确 eccDNA 与 gDNA 的配对关系

#### Scenario: samplesheet 仅含 eccDNA 样本（不启用 control）

- **WHEN** samplesheet 中所有行的 `data_type` 均为 `eccDNA`
- **THEN** 流程 SHALL 不启用 gDNA control
- **AND** ECCSPLORER 模块 SHALL 仅接收 eccDNA reads + 参考基因组 + 可选数据库
- **AND** ECCsplorer 命令 SHALL 不包含 dataset B 参数（保持现有行为）

#### Scenario: samplesheet 无 data_type 列（向后兼容）

- **WHEN** samplesheet 不包含 `data_type` 列（向后兼容旧格式）
- **THEN** 流程 SHALL 视所有样本为 `eccDNA`，不启用 control

### Requirement: `check_samplesheet.py` 兼容现有 `data_type` 列格式 + 新增 `group` 列

`check_samplesheet.py` SHALL 同时识别 `data_type` 和 `datatype` 两种列名，接受 camelCase（`eccDNA`/`gDNA`）或小写（`eccdna`/`gdna`）任一写法的值，并识别可选的 `group` 分组列。

#### Scenario: 输入 samplesheet 使用 `data_type` 列名与 camelCase 值

- **WHEN** 输入 samplesheet 头部包含 `data_type` 列，某行值为 `eccDNA`
- **THEN** 脚本 SHALL 接受该行，将值归一化为小写 `eccdna` 存入 meta
- **AND** 输出 samplesheet SHALL 保留原列名 `data_type`

#### Scenario: 输入 samplesheet 使用 `datatype` 列名与小写值（旧格式）

- **WHEN** 输入 samplesheet 头部包含 `datatype` 列，某行值为 `gdna`
- **THEN** 脚本 SHALL 接受该行，将值归一化为小写 `gdna` 存入 meta
- **AND** 输出 samplesheet SHALL 保留原列名 `datatype`

#### Scenario: 输入 samplesheet 包含 `group` 列

- **WHEN** 输入 samplesheet 头部包含 `group` 列，某行值为 `A`
- **THEN** 脚本 SHALL 接受该行，将 `group` 值存入 meta
- **AND** 输出 samplesheet SHALL 保留 `group` 列

#### Scenario: 输入 samplesheet 不包含 `group` 列

- **WHEN** 输入 samplesheet 头部不包含 `group` 列
- **THEN** 脚本 SHALL 视所有样本的 group 为空字符串
- **AND** 输出 samplesheet SHALL 不包含 `group` 列

### Requirement: gDNA 测试数据与对比测试

系统 SHALL 提供 gDNA 测试数据与两个 samplesheet（with_gdna / 无 gDNA），用于对比验证 control 对 eccDNA 检测结果的影响。

#### Scenario: 对比测试

- **WHEN** 用户分别运行 `test_local`（无 control）和 `test_local_gdna`（有 control）profile
- **THEN** 两次运行 SHALL 使用相同的 eccDNA 测试数据（circdna_1/2/3）
- **AND** 输出 SHALL 分别保存到 `circdna.nf/results/test_local/` 和 `circdna.nf/results/test_local_gdna/`
- **AND** 用户可对比两次运行的 `*_candidates.bed` 数量差异验证 control 效果

### Requirement: 统一输出目录规范

所有 test profile 的 `--outdir` SHALL 统一设置为 `circdna.nf/results/<profile_name>/`，不再使用 `/tmp/` 或 `circdna_test_local_*` 等旧式输出目录。

#### Scenario: test_local 输出

- **WHEN** 用户运行 `-profile test_local,docker`
- **THEN** 输出 SHALL 保存到 `circdna.nf/results/test_local/`

#### Scenario: test_local_gdna 输出

- **WHEN** 用户运行 `-profile test_local_gdna,docker`
- **THEN** 输出 SHALL 保存到 `circdna.nf/results/test_local_gdna/`

## MODIFIED Requirements

### Requirement: ECCSPLORER 模块输入

ECCSPLORER 模块的 input 块 SHALL 支持可选的 gDNA control FASTQ 输入，当未提供 control 时保持现有行为。control 完全由 subworkflow 根据 samplesheet 的 `data_type` + `group` 列路由，不通过命令行参数传入。

#### Scenario: 有 control 输入

- **WHEN** samplesheet 中存在 gDNA 样本，subworkflow 将匹配的 gDNA control R1/R2 传入
- **THEN** ECCSPLORER 模块 SHALL 接收 control R1/R2 作为额外输入
- **AND** script 块 SHALL 将 control 作为 dataset B（位置参数 3、4）传递给 ECCsplorer

#### Scenario: 无 control 输入（向后兼容）

- **WHEN** samplesheet 中无 gDNA 样本，subworkflow 未传入 control
- **THEN** ECCSPLORER 模块 SHALL 仅接收 eccDNA reads + fasta（保持现有行为）
- **AND** script 块 SHALL 不传递 dataset B 参数

### Requirement: eccdna_mode subworkflow

eccdna_mode subworkflow SHALL 根据 samplesheet 的 `data_type` 列分离 eccDNA 和 gDNA 样本，通过 `group` 列建立配对关系，将 gDNA 样本的 reads 作为 control 通道路由到 ECCSPLORER 模块。

#### Scenario: samplesheet 含 gDNA 样本且有 `group` 列

- **WHEN** samplesheet 中存在 data_type=gDNA 的行，且包含 `group` 列
- **THEN** subworkflow SHALL 按 `group` 列值将 gDNA reads 与同 group 的 eccDNA 样本配对
- **AND** 配对后的 control FASTQ SHALL 与 eccDNA reads combine 后传入 ECCSPLORER
- **AND** 无 gDNA 配对的 group 下的 eccDNA 样本 SHALL 走无 control 路径

#### Scenario: samplesheet 含 gDNA 样本但无 `group` 列（仅 1 个 gDNA）

- **WHEN** samplesheet 中存在 1 个 gDNA 样本，无 `group` 列
- **THEN** subworkflow SHALL 让所有 eccDNA 样本共用该 gDNA control

#### Scenario: samplesheet 不含 gDNA 样本

- **WHEN** samplesheet 中所有样本 data_type=eccDNA
- **THEN** subworkflow SHALL 不传入 control，ECCSPLORER 仅接收 eccDNA reads + fasta
