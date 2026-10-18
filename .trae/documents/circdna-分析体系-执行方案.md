# circdna.nf 三代 eccDNA 分析体系完善计划

> 本文档在 `circdna-分析体系-蓝图设计.md` 基础上，针对 `circdna.nf` 和 `circdnalr.nf` 两个流程给出可执行的修改方案。
>
> **范围说明**：本计划仅覆盖 `circdna.nf` 与 `circdnalr.nf`。`circrna.nf` 已有其他任务在完善，本次不修改。

## 一、容器与模块构建策略（强制约定）

### 1.1 容器来源优先级

1. **优先**：nf-core 标准模块自带的 conda/docker 环境
2. **次选**：`quay.io/biocontainers/<tool>` 上的 biocontainer 镜像
3. **自构建**：若上述均无，先在用户自己的频道构建：
   - conda 包：上传至 `anaconda.org/yangmingsi/<tool>`
   - docker 镜像：上传至 `quay.io/bioinfortools/<tool>:<version>`

### 1.2 模块来源优先级

1. **优先**：`nf-core modules install <tool>` 安装到 `circdna.nf/modules/nf-core/`
2. **次选**：若 nf-core 无对应模块，先在 **bio.nf** 仓库新建分支构建模块，验证通过后复制到 `circdna.nf/modules/local/`
3. **subworkflow 同样遵循**：先在 bio.nf 验证，再复制到 `circdna.nf/subworkflows/local/`

### 1.3 双分支策略

| 改动类型                                                                   | 目标分支    | 版本号变更      |
| -------------------------------------------------------------------------- | ----------- | --------------- |
| 三代相关（NanoPlot、ecc_finder、SV、CIDER-seq、长读 annotation、统一输出） | `circdnalr` | v3.1.0 → v3.2.0 |
| 二代相关（fastp、二代统一输出）                                            | `master`    | v3.2.0 → v3.2.1 |

## 二、当前状态分析

### 2.1 circdna.nf（v3.1.0）

**已具备**：

- 二代流程：`FastQC → TrimGalore → BWA → BAM_PREPROCESSING → Circle-Map/Circle_Finder/CircExplorer2/AmpliconArchitect/Unicycler`
- 三代流程（`circdnalr` 分支）：`protocol=short_read|pacbio|ont`，`entrypoint=cleaned_fastq|raw_fastq|subreads|hifi_bam`
- 三引擎：CReSIL（mapping）、FLED（mapping）、Flye（assembly）
- 预处理：PacBio（PBCCS + LIMA）、ONT（CHOPPER + PYCHOPPER）
- 结果过滤：`LONG_READ_FILTERING`（min_read_support + blacklist + repeats，通过 `FILTER_ECCDNA_BY_SUPPORT` 模块实现）
- samplesheet 扩展：`schema_input.json` 已支持 `input_bam`、`entrypoint` 列

**缺失项与实际进度（2026-07-24 代码库验证，circdnalr 分支未提交改动）**：

| 功能                  | 计划要求                                   | 实际状态  | 说明                                                                                                                                                                              |
| --------------------- | ------------------------------------------ | --------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| 长读 QC (NanoPlot)    | 集成到 LONG_READ_PREPROCESSING             | ✅ 已实现 | `modules/nf-core/nanoplot/` 已安装；`long_read_preprocessing/main.nf` 已含 NANOPLOT 调用 + `skip_long_read_qc` 开关；`conf/modules.config` 已有 publishDir                        |
| 第四引擎 ecc_finder   | map + asm 子流程                           | ✅ 已实现 | `modules/local/ecc_finder/{map_ont,map_sr,asm_ont,asm_sr}/` 已从 bio.nf 迁移；`subworkflows/local/eccfinder_pipeline/main.nf` 已创建；`eccfinder_mode` 参数已加入 nextflow.config |
| ecc_finder 集成主流程 | workflows/circdna.nf 调用                  | ✅ 已实现 | `long_read_identifier` 默认值已含 `eccfinder`                                                                                                                                     |
| SV 检测               | Sniffles2/cuteSV/SVIM                      | ❌ 待实现 | Phase 3                                                                                                                                                                           |
| CIDER-seq 处理        | DeConcat → TideHunter → CircleSeeker       | ❌ 待实现 | Phase 4                                                                                                                                                                           |
| eccDNA annotation     | Gene/TE/Repeat/CDS/Intergenic              | ❌ 待实现 | Phase 5                                                                                                                                                                           |
| 统一输出              | Unified eccDNA BED/TSV/FASTA               | ❌ 待实现 | Phase 5                                                                                                                                                                           |
| Metadata 分层系统     | Level A-E + 按物种 TSV                     | ❌ 待实现 | Phase 7                                                                                                                                                                           |
| 二代 fastp 替代       | 可选 fastp 替代 TrimGalore                 | ❌ 待实现 | Phase 8（master 分支）                                                                                                                                                            |
| Stub 测试验证         | test_pacbio_lr / test_nanopore_lr          | ❌ 待验证 | 配置文件已建，需跑 stub 测试                                                                                                                                                      |
| CHANGELOG + 版本号    | v3.2.0                                     | ⚠️ 部分   | `manifest.version` 已改为 3.2.0，但 CHANGELOG.md 未更新                                                                                                                           |
| 临时文件清理          | 删除 genome.txt/program.txt/scripts/err 等 | ❌ 待清理 | 工作目录残留临时文件                                                                                                                                                              |

> **结论**：Phase 0（Bug 修复）、Phase 1（资产迁移）、Phase 2（ecc_finder 引擎）、Phase 6（NanoPlot QC）的核心代码**均已实现**（未提交）。剩余工作集中在 Phase 3-5（SV/CIDER-seq/Annotation/Unified Output）、Phase 7（Metadata）、Phase 8（二代 master）、以及收尾工作（stub 测试、CHANGELOG、circdnalr.nf 删除、临时文件清理）。

**已知问题状态（2026-07-24 验证）**：

> 经代码库实际验证，v3.1.0 中的 5 个已知问题已有 4 个在当前代码中修复，仅剩 1 个待处理。

1. ✅ **[已修复] FLED_PIPELINE 输入不匹配 + LONG_READ_MAPPING 未 emit bai**：`long_read_mapping/main.nf` 现已调用 `SAMTOOLS_INDEX`，并通过 `join(by: [0])` 将 sorted_bam 与 bai 合并为 3 元组 `[val(meta), bam, bai]` emit。`fled_pipeline/main.nf` 的 `take` 声明 2 个参数（`mapped_reads` + `genome_fasta`），其中 `mapped_reads` 内部为 3 元组，在 FLED 调用时通过 `.map { meta, bam, bai -> ... }` 正确拆分。
2. ✅ **[已修复] FLYE_PIPELINE ext.args 设置方式不生效**：`flye_pipeline/main.nf` 中已无 `.ext.args` 赋值。`ext.args` 已正确在 `conf/modules.config` 中通过 `withName: 'FLYE' { ext.args = { meta.platform == "pacbio" ? "--pacbio-hifi" : "--nano-hq" } }` 设置。
3. ✅ **[已修复] LONG_READ_FILTERING 未实现 min_read_support 过滤**：`long_read_filtering/main.nf` 现已通过 `FILTER_ECCDNA_BY_SUPPORT` 模块（`modules/local/filter_eccdna_by_support/`）+ `bin/filter_by_read_support.py` 实现基于 `params.min_read_support` 的 read 支持数过滤，随后再进行 blacklist + repeats 过滤。
4. ❌ **[待修复] 长读测试配置缺失**：`conf/` 下无 `test_pacbio_lr.config` 或 `test_nanopore_lr.config`，现有测试配置（`test.config`、`test_AA.config` 等）均针对二代数据。需新建长读专用测试配置。
5. ✅ **[已修复] 主流程长读分支 minimap2 preset 切换不够灵活**：`conf/modules.config` 中 `MINIMAP2_ALIGN` 的 `ext.args` 已改为 `{ meta.platform == "pacbio" ? "-x map-hifi" : "-x map-ont" }`，使用 `meta.platform` 驱动而非 `params.protocol`。

### 2.2 circdnalr.nf（v1.0.0dev）

**可保留资产**（迁入 `circdna.nf` 的 `circdnalr` 分支）：

| 资产                          | 路径                                                                                                                             | 保留理由                                                                                                                                      | 迁入目标                                                                                     |
| ----------------------------- | -------------------------------------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------- |
| **NanoPlot 模块**             | `circdnalr.nf/modules/nf-core/nanoplot/`                                                                                         | 长读 QC 必需，circdna.nf 缺失                                                                                                                 | `circdna.nf/modules/nf-core/nanoplot/`（通过 `nf-core modules install nanoplot` 安装最新版） |
| **平台自动切换设计**          | `circdnalr.nf/conf/modules.config` 中 `MINIMAP2_ALIGN.ext.args = { meta.platform == 'pacbio' ? '-ax map-hifi' : '-ax map-ont' }` | **核心设计**：通过 meta.platform 驱动 minimap2 preset 自动切换。circdna.nf 已采纳此设计（`-x map-hifi`/`-x map-ont`），仅参考 publishDir 模式 | ✅ 已整合到 circdna.nf 的 `conf/modules.config`，保持 `meta.platform` 驱动                   |
| **多识别器并行机制**          | `circdnalr.nf/workflows/circdna_lr.nf` 中 `circle_identifier.tokenize(',').collect { it.trim().toLowerCase() }`                  | 优雅的逗号分隔多选设计，可同时跑 CReSIL+ecc_finder_map+ecc_finder_asm                                                                         | 参考，circdna.nf 已有类似机制（`long_read_identifier`），确保一致性                          |
| **test_pacbio.config 思路**   | `circdnalr.nf/conf/test_pacbio.config`                                                                                           | 测试样本来源（PRJNA896947 CHO-K1, ~33MB 20k reads）、资源限制（2 CPU/6GB/6h）                                                                 | `circdna.nf/conf/test_pacbio_lr.config`，需将 GitHub URL 数据本地化                          |
| **test_nanopore.config 思路** | `circdnalr.nf/conf/test_nanopore.config`                                                                                         | 测试样本来源（PRJNA1251131, ~9MB 1.5k reads）、资源限制                                                                                       | `circdna.nf/conf/test_nanopore_lr.config`，需将 GitHub URL 数据本地化                        |
| **长读模块 publishDir 模式**  | `circdnalr.nf/conf/modules.config`                                                                                               | NanoPlot/minimap2/samtools/ecc_finder 的完整输出目录结构（qc/alignment/eccdna/multiqc）                                                       | 参考，整合到 circdna.nf 的 `conf/modules.config`                                             |
| **WorkflowCircdnalr.groovy**  | `circdnalr.nf/lib/WorkflowCircdnalr.groovy`                                                                                      | MultiQC summary 生成逻辑、引用文本                                                                                                            | 参考迁移到 circdna.nf 的 `lib/WorkflowCircdna.groovy`                                        |

> **ecc_finder 模块来源决策**：circdnalr.nf 的 `modules/local/eccfinder/{map,asm}.nf` 为简化版（单文件、无 meta.yml、用 `map-sr`/`asm-sr` 统一命令通过 `--data-type` 切换平台）。**bio.nf 已有更完整的 ecc_finder 模块**（`map_ont`/`map_sr`/`asm_ont`/`asm_sr` 四个子模块，含完整 meta.yml + environment.yml + tests + testdata），按模块来源优先级策略，**应从 bio.nf 迁移完整版**，放弃 circdnalr.nf 的简化版。

> **CReSIL 模块说明**：circdnalr.nf 用 CReSIL 1.1.1（旧版单 `detect` 命令），bio.nf 用 CReSIL 1.2.1（新版拆为 trim/identify/identify_wgls/annotate/visualize 5 个模块）。circdna.nf v3.1.0 已集成 bio.nf 的完整版本，无需从 circdnalr.nf 迁移。

> **samplesheet schema 不一致问题**：circdnalr.nf 的 `assets/schema_input.json` 定义格式为 `sample,fastq,platform`（单 fastq + 平台列），但 `samplesheets/circdnalr_*.csv` 实际生成的是 `sample,fastq_1,fastq_2`（双 fastq 列、无 platform 列、第二列留空），二者不匹配。迁移时需统一为 circdna.nf 的 `sample,fastq_1,fastq_2` 格式（已在 v3.1.0 中支持长读单端）。

**直接删除/放弃**：

- `main.nf`、`workflows/circdna_lr.nf`、`nextflow.config`、`nextflow_schema.json` —— 已被 circdna.nf v3.1.0 取代
- `lib/Utils.groovy`、`conf/base.config` —— 已被 circdna.nf v3.1.0 取代
- `samplesheets/` —— 已在 circdna.nf/samplesheets/circdnalr\_\*.csv 中存在，且路径写死需重新生成
- `CHANGE&FIX/20260701.md` —— 历史记录，无技术价值
- `CRESIL_DETECT` 简化版（`modules/local/cresil/detect.nf`，v1.1.1）—— circdna.nf 已有更完整的 CReSIL 1.2.1（trim+identify+identify_wgls+annotate+visualize）
- `ECCFINDER_MAP`/`ECCFINDER_ASM` 简化版（`modules/local/eccfinder/{map,asm}.nf`）—— bio.nf 有更完整的 4 子模块版本
- `assets/email_template.*`、`assets/multiqc_config.yml` —— nf-core 通用模板，circdna.nf 已有
- `tree.txt`、`acc_list.txt` —— 辅助文件，无迁移价值

**结论（2026-07-24 最终评估）**：

经逐文件核查，circdnalr.nf 的**所有可保留资产均已迁移**到 circdna.nf 的 circdnalr 分支：

| circdnalr.nf 资产                                       | 迁移状态      | circdna.nf 对应位置                                                                            |
| ------------------------------------------------------- | ------------- | ---------------------------------------------------------------------------------------------- |
| NanoPlot 模块 (`modules/nf-core/nanoplot/`)             | ✅ 已迁移     | `circdna.nf/modules/nf-core/nanoplot/`（含 tests）                                             |
| 平台自动切换设计 (`conf/modules.config` MINIMAP2_ALIGN) | ✅ 已采纳     | `circdna.nf/conf/modules.config`（`-x map-hifi`/`-x map-ont`）                                 |
| 长读 publishDir 模式                                    | ✅ 已整合     | `circdna.nf/conf/modules.config`（long_read/ 目录结构）                                        |
| 多识别器并行机制 (`workflows/circdna_lr.nf`)            | ✅ 已采纳     | `circdna.nf/workflows/circdna.nf`（`long_read_identifier` 逗号分隔）                           |
| test_pacbio/test_nanopore 配置思路                      | ✅ 已参考重建 | `circdna.nf/conf/test_pacbio_lr.config` + `test_nanopore_lr.config`（数据本地化）              |
| WorkflowCircdnalr.groovy                                | ❌ 无迁移价值 | 仅含 NanoPlot/Minimap2/Samtools/MultiQC 4 条引用，circdna.nf 已有自己的 WorkflowCircdna.groovy |
| ecc_finder 简化版 (`modules/local/eccfinder/`)          | ❌ 已放弃     | 改用 bio.nf 完整版 4 子模块                                                                    |
| CReSIL v1.1.1 简化版 (`modules/local/cresil/detect.nf`) | ❌ 已放弃     | circdna.nf 已有 CReSIL v1.2.1 完整版                                                           |
| samplesheets/ (路径写死)                                | ❌ 无迁移价值 | circdna.nf/samplesheets/ 已重新生成                                                            |
| assets/、docs/、.github/                                | ❌ 无迁移价值 | nf-core 通用模板，circdna.nf 已有                                                              |
| lib/Utils.groovy、conf/base.config                      | ❌ 已被取代   | circdna.nf 已有更完整版本                                                                      |

> **circdnalr.nf 已无任何剩余保留价值。所有有用资产均已迁移完毕。**
>
> **处理决策**：直接删除整个 `circdnalr.nf/` 目录（依赖 GitHub 版本控制保留历史，无需本地归档）。删除时机：在 circdna.nf circdnalr 分支的当前未提交改动提交并验证 stub 测试通过后执行。

## 三、修改方案

### 3.1 circdna.nf（circdnalr 分支，v3.2.0）

#### 3.1.1 资产迁移

##### A. ecc_finder 模块迁移（从 bio.nf 迁移完整版）

**来源**：`bio.nf/modules/ecc_finder/{map_ont,map_sr,asm_ont,asm_sr}/`（每个子模块含 `main.nf` + `meta.yml` + `environment.yml` + `tests/` + `testdata/`）
**目标**：`circdna.nf/modules/local/ecc_finder/{map_ont,map_sr,asm_ont,asm_sr}/`

> **为什么不从 circdnalr.nf 迁移**：circdnalr.nf 的 ecc_finder 为简化版（单文件 `map.nf`/`asm.nf`，无 meta.yml，用 `map-sr`/`asm-sr` 统一命令通过 `--data-type` 切换平台，无测试）。bio.nf 的版本按平台拆分为 4 个独立子模块（`map_ont`/`map_sr`/`asm_ont`/`asm_sr`），更符合 ecc_finder 原生设计，且自带完整测试和 testdata。按模块来源优先级策略，优先使用 bio.nf 完整版。

**迁移步骤**：

1. 从 bio.nf 的 `ecc_finder` 分支复制 4 个子模块到 `circdna.nf/modules/local/ecc_finder/`：
   ```
   cp -r bio.nf/modules/ecc_finder/map_ont  circdna.nf/modules/local/ecc_finder/map_ont/
   cp -r bio.nf/modules/ecc_finder/map_sr   circdna.nf/modules/local/ecc_finder/map_sr/
   cp -r bio.nf/modules/ecc_finder/asm_ont  circdna.nf/modules/local/ecc_finder/asm_ont/
   cp -r bio.nf/modules/ecc_finder/asm_sr   circdna.nf/modules/local/ecc_finder/asm_sr/
   ```
2. 保留各子模块的 `main.nf` + `meta.yml` + `environment.yml`，可选保留 `tests/` 和 `testdata/`（用于模块级测试）
3. 容器配置：bio.nf 版本使用 `quay.io/bioinfortools/ecc_finder:1.0.0`（用户自构建频道），已满足策略
4. conda 环境：`bioconda::ecc_finder=1.0.0`
5. 输入输出签名（以 `map_ont` 为例）：
   - 输入：`tuple val(meta), path(idx), path(query), path(ref)` — 注意 bio.nf 版本需要 minimap2 索引
   - 输出：`tuple val(meta), path("*.csv"), path("*.fasta")` + topic versions 通道

**模块与平台对应关系**：

| 子模块               | 命令                    | 适用平台    | 用途                       |
| -------------------- | ----------------------- | ----------- | -------------------------- |
| `ECC_FINDER_MAP_ONT` | `ecc_finder.py map-ont` | ONT         | mapping-based eccDNA 检测  |
| `ECC_FINDER_MAP_SR`  | `ecc_finder.py map-sr`  | PacBio HiFi | mapping-based eccDNA 检测  |
| `ECC_FINDER_ASM_ONT` | `ecc_finder.py asm-ont` | ONT         | assembly-based eccDNA 检测 |
| `ECC_FINDER_ASM_SR`  | `ecc_finder.py asm-sr`  | PacBio HiFi | assembly-based eccDNA 检测 |

**新增** `modules/local/ecc_finder/environment.yml`（所有子模块共用）：

```yaml
name: ecc_finder
channels:
  - conda-forge
  - bioconda
dependencies:
  - bioconda::ecc_finder=1.0.0
```

##### B. NanoPlot 模块安装

**方式**：通过 `nf-core modules install nanoplot` 直接安装到 `circdna.nf/modules/nf-core/nanoplot/`

**验证**：`nf-core modules list local` 应显示 nanoplot

##### C. 测试配置迁移（Phase 0 Bug 4 修复）

> 数据来源参考 circdnalr.nf 的 `conf/test_pacbio.config` 和 `conf/test_nanopore.config`。
>
> - PacBio 测试数据：PRJNA896947 CHO-K1（~33MB, 20k reads），参考基因组 chr1_subset.fa
> - Nanopore 测试数据：PRJNA1251131（~9MB, 1.5k reads），参考基因组 chr22_subset.fa

**新建** `circdna.nf/conf/test_pacbio_lr.config`：

```groovy
params {
    config_profile_name        = 'Test PacBio Long-Read profile'
    config_profile_description = 'Minimal PacBio eccDNA dataset for circdna.nf long-read pipeline'

    max_cpus   = 2
    max_memory = '6.GB'
    max_time   = '6.h'

    // Input data - PacBio HiFi eccDNA (PRJNA896947 CHO-K1, ~33MB, 20k reads)
    input     = './samplesheets/test_pacbio_lr.csv'
    fasta     = './testdatasets/reference/chr1_subset.fa'

    protocol              = 'pacbio'
    entrypoint            = 'cleaned_fastq'
    long_read_identifier  = 'cresil,fled,flye,eccfinder'
    circle_identifier      = ''
    input_format           = 'FASTQ'
    igenomes_ignore        = true
}
```

**新建** `circdna.nf/conf/test_nanopore_lr.config`：

```groovy
params {
    config_profile_name        = 'Test Nanopore Long-Read profile'
    config_profile_description = 'Minimal Nanopore eccDNA dataset for circdna.nf long-read pipeline'

    max_cpus   = 2
    max_memory = '6.GB'
    max_time   = '6.h'

    // Input data - Nanopore eccDNA (PRJNA1251131, ~9MB, 1.5k reads)
    input     = './samplesheets/test_nanopore_lr.csv'
    fasta     = './testdatasets/reference/chr22_subset.fa'

    protocol              = 'ont'
    entrypoint            = 'cleaned_fastq'
    long_read_identifier  = 'cresil,fled,flye,eccfinder'
    circle_identifier      = ''
    input_format           = 'FASTQ'
    igenomes_ignore        = true
}
```

**注册 profiles**：在 `nextflow.config` 的 `profiles` 块中添加：

```groovy
test_pacbio_lr     { includeConfig 'conf/test_pacbio_lr.config' }
test_nanopore_lr   { includeConfig 'conf/test_nanopore_lr.config' }
```

**测试数据本地化**：需将 circdnalr.nf 的 GitHub URL 数据下载到 `circdna.nf/testdatasets/` 和 `circdna.nf/samplesheets/` 下，避免依赖外部仓库。

#### 3.1.2 新增 ecc_finder 子流程（第四引擎）

**新建文件**：`circdna.nf/subworkflows/local/eccfinder_pipeline/main.nf`

> **设计说明**：bio.nf 的 ecc_finder 按平台拆分为 4 个子模块（`map_ont`/`map_sr`/`asm_ont`/`asm_sr`），子流程根据 `params.protocol` 自动选择对应平台模块。`eccfinder_mode` 控制运行 map（mapping-based）、asm（assembly-based）还是 both。

```groovy
include { ECC_FINDER_MAP_ONT } from '../../../modules/local/ecc_finder/map_ont/main'
include { ECC_FINDER_MAP_SR  } from '../../../modules/local/ecc_finder/map_sr/main'
include { ECC_FINDER_ASM_ONT } from '../../../modules/local/ecc_finder/asm_ont/main'
include { ECC_FINDER_ASM_SR  } from '../../../modules/local/ecc_finder/asm_sr/main'

workflow ECCFINDER_PIPELINE {
    take:
    reads           // [ val(meta), fastq ]
    genome_fasta    // file

    main:
    def mode = params.eccfinder_mode   // map | asm | both
    def all_candidates = channel.empty()

    // mapping-based 检测
    if (mode == 'map' || mode == 'both') {
        if (params.protocol == 'ont') {
            ECC_FINDER_MAP_ONT ( reads, genome_fasta )
                .csv
                .set { map_candidates }
        } else {
            ECC_FINDER_MAP_SR ( reads, genome_fasta )
                .csv
                .set { map_candidates }
        }
        all_candidates = all_candidates.mix(map_candidates)
    }

    // assembly-based 检测
    if (mode == 'asm' || mode == 'both') {
        if (params.protocol == 'ont') {
            ECC_FINDER_ASM_ONT ( reads, genome_fasta )
                .csv
                .set { asm_candidates }
        } else {
            ECC_FINDER_ASM_SR ( reads, genome_fasta )
                .csv
                .set { asm_candidates }
        }
        all_candidates = all_candidates.mix(asm_candidates)
    }

    emit:
    eccdna_candidates    // [ val(meta), csv ]
    versions             // topic channel
}
```

**新增参数**（`nextflow.config` + `nextflow_schema.json`）：

```groovy
// ecc_finder options
eccfinder_mode               = 'map'   // map | asm | both
eccfinder_min_length         = 100
```

**主流程集成**（`workflows/circdna.nf`）：

在 `if (run_flye) { ... }` 之后增加：

```groovy
def run_eccfinder = ("eccfinder" in lr_branch)
if (run_eccfinder) {
    ECCFINDER_PIPELINE (
        ch_preprocessed_fastq,
        ch_fasta
    )
    .eccdna_candidates
    .set { ch_eccfinder_candidates }

    LONG_READ_FILTERING ( ch_eccfinder_candidates )
    .filtered_candidates
    .set { ch_filtered_eccfinder }
}
```

更新 `long_read_identifier` 默认值为 `cresil,fled,flye,eccfinder`，并在 schema 中扩展 enum。

#### 3.1.3 新增 SV 检测子流程（WGS 模式）

**适用数据**：PacBio HiFi WGS（ERR11838731/SRR30359583）、ONT WGS（ERR12723706）

**模块来源策略**：

1. 尝试 `nf-core modules install sniffles2` / `nf-core modules install cutesv` / `nf-core modules install svim`
2. 若 nf-core 无对应模块，按 1.2 节策略在 bio.nf 新建分支构建后复制到 `circdna.nf/modules/local/`
3. 容器优先 `quay.io/biocontainers/<tool>`，无则自构建上传至 `quay.io/bioinfortools/<tool>` 和 `anaconda.org/yangmingsi/<tool>`

**新建子流程**：`circdna.nf/subworkflows/local/sv_detection_pipeline/main.nf`

```groovy
include { SNIFFLES2 } from '../../../modules/nf-core/sniffles2/main'
include { CUTESV    } from '../../../modules/nf-core/cutesv/main'
include { SVIM      } from '../../../modules/nf-core/svim/main'

workflow SV_DETECTION_PIPELINE {
    take:
    bam_sorted      // [ val(meta), bam, bai ]
    genome_fasta    // file

    main:
    def sv_tools = params.sv_identifier.split(',')
    if ('sniffles2' in sv_tools) {
        SNIFFLES2 ( bam_sorted, genome_fasta )
    }
    if ('cutesv' in sv_tools) {
        CUTESV ( bam_sorted, genome_fasta )
    }
    if ('svim' in sv_tools) {
        SVIM ( bam_sorted, genome_fasta )
    }

    emit:
    sv_calls        // [ val(meta), vcf ]
    versions
}
```

**新增参数**：

```groovy
// SV detection options
sv_identifier                = 'sniffles2'  // sniffles2|cutesv|svim (comma-separated)
run_sv_for_wgs               = false
min_sv_length                = 100
```

**主流程集成**：

在 `LONG_READ_MAPPING` 之后，根据 metadata 中 `analysis_role == 'reference_genome'` 或 `experimental_type == 'WGS_reference'` 触发：

```groovy
if (params.run_sv_for_wgs && params.protocol in ["pacbio", "ont"]) {
    SV_DETECTION_PIPELINE ( ch_long_read_mapped, ch_fasta )
}
```

#### 3.1.4 新增 CIDER-seq 专门处理子流程

**适用数据**：SRR16958693（Amaranthus_CIDER）、SRR26069818（Alopecurus_RCA）

**模块来源策略**：

- `deconcat`、`tidehunter`、`circleseeker` 在 nf-core 中均无现成模块
- 按 1.2 节策略，先在 bio.nf 新建分支 `ciderseq`，分别构建 `modules/deconcat/`、`modules/tidehunter/`、`modules/circleseeker/`
- 容器检查顺序：
  1. `quay.io/biocontainers/deconcat`、`quay.io/biocontainers/tidehunter`、`quay.io/biocontainers/circleseeker`
  2. 若无，构建 docker 镜像推送至 `quay.io/bioinfortools/<tool>:<version>`
  3. conda 包推送至 `anaconda.org/yangmingsi/<tool>`
- 验证通过后复制到 `circdna.nf/modules/local/{deconcat,tidehunter,circleseeker}/`

**新建子流程**：`circdna.nf/subworkflows/local/ciderseq_pipeline/main.nf`

```groovy
include { DECONCAT     } from '../../../modules/local/deconcat/main'
include { TIDEHUNTER   } from '../../../modules/local/tidehunter/main'
include { CIRCLESEEKER } from '../../../modules/local/circleseeker/main'

workflow CIDERSEQ_PIPELINE {
    take:
    reads           // [ val(meta), fastq ]
    genome_fasta    // file

    main:
    DECONCAT ( reads )
        .deconcat_reads
        .set { deconcat_reads }

    TIDEHUNTER ( deconcat_reads )
        .tandem_cleaned_reads
        .set { cleaned_reads }

    CIRCLESEEKER ( cleaned_reads, genome_fasta )
        .circular_candidates
        .set { eccdna_candidates }

    emit:
    eccdna_candidates
    versions
}
```

**新增参数**：

```groovy
// CIDER-seq options
ciderseq_enabled             = false
deconcat_min_length          = 200
tidehunter_tid               = 2
circleseeker_min_circle_length = 100
```

**主流程集成**：

在 `LONG_READ_PREPROCESSING` 之后，根据 `--ciderseq_enabled` 或 metadata 中 `CIDER_seq == 'yes'` 触发：

```groovy
if (params.ciderseq_enabled) {
    CIDERSEQ_PIPELINE ( ch_preprocessed_fastq, ch_fasta )
}
```

#### 3.1.5 新增 NanoPlot 长读 QC 步骤

**修改文件**：`circdna.nf/subworkflows/local/long_read_preprocessing/main.nf`

在 `LONG_READ_PREPROCESSING` 开头增加 NanoPlot QC：

```groovy
include { NANOPLOT } from '../../../modules/nf-core/nanoplot/main'

workflow LONG_READ_PREPROCESSING {
    take:
    reads

    main:
    // QC with NanoPlot (works for both PacBio and ONT)
    if (!params.skip_long_read_qc) {
        NANOPLOT ( reads.map { meta, fastq, bam, ep -> [meta, fastq ?: bam] } )
        ch_nanoplot_versions = NANOPLOT.out.versions
    }

    // 既有 PacBio/ONT 预处理逻辑保持不变
    if (params.protocol == "pacbio") { /* 既有逻辑 */ }
    else if (params.protocol == "ont") { /* 既有逻辑 */ }

    emit:
    preprocessed_fastq
    nanoplot_qc = NANOPLOT.out.multiqc_files
    versions
}
```

**新增配置**（`conf/modules.config`）：

```groovy
withName: 'NANOPLOT' {
    ext.args = ''
    publishDir = [
        path: { "${params.outdir}/qc/nanoplot/${meta.id}" },
        mode: params.publish_dir_mode,
        saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    ]
}
```

#### 3.1.6 新增 eccDNA Annotation 子流程

**模块来源策略**：

- `bedtools/intersect` 已存在于 `circdna.nf/modules/local/bedtools/intersect/`
- annotation 整合脚本为自定义 Python，新建 `modules/local/annotate_eccdna/`
- 按 1.2 节策略，先在 bio.nf 新建 `modules/annotate_eccdna/`，验证后复制

**新建脚本**：`bin/annotate_eccdna.py`

```python
#!/usr/bin/env python3
"""Annotate eccDNA candidates with Gene/TE/Repeat/CDS/Intergenic regions."""
import argparse
import pandas as pd
import pybedtools

def annotate(input_bed, gene_bed, te_bed, repeat_bed, cds_bed, output_tsv):
    ecc = pybedtools.BedTool(input_bed)
    # ... 逐个 intersect 并合并注释列
    # 输出统一字段：eccDNA_ID, chromosome, start, end, length, gene, TE, repeat, CDS, region_type
    pass

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', required=True)
    parser.add_argument('--gene-bed', required=True)
    parser.add_argument('--te-bed', required=True)
    parser.add_argument('--repeat-bed', required=True)
    parser.add_argument('--cds-bed', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    annotate(args.input, args.gene_bed, args.te_bed, args.repeat_bed, args.cds_bed, args.output)
```

**新建子流程**：`circdna.nf/subworkflows/local/eccdna_annotation/main.nf`

```groovy
include { ANNOTATE_ECCDNA } from '../../../modules/local/annotate_eccdna/main'

workflow ECCDNA_ANNOTATION {
    take:
    eccdna_candidates    // [ val(meta), bed ]
    gene_bed
    te_bed
    repeat_bed
    cds_bed

    main:
    ANNOTATE_ECCDNA (
        eccdna_candidates,
        gene_bed,
        te_bed,
        repeat_bed,
        cds_bed
    )
    .annotated_bed
    .set { annotated_eccdna }

    emit:
    annotated_eccdna
    versions
}
```

**新增参数**：

```groovy
// Annotation options
gene_bed                     = null
te_bed                       = null
cds_bed                      = null
// 注：repeat_bed 已存在，用于过滤；此处新增的 te_bed 等用于注释
```

#### 3.1.7 新增统一输出接口（Unified eccDNA Table）

**新建脚本**：`bin/unify_eccdna_output.py`

输入：CReSIL/FLED/Flye/ecc_finder/CircleSeeker 各引擎的原始输出
输出：统一 TSV/BED/FASTA，字段遵循计划第 13 节：

```
eccDNA_ID | species | sample_id | chromosome | start | end | length | strand |
support_reads | support_method | platform | technology | circularity_score |
breakpoint_confidence | gene | TE | repeat | GC | copy_number | abundance
```

**新建子流程**：`circdna.nf/subworkflows/local/unified_output/main.nf`

```groovy
workflow UNIFIED_OUTPUT {
    take:
    cresil_results
    fled_results
    flye_results
    eccfinder_results
    ciderseq_results
    annotated_eccdna

    main:
    // 调用 unify_eccdna_output.py 合并所有引擎结果
    // 输出: ${meta.id}_unified_eccdna.tsv / .bed / .fa

    emit:
    unified_tsv
    unified_bed
    unified_fasta
}
```

**主流程集成**（在 MultiQC 之前）：

```groovy
UNIFIED_OUTPUT (
    run_cresil ? ch_filtered_cresil : channel.empty(),
    run_fled ? ch_filtered_fled : channel.empty(),
    run_flye ? ch_flye_contigs : channel.empty(),
    run_eccfinder ? ch_filtered_eccfinder : channel.empty(),
    run_ciderseq ? ch_ciderseq_candidates : channel.empty(),
    ch_annotated_eccdna
)
```

#### 3.1.8 主流程 `workflows/circdna.nf` 整体调整

**长读分支扩展**（`if (params.protocol in ["pacbio", "ont"])` 块）：

```
INPUT_CHECK → LONG_READ_PREPROCESSING (含 NanoPlot QC)
            ├─ CIDER-seq 分支（条件触发）
            │   └─ DECONCAT → TIDEHUNTER → CIRCLESEEKER
            ├─ CReSIL 分支（已有）
            │   └─ CRESIL_PIPELINE → LONG_READ_FILTERING
            ├─ FLED 分支（已有）
            │   └─ LONG_READ_MAPPING → FLED_PIPELINE → LONG_READ_FILTERING
            ├─ Flye 分支（已有）
            │   └─ FLYE_PIPELINE
            ├─ ecc_finder 分支（新增）
            │   └─ ECCFINDER_PIPELINE → LONG_READ_FILTERING
            └─ SV 检测分支（WGS 模式条件触发）
                └─ SV_DETECTION_PIPELINE

ECCDNA_ANNOTATION (gene/te/repeat/cds) → UNIFIED_OUTPUT → MultiQC
```

**新增 includes**（文件顶部）：

```groovy
include { ECCFINDER_PIPELINE } from '../subworkflows/local/eccfinder_pipeline/main'
include { SV_DETECTION_PIPELINE } from '../subworkflows/local/sv_detection_pipeline/main'
include { CIDERSEQ_PIPELINE } from '../subworkflows/local/ciderseq_pipeline/main'
include { ECCDNA_ANNOTATION } from '../subworkflows/local/eccdna_annotation/main'
include { UNIFIED_OUTPUT } from '../subworkflows/local/unified_output/main'
```

#### 3.1.9 Metadata 分层系统

**新建目录**：`circdna.nf/metadata/eccDNA/`

按计划第 4 节，每个物种一个 TSV：

- `Oryza_sativa.tsv`、`Triticum_aestivum.tsv`、`Arabidopsis_thaliana.tsv` 等

**字段**：遵循计划第 4 节 eccDNA Metadata 字段表（含 `eccdna_enrichment`、`exonuclease_treatment`、`plasmid_safe`、`RCA`、`phi29`、`debranching`、`T7_treatment`、`linearization`、`CIDER_seq`、`mobilome_seq`、`circle_seq`、`eccDNA_confidence` 等）。

**关联 samplesheet**：在 `assets/schema_input.json` 增加可选列 `metadata_tsv`，指向物种 metadata 文件。

#### 3.1.10 配置文件更新

**`nextflow.config`** 新增 params：

```groovy
// ecc_finder options
eccfinder_mode               = 'map'
eccfinder_min_length         = 100

// SV detection options
sv_identifier                = 'sniffles2'
run_sv_for_wgs               = false
min_sv_length                = 100

// CIDER-seq options
ciderseq_enabled             = false
deconcat_min_length          = 200
tidehunter_tid               = 2
circleseeker_min_circle_length = 100

// Annotation options
gene_bed                     = null
te_bed                       = null
cds_bed                      = null

// Unified output
output_unified_eccdna        = true
```

**`nextflow_schema.json`** 新增对应 schema 定义。

**`conf/base.config`** 资源分配：

```groovy
withName: 'NANOPLOT' { label 'process_low' }
withName: 'ECCFINDER_MAP' { label 'process_high' }
withName: 'ECCFINDER_ASM' { label 'process_high' }
withName: 'SNIFFLES2' { label 'process_medium' }
withName: 'CUTESV' { label 'process_medium' }
withName: 'SVIM' { label 'process_medium' }
withName: 'DECONCAT' { label 'process_medium' }
withName: 'TIDEHUNTER' { label 'process_medium' }
withName: 'CIRCLESEEKER' { label 'process_high' }
withName: 'ANNOTATE_ECCDNA' { label 'process_low' }
```

**`conf/modules.config`** 新增 publishDir 配置：

```groovy
withName: 'ECCFINDER_MAP' {
    publishDir = [
        path: { "${params.outdir}/eccdna/eccfinder_map/${meta.id}" },
        mode: params.publish_dir_mode
    ]
}
withName: 'ECCFINDER_ASM' {
    publishDir = [
        path: { "${params.outdir}/eccdna/eccfinder_asm/${meta.id}" },
        mode: params.publish_dir_mode
    ]
}
withName: 'SNIFFLES2' {
    publishDir = [
        path: { "${params.outdir}/sv/sniffles2/${meta.id}" },
        mode: params.publish_dir_mode
    ]
}
// ... (cutesv, svim, deconcat, tidehunter, circleseeker 同上)
```

#### 3.1.11 版本与变更记录

- `CHANGELOG.md` 新增 `v3.2.0` 条目，记录所有 circdnalr 分支改动
- `CHANGES&FIX/20260724.md` 详细记录本次完善
- `manifest.version` 改为 `3.2.0`

#### 3.1.12 现有 bug 修复状态（Phase 0）

> 2026-07-24 代码库验证结果：5 个已知问题中 4 个已在当前代码中修复，仅 Bug 4 待处理。

**Bug 1: ✅ 已修复 — FLED_PIPELINE 输入不匹配 + LONG_READ_MAPPING 未 emit bai**

- **当前状态**：`long_read_mapping/main.nf` 已调用 `SAMTOOLS_INDEX`，通过 `join(by: [0])` 合并 sorted_bam + bai，emit 3 元组 `[val(meta), bam, bai]`。`fled_pipeline/main.nf` 的 `take` 接收 `mapped_reads`（内部 3 元组）+ `genome_fasta`，在 FLED 调用时通过 `.map { meta, bam, bai -> ... }` 正确拆分。
- **无需额外操作**。

**Bug 2: ✅ 已修复 — FLYE_PIPELINE ext.args 设置方式不生效**

- **当前状态**：`flye_pipeline/main.nf` 中无 `.ext.args` 赋值。`ext.args` 已在 `conf/modules.config` 中通过 `withName: 'FLYE' { ext.args = { meta.platform == "pacbio" ? "--pacbio-hifi" : "--nano-hq" } }` 正确设置。
- **无需额外操作**。

**Bug 3: ✅ 已修复 — LONG_READ_FILTERING 未实现 min_read_support 过滤**

- **当前状态**：`long_read_filtering/main.nf` 已通过 `FILTER_ECCDNA_BY_SUPPORT` 模块（`modules/local/filter_eccdna_by_support/`）+ `bin/filter_by_read_support.py` 实现基于 `params.min_read_support` 的 read 支持数过滤，流程为：先按 read support 过滤 → 再按 blacklist 过滤 → 再按 repeats 过滤。
- **无需额外操作**。

**Bug 4: ❌ 待修复 — 长读测试配置缺失**

- **问题**：`conf/` 下无 `test_pacbio_lr.config` 或 `test_nanopore_lr.config`，现有测试配置（`test.config`、`test_AA.config` 等）均针对二代数据
- **修复**：新建 `conf/test_pacbio_lr.config` 和 `conf/test_nanopore_lr.config`（见 3.1.1 C 节），并在 `nextflow.config` 的 profiles 中注册
- **验证**：`nextflow run . -stub -profile test_pacbio_lr` 和 `nextflow run . -stub -profile test_nanopore_lr` 成功

**Bug 5: ✅ 已修复 — minimap2 preset 切换为 meta.platform 驱动**

- **当前状态**：`conf/modules.config` 中 `MINIMAP2_ALIGN` 的 `ext.args` 已改为 `{ meta.platform == "pacbio" ? "-x map-hifi" : "-x map-ont" }`，使用 `meta.platform` 驱动。
- **注意**：circdnalr.nf 中使用 `-ax map-hifi`（含 `a` 前缀，输出 SAM 格式），circdna.nf 中使用 `-x map-hifi`（仅 preset）。两者均可工作，但 `-ax` 会显式指定输出格式为 SAM。保持 circdna.nf 现有 `-x` 写法即可（minimap2 默认输出 SAM）。
- **无需额外操作**。

### 3.2 circdna.nf（master 分支，v3.2.1）

> 二代流程改动较小，仅做必要完善。

#### 3.2.1 增加可选 fastp 替代 TrimGalore

**模块来源**：`nf-core modules install fastp`

**新增参数**：`trimmer`：enum `["trimgalore", "fastp"]`，default `trimgalore`

**修改** `workflows/circdna.nf`：

```groovy
if (params.trimmer == "fastp") {
    FASTP ( ch_cat_fastq )
} else {
    TRIMGALORE ( ch_cat_fastq )
}
```

#### 3.2.2 二代流程输出接口标准化

复用 `UNIFIED_OUTPUT` 子流程（与三代共用），将 Circle-Map/Circle_Finder/CircExplorer2 输出转换为统一 BED/TSV。

#### 3.2.3 版本与变更记录

- `CHANGELOG.md` 新增 `v3.2.1` 条目
- `CHANGES&FIX/20260724_master.md` 记录二代改动

### 3.3 circdnalr.nf 处理

**步骤 1**：将 3.1.1 中列出的资产复制到 `circdna.nf` 的 `circdnalr` 分支。

**步骤 2**：完成资产迁移并验证后，将 `circdnalr.nf/` 整体归档为 `circdnalr.nf.archive/` 或直接删除。

**步骤 3**：在 `circdna.nf/CHANGELOG.md` 中明确记录"合并 circdnalr.nf 资产"的变更条目。

## 四、实施顺序与依赖

> **当前进度总结（2026-07-24）**：Phase 0/1/2/6 的核心代码已实现（未提交）。当前应先完成**收尾验证**（stub 测试 + 提交），再推进 Phase 3-5/7/8。

```
Phase 0: 现有 Bug 修复（circdnalr 分支）✅ 代码已完成
├─ Bug 1: ✅ 已修复（long_read_mapping 已调用 SAMTOOLS_INDEX，emit 3 元组）
├─ Bug 2: ✅ 已修复（flye_pipeline 无 .ext.args，在 conf/modules.config 中设置）
├─ Bug 3: ✅ 已修复（long_read_filtering 已实现 min_read_support 过滤）
├─ Bug 4: ✅ 代码已完成 — test_pacbio_lr.config / test_nanopore_lr.config 已建并注册 profile
│   └─ ❌ 待验证: nextflow run . -stub -profile test_pacbio_lr / test_nanopore_lr
└─ Bug 5: ✅ 已修复（MINIMAP2_ALIGN ext.args 已使用 meta.platform 驱动）

Phase 1: 资产迁移（circdnalr 分支）✅ 已完成
├─ ✅ 从 bio.nf 复制 ecc_finder 4 子模块（map_ont/map_sr/asm_ont/asm_sr）
├─ ✅ 安装 NanoPlot 模块（modules/nf-core/nanoplot/）
├─ ✅ 新建测试配置 test_pacbio_lr.config / test_nanopore_lr.config（本地化测试数据）
├─ ✅ 整合 circdnalr.nf 的 conf/modules.config 中的长读 publishDir 模式
├─ ✅ 整合 circdnalr.nf 的平台自动切换设计（meta.platform 驱动 minimap2）
└─ ❌ 待执行：删除 circdnalr.nf（stub 测试通过后）

Phase 2: 引擎扩展（circdnalr 分支）✅ 已完成
├─ ✅ 新建 ECCFINDER_PIPELINE 子流程（subworkflows/local/eccfinder_pipeline/main.nf）
├─ ✅ 集成到主流程（workflows/circdna.nf，long_read_identifier 默认含 eccfinder）
├─ ✅ 新增 eccfinder_mode 参数（nextflow.config）
├─ ⚠️ 待补充：nextflow_schema.json 中 eccfinder_mode 的 schema 定义
└─ ❌ 待验证：Stub 测试

Phase 3: SV 检测（circdnalr 分支）❌ 待实现
├─ 检查 nf-core 是否有 sniffles2/cutesv/svim 模块
│   ├─ 有：nf-core modules install
│   └─ 无：在 bio.nf 新建分支构建 → 复制到 modules/local/
├─ 检查 quay.io/biocontainers 是否有对应镜像
│   └─ 无：构建并推送至 quay.io/bioinfortools + anaconda.org/yangmingsi
├─ 新建 SV_DETECTION_PIPELINE 子流程
└─ 集成到主流程

Phase 4: CIDER-seq（circdnalr 分支）❌ 待实现
├─ 在 bio.nf 新建分支 ciderseq
├─ 构建 deconcat/tidehunter/circleseeker 三个模块
├─ 检查 biocontainers；无则自构建上传至 quay.io/bioinfortools + anaconda.org/yangmingsi
├─ 复制到 circdna.nf/modules/local/
├─ 新建 CIDERSEQ_PIPELINE 子流程
└─ 集成到主流程

Phase 5: Annotation + 统一输出（circdnalr 分支）❌ 待实现
├─ 在 bio.nf 新建 annotate_eccdna 模块
├─ 复制到 circdna.nf/modules/local/
├─ 编写 bin/unify_eccdna_output.py + bin/annotate_eccdna.py
├─ 新建 ECCDNA_ANNOTATION + UNIFIED_OUTPUT 子流程
└─ 集成到主流程

Phase 6: NanoPlot QC 集成（circdnalr 分支）✅ 已完成
├─ ✅ 修改 long_read_preprocessing/main.nf（已含 NANOPLOT 调用 + skip_long_read_qc）
├─ ✅ 添加 NANOPLOT publishDir 配置（conf/modules.config）
└─ ❌ 待验证：Stub 测试

Phase 7: Metadata 分层系统（circdnalr 分支）❌ 待实现
├─ 新建 metadata/eccDNA/ 目录
├─ 按物种生成 TSV（基于 plan/metadata.csv 的 32 条 eccDNA 记录）
└─ 扩展 schema_input.json

Phase 8: 二代完善（master 分支）❌ 待实现
├─ 安装 fastp 模块
├─ 增加 trimmer 参数
├─ 二代流程接入 UNIFIED_OUTPUT
└─ 版本升至 v3.2.1

收尾工作（circdnalr 分支提交前）❌ 待执行
├─ Stub 测试验证：nextflow run . -stub -profile test_pacbio_lr / test_nanopore_lr
├─ 更新 nextflow_schema.json：补全 eccfinder_mode 等 schema 定义
├─ 更新 conf/base.config：新增模块资源分配
├─ 更新 CHANGELOG.md：v3.2.0 条目
├─ 新建 CHANGES&FIX/20260724.md
├─ 清理临时文件：genome.txt、program.txt、scripts/err、scripts/errors.txt、scripts/reslts 等
├─ 提交 circdnalr 分支改动
└─ 删除 circdnalr.nf/ 目录
```

### 4.1 数据集 → Phase 映射（基于 plan/metadata.csv）

> metadata.csv 共 34 条记录（32 条 eccDNA + 2 条 circRNA）。按 `analysis_role` 和 `experimental_type` 分层映射到各 Phase 的验证数据：

| Phase                       | 数据集 (run_id)         | 物种        | 平台        | experimental_type     | analysis_role    | 用途                            |
| --------------------------- | ----------------------- | ----------- | ----------- | --------------------- | ---------------- | ------------------------------- |
| **Phase 0/1/2 收尾**        | test_pacbio_lr (本地)   | —           | PacBio      | —                     | —                | Stub 测试验证                   |
| **Phase 0/1/2 收尾**        | test_nanopore_lr (本地) | —           | ONT         | —                     | —                | Stub 测试验证                   |
| **Phase 2 ecc_finder**      | ERR6326020              | Arabidopsis | ONT         | eccDNA_enriched       | primary_eccDNA   | ecc_finder benchmark 原始数据集 |
| **Phase 2 ecc_finder**      | SRR24335762             | Arabidopsis | ONT         | mobilome_seq          | primary_eccDNA   | mobilome-seq eccDNA             |
| **Phase 3 SV (PacBio WGS)** | ERR11838731             | Oryza       | PacBio HiFi | WGS_reference         | reference_genome | HiFi WGS SV 验证                |
| **Phase 3 SV (PacBio WGS)** | SRR30359583             | Amaranthus  | PacBio HiFi | WGS_reference         | reference_genome | HiFi WGS SV 验证                |
| **Phase 3 SV (ONT WGS)**    | ERR12723706             | Triticum    | ONT         | WGS_reference         | reference_genome | ONT WGS SV 验证                 |
| **Phase 3 SV (ONT WGS)**    | SRR24335749             | Arabidopsis | ONT         | WGS_reference         | reference_genome | ONT WGS SV 验证                 |
| **Phase 4 CIDER-seq**       | SRR16958693             | Amaranthus  | PacBio      | CIDER-seq             | primary_eccDNA   | CIDER-seq 主测试                |
| **Phase 4 CIDER-seq (RCA)** | SRR26069818             | Alopecurus  | PacBio      | eccDNA_enriched (RCA) | primary_eccDNA   | RCA-PacBio 测试                 |
| **Phase 5 Annotation**      | 任意 primary_eccDNA     | —           | —           | —                     | —                | 提供基因/TE/CDS BED 注释        |
| **Phase 7 Metadata**        | 全部 32 条 eccDNA       | 10 物种     | —           | —                     | —                | 按物种生成 TSV                  |
| **Phase 8 二代 (master)**   | ERR11535564             | Artemisia   | Illumina    | eccDNA_enriched       | primary_eccDNA   | Illumina eccDNA 测试            |
| **Phase 8 二代 (master)**   | CRR3168890              | Cynodon     | Illumina    | eccDNA_enriched       | primary_eccDNA   | Illumina eccDNA 测试            |
| **Phase 8 二代 (master)**   | ERR6004146              | Beta        | Illumina    | eccDNA_enriched       | primary_eccDNA   | circSeq benchmark               |

**eccDNA 数据分层（对应计划第 5 节 Level A-E）**：

| Level | 条件                          | 数据集数量 | 代表 run_id                                                                              |
| ----- | ----------------------------- | ---------- | ---------------------------------------------------------------------------------------- |
| A     | Enriched eccDNA + Illumina    | 9          | ERR10889836, ERR6326000, ERR11535564, CRR3168890, CRR1082975, ERR6004146 等              |
| B     | Enriched eccDNA + PacBio      | 2          | SRR16958693 (CIDER), SRR26069818 (RCA)                                                   |
| C     | Enriched eccDNA + ONT         | 7          | ERR12724336, ERR6326020, SRR24335762, SRR24693334, SRR36603439, SRR28004411, SRR31773424 |
| D     | WGS / genome-skimming         | 10         | ERR11838731, SRR30359583, ERR12723706, ERR11535563 等                                    |
| E     | RNA-seq（不进入 eccDNA 流程） | 2          | ERR10889820, ERR12724677（归 circrna.nf）                                                |

> **A+B+C = 18 条核心 eccDNA 数据**，是 circdna.nf 的主要处理对象。D 类 10 条作为 WGS 参考用于 SV 验证和基因组背景。

## 五、假设与决策

| 决策项                   | 选择                                                                                     | 理由                                                                                                                                 |
| ------------------------ | ---------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------ |
| 分支策略                 | 三代 → circdnalr / 二代 → master                                                         | 与用户既有工作流一致，避免分支冲突                                                                                                   |
| 版本策略                 | 每次合并升版本号                                                                         | 遵循 Semantic Versioning，便于追踪                                                                                                   |
| 容器来源                 | nf-core > quay.io biocontainers > 自构建 quay.io/bioinfortools + anaconda.org/yangmingsi | 遵循用户约定，降低维护成本                                                                                                           |
| 模块来源                 | nf-core > bio.nf 新建后复制 > 直接 local                                                 | 符合用户 bio.nf 中转验证工作流                                                                                                       |
| **ecc_finder 模块来源**  | **从 bio.nf 迁移完整版（4 子模块）**                                                     | bio.nf 有 map_ont/map_sr/asm_ont/asm_sr 完整模块（含测试+testdata），优于 circdnalr.nf 的简化版（无 meta.yml、无测试）               |
| ecc_finder 引擎          | 新增为第四引擎                                                                           | ONT eccDNA-enriched 数据友好，补全三引擎不足                                                                                         |
| **Bug 修复优先级**       | **Phase 0 仅剩 Bug 4（测试配置）**                                                       | v3.1.0 的 5 个已知问题中 4 个（FLED 输入、FLYE ext.args、min_read_support、minimap2 preset）已在当前代码中修复，仅需新建长读测试配置 |
| **minimap2 preset 切换** | **改为 meta.platform 驱动**                                                              | 参考 circdnalr.nf 设计，比 params.protocol 更灵活，支持混合平台 samplesheet                                                          |
| SV 检测工具              | Sniffles2 + cuteSV + SVIM                                                                | 三工具 consensus 提高置信度                                                                                                          |
| CIDER-seq 处理           | 新建独立子流程                                                                           | 实验方法差异显著，需专门处理链路                                                                                                     |
| NanoPlot 集成位置        | LONG_READ_PREPROCESSING 内                                                               | 与既有预处理链路协同                                                                                                                 |
| 统一输出格式             | TSV + BED + FASTA 三件套                                                                 | 适配下游不同分析需求                                                                                                                 |
| Metadata 分层            | 按物种组织 TSV                                                                           | 与计划第 4 节一致，便于多物种管理                                                                                                    |
| circdnalr.nf 处理        | 资产迁移后删除/归档                                                                      | 避免维护两套相似流程                                                                                                                 |

## 六、验证步骤

### 6.1 circdna.nf（circdnalr 分支）

0. **Bug 修复验证（Phase 0）**：

   - **Bug 1** ✅：已验证 `long_read_mapping/main.nf` 调用 `SAMTOOLS_INDEX` 并 emit 3 元组；`fled_pipeline/main.nf` 正确接收 3 元组
   - **Bug 2** ✅：已验证 `flye_pipeline/main.nf` 无 `.ext.args` 赋值；`conf/modules.config` 中 `withName: 'FLYE'` 正确设置 `meta.platform` 驱动的 ext.args
   - **Bug 3** ✅：已验证 `long_read_filtering/main.nf` 调用 `FILTER_ECCDNA_BY_SUPPORT` 模块；`bin/filter_by_read_support.py` 存在
   - **Bug 4** ❌：待验证 — 新建配置后 `nextflow run . -stub -profile test_pacbio_lr` 和 `nextflow run . -stub -profile test_nanopore_lr` 成功
   - **Bug 5** ✅：已验证 `conf/modules.config` 中 `MINIMAP2_ALIGN` 使用 `meta.platform` 驱动

1. **资产迁移验证**：

   - `nf-core modules list` 确认 NanoPlot 已安装
   - `ls modules/local/ecc_finder/` 确认 4 个子模块（map_ont/map_sr/asm_ont/asm_sr）存在且含 main.nf + meta.yml + environment.yml
   - `nextflow run . -entry CIRCDNA -stub -profile test_pacbio_lr` 验证加载

2. **ecc_finder 验证**：

   - `--protocol ont --long_read_identifier eccfinder --eccfinder_mode map --entrypoint cleaned_fastq`
   - 使用 SRR24335762（Arabidopsis ONT eccDNA）小样本测试
   - `--protocol pacbio --long_read_identifier eccfinder --eccfinder_mode both --entrypoint cleaned_fastq`
   - 使用 SRR24335749（Arabidopsis PacBio WGS）小样本测试

3. **SV 检测验证**：

   - `--protocol pacbio --run_sv_for_wgs --sv_identifier sniffles2,cutesv,svim`
   - 使用 ERR11838731（Oryza PacBio HiFi WGS）小样本测试

4. **CIDER-seq 验证**：

   - `--protocol pacbio --ciderseq_enabled`
   - 使用 SRR16958693（Amaranthus CIDER-seq）小样本测试

5. **Annotation + 统一输出验证**：

   - 提供 `--gene_bed --te_bed --cds_bed`
   - 检查 `unified_eccdna.tsv` 字段完整性（19 个字段）

6. **NanoPlot QC 验证**：

   - `--protocol ont --entrypoint cleaned_fastq`
   - 检查 `qc/nanoplot/<sample>/` 目录生成

7. **回归测试**：
   - `--protocol short_read --circle_identifier circle_map_realign` 二代流程不变
   - `--protocol pacbio --long_read_identifier cresil,fled,flye` 三引擎回归

### 6.2 circdna.nf（master 分支）

1. **fastp 验证**：

   - `--trimmer fastp` 二代流程测试
   - `--trimmer trimgalore` 回归测试

2. **统一输出验证**：
   - 检查二代流程的 `unified_eccdna.tsv` 字段

### 6.3 整体一致性验证

- `nf-core lint circdna.nf/` 无致命错误
- `unified_eccdna.tsv` 字段名与下游 Snakemake 整合层期望一致

## 七、修改文件清单

### 7.1 circdna.nf（circdnalr 分支，v3.2.0）

**新增文件**：

- `modules/local/ecc_finder/{map_ont,map_sr,asm_ont,asm_sr}/{main.nf,meta.yml,environment.yml}`（来自 bio.nf，4 个子模块）
- `modules/local/deconcat/{main.nf,meta.yml,environment.yml}`（来自 bio.nf）
- `modules/local/tidehunter/{main.nf,meta.yml,environment.yml}`（来自 bio.nf）
- `modules/local/circleseeker/{main.nf,meta.yml,environment.yml}`（来自 bio.nf）
- `modules/local/annotate_eccdna/{main.nf,meta.yml,environment.yml,templates/annotate_eccdna.py}`（来自 bio.nf）
- `modules/nf-core/nanoplot/`（通过 nf-core 安装）
- `modules/nf-core/sniffles2/`（通过 nf-core 安装或来自 bio.nf）
- `modules/nf-core/cutesv/`（同上）
- `modules/nf-core/svim/`（同上）
- `subworkflows/local/eccfinder_pipeline/{main.nf,meta.yml}`
- `subworkflows/local/sv_detection_pipeline/{main.nf,meta.yml}`
- `subworkflows/local/ciderseq_pipeline/{main.nf,meta.yml}`
- `subworkflows/local/eccdna_annotation/{main.nf,meta.yml}`
- `subworkflows/local/unified_output/{main.nf,meta.yml}`
- `bin/unify_eccdna_output.py`
- `bin/annotate_eccdna.py`
- `conf/test_pacbio_lr.config`
- `conf/test_nanopore_lr.config`
- `metadata/eccDNA/*.tsv`（每个物种一个）

**修改文件**：

- `main.nf`（help 信息）
- `nextflow.config`（新增 params + 注册 test_pacbio_lr / test_nanopore_lr profiles）
- `nextflow_schema.json`（新增 schema）
- `workflows/circdna.nf`（新增 includes、长读分支扩展）
- `subworkflows/local/long_read_preprocessing/main.nf`（集成 NanoPlot）
- `conf/base.config`（资源分配）
- `conf/modules.config`（新增 publishDir 配置）
- `assets/schema_input.json`（增加 metadata_tsv 列）
- `CHANGELOG.md`（v3.2.0 条目）
- `CHANGES&FIX/20260724.md`
- `nextflow.config` 中 `manifest.version`（→ 3.2.0）
- `modules.json`（更新 nf-core 模块版本）

> **注**：Phase 0 中 Bug 1/2/3/5 已在当前代码中修复，无需修改 `long_read_mapping/main.nf`、`fled_pipeline/main.nf`、`flye_pipeline/main.nf`、`long_read_filtering/main.nf`、`conf/modules.config`（MINIMAP2_ALIGN 部分）。仅 Bug 4 需新建测试配置文件。

### 7.2 circdna.nf（master 分支，v3.2.1）

**新增文件**：

- `modules/nf-core/fastp/`（通过 nf-core 安装）
- `CHANGES&FIX/20260724_master.md`

**修改文件**：

- `workflows/circdna.nf`（增加 trimmer 分支、统一输出）
- `nextflow.config`（新增 `trimmer` 参数）
- `nextflow_schema.json`（新增 `trimmer` schema）
- `CHANGELOG.md`（v3.2.1 条目）
- `manifest.version`（→ 3.2.1）

### 7.3 circdnalr.nf 处理

- 复制资产到 `circdna.nf` 的 `circdnalr` 分支
- 归档为 `circdnalr.nf.archive/` 或删除

### 7.4 bio.nf 新建分支与模块（外部依赖）

| bio.nf 分支 | 模块                       | 目的                   |
| ----------- | -------------------------- | ---------------------- |
| `ciderseq`  | `modules/deconcat/`        | PacBio concatemer 去除 |
| `ciderseq`  | `modules/tidehunter/`      | 串联重复检测           |
| `ciderseq`  | `modules/circleseeker/`    | 环状分子识别           |
| `annotate`  | `modules/annotate_eccdna/` | eccDNA 注释整合        |

**容器构建清单**（仅在 biocontainer 缺失时）：

| 工具         | quay.io 镜像                                   | anaconda 包                                      |
| ------------ | ---------------------------------------------- | ------------------------------------------------ |
| deconcat     | `quay.io/bioinfortools/deconcat:<version>`     | `anaconda.org/yangmingsi/deconvert:<version>`    |
| tidehunter   | `quay.io/bioinfortools/tidehunter:<version>`   | `anaconda.org/yangmingsi/tidehunter:<version>`   |
| circleseeker | `quay.io/bioinfortools/circleseeker:<version>` | `anaconda.org/yangmingsi/circleseeker:<version>` |

## 八、风险与回退

| 风险                                        | 回退方案                                                        |
| ------------------------------------------- | --------------------------------------------------------------- |
| ecc_finder 容器构建失败                     | 第一版仅启用 map 模式，asm 后续迭代                             |
| SV 工具 nf-core 模块不存在                  | 在 bio.nf 新建分支后复制                                        |
| CIDER-seq 工具无 biocontainer               | 自构建上传至 quay.io/bioinfortools + anaconda.org/yangmingsi    |
| 统一输出脚本字段映射错误                    | 提供 `--output_unified_eccdna false` 关闭                       |
| master/circdnalr 分支合并冲突               | 优先保留 circdnalr 分支的长读改动，master 分支仅保留二代改动    |
| deconcat/tidehunter/circleseeker 工具不可用 | 第一版先交付 CReSIL+FLED+Flye+ecc_finder+SV，CIDER-seq 后续迭代 |

## 九、后续 Snakemake 整合层（Phase 7 之后）

按计划第 14-16 节，最终在 Snakemake 中消费 `unified_eccdna.tsv`：

```
integration/
├── Snakefile
├── config/
│   ├── config.yaml
│   └── eccDNA.yaml
├── rules/
│   ├── eccDNA.smk
│   ├── annotation.smk
│   ├── overlap.smk
│   └── integration.smk
├── scripts/
│   ├── map_circRNA_eccDNA.py
│   ├── annotate_TE.py
│   └── integrate_results.py
└── results/
```

输入：

- `circdna.nf` → `unified_eccdna.tsv` + `unified_eccdna.bed` + `unified_eccdna.fa`
- circrna.nf 输出（由其他任务负责）
- 基因组注释、TE 注释、表达矩阵

输出：`Circular DNA-RNA atlas`（按计划第 18 节的三类关联）。

---

**计划生成时间**：2026-07-24
**分支策略**：circdnalr 分支（三代）+ master 分支（二代）
**版本目标**：circdna.nf v3.2.0/v3.2.1
**范围**：仅 circdna.nf 与 circdnalr.nf；circrna.nf 由其他任务处理
