# ECCSPLORER / ecc_finder 原子化拆分方案（slim 版本）

> 分析日期: 2026-08-07 | 分支: ECCsplorer | 状态: 执行阶段
>
> **核心目标**：将 ECCsplorer 和 ecc_finder 的外部工具调用替换为 nf-core 标准模块，仅保留独有的结果处理逻辑。
> **策略**：旧模块保留不动，新 slim 版本与原版并行共存，通过参数 `circle_identifier` 区分。
> **构建顺序**：先在 bio.nf 构建模块 → 再接入 circdna.nf → 最后 test_local 测试。
> **镜像策略**：优先使用 quay.io/biocontainers 和用户频道 bioinfortools 的现有镜像，无需从头构建。

---

## 1. 源码级拆解分析

### 1.1 ECCsplorer Map 模式 — 16 步内部调用链

基于 [ECCsplorer.py](file:///Users/siyangming/nextflow_nf_core/ECCsplorer/ECCsplorer.py)、[eccMapper.py](file:///Users/siyangming/nextflow_nf_core/ECCsplorer/lib/eccMapper.py)、[eccPrepare.py](file:///Users/siyangming/nextflow_nf_core/ECCsplorer/lib/eccPrepare.py)、[config.py](file:///Users/siyangming/nextflow_nf_core/ECCsplorer/lib/config.py) 源码：

| 步骤 | 操作                           | 原工具                                   | slim 版替代                                              | 属性         |
| ---- | ------------------------------ | ---------------------------------------- | -------------------------------------------------------- | ------------ |
| 1    | FASTQ adapter trimming         | Trimmomatic                              | 跳过（复用 BAM_PREPROCESSING 上游 QC）                   | **跳过**     |
| 2    | FASTQ→FASTA 格式转换           | seqtk seq -A                             | `nf-core seqtk/seq` ✓                                    | **nf-core**  |
| 3    | segemehl 索引构建              | segemehl.x -x -d                         | `nf-core segemehl/index` ✓                               | **nf-core**  |
| 4    | segemehl split-read 比对       | segemehl.x --splits --MEOP               | `nf-core segemehl/align` ✓                               | **nf-core**  |
| 5    | Split-read 检测 + 区域合并     | haarz.x split + bedtools merge           | `bio.nf/modules/haarz/` + `nf-core bedtools/merge`       | **混合**     |
| 6    | SAM→BED 转换                   | samtools→sort→bamtobed                   | `nf-core samtools/view` + `sort` + `bedtools/bamtobed` ✓ | **nf-core**  |
| 7    | Discordant read 提取           | samtools view -f/-G                      | `nf-core samtools/view` + `bedtools/bamtobed` ✓          | **nf-core**  |
| 8    | Discordant read 区域检测       | bedtools genomecov + merge               | `nf-core bedtools/genomecov` + `bedtools/merge` ✓        | **nf-core**  |
| 9    | 窗口覆盖度 + **峰值检测**      | bedtools coverage + **scipy find_peaks** | `nf-core bedtools/coverage` + **自定义 Python**          | **专有逻辑** |
| 10   | 候选区域提取（SR∩DR∩coverage） | bedtools intersect                       | `nf-core bedtools/intersect` + **自定义 Python**         | **专有逻辑** |
| 11   | 候选序列提取                   | bedtools getfasta                        | `nf-core bedtools/getfasta` ✓                            | **nf-core**  |
| 12   | BLAST 注释                     | blastn                                   | `nf-core blast/blastn` ✓                                 | **nf-core**  |
| 13   | 候选区 per-base 覆盖度         | bedtools coverage -d                     | `nf-core bedtools/coverage` ✓                            | **nf-core**  |
| 14   | RPM 归一化 + Fold Enrichment   | R/pyRserve                               | **自定义 R 脚本**                                        | **专有逻辑** |
| 15   | Manhattan plot + 候选可视化    | R (ggplot2)                              | **自定义 R 脚本**                                        | **专有逻辑** |
| 16   | HTML 报告生成                  | eccHTML_templates.py                     | **自定义 Python 脚本**                                   | **专有逻辑** |

### 1.2 ecc_finder MAP_SR_slim — 简化步骤（复用 BAM_PREPROCESSING 输出）

| 步骤 | 操作                             | 原工具                | slim 版替代                                            | 属性         |
| ---- | -------------------------------- | --------------------- | ------------------------------------------------------ | ------------ |
| 1    | Read QC                          | fastp                 | **跳过**（BAM_PREPROCESSING 上游已 QC）                | **跳过**     |
| 2    | BWA mem 比对                     | bwa mem               | **跳过**（复用 `BAM_PREPROCESSING` 产出的 sorted BAM） | **核心优化** |
| 3    | SAM name-sort                    | pysam.sort -n         | `nf-core samtools/sort -n` ✓                           | **nf-core**  |
| 4    | 富集位点检测                     | Genrich               | `bio.nf/modules/genrich/`（bioconda 现有）             | **bio.nf**   |
| 5    | Split-read 检测                  | tidehunter            | `bio.nf/modules/tidehunter/`（bioconda 现有）          | **bio.nf**   |
| 6    | 候选合并 + 边界精修 + 置信度评分 | pybedtools + 内部算法 | **自定义 Python 脚本**                                 | **专有逻辑** |

### 1.3 ecc_finder ASM_SR_slim — 简化步骤

| 步骤 | 操作            | 原工具          | slim 版替代            | 属性         |
| ---- | --------------- | --------------- | ---------------------- | ------------ |
| 1    | Read QC         | fastp           | **跳过**               | **跳过**     |
| 2    | De novo 组装    | unicycler -1 -2 | `nf-core unicycler` ✓  | **nf-core**  |
| 3    | eccDNA 候选识别 | 自定义算法      | **自定义 Python 脚本** | **专有逻辑** |

---

## 2. 容器与 Conda 环境完整矩阵

> **策略**：
>
> - 标准工具(genrich/tidehunter/haarz/segemehl/unicycler/blast) → 使用 quay.io/biocontainers 现有镜像
> - 专有逻辑(ECCsplorer_slim/ecc_finder_slim) → **构建最小化镜像**，仅含 Python/R 脚本依赖，不含任何外部工具（bwa/segemehl/unicycler 等已由 nf-core 模块替代）
> - Conda 包推送至 `yangmingsi` 频道，Docker 镜像推送至 `quay.io/bioinfortools`

### 2.1 标准工具模块（全部已有镜像，无需构建）

| 模块               | bioconda 包                  | Docker 镜像 (quay.io)                                         | Apptainer/Singularity (galaxy depot)                                       |
| ------------------ | ---------------------------- | ------------------------------------------------------------- | -------------------------------------------------------------------------- |
| `genrich`          | `bioconda::genrich=0.6.1`    | `quay.io/biocontainers/genrich:0.6.1--h577a1d6_5` (~6 MB)     | `https://depot.galaxyproject.org/singularity/genrich:0.6.1--h577a1d6_5`    |
| `tidehunter`       | `bioconda::tidehunter=1.5.6` | `quay.io/biocontainers/tidehunter:1.5.6--h7f5d12c_0` (~15 MB) | `https://depot.galaxyproject.org/singularity/tidehunter:1.5.6--h7f5d12c_0` |
| `haarz` (segemehl) | `bioconda::segemehl=0.3.4`   | `quay.io/biocontainers/segemehl:0.3.4--hc2ea5fd_5` (~50 MB)   | `https://depot.galaxyproject.org/singularity/segemehl:0.3.4--hc2ea5fd_5`   |
| `segemehl/index`   | `bioconda::segemehl=0.3.4`   | 同上（已有 nf-core 模块）                                     | 同上                                                                       |
| `segemehl/align`   | `bioconda::segemehl=0.3.4`   | 同上（已有 nf-core 模块）                                     | 同上                                                                       |
| `unicycler`        | `bioconda::unicycler=0.5.1`  | `community.wave.seqera.io/library/unicycler:0.5.1`            | seqera container registry                                                  |
| `blast/blastn`     | `bioconda::blast=2.16.0`     | nf-core 标准镜像                                              | 同上                                                                       |

### 2.2 专有逻辑模块 — 需构建最小化 Slim 镜像

> **与原版的本质区别**：slim 镜像**不含任何外部生物信息学工具**（bwa, segemehl, unicycler, bedtools, samtools, blast 等），这些已全部由 nf-core 标准模块接管。仅保留 Python 数值计算库 + R 绑图库 + 自定义脚本。

#### 2.2.1 ECCsplorer_slim 镜像

| 属性         | 值                                                     |
| ------------ | ------------------------------------------------------ |
| **镜像名**   | `quay.io/bioinfortools/eccsplorer_slim:1.0.0`          |
| **Conda 包** | `yangmingsi::eccsplorer_slim=1.0.0`                    |
| **预估大小** | ~400 MB（vs 原版 ~1.5 GB，瘦身 73%）                   |
| **依赖**     | Python 3.9+, numpy, scipy, biopython, R 4.x, r-ggplot2 |

**只含的依赖（最小化）**：

```yaml
# environment.yml
name: eccsplorer_slim
channels:
  - conda-forge
  - bioconda
dependencies:
  - python>=3.9
  - numpy
  - scipy
  - biopython
  - r-base>=4.0
  - r-ggplot2
```

**不含的依赖（已由 nf-core 模块替代）**：
~~segemehl, haarz, samtools, bedtools, blast, trimmomatic, seqtk, Rserve~~

#### 2.2.2 ecc_finder_slim 镜像

| 属性         | 值                                                           |
| ------------ | ------------------------------------------------------------ |
| **镜像名**   | `quay.io/bioinfortools/ecc_finder_slim:1.0.0`                |
| **Conda 包** | `yangmingsi::ecc_finder_slim=1.0.0`                          |
| **预估大小** | ~250 MB（vs 原版 ~2 GB，瘦身 87%）                           |
| **依赖**     | Python 3.9+, numpy, pandas, matplotlib, pybedtools, bedtools |

**只含的依赖（最小化）**：

```yaml
# environment.yml
name: ecc_finder_slim
channels:
  - conda-forge
  - bioconda
dependencies:
  - python>=3.9
  - numpy
  - pandas
  - matplotlib
  - pybedtools
  - bedtools
```

**不含的依赖（已由 nf-core 模块替代）**：
~~bwa, minimap2, unicycler, fastp, seqtk, cd-hit, tidehunter, genrich, samtools~~

> **注**：ecc_finder_slim 仍需要 `pybedtools` + `bedtools` 的原因是 `merge_score` 脚本在内部使用 pybedtools 做 BED 操作。如果后续将这些操作也改为 nf-core bedtools 模块调用，可进一步瘦身至 ~100 MB。

### 2.3 镜像大小对比

| 镜像                  | 大小        | 瘦身比   | 来源               |
| --------------------- | ----------- | -------- | ------------------ |
| genrich               | ~6 MB       | —        | biocontainers 现有 |
| tidehunter            | ~15 MB      | —        | biocontainers 现有 |
| segemehl (haarz)      | ~50 MB      | —        | biocontainers 现有 |
| **eccsplorer (原版)** | ~1.5 GB     | 基线     | bioinfortools 现有 |
| **eccsplorer_slim**   | **~400 MB** | **↓73%** | **需新建**         |
| **ecc_finder (原版)** | ~2 GB       | 基线     | bioinfortools 现有 |
| **ecc_finder_slim**   | **~250 MB** | **↓87%** | **需新建**         |

---

## 2.4 Slim 镜像构建规格

### 2.4.1 ECCsplorer_slim Conda Recipe

目录：`bio.nf/modules/eccsplorer_slim/conda-recipe/`

```yaml
# meta.yaml
{% set name = "eccsplorer_slim" %}
{% set version = "1.0.0" %}

package:
  name: {{ name }}
  version: {{ version }}

source:
  path: ../bin  # 指向自定义 Python/R 脚本目录

build:
  number: 0
  noarch: python
  script: "{{ PYTHON }} -m pip install . --no-deps -vv"

requirements:
  host:
    - python>=3.9
    - pip
  run:
    - python>=3.9
    - numpy
    - scipy
    - biopython
    - r-base>=4.0
    - r-ggplot2

test:
  imports:
    - numpy
    - scipy
    - Bio
  commands:
    - python -c "import numpy; import scipy; from Bio import SeqIO; print('OK')"
    - Rscript -e "library(ggplot2); print('OK')"

about:
  home: https://quay.io/bioinfortools/eccsplorer_slim
  license: GPL-3.0
  summary: Minimal ECCsplorer analysis scripts (no external bioinfo tools)
```

```bash
# build.sh
#!/bin/bash
mkdir -p ${PREFIX}/bin
cp -r ${SRC_DIR}/bin/* ${PREFIX}/bin/
chmod +x ${PREFIX}/bin/*.py ${PREFIX}/bin/*.R
```

### 2.4.2 ECCsplorer_slim Dockerfile

```dockerfile
# Dockerfile
FROM condaforge/mambaforge:latest

# Install minimal dependencies via conda
COPY environment.yml /tmp/environment.yml
RUN mamba env create -f /tmp/environment.yml && \
    mamba clean -afy

# Copy custom scripts
COPY bin/ /opt/eccsplorer_slim/bin/

ENV PATH="/opt/conda/envs/eccsplorer_slim/bin:/opt/eccsplorer_slim/bin:${PATH}"
ENV CONDA_DEFAULT_ENV=eccsplorer_slim

# Activate conda env on run
SHELL ["conda", "run", "-n", "eccsplorer_slim", "/bin/bash", "-c"]
```

### 2.4.3 ecc_finder_slim Conda Recipe

目录：`bio.nf/modules/ecc_finder_slim/conda-recipe/`

```yaml
# meta.yaml
{% set name = "ecc_finder_slim" %}
{% set version = "1.0.0" %}

package:
  name: {{ name }}
  version: {{ version }}

source:
  path: ../bin

build:
  number: 0
  noarch: python
  script: "{{ PYTHON }} -m pip install . --no-deps -vv"

requirements:
  host:
    - python>=3.9
    - pip
  run:
    - python>=3.9
    - numpy
    - pandas
    - matplotlib
    - pybedtools
    - bedtools

test:
  imports:
    - numpy
    - pandas
    - matplotlib
    - pybedtools
  commands:
    - python -c "import numpy, pandas, matplotlib, pybedtools; print('OK')"
    - bedtools --version

about:
  home: https://quay.io/bioinfortools/ecc_finder_slim
  license: MIT
  summary: Minimal ecc_finder analysis scripts (no external bioinfo tools)
```

### 2.4.4 ecc_finder_slim Dockerfile

```dockerfile
# Dockerfile
FROM condaforge/mambaforge:latest

COPY environment.yml /tmp/environment.yml
RUN mamba env create -f /tmp/environment.yml && \
    mamba clean -afy

COPY bin/ /opt/ecc_finder_slim/bin/

ENV PATH="/opt/conda/envs/ecc_finder_slim/bin:/opt/ecc_finder_slim/bin:${PATH}"
ENV CONDA_DEFAULT_ENV=ecc_finder_slim

SHELL ["conda", "run", "-n", "ecc_finder_slim", "/bin/bash", "-c"]
```

### 2.4.5 构建与推送命令

```bash
# === ECCsplorer_slim ===
# 1. 构建 conda 包
cd bio.nf/modules/eccsplorer_slim/
conda build conda-recipe/

# 2. 本地验证
conda create -n test_slim --use-local eccsplorer_slim=1.0.0
conda activate test_slim
python -c "import numpy, scipy; from Bio import SeqIO; print('eccsplorer_slim OK')"

# 3. 推送至 yangmingsi 频道
anaconda login
anaconda upload /path/to/conda-bld/noarch/eccsplorer_slim-1.0.0-py_0.conda

# 4. 构建 Docker 镜像
docker build -t quay.io/bioinfortools/eccsplorer_slim:1.0.0 .

# 5. 验证 Docker 镜像
docker run --rm quay.io/bioinfortools/eccsplorer_slim:1.0.0 \
    python -c "import numpy, scipy; from Bio import SeqIO; print('OK')"

# 6. 推送至 quay.io
docker login quay.io -u bioinfortools
docker push quay.io/bioinfortools/eccsplorer_slim:1.0.0

# === ecc_finder_slim ===
# 同上流程，替换路径和镜像名
cd bio.nf/modules/ecc_finder_slim/
conda build conda-recipe/
anaconda upload /path/to/conda-bld/noarch/ecc_finder_slim-1.0.0-py_0.conda
docker build -t quay.io/bioinfortools/ecc_finder_slim:1.0.0 .
docker push quay.io/bioinfortools/ecc_finder_slim:1.0.0
```

### 2.4.6 模块 main.nf 中的容器引用（更新为 slim 镜像）

ECCsplorer 专有逻辑模块的容器指令：

```nextflow
container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/eccsplorer_slim:1.0.0--0' :
    'quay.io/bioinfortools/eccsplorer_slim:1.0.0' }"
```

ecc_finder 专有逻辑模块的容器指令：

```nextflow
container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/ecc_finder_slim:1.0.0--0' :
    'quay.io/bioinfortools/ecc_finder_slim:1.0.0' }"
```

---

## 3. 新旧共存架构

```
circdna.nf/
├── modules/
│   ├── nf-core/                       # 现有 nf-core 模块（无需新增）
│   │   ├── blast/blastn/              # nf-core modules install（如需）
│   │   └── ... (existing)
│   └── local/
│       ├── eccsplorer/                # ★ 保留：原版 ECCsplorer（黑盒调用）
│       └── ecc_finder/                # ★ 保留：原版 ecc_finder（黑盒调用）
│
├── subworkflows/local/
│   ├── eccsplorer_pipeline/           # ★ 保留：原版 ECCSPLORER_PIPELINE
│   ├── ecc_finder_pipeline/           # ★ 保留：原版 ECC_FINDER_PIPELINE
│   ├── eccsplorer_slim_pipeline/      # ★ 新增
│   └── ecc_finder_slim_pipeline/      # ★ 新增

bio.nf/
├── modules/
│   ├── genrich/                       # ★ 新增 → quay.io/biocontainers/genrich
│   ├── tidehunter/                    # ★ 新增 → quay.io/biocontainers/tidehunter
│   ├── haarz/                         # ★ 新增 → quay.io/biocontainers/segemehl (复用)
│   ├── eccsplorer_slim/              # ★ 新增 → quay.io/bioinfortools/eccsplorer_slim (需构建)
│   │   ├── conda-recipe/             #   meta.yaml + build.sh
│   │   ├── Dockerfile
│   │   ├── environment.yml
│   │   ├── bin/                       #   自定义 Python/R 脚本
│   │   ├── peak_detect/
│   │   ├── candidate_extract/
│   │   ├── coverage_profile/
│   │   ├── normalize/
│   │   ├── visualize/
│   │   └── html_report/
│   └── ecc_finder_slim/              # ★ 新增 → quay.io/bioinfortools/ecc_finder_slim (需构建)
│       ├── conda-recipe/             #   meta.yaml + build.sh
│       ├── Dockerfile
│       ├── environment.yml
│       ├── bin/                       #   自定义 Python 脚本
│       ├── merge_score/
│       └── asm_filter/
```

### 参数区分

在 [test_local.config](file:///Users/siyangming/nextflow_nf_core/circdna.nf/conf/test_local.config#L26-27) 中：

```groovy
// 原版（黑盒模式）— 保留
// circle_identifier   = 'eccsplorer_map,ecc_finder_map_sr,ecc_finder_asm_sr'

// slim 版（原子化模式）— 新增，用于测试新的轻量版流程
circle_identifier   = 'eccsplorer_map_slim,ecc_finder_map_sr_slim,ecc_finder_asm_sr_slim'
```

---

## 4. 数据流（slim 版）

### 4.1 ECCSPLORER_slim 数据流

```
SEGEMEHL_INDEX [nf-core] → .idx
     ↓
SEGEMEHL_ALIGN [nf-core] → .sam + .sngl.bed
     │
     ├─→ HAARZ [bio.nf, biocontainers/segemehl] → SR.bed
     ├─→ SAMTOOLS_VIEW→SORT→BAMTOBED [nf-core] → all.bed
     ├─→ SAMTOOLS_VIEW(-f/-G)→BAMTOBED [nf-core] → DR.bed
     │
     ├─→ BEDTOOLS_GENOMECOV→MERGE [nf-core] → DR_regions.bed
     ├─→ BEDTOOLS_MAKEWINDOWS [nf-core] → window.bed
     ├─→ BEDTOOLS_COVERAGE [nf-core] → coverage.tsv
     │
     ├─→ ECCSPLORER_PEAK_DETECT [bio.nf, bioinfortools/eccsplorer] → peak_regions.bed
     ├─→ ECCSPLORER_CANDIDATE_EXTRACT [bio.nf, bioinfortools/eccsplorer] → candidates.bed
     ├─→ BEDTOOLS_GETFASTA [nf-core] → sequences.fasta
     ├─→ BLAST_BLASTN [nf-core] → annotation
     ├─→ ECCSPLORER_COVERAGE_PROFILE [bio.nf, bioinfortools/eccsplorer]
     ├─→ ECCSPLORER_NORMALIZE [bio.nf, bioinfortools/eccsplorer] → RPM + fold enrichment
     ├─→ ECCSPLORER_VISUALIZE [bio.nf, bioinfortools/eccsplorer] → plots
     └─→ ECCSPLORER_HTML_REPORT [bio.nf, bioinfortools/eccsplorer] → HTML
```

### 4.2 ECC_FINDER_slim 数据流

```
BAM_PREPROCESSING.out.bam_sorted  (已比对好的 BAM — 跳过 fastp + bwa!)
     │
  SAMTOOLS_SORT (-n) [nf-core] → name-sorted BAM
     │
     ├─→ GENRICH [bio.nf, biocontainers/genrich] → enrichment_sites.bed
     ├─→ TIDEHUNTER [bio.nf, biocontainers/tidehunter] → split_reads.bed
     │
     └─→ ECC_FINDER_MERGE_SCORE [bio.nf, bioinfortools/ecc_finder] → candidates.csv + candidates.fasta

Raw FASTQ
     │
  UNICYCLER [nf-core] → assembly.fasta
     │
  ECC_FINDER_ASM_FILTER [bio.nf, bioinfortools/ecc_finder] → eccDNA_asm.fasta
```

---

## 5. 实施步骤

### 阶段 0：构建 Slim 镜像（先决条件）

> 在执行任何模块构建之前，需先完成两个 slim 镜像的 conda 包 + Docker 镜像构建。

| 步骤 | 操作                                          | 详细规格                                                                                                                       |
| ---- | --------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------ |
| 0.1  | 编写 ECCsplorer_slim 专有 Python/R 脚本       | 从 [ECCsplorer 源码](file:///Users/siyangming/nextflow_nf_core/ECCsplorer/) 提取 → [§2.4.1](#241-eccsplorer_slim-conda-recipe) |
| 0.2  | 构建 `eccsplorer_slim` conda 包 + Docker 镜像 | [§2.4.5](#245-构建与推送命令)                                                                                                  |
| 0.3  | 编写 ecc_finder_slim 专有 Python 脚本         | 从 [ecc_finder 源码](https://github.com/njaupan/ecc_finder) 提取 → [§2.4.3](#243-ecc_finder_slim-conda-recipe)                 |
| 0.4  | 构建 `ecc_finder_slim` conda 包 + Docker 镜像 | [§2.4.5](#245-构建与推送命令)                                                                                                  |

### 阶段 1：bio.nf 模块构建（不涉及 circdna.nf）

#### 1.1 构建标准工具模块（全部使用现有 bioconda 镜像）

| 步骤 | 模块           | 位置                                       | environment.yml                          | container (3引擎支持)                                                                                                                                |
| ---- | -------------- | ------------------------------------------ | ---------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------- |
| 1.1a | `GENRICH`      | `bio.nf/modules/genrich/`                  | `bioconda::genrich=0.6.1`                | docker: `quay.io/biocontainers/genrich:0.6.1--h577a1d6_5` / apptainer: `https://depot.galaxyproject.org/singularity/genrich:0.6.1--h577a1d6_5`       |
| 1.1b | `TIDEHUNTER`   | `bio.nf/modules/tidehunter/`               | `bioconda::tidehunter=1.5.6`             | docker: `quay.io/biocontainers/tidehunter:1.5.6--h7f5d12c_0` / apptainer: `https://depot.galaxyproject.org/singularity/tidehunter:1.5.6--h7f5d12c_0` |
| 1.1c | `HAARZ`        | `bio.nf/modules/haarz/`                    | `bioconda::segemehl=0.3.4`（haarz 内置） | docker: `quay.io/biocontainers/segemehl:0.3.4--hc2ea5fd_5` / apptainer: `https://depot.galaxyproject.org/singularity/segemehl:0.3.4--hc2ea5fd_5`     |
| 1.1d | `blast/blastn` | `circdna.nf/modules/nf-core/blast/blastn/` | `nf-core modules install blast/blastn`   | nf-core 标准镜像                                                                                                                                     |

每个模块的标准模板（以 genrich 为例）：

```nextflow
// bio.nf/modules/genrich/main.nf
process GENRICH {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/genrich:0.6.1--h577a1d6_5' :
        'quay.io/biocontainers/genrich:0.6.1--h577a1d6_5' }"

    input:
    tuple val(meta), path(bam)       // name-sorted BAM

    output:
    tuple val(meta), path("*.bed"), emit: bed
    path "versions.yml"            , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    Genrich -t ${bam} -o ${prefix}.site -v $args
    cut -f1-3 ${prefix}.site > ${prefix}.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        genrich: \$(Genrich 2>&1 | head -1 | sed 's/.*Genrich v//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bed
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        genrich: 0.6.1
    END_VERSIONS
    """
}
```

```yaml
# bio.nf/modules/genrich/environment.yml
name: genrich
channels:
  - conda-forge
  - bioconda
dependencies:
  - bioconda::genrich=0.6.1
```

#### 1.2 构建 ECCsplorer 专有逻辑模块（bio.nf，复用 bioinfortools/eccsplorer 镜像）

所有 6 个模块使用同一个容器 `quay.io/bioinfortools/eccsplorer:2022.01.1.1`：

| 步骤 | Process 名                     | 位置                                                | 输入 → 输出                                                  | 核心逻辑来源                                                                                                                               |
| ---- | ------------------------------ | --------------------------------------------------- | ------------------------------------------------------------ | ------------------------------------------------------------------------------------------------------------------------------------------ |
| 1.2a | `ECCSPLORER_PEAK_DETECT`       | `bio.nf/modules/eccsplorer_slim/peak_detect/`       | `coverage.tsv` → `peak_regions.bed`                          | [eccMapper.py L44-107](file:///Users/siyangming/nextflow_nf_core/ECCsplorer/lib/eccMapper.py#L44-L107) `get_rough_coverage` + `find_peaks` |
| 1.2b | `ECCSPLORER_CANDIDATE_EXTRACT` | `bio.nf/modules/eccsplorer_slim/candidate_extract/` | `SR.bed + peak_all.bed + peak_DR.bed` → `candidates.bed`     | [eccMapper.py L450-513](file:///Users/siyangming/nextflow_nf_core/ECCsplorer/lib/eccMapper.py#L450-L513) `extract_candidate_regions`       |
| 1.2c | `ECCSPLORER_COVERAGE_PROFILE`  | `bio.nf/modules/eccsplorer_slim/coverage_profile/`  | `candidates.bed + all.bed + ref.fa` → per-candidate coverage | [eccMapper.py L146-240](file:///Users/siyangming/nextflow_nf_core/ECCsplorer/lib/eccMapper.py#L146-L240) `analyze_candidate_region`        |
| 1.2d | `ECCSPLORER_NORMALIZE`         | `bio.nf/modules/eccsplorer_slim/normalize/`         | 原始覆盖度 + mapped bases → 归一化 RPM + fold enrichment     | [eccDNA_Rcodes.py](file:///Users/siyangming/nextflow_nf_core/ECCsplorer/lib/eccDNA_Rcodes.py) R 函数                                       |
| 1.2e | `ECCSPLORER_VISUALIZE`         | `bio.nf/modules/eccsplorer_slim/visualize/`         | 归一化数据 → PNG plots                                       | [eccDNA_Rcodes.py](file:///Users/siyangming/nextflow_nf_core/ECCsplorer/lib/eccDNA_Rcodes.py) R 函数                                       |
| 1.2f | `ECCSPLORER_HTML_REPORT`       | `bio.nf/modules/eccsplorer_slim/html_report/`       | 候选数据 → HTML report                                       | [eccHTML_templates.py](file:///Users/siyangming/nextflow_nf_core/ECCsplorer/lib/eccHTML_templates.py) HTML 模板                            |

#### 1.3 构建 ecc_finder 专有逻辑模块（bio.nf，复用 bioinfortools/ecc_finder 镜像）

| 步骤 | Process 名               | 位置                                          | 输入 → 输出                                                                       | 核心逻辑来源                                                                                    |
| ---- | ------------------------ | --------------------------------------------- | --------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------- |
| 1.3a | `ECC_FINDER_MERGE_SCORE` | `bio.nf/modules/ecc_finder_slim/merge_score/` | `enrichment.bed + split_reads.bed + ref.fa` → `candidates.csv + candidates.fasta` | [ecc_finder map-sr.py](https://github.com/njaupan/ecc_finder/blob/main/map-sr.py) 合并+评分逻辑 |
| 1.3b | `ECC_FINDER_ASM_FILTER`  | `bio.nf/modules/ecc_finder_slim/asm_filter/`  | `assembly.fasta` → `eccDNA_asm.fasta`                                             | [ecc_finder asm-sr.py](https://github.com/njaupan/ecc_finder/blob/main/asm-sr.py) 过滤逻辑      |

---

### 阶段 2：circdna.nf slim 子工作流构建

#### 2.1 创建 `eccsplorer_slim_pipeline`

文件：`circdna.nf/subworkflows/local/eccsplorer_slim_pipeline/main.nf`

```nextflow
//
// ECCSPLORER_slim — 原子化 eccDNA 检测子工作流
// 使用 nf-core 标准模块 + bio.nf 自建模块，替代原版 ECCsplorer 黑盒调用
//

include { SEGEMEHL_INDEX    } from '../../../modules/nf-core/segemehl/index/main'
include { SEGEMEHL_ALIGN    } from '../../../modules/nf-core/segemehl/align/main'
include { HAARZ             } from '../../../../bio.nf/modules/haarz/main'
include { SAMTOOLS_VIEW as SAMTOOLS_VIEW_SAM2BAM } from '../../../modules/nf-core/samtools/view/main'
include { SAMTOOLS_SORT     } from '../../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX    } from '../../../modules/nf-core/samtools/index/main'
include { BEDTOOLS_BAMTOBED } from '../../../modules/nf-core/bedtools/bamtobed/main'
include { BEDTOOLS_GENOMECOV } from '../../../modules/nf-core/bedtools/genomecov/main'
include { BEDTOOLS_MERGE    } from '../../../modules/nf-core/bedtools/merge/main'
include { BEDTOOLS_MAKEWINDOWS } from '../../../modules/nf-core/bedtools/makewindows/main'
include { BEDTOOLS_COVERAGE } from '../../../modules/nf-core/bedtools/coverage/main'
include { BEDTOOLS_GETFASTA } from '../../../modules/nf-core/bedtools/getfasta/main'
include { ECCSPLORER_PEAK_DETECT } from '../../../../bio.nf/modules/eccsplorer_slim/peak_detect/main'
include { ECCSPLORER_CANDIDATE_EXTRACT } from '../../../../bio.nf/modules/eccsplorer_slim/candidate_extract/main'
include { ECCSPLORER_COVERAGE_PROFILE } from '../../../../bio.nf/modules/eccsplorer_slim/coverage_profile/main'
include { ECCSPLORER_NORMALIZE } from '../../../../bio.nf/modules/eccsplorer_slim/normalize/main'
include { ECCSPLORER_VISUALIZE } from '../../../../bio.nf/modules/eccsplorer_slim/visualize/main'
include { ECCSPLORER_HTML_REPORT } from '../../../../bio.nf/modules/eccsplorer_slim/html_report/main'

workflow ECCSPLORER_SLIM_PIPELINE {
    take:
    reads           // channel: [meta, [r1, r2]]
    fasta_meta      // channel: [meta, ref_fasta]

    main:
    ch_versions = channel.empty()

    // Step 0: segemehl index + align
    SEGEMEHL_INDEX ( fasta_meta )
    SEGEMEHL_ALIGN ( reads, SEGEMEHL_INDEX.out.idx, fasta_meta )
    ch_versions = ch_versions.mix(SEGEMEHL_INDEX.out.versions, SEGEMEHL_ALIGN.out.versions)

    // Step 1: Split-read detection via haarz
    HAARZ ( SEGEMEHL_ALIGN.out.sngl_bed )
    ch_versions = ch_versions.mix(HAARZ.out.versions)

    // Step 2: SAM → sorted BAM → BED (all alignments)
    SAMTOOLS_VIEW_SAM2BAM ( SEGEMEHL_ALIGN.out.sam, fasta_meta, [] )
    SAMTOOLS_SORT ( SAMTOOLS_VIEW_SAM2BAM.out.bam, fasta_meta, 'coordinate' )
    SAMTOOLS_INDEX ( SAMTOOLS_SORT.out.bam )
    BEDTOOLS_BAMTOBED ( SAMTOOLS_SORT.out.bam.join(SAMTOOLS_INDEX.out.bai) )
    ch_versions = ch_versions.mix(SAMTOOLS_VIEW_SAM2BAM.out.versions, SAMTOOLS_SORT.out.versions_samtools, BEDTOOLS_BAMTOBED.out.versions)

    // Step 3: Discordant read extraction (SAM flags 2, 83, 163)
    // ... (multiple SAMTOOLS_VIEW with flags -G 2, -f 83, -f 163 → BEDTOOLS_BAMTOBED → concat)
    // Simplified: emit DR.bed

    // Step 4: Coverage + peak detection → candidate extraction
    // BEDTOOLS_GENOMECOV → BEDTOOLS_MERGE → DR_regions.bed
    // BEDTOOLS_MAKEWINDOWS → window.bed
    // BEDTOOLS_COVERAGE → coverage.tsv
    // ECCSPLORER_PEAK_DETECT → peak_regions.bed
    // ECCSPLORER_CANDIDATE_EXTRACT → candidates.bed
    // BEDTOOLS_GETFASTA → sequences.fasta

    // Step 5: Analysis + visualization
    // ECCSPLORER_COVERAGE_PROFILE → per-candidate coverage
    // ECCSPLORER_NORMALIZE → RPM + fold enrichment
    // ECCSPLORER_VISUALIZE → plots
    // ECCSPLORER_HTML_REPORT → HTML

    emit:
    candidates_bed
    versions = ch_versions
}
```

#### 2.2 创建 `ecc_finder_slim_pipeline`

文件：`circdna.nf/subworkflows/local/ecc_finder_slim_pipeline/main.nf`

```nextflow
//
// ECC_FINDER_slim — 原子化 eccDNA 检测子工作流 (MAP_SR + ASM_SR)
// 复用 BAM_PREPROCESSING 产出的 sorted BAM，跳过 BWA 重复比对
//

include { SAMTOOLS_SORT as SAMTOOLS_SORT_NAME } from '../../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX } from '../../../modules/nf-core/samtools/index/main'
include { GENRICH          } from '../../../../bio.nf/modules/genrich/main'
include { TIDEHUNTER       } from '../../../../bio.nf/modules/tidehunter/main'
include { ECC_FINDER_MERGE_SCORE } from '../../../../bio.nf/modules/ecc_finder_slim/merge_score/main'
include { UNICYCLER        } from '../../../modules/nf-core/unicycler/main'
include { ECC_FINDER_ASM_FILTER } from '../../../../bio.nf/modules/ecc_finder_slim/asm_filter/main'

workflow ECC_FINDER_SLIM_PIPELINE {
    take:
    bam_sorted      // channel: [meta, bam] (from BAM_PREPROCESSING)
    bam_sorted_bai  // channel: [meta, bai]
    fasta_meta      // channel: [meta, ref_fasta]
    reads           // channel: [meta, [r1, r2]] (ASM_SR needs FASTQ)
    run_map_sr      // boolean
    run_asm_sr      // boolean

    main:
    ch_versions = channel.empty()

    if (run_map_sr) {
        // Name-sort BAM for Genrich (requires name-sorted input)
        SAMTOOLS_SORT_NAME ( bam_sorted, fasta_meta, 'name' )
        SAMTOOLS_INDEX ( SAMTOOLS_SORT_NAME.out.bam )

        ch_name_sorted = SAMTOOLS_SORT_NAME.out.bam.join(SAMTOOLS_INDEX.out.bai)
        ch_versions = ch_versions.mix(SAMTOOLS_SORT_NAME.out.versions_samtools)

        // Genrich peak calling
        GENRICH ( ch_name_sorted.map { meta, bam, bai -> [meta, bam] } )
        ch_versions = ch_versions.mix(GENRICH.out.versions)

        // TideHunter split-read detection
        TIDEHUNTER (
            ch_name_sorted.map { meta, bam, bai -> [meta, bam] },
            fasta_meta
        )
        ch_versions = ch_versions.mix(TIDEHUNTER.out.versions)

        // Merge + score
        ECC_FINDER_MERGE_SCORE (
            GENRICH.out.bed,
            TIDEHUNTER.out.bed,
            fasta_meta
        )
        ch_versions = ch_versions.mix(ECC_FINDER_MERGE_SCORE.out.versions)
    }

    if (run_asm_sr) {
        ch_shortreads = reads.map { meta, r ->
            def rlist = r instanceof List ? r : [r]
            [meta, rlist, []]
        }
        UNICYCLER ( ch_shortreads )
        ch_versions = ch_versions.mix(UNICYCLER.out.versions)

        ECC_FINDER_ASM_FILTER ( UNICYCLER.out.scaffolds )
        ch_versions = ch_versions.mix(ECC_FINDER_ASM_FILTER.out.versions)
    }

    emit:
    map_csv     = run_map_sr ? ECC_FINDER_MERGE_SCORE.out.csv : channel.empty()
    map_fasta   = run_map_sr ? ECC_FINDER_MERGE_SCORE.out.fasta : channel.empty()
    asm_fasta   = run_asm_sr ? ECC_FINDER_ASM_FILTER.out.fasta : channel.empty()
    versions    = ch_versions
}
```

#### 2.3 更新 `workflows/circdna.nf`：接入 slim 子工作流

在现有 `ECCSPLORER_PIPELINE` 调用之后，添加 slim 版本并行调用：

```nextflow
// Slim 版本控制（在现有 ECCSPLORER_PIPELINE 调用之后添加）
def use_eccsplorer_slim  = params.circle_identifier.contains('eccsplorer_map_slim')
def use_ecc_finder_map_sr_slim = params.circle_identifier.contains('ecc_finder_map_sr_slim')
def use_ecc_finder_asm_sr_slim = params.circle_identifier.contains('ecc_finder_asm_sr_slim')

if (use_eccsplorer_slim) {
    ECCSPLORER_SLIM_PIPELINE ( reads, fasta_meta )
}

if (use_ecc_finder_map_sr_slim || use_ecc_finder_asm_sr_slim) {
    ECC_FINDER_SLIM_PIPELINE (
        BAM_PREPROCESSING.out.bam_sorted,
        BAM_PREPROCESSING.out.bam_sorted_bai,
        fasta_meta,
        reads,
        use_ecc_finder_map_sr_slim,
        use_ecc_finder_asm_sr_slim
    )
}
```

#### 2.4 更新 `conf/modules.config`

```groovy
// ECCSPLORER_slim 配置
withName: 'ECCSPLORER_PEAK_DETECT' {
    publishDir = [
        path: { "${params.outdir}/eccsplorer_slim/peak_detect" },
        mode: params.publish_dir_mode,
        saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    ]
}
withName: 'ECCSPLORER_CANDIDATE_EXTRACT' {
    publishDir = [
        path: { "${params.outdir}/eccsplorer_slim/candidates" },
        mode: params.publish_dir_mode,
        saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    ]
}
// ...

// ECC_FINDER_slim 配置
withName: 'GENRICH' {
    ext.args = '-v'
    publishDir = [
        path: { "${params.outdir}/ecc_finder_slim/genrich" },
        mode: params.publish_dir_mode,
        saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    ]
}
withName: 'TIDEHUNTER' {
    publishDir = [
        path: { "${params.outdir}/ecc_finder_slim/tidehunter" },
        mode: params.publish_dir_mode,
        saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    ]
}
// ...
```

---

### 阶段 3：test_local 配置更新

[test_local.config](file:///Users/siyangming/nextflow_nf_core/circdna.nf/conf/test_local.config#L26-27) 修改：

```groovy
// 原版（黑盒模式）— 保留，注释掉
// circle_identifier   = 'eccsplorer_map,ecc_finder_map_sr,ecc_finder_asm_sr'

// slim 版（原子化模式）— 新增，用于测试新的轻量版流程
circle_identifier   = 'eccsplorer_map_slim,ecc_finder_map_sr_slim,ecc_finder_asm_sr_slim'
```

---

### 阶段 4：测试与验证

```bash
# 1. bio.nf 模块 stub 测试（验证语法）
cd /Users/siyangming/nextflow_nf_core/bio.nf/modules/genrich && nf-test test
cd /Users/siyangming/nextflow_nf_core/bio.nf/modules/tidehunter && nf-test test
cd /Users/siyangming/nextflow_nf_core/bio.nf/modules/haarz && nf-test test
cd /Users/siyangming/nextflow_nf_core/bio.nf/modules/eccsplorer_slim/peak_detect && nf-test test
# ... (each module)

# 2. circdna.nf stub 测试（语法验证）
cd /Users/siyangming/nextflow_nf_core/circdna.nf
nextflow run main.nf -profile test_local -stub

# 3. circdna.nf 真实数据测试（功能验证，slim 版）
nextflow run main.nf -profile test_local \
    --circle_identifier 'eccsplorer_map_slim,ecc_finder_map_sr_slim,ecc_finder_asm_sr_slim'

# 4. 与原版对比验证（切换回原版）
nextflow run main.nf -profile test_local \
    --circle_identifier 'eccsplorer_map,ecc_finder_map_sr,ecc_finder_asm_sr'
```

---

## 6. 关键设计决策

| 决策                       | 结论                                                                                                           |
| -------------------------- | -------------------------------------------------------------------------------------------------------------- |
| 旧模块保留还是删除？       | **保留**，新旧并行共存                                                                                         |
| 简化版命名                 | `ECCSPLORER_slim` / `ECC_FINDER_MAP_SR_slim` / `ECC_FINDER_ASM_SR_slim`                                        |
| 参数命名                   | `circle_identifier = 'eccsplorer_map_slim,...'`                                                                |
| 构建顺序                   | **阶段0** 构建 slim 镜像 → **阶段1** bio.nf 模块 → **阶段2** circdna.nf 接入                                   |
| segemehl 比对 vs BWA       | **保留 segemehl**（split-read 感知）                                                                           |
| ecc_finder MAP_SR BWA 比对 | **跳过**（复用 BAM_PREPROCESSING BWA BAM）                                                                     |
| ECCsplorer_slim 镜像       | **需新建** `quay.io/bioinfortools/eccsplorer_slim:1.0.0` (~400 MB, ↓73%), conda: `yangmingsi::eccsplorer_slim` |
| ecc_finder_slim 镜像       | **需新建** `quay.io/bioinfortools/ecc_finder_slim:1.0.0` (~250 MB, ↓87%), conda: `yangmingsi::ecc_finder_slim` |
| genrich 镜像               | `quay.io/biocontainers/genrich:0.6.1--h577a1d6_5`（~6 MB，现有）                                               |
| tidehunter 镜像            | `quay.io/biocontainers/tidehunter:1.5.6--h7f5d12c_0`（~15 MB，现有）                                           |
| haarz 镜像                 | 复用 `quay.io/biocontainers/segemehl:0.3.4--hc2ea5fd_5`（~50 MB，现有，haarz 内置）                            |

---

## 7. 文件变更清单

### 新增文件（bio.nf）

```
bio.nf/modules/genrich/main.nf                    # Genrich 模块 → biocontainers/genrich
bio.nf/modules/genrich/environment.yml
bio.nf/modules/genrich/meta.yml
bio.nf/modules/genrich/tests/main.nf.test
bio.nf/modules/genrich/tests/nextflow.config

bio.nf/modules/tidehunter/main.nf                 # TideHunter 模块 → biocontainers/tidehunter
bio.nf/modules/tidehunter/environment.yml
bio.nf/modules/tidehunter/meta.yml
bio.nf/modules/tidehunter/tests/main.nf.test
bio.nf/modules/tidehunter/tests/nextflow.config

bio.nf/modules/haarz/main.nf                      # haarz 模块 → biocontainers/segemehl (复用)
bio.nf/modules/haarz/environment.yml
bio.nf/modules/haarz/meta.yml

bio.nf/modules/eccsplorer_slim/conda-recipe/     # Conda recipe → yangmingsi::eccsplorer_slim
bio.nf/modules/eccsplorer_slim/Dockerfile         # Docker → quay.io/bioinfortools/eccsplorer_slim
bio.nf/modules/eccsplorer_slim/environment.yml
bio.nf/modules/eccsplorer_slim/bin/               # 自定义 Python/R 脚本
bio.nf/modules/eccsplorer_slim/peak_detect/main.nf
bio.nf/modules/eccsplorer_slim/candidate_extract/main.nf
bio.nf/modules/eccsplorer_slim/coverage_profile/main.nf
bio.nf/modules/eccsplorer_slim/normalize/main.nf
bio.nf/modules/eccsplorer_slim/visualize/main.nf
bio.nf/modules/eccsplorer_slim/html_report/main.nf

bio.nf/modules/ecc_finder_slim/conda-recipe/     # Conda recipe → yangmingsi::ecc_finder_slim
bio.nf/modules/ecc_finder_slim/Dockerfile         # Docker → quay.io/bioinfortools/ecc_finder_slim
bio.nf/modules/ecc_finder_slim/environment.yml
bio.nf/modules/ecc_finder_slim/bin/               # 自定义 Python 脚本
bio.nf/modules/ecc_finder_slim/merge_score/main.nf
bio.nf/modules/ecc_finder_slim/asm_filter/main.nf
```

### 新增文件（circdna.nf）

```
circdna.nf/subworkflows/local/eccsplorer_slim_pipeline/main.nf
circdna.nf/subworkflows/local/eccsplorer_slim_pipeline/meta.yml
circdna.nf/subworkflows/local/ecc_finder_slim_pipeline/main.nf
circdna.nf/subworkflows/local/ecc_finder_slim_pipeline/meta.yml
```

### 修改文件

```
circdna.nf/workflows/circdna.nf       # 添加 slim 子工作流 include + 调用
circdna.nf/conf/modules.config        # 添加 slim process 资源/publishDir 配置
circdna.nf/conf/test_local.config     # 切换到 slim 版 circle_identifier
circdna.nf/CHANGELOG.md               # 添加优化记录
```

### 不修改（保留）

```
circdna.nf/modules/local/eccsplorer/              # 原版保留
circdna.nf/modules/local/ecc_finder/              # 原版保留
circdna.nf/subworkflows/local/eccsplorer_pipeline/ # 原版保留
circdna.nf/subworkflows/local/ecc_finder_pipeline/ # 原版保留
```

---

_本计划基于 2026-08-07 对 ECCsplorer 源码、ecc_finder 源码的深度阅读，以及 bioconda/quay.io 容器可用性的实际调研。所有工具均有现有镜像，无需从头构建。_
