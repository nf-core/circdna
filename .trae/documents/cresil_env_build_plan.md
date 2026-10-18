# CReSIL Conda 环境和 Docker 环境构建计划

## 一、项目调研结论

### 1.1 项目概况
- **项目名称**: CReSIL (Accurate Identification of Extrachromosomal Circular DNA from Long-read Sequences)
- **版本**: V1.2.0
- **源码位置**: `/Users/siyangming/nextflow_nf_core/bio.nf/cresil-master`
- **项目类型**: Python 包，提供 CLI 工具（cresil trim / identify / identify_wgls / annotate / visualize）

### 1.2 现有配置
- 已存在 [environment.yml](file:///Users/siyangming/nextflow_nf_core/bio.nf/cresil-master/environment.yml)，包含完整的 conda 依赖
- 已存在 [setup.py](file:///Users/siyangming/nextflow_nf_core/bio.nf/cresil-master/setup.py)，用于 pip 安装 cresil 包
- 不存在 Dockerfile

### 1.3 依赖分析（来自 environment.yml）
| 类别 | 依赖项 |
|------|--------|
| Python 基础 | python=3.8.15, pip=23.3.2 |
| 生物信息学工具 | minimap2=2.26, samtools=1.18, bedtools=2.31.1, bcftools=1.17, tabix=1.11, bioawk=1.0, mosdepth=0.3.3 |
| Python 包 | biopython=1.81, mappy=2.26, pysam=0.21.0, pybedtools=0.9.1, pandas=2.0.3, matplotlib=3.7.3, networkx=3.1, tqdm=4.66.1, python-intervaltree=3.1.0, python-graphviz=0.20.1 |
| 其他工具 | openjdk=17.0.3, igvtools=2.14.1, nanofilt=2.8.0, medaka=1.11.3, parallel=20231122, jedi=0.19.1, ipykernel=6.26.0 |

### 1.4 bio.nf 项目规范参考
- 模块目录结构: `modules/<tool_name>/` 下包含 `environment.yml`、`meta.yml`、`main.nf`
- Conda channels: conda-forge, bioconda, defaults
- 容器策略: 优先使用 nf-core 提供的容器，否则使用 Galaxy Singularity/quay.io 容器

## 二、需要创建/修改的文件

### 2.1 Conda 环境相关
- 验证现有 `cresil-master/environment.yml` 是否可用
- 可选：创建独立的 `modules/cresil/environment.yml`（按 bio.nf 模块规范）

### 2.2 Docker 环境相关
- 新建 `cresil-master/Dockerfile` — 基于 conda 环境构建 Docker 镜像
- 新建 `cresil-master/docker-compose.yml`（可选，方便本地测试）

### 2.3 模块整合（可选）
- 新建 `modules/cresil/` 目录结构（main.nf, meta.yml, environment.yml）

## 三、实施步骤

### 步骤 1: 验证并优化 Conda 环境配置
1. 检查 `environment.yml` 中所有包的版本可用性
2. 确认 conda channels 顺序（conda-forge > bioconda，确保依赖解析正确）
3. 在本地使用 `conda env create -f environment.yml` 测试环境创建

### 步骤 2: 构建 Docker 镜像（Dockerfile）
1. **基础镜像选择**: 使用 `mambaorg/micromamba` 或 `continuumio/miniconda3` 作为基础镜像
   - 推荐 `mambaorg/micromamba`：体积更小，依赖解析更快
2. **Dockerfile 结构**:
   ```
   FROM mambaorg/micromamba:latest
   COPY environment.yml /tmp/environment.yml
   RUN micromamba create -f /tmp/environment.yml && micromamba clean -a -y
   COPY . /app
   WORKDIR /app
   RUN pip install .
   ENV PATH /opt/conda/envs/cresil/bin:$PATH
   CMD ["cresil", "--help"]
   ```
3. 优化镜像大小：多阶段构建，清理 conda 缓存和 pip 缓存

### 步骤 3: 构建并测试 Docker 镜像
1. 执行 `docker build -t cresil:1.2.0 .`
2. 运行 `docker run --rm cresil:1.2.0 cresil --help` 验证安装
3. 使用 example 数据进行基本功能测试

### 步骤 4: （可选）按 bio.nf 模块规范整合
1. 创建 `modules/cresil/` 目录
2. 复制/创建 `environment.yml`
3. 创建 `meta.yml` 描述工具信息
4. 创建 `main.nf` 基础流程模板

## 四、潜在依赖与注意事项

### 4.1 依赖风险
- **medaka=1.11.3**: 体积较大，依赖复杂（含 tensorflow），可能导致镜像体积庞大（>5GB）
- **openjdk=17.0.3**: igvtools 依赖 Java
- **Python 3.8.15**: 较旧版本，需确认所有依赖兼容性
- **ipykernel / jedi**: 可能非运行必需，可考虑移除以减小镜像体积

### 4.2 构建注意事项
- Conda 环境解析可能较慢，建议使用 mamba 替代 conda
- Docker 镜像构建时需注意网络问题（conda 下载慢）
- 若 medaka 不是必需依赖，可考虑分层镜像或提供轻量版

### 4.3 版本兼容性
- 所有依赖版本已在 environment.yml 中锁定
- 需注意 bioconda 和 conda-forge 渠道的包版本可用性

## 五、风险处理

| 风险 | 影响 | 应对策略 |
|------|------|----------|
| conda 依赖解析失败 | 环境构建失败 | 使用 mamba 加速解析；调整 channels 优先级；锁定更精确的版本 |
| Docker 镜像体积过大 | 部署/传输慢 | 多阶段构建；清理缓存；移除非必要依赖（如 ipykernel） |
| medaka 安装失败 | 部分功能不可用 | 确认 medaka 版本与 python 3.8 兼容性；考虑使用单独的 medaka 容器 |
| 构建网络问题 | 下载超时 | 使用国内 conda 镜像源；提前下载依赖包 |

## 六、验证标准

### Conda 环境验证
- [ ] `conda env create -f environment.yml` 成功执行
- [ ] `conda activate cresil && pip install .` 成功
- [ ] `cresil --help` 正常输出帮助信息
- [ ] `cresil trim --help` 等子命令可用

### Docker 环境验证
- [ ] `docker build` 成功完成
- [ ] `docker run cresil:1.2.0 cresil --help` 正常输出
- [ ] 镜像体积控制在合理范围（建议 < 10GB）
