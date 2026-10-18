# ecc_finder Nextflow 模块构建 - 实现计划

## [x] Task 1: 创建 conda 环境配置文件
- **Priority**: high
- **Depends On**: None
- **Description**: 
  - 基于原始 ecc_finder.yaml 创建精简的 environment.yml，只保留必要依赖
  - 包含：python、pysam、numpy、pandas、matplotlib、pybedtools、minimap2、bwa、samtools、bedtools、fastp、seqtk、cd-hit、unicycler、tidehunter、genrich、requests
- **Acceptance Criteria Addressed**: AC-1
- **Notes**: 已完成，ecc_finder conda 环境已存在并验证通过

## [x] Task 2: 创建 Dockerfile
- **Priority**: high
- **Depends On**: Task 1
- **Description**: 
  - 创建 Dockerfile，基于 continuumio/miniconda3
  - 安装 conda 环境依赖，通过 pip install . 安装 ecc_finder
  - tidehunter 和 genrich 通过 bioconda 安装（非 PyPI）
- **Notes**: 已完成并修复 pip 环境问题

## [x] Task 3-6: 创建四个 Nextflow 模块
- **Priority**: high
- **Depends On**: Task 1
- **Description**: 
  - 创建 map-ont、map-sr、asm-ont、asm-sr 模块
  - 每个模块包含 main.nf、meta.yml、environment.yml
  - 测试文件 main.nf.test 和 nextflow.config 已创建
- **Notes**: 全部完成

## [x] Task 7: 修复模块 ecc_finder.py 调用路径问题并创建缺失测试文件
- **Priority**: high
- **Depends On**: Tasks 3-6
- **Description**: 
  - 修改 Dockerfile 使用 `pip install .` 安装 ecc_finder 到 PATH
  - 模块脚本改为直接调用 `ecc_finder.py`
  - 创建所有缺失的测试文件
- **Notes**: 已完成

## [x] Task 8: 创建 conda-recipe 并构建 conda 包
- **Priority**: high
- **Depends On**: Task 7
- **Description**: 
  - 创建 `ecc_finder/conda-recipe/meta.yaml` 和 `build.sh`
  - 使用 conda_build 环境构建 osx-arm64 平台的 conda 包
  - 固定 Python 版本 <3.13 避免兼容性问题
  - 添加 C/C++ 编译器和 wheel 到构建依赖
- **Acceptance Criteria Addressed**: AC-1
- **Notes**: 包文件: ecc_finder-1.0.0-py312h70deae4_0.conda (210KB)

## [x] Task 9: 上传 conda 包到 anaconda.org
- **Priority**: high
- **Depends On**: Task 8
- **Description**: 
  - 使用 anaconda 命令上传构建好的 conda 包
  - 上传到 anaconda.org/YangmingSi/ecc_finder
- **Acceptance Criteria Addressed**: AC-1
- **Notes**: 包已公开可用，安装命令: `conda install --channel YangmingSi ecc_finder`

## [x] Task 10: 构建 Docker 镜像并推送到 quay.io
- **Priority**: high
- **Depends On**: Task 7
- **Description**: 
  - 使用 `--platform linux/amd64` 构建 Docker 镜像
  - 镜像标签：`quay.io/bioinfortools/ecc_finder:1.0.0`
  - 推送镜像到 quay.io
  - 验证容器中 ecc_finder.py 可运行
- **Acceptance Criteria Addressed**: AC-2
- **Notes**: 镜像 4.76GB，验证 ecc_finder.py --version 输出 v1.0.0，TideHunter 和 Genrich 均可用

## [x] Task 11: 运行模块测试并验证
- **Priority**: high
- **Depends On**: Task 10
- **Description**: 
  - 运行所有四个模块的 stub 测试
  - 修复 map_sr 文件名冲突（创建 BWA 索引）
  - 更新所有测试快照
- **Acceptance Criteria Addressed**: AC-3, AC-4, AC-5, AC-6, AC-7
- **Notes**: 四个模块 stub 测试全部通过

## [x] Task 12: 更新模块版本引用
- **Priority**: medium
- **Depends On**: Task 9, Task 10
- **Description**: 
  - 所有模块 container 引用指向 quay.io/bioinfortools/ecc_finder:1.0.0
  - 所有模块 environment.yml 包含完整依赖（含 tidehunter、genrich）
  - meta.yml 文件信息正确
- **Acceptance Criteria Addressed**: AC-7
- **Notes**: 已验证所有引用正确

# Task Dependencies
- Task 8 depends on Task 7
- Task 9 depends on Task 8
- Task 10 depends on Task 7
- Task 11 depends on Task 10
- Task 12 depends on Task 9 and Task 10
