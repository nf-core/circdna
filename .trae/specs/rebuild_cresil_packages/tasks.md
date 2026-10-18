# Rebuild and Publish CReSIL Docker/Conda Packages - The Implementation Plan

## [x] Task 1: 更新所有版本号到 1.2.1
- **Priority**: high
- **Depends On**: None
- **Description**: 
  - 将 CReSIL 所有版本号从 1.2.0 统一更新到 1.2.1
  - 更新文件：
    - `cresil/cresil/__init__.py`: `__version__`
    - `cresil/Dockerfile`: `LABEL version`
    - `cresil/conda-recipe/meta.yaml`: `version`
    - `cresil/README.md`: 版本引用
- **Acceptance Criteria Addressed**: AC-1
- **Test Requirements**:
  - `programmatic` TR-1.1: `cresil/__init__.py` 中 `__version__` 为 `'1.2.1'`
  - `programmatic` TR-1.2: `Dockerfile` 中 `LABEL version` 为 `"1.2.1"`
  - `programmatic` TR-1.3: `conda-recipe/meta.yaml` 中 version 为 `"1.2.1"`
  - `programmatic` TR-1.4: `README.md` 中版本号引用更新为 1.2.1
- **Notes**: 统一所有版本号，便于 GitHub 备份

## [x] Task 2: 构建并推送 Docker 镜像
- **Priority**: high
- **Depends On**: Task 1
- **Description**: 
  - 使用 `--platform linux/amd64` 构建 Docker 镜像
  - 镜像标签：`quay.io/bioinfortools/cresil:1.2.1`
  - 构建目录：`/Users/siyangming/nextflow_nf_core/bio.nf/cresil/`
  - 推送镜像到 quay.io
- **Acceptance Criteria Addressed**: AC-2, AC-3, AC-7
- **Test Requirements**:
  - `programmatic` TR-2.1: 镜像构建成功，docker images 可见
  - `programmatic` TR-2.2: 容器中 `cresil --version` 输出 `cresil 1.2.1`
  - `programmatic` TR-2.3: 镜像成功推送到 quay.io
  - `programmatic` TR-2.4: 使用测试数据运行 identify_wgls 不再出现 KeyError
- **Notes**: 基础镜像为 `quay.io/bioinfortools/cresil:1.2.0`，修复了 strand 比较 bug（1/-1 → '+/-'）和 bedtools 排序问题

## [x] Task 3: 构建并发布 Conda 包
- **Priority**: high
- **Depends On**: Task 1
- **Description**: 
  - 使用 conda build 构建 linux-64 平台的 Conda 包 1.2.1
  - 构建环境：Docker（quay.io/bioinfortools/cresil:1.2.1 容器内安装 conda-build）
  - 发布到 anaconda.org/yangmingsi/cresil
- **Acceptance Criteria Addressed**: AC-4, AC-5
- **Test Requirements**:
  - `programmatic` TR-3.1: Conda 包构建成功，生成 `cresil-1.2.1-py38_0.conda`
  - `programmatic` TR-3.2: 包文件中包含修复后的 identify_wgls.py
  - `programmatic` TR-3.3: 包成功发布到 Anaconda，版本为 1.2.1
- **Notes**: 构建耗时约 25 分钟（主要是依赖解析）

## [x] Task 4: 更新 Nextflow 模块版本引用
- **Priority**: medium
- **Depends On**: Task 2, Task 3
- **Description**: 
  - 更新所有 CReSIL Nextflow 模块的 container 版本到 1.2.1
  - 更新 environment.yml 的 conda 版本引用
  - 更新 meta.yml 中的版本信息
  - 涉及模块：trim, identify, identify_wgls, annotate, visualize
- **Acceptance Criteria Addressed**: AC-6
- **Test Requirements**:
  - `programmatic` TR-4.1: 所有模块 main.nf 中的 container 标签更新为 1.2.1
  - `programmatic` TR-4.2: Singularity 镜像路径同步更新
  - `programmatic` TR-4.3: 所有模块 stub 测试通过
- **Notes**: 5 个模块都需要更新

# Task Dependencies
- Task 2 depends on Task 1
- Task 3 depends on Task 1
- Task 4 depends on Task 2 and Task 3
- Task 2 和 Task 3 可以并行执行