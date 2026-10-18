# Tasks

- [x] Task 1: 创建 ecc_finder 和 ECCsplorer 的 `slim` 分支

  - [x] 1.1 ecc_finder 本地仓库创建 `slim` 分支（从 main 分出）
  - [x] 1.2 ECCsplorer 本地仓库创建 `slim` 分支（从 master 分出）
  - [x] 1.3 确认 master/main 分支 Dockerfile 和 conda-recipe 未被修改

- [x] Task 2: 从 ecc_finder 源码提取 slim 专有脚本（代码+构建完成）

  - [x] 2.1 merge_score.py 独立脚本
  - [x] 2.2 asm_filter.py 独立脚本
  - [x] 2.3 最小化 environment.yml
  - [x] 2.4 conda-recipe/meta.yaml + build.sh
  - [x] 2.5 Dockerfile（slim版）
  - [x] 2.6 conda 包已构建（`/opt/homebrew/Caskroom/miniforge/base/envs/conda_build/conda-bld/noarch/ecc_finder_slim-1.0.0-0.conda`）
  - [ ] 2.7 Docker 构建/推送（网络不可达，待后续执行：`cd ecc_finder && docker build -t quay.io/bioinfortools/ecc_finder_slim:1.0.0 . && docker push`）
  - [ ] 2.8 anaconda 推送（需交互式登录：`anaconda login && anaconda upload ...`）
  - [x] 2.9 bio.nf 模块 merge_score/main.nf 和 asm_filter/main.nf 已创建

- [x] Task 3: 从 ECCsplorer 源码提取 6 个 slim 专有脚本（代码+构建完成）

  - [x] 3.1-3.6 6 个独立脚本全部创建
  - [x] 3.7-3.9 构建文件全部创建
  - [x] 3.10 conda 包已构建（`/opt/homebrew/Caskroom/miniforge/base/envs/conda_build/conda-bld/noarch/eccsplorer_slim-1.0.0-0.conda`）
  - [ ] 3.10b Docker 构建/推送（网络不可达）
  - [x] 3.11 bio.nf 6 个模块 main.nf 已创建

- [x] Task 4: 更新 circdna.nf slim 子工作流接入专有逻辑模块

  - [x] 4.1 ecc_finder_slim_pipeline 接入 MERGE_SCORE 和 ASM_FILTER
  - [x] 4.2 eccsplorer_slim_pipeline 接入 6 个 ECCsplorer 专有逻辑模块
  - [x] 4.3 modules.config 添加新模块 publishDir 配置

- [x] Task 5: 验证
  - [x] 5.1 Stub 编译+channel 连接通过（16 个 slim 进程全部成功编译提交）
  - [ ] 5.2 Docker profile 完整运行（需 Docker 环境+镜像构建后）
