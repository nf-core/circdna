# Tasks

## 阶段一：调研与构建 Conda 包和 Docker 镜像

- [x] Task 1: 调研 ECCsplorer 可用性

  - [x] SubTask 1.1: 检查 `quay.io/biocontainers/eccsplorer` 是否存在 → 不存在
  - [x] SubTask 1.2: 检查 bioconda 是否有 `eccsplorer` 包 → 不存在（bioconda-recipes 仓库搜索无结果）
  - [x] SubTask 1.3: 记录调研结果 → 需要构建自定义 conda 包和 Docker 镜像

- [x] Task 2: 构建 ECCsplorer conda 包（若 Task 1 确认 bioconda 不存在）

  - [x] SubTask 2.1: 创建 `conda-recipe/meta.yaml`，声明包名 `eccsplorer`、版本 `2022.01.1.1`、Python 3.7、numpy、biopython、scipy、pyRserve、R (ggplot2/ggrepel/gridExtra/dplyr) 等依赖
  - [x] SubTask 2.2: 系统级依赖（blast+、segemehl、samtools≥1.9、bedtools≥2.28.0、RepeatExplorer2、Trimmomatic）作为 run requirements 声明（已有 bioconda 包则引用，否则单独打包）
  - [x] SubTask 2.3: 编写 `conda-recipe/build.sh`，从 GitHub 克隆 ECCsplorer 源码并配置 `lib/config.py`
  - [ ] SubTask 2.4: 本地构建验证：`conda build conda-recipe/` + `conda create -n eccsplorer-test --use-local eccsplorer`（需用户手动执行）
  - [ ] SubTask 2.5: 使用 testdata 验证 conda 包可运行 ECCsplorer.py（需用户手动执行）
  - [ ] SubTask 2.6: 登录 anaconda.org，推送至用户频道：`anaconda upload` → `https://anaconda.org/siyangming/eccsplorer`（需用户手动执行）

- [x] Task 3: 构建 ECCsplorer Docker 镜像（若 Task 1 确认 biocontainers 不存在）
  - [x] SubTask 3.1: 编写 Dockerfile，优先基于 conda 包（`conda install -c siyangming eccsplorer=2022.01.1.1`），否则直接从源码安装所有依赖
  - [x] SubTask 3.2: 安装系统级依赖（blast+、segemehl、samtools、bedtools、RepeatExplorer2、R 包、Trimmomatic、seqtk）
  - [x] SubTask 3.3: 克隆 ECCsplorer 仓库至镜像内，配置 lib/config.py 中的第三方工具路径
  - [ ] SubTask 3.4: 本地构建并测试镜像：`docker build -t eccsplorer:2022.01.1.1 .`（需用户手动执行）
  - [ ] SubTask 3.5: 使用 testdata 验证镜像可运行 ECCsplorer.py（需用户手动执行）
  - [ ] SubTask 3.6: 推送镜像至 `quay.io/siyangming/eccsplorer:2022.01.1.1`（需用户手动执行）

## 阶段二：在 bio.nf 构建 eccsplorer 模块

- [x] Task 4: 创建 bio.nf eccsplorer 分支与模块结构

  - [ ] SubTask 4.1: 在 `bio.nf/` 创建分支 `eccsplorer`（需用户手动执行 git 操作）
  - [x] SubTask 4.2: 创建 `bio.nf/modules/eccsplorer/` 目录结构（含 tests/、testdata/）
  - [x] SubTask 4.3: 创建 `environment.yml`，channels 列表含 `siyangming`（自定义频道），依赖 `- siyangming::eccsplorer=2022.01.1.1`

- [x] Task 5: 编写 main.nf 模块逻辑

  - [x] SubTask 5.1: 定义 process ECCSPLORER，label `process_high`
  - [x] SubTask 5.2: 定义输入通道：`tuple val(meta), path(reads_r1), path(reads_r2), path(fasta)`
  - [x] SubTask 5.3: 定义输出通道：candidates_bed、junction_reads、versions
  - [x] SubTask 5.4: 编写 script 块调用 `python ECCsplorer.py`，支持 mapping 模式（comparative 通过 ext.args 扩展）
  - [x] SubTask 5.5: 编写 stub 块生成空输出文件
  - [x] SubTask 5.6: 配置 container 指向 `quay.io/siyangming/eccsplorer:2022.01.1.1`
  - [x] SubTask 5.7: 支持 `task.ext.args` 传递额外参数

- [x] Task 6: 编写 meta.yml 和测试
  - [x] SubTask 6.1: 编写 meta.yml（name、description、keywords、tools、input、output、authors）
  - [x] SubTask 6.2: 编写 tests/main.nf.test（stub 模式测试用例）
  - [x] SubTask 6.3: 编写 tests/nextflow.config（docker.enabled、`--platform linux/amd64`）
  - [ ] SubTask 6.4: 运行 `nf-test test` 验证 stub 模式通过（需 Docker 和镜像构建后执行）

## 阶段三：拷贝模块至 circdna.nf 并接入

- [x] Task 7: 拷贝模块并替换占位实现

  - [x] SubTask 7.1: 将 `bio.nf/modules/eccsplorer/` 拷贝至 `circdna.nf/modules/local/eccsplorer/`（覆盖占位 main.nf）
  - [x] SubTask 7.2: 调整 include 路径适配 circdna.nf 目录结构
  - [x] SubTask 7.3: 确保 environment.yml 同步拷贝并指向 `siyangming::eccsplorer`

- [x] Task 8: 调整 eccdna_mode 子工作流

  - [x] SubTask 8.1: 修改 `subworkflows/local/eccdna_mode/main.nf` 的 ECCSPLORER 调用，输入从 `ch_bam_bai` 改为 `reads`（FASTQ）
  - [x] SubTask 8.2: 传递 `fasta` 通道至 ECCSPLORER（已有 fasta_meta 通道可复用）
  - [x] SubTask 8.3: 添加 control reads 通道支持（可选，通过参数控制）
  - [x] SubTask 8.4: 验证 emit 通道 candidates_bed 和 junction_reads 正确连接

- [x] Task 9: 更新配置文件
  - [x] SubTask 9.1: 在 `conf/modules.config` 添加 ECCSPLORER 资源配置（label process_high，内存 16GB+）
  - [x] SubTask 9.2: 在 `nextflow.config` 添加 ECCsplorer 参数（`eccsplorer_trim_reads`）
  - [x] SubTask 9.3: 在 `conf/test_local.config` 添加 ECCsplorer 测试参数（使用默认值，无需额外配置）

## 阶段四：版本更新与文档

- [x] Task 10: 更新版本号与 CHANGELOG

  - [x] SubTask 10.1: 在 `CHANGELOG.md` 顶部新增 `## v4.1.0 - [2026-08-02]` 条目
  - [x] SubTask 10.2: 更新 `nextflow.config` manifest `version = '4.1.0'`
  - [ ] SubTask 10.3: 运行 `nf-core schema build` 重新生成 nextflow_schema.json（需用户手动执行）
  - [x] SubTask 10.4: 更新 CHANGELOG 描述 ECCsplorer 真实检测实现、输入接口变更、conda 包与 Docker 镜像发布

- [x] Task 11: 更新 AGENTS.md 工作流程规则
  - [x] SubTask 11.1: 在 AGENTS.md 新增 "12. 第三方模块构建工作流程" 章节
  - [x] SubTask 11.2: 编写规则：nf-core 未提供软件 → bio.nf 构建模块 → 拷贝至目标流程
  - [x] SubTask 11.3: 编写规则：bioconda 不可用 → 构建自定义 conda 包 → 推送至用户 anaconda.org 频道
  - [x] SubTask 11.4: 编写规则：biocontainers 不可用 → 构建自定义 Docker 镜像 → 推送至用户 quay.io 频道
  - [x] SubTask 11.5: 参考本次 ECCsplorer 实现作为示例

## 阶段五：验证与提交

- [x] Task 12: 本地验证

  - [x] SubTask 12.1: `nextflow run main.nf -profile test_local,docker --mode eccdna -stub-run` 验证 eccdna 模式（使用临时 container override，因 quay.io/siyangming/eccsplorer:2022.01.1.1 镜像尚未构建）
  - [x] SubTask 12.2: 确认 ECCSPLORER stub 模式正常工作（3/3 样本通过，输出 circdna*{1,2,3}\_candidates.bed 和 circdna*{1,2,3}\_junction_reads.txt）
  - [x] SubTask 12.3: 确认下游通道（eccsplorer_bed）正确连接（Pipeline completed successfully, 52 tasks succeeded）

- [ ] Task 13: 提交至 GitHub
  - [ ] SubTask 13.1: 提交 conda-recipe 至 ECCsplorer fork 或独立仓库
  - [ ] SubTask 13.2: 提交 bio.nf eccsplorer 分支至 GitHub
  - [ ] SubTask 13.3: 提交 circdna.nf v4.1.0 至 GitHub（含 tag）
  - [ ] SubTask 13.4: 提交 AGENTS.md 更新

# Task Dependencies

- Task 2 与 Task 3 依赖 Task 1（确认 conda 包/镜像不存在才构建）
- Task 3（Docker 镜像）可基于 Task 2（conda 包），亦可并行构建
- Task 4-6 可并行（模块结构、main.nf、meta.yml 独立编写）
- Task 7 依赖 Task 5 + Task 6（需要完整模块才能拷贝）
- Task 8 依赖 Task 7（模块就位后才能调整工作流）
- Task 9 依赖 Task 8（工作流调整后才能更新配置）
- Task 10 依赖 Task 9（所有代码变更完成后才更新版本）
- Task 11 可与 Task 4-10 并行（AGENTS.md 独立更新）
- Task 12 依赖 Task 10
- Task 13 依赖 Task 12

# 并行化说明

- 阶段一（Task 1-3）与阶段四的 Task 11（AGENTS.md）可并行
- 阶段二的 Task 4、5、6 可并行编写不同文件
- 阶段三必须串行执行（模块拷贝 → 工作流调整 → 配置更新）
