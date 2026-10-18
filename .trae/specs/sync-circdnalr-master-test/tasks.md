# Tasks

- [x] Task 1: 在 circdnalr 分支合并 master 并解决冲突

  - [x] 1.1 切换到 circdnalr 分支，拉取最新 origin/master
  - [x] 1.2 执行 `git merge master`，列出冲突清单（38 处冲突）
  - [x] 1.3 按类别逐一解决冲突：
    - 核心工作流 `workflows/circdna.nf`（同时保留短读长 + 长读长分支）
    - `nextflow.config` / `nextflow_schema.json`（合并参数：长读长 protocol/long_read_identifier + master 短读长参数，版本 v4.1.0）
    - `conf/modules.config`（合并长读长与短读长模块配置，自动合并无冲突）
    - `bin/check_samplesheet.py` / `subworkflows/local/input_check/main.nf`（合并两套校验/解析逻辑）
    - `assets/schema_input.json`、`modules/local/samplesheet_check/main.nf`
    - samplesheets：circdna*\* 系列以 master 为准（含 data_type 列、PlanteccDNADB 路径），circdnalr*\* 长读系列用 master 路径；保留 master 的 ONT/PacBio 测试数据文件（ont/pacbio 目录）；circdna_tgs_clean.csv 保留 HEAD 并修正路径
    - `.gitignore`、`CHANGELOG.md`、`SERVER_RUN_GUIDE.md` 合并两侧内容
  - [x] 1.4 验证合并结果：`git status` 无冲突标记；`nextflow config -profile test_local` 可解析；`python -m py_compile bin/check_samplesheet.py` 通过
  - [x] 1.5 提交合并（commit 33445d1）

- [x] Task 2: 准备长读长测试配置与样本表

  - [x] 2.1 更新 `conf/test_nanopore_lr.config`：`long_read_identifier = 'cresil,fled,flye,eccfinder'`，输入指向服务器内置 `ont_eccdna_smoke.fastq.gz`，fasta 指向拟南芥参考基因组（PublicDB）
  - [x] 2.2 更新 `conf/test_pacbio_lr.config`：`long_read_identifier = 'cresil,fled,flye,eccfinder'`，输入指向服务器内置 `pacbio_eccdna_smoke.fastq.gz`，fasta 指向水稻参考基因组（PublicDB）
  - [x] 2.3 创建 `samplesheets/test_ont_lr.csv` 与更新 `samplesheets/test_pacbio_lr.csv`（列：sample,fastq_1,fastq_2），路径使用服务器绝对路径
  - [x] 2.4 本地运行 `nextflow config -profile test_nanopore_lr` 与 `test_pacbio_lr` 验证参数解析；修复 check_samplesheet.py 长读输出分支顺序 bug（commit 07664d8）
  - [x] 2.5 提交本任务改动（commit 07664d8）

- [x] Task 3: 推送合并结果并同步到服务器

  - [x] 3.1 推送 `circdnalr` 到 `origin/circdnalr`（afb473c..e685e94）
  - [x] 3.2 SSH 到服务器，在 `/data1/users/siyangming/PlanteccDNADB/circdna.nf` 执行 `git fetch origin` + `git checkout circdnalr`（服务器工作区原本干净）
  - [x] 3.3 核对服务器上测试数据存在：`testdatasets/ont/ont_eccdna_smoke.fastq.gz`（6.5M）、`testdatasets/pacbio/pacbio_eccdna_smoke.fastq.gz`（25M）、参考基因组（拟南芥 32M、水稻 99M）
  - [x] 3.4 验证服务器 nextflow 环境可用（conda env `nextflow` → nextflow 26.04.6；server.config docker.enabled=true）

- [ ] Task 4: 服务器运行 ONT 长读长测试（四引擎）

  - [ ] 4.1 运行 `nextflow run main.nf -profile server,test_nanopore_lr --outdir ...`（或等价命令，含 -with-report/-with-trace/-resume）
  - [ ] 4.2 检查 CReSIL、FLED、FLYE、ECCFINDER 四引擎均产生输出
  - [ ] 4.3 记录退出状态、trace/report 概要、主要输出文件

- [ ] Task 5: 服务器运行 PacBio 长读长测试（四引擎）

  - [ ] 5.1 运行 `nextflow run main.nf -profile server,test_pacbio_lr --outdir ...`
  - [ ] 5.2 检查 CReSIL、FLED、FLYE、ECCFINDER 四引擎均产生输出
  - [ ] 5.3 记录退出状态、trace/report 概要、主要输出文件

- [ ] Task 6: 测试结果验证与记录
  - [ ] 6.1 汇总 ONT/PacBio 测试结果（命令、状态、输出路径、产出文件）
  - [ ] 6.2 记录至 `.trae/changes/` 或本 spec 目录（不创建冗余文档，遵循 AGENTS.md 12 节归属规则）
  - [ ] 6.3 若测试失败，定位原因并修复（可复用 `-resume` 快速迭代），直到通过

# Task Dependencies

- [Task 2] depends on [Task 1]（需先合并得到长读长模块 + master 测试数据）
- [Task 3] depends on [Task 1] and [Task 2]（先提交合并与配置再推送）
- [Task 4] and [Task 5] depend on [Task 3]，两者可并行
- [Task 6] depends on [Task 4] and [Task 5]
