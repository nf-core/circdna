# Tasks

- [x] Task 1: 将 master 提交合并到当前分支（ECCsplorer），冲突以 master 为准

  - [x] SubTask 1.1: `git stash -u`（暂存当前 80 项未提交修改，含未跟踪文件）
  - [x] SubTask 1.2: `git merge master`；冲突文件以 master 为准（`git checkout --theirs <conflicts>` + `git add`）
  - [x] SubTask 1.3: 提交 merge；`git stash pop` 恢复 slim 独有修改（pop 冲突以 master 为准）
  - [x] SubTask 1.4: 确认 merge 后 master 提交已并入、数据文件与 master 版 server.config 就位、slim 独有文件保留

- [x] Task 2: 测试数据就位确认（相对路径）

  - [x] SubTask 2.1: 确认 merge 后 `testdatasets/ont/ont_eccdna_smoke/regular.fastq.gz`、`testdatasets/pacbio/pacbio_eccdna_smoke/regular.fastq.gz` 存在且可读
  - [x] SubTask 2.2: 将 `samplesheets/test_ont.csv`/`test_pacbio.csv` 中 fastq 路径改写为**相对路径**（`testdatasets/ont/ont_eccdna_smoke.fastq.gz` 等），确认 `platform`/`protocol` 列（ont|pacbio / long_read）
  - [x] SubTask 2.3: 确认 `testdatasets/extract_test_data.sh` 就位

- [x] Task 3: 制作 control 配对 test\_ samplesheet

  - [x] SubTask 3.1: 基于 `samplesheets/test_real_integrated.csv`（6 样本 gdna/circdna × 3）补 `pair` 列（gdna_1/circdna_1=p1，gdna_2/circdna_2=p2，gdna_3/circdna_3=p3）
  - [x] SubTask 3.2: 生成 `samplesheets/test_eccsplorer_pair.csv`，列 `sample,fastq_1,fastq_2,datatype,pair`
  - [x] SubTask 3.3: 用 `bin/check_samplesheet.py` 验证新 samplesheet 解析通过

- [x] Task 4: 参考基因组入口确认（master server.config，不下载）

  - [x] SubTask 4.1: 确认 master 版 `conf/server.config` 就位（`fasta_base_path=/data1/users/siyangming/PublicDB/reference` + `genomes` 块含 Arabidopsis_thaliana/Oryza_sativa 等）
  - [x] SubTask 4.2: 记录服务器运行命令模板（`-profile server` + `--genome <species>` + `--input samplesheets/test_*.csv` + `--circle_identifier ...`）到文档

- [x] Task 5: ECCsplorer 分支冻结修改保留价值评估（merge 后）

  - [x] SubTask 5.1: 分析 merge 后 `git status` 剩余修改，分类：slim 链新增（ECCsplorer 独有）/ 冗余废弃 / 配置修改
  - [x] SubTask 5.2: 确认冲突文件已按 master 版处理（server.config/modules.config/nextflow.config/check_samplesheet.py 等）
  - [x] SubTask 5.3: 产出保留/合并/丢弃清单（写入 .trae/documents 报告），呈现用户确认

- [x] Task 6: 实际合并与清理（清单确认后执行）

  - [x] SubTask 6.1: 提交有价值修改（slim 模块链、vendored 模块、新子工作流、params、test\_ samplesheet）到 ECCsplorer 分支
  - [x] SubTask 6.2: 清理无价值/冗余文件（确认后删除）
  - [x] SubTask 6.3: 提交后主流水线 slim 标识符 stub 编译通过

- [x] Task 7: 交付说明
  - [x] SubTask 7.1: 汇总交付物（test\_ 文件清单 + 服务器运行命令模板 + 数据相对路径说明）写入文档
  - [x] SubTask 7.2: 主流水线回归（3 个 slim circle_identifier + all_slim）stub 编译通过

# Task Dependencies

- Task 1 先行（merge master 是后续构建的基础）
- Task 2 依赖 Task 1（数据随 merge 就位）
- Task 3 独立于 Task 1/2（短读 samplesheet），可并行
- Task 4 依赖 Task 1（server.config 随 merge 就位）
- Task 5 依赖 Task 1（merge 后评估），可与 Task 2/3/4 并行
- Task 6 依赖 Task 5（清单确认）
- Task 7 依赖 Task 2/3/4/6
