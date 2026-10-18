# Tasks

- [x] Task 1: 修复 TIDEHUNTER 模块 CLI

  - [x] 1.1 main.nf：input 改为 reads（fasta/fastq），output 改为 `${prefix}.fasta`，命令 `TideHunter ${reads} > ${prefix}.fasta`
  - [x] 1.2 更新 stub 与 versions
  - [x] 1.3 更新 meta.yml
  - [x] 1.4 stub 编译验证

- [x] Task 2: 提取 clu_prepare 脚本与模块

  - [x] 2.1 从 eccPrepare.py prepare_for_clustering 提取 `clu_prepare.py`（读长统一/子采样/前缀拼接 → REPEATEXPLORER_READY.fa）
  - [x] 2.2 加入 eccsplorer_slim 镜像并重建推送
  - [x] 2.3 新建 `modules/local/eccsplorer_slim/clu_prepare/main.nf`（FASTQ 对 + taxon → REPEATEXPLORER_READY.fa）

- [x] Task 3: 新建 REPEATEXPLORER2 模块（bio.nf）

  - [x] 3.1 新建 `bio.nf/modules/repeatexplorer2/main.nf`（seqclust 聚类，镜像 docker.1ms.run/kavonrtep/repeatexplorer:2.3.8，命令 /repex_tarean/seqclust）
  - [x] 3.2 输入 REPEATEXPLORER_READY.fa + taxon，输出 seqclust 聚类目录
  - [x] 3.3 新建 meta.yml / environment.yml
  - [x] 3.4 stub 编译验证

- [x] Task 4: 提取 clu_candidates 脚本与模块

  - [x] 4.1 提取 `clu_candidates.py`（解析 seqclust 输出 → cluster_candidates.csv / comparative_cluster_table.csv）
  - [x] 4.2 加入 eccsplorer_slim 镜像并重建推送
  - [x] 4.3 新建 `modules/local/eccsplorer_slim/clu_candidates/main.nf`

- [x] Task 5: 新建 eccsplorer_clu_slim 子工作流并接入

  - [x] 5.1 新建 `subworkflows/local/eccsplorer_clu_slim/main.nf`（clu_prepare → repex_tarean → clu_candidates）
  - [x] 5.2 workflows/circdna.nf circle_identifier 新增 `eccsplorer_clu_slim`（自动启用 map）
  - [x] 5.3 conf/modules.config 添加新模块 publishDir
  - [x] 5.4 stub 编译验证

- [x] Task 6: 独立测试模块

  - [x] 6.1 生成模拟串联重复长读 `tests/tidehunter/data/sim_tandem.fa`
  - [x] 6.2 创建 `tests/tidehunter/test.nf`，运行验证 consensus fasta 非空含重复单元
  - [x] 6.3 创建 `tests/eccsplorer_clu/test.nf`（circdna_1 + gdna_1），运行验证 cluster 输出非空
  - [x] 6.4 输出 publish 到 results/slim_run/tidehunter_test/ 与 results/slim_run/eccsplorer_clu_test/

- [x] Task 7: 回归验证
  - [x] 7.1 主流水线 -resume（3 个 slim circle_identifier）仍成功
  - [x] 7.2 更新 modularize-eccsplorer-modes（Task 14/15）与 implement-eccsplorer-clu-all（Task 5）勾选

# Task Dependencies

- [Task 3] 依赖 [Task 2]；[Task 5] 依赖 [Task 2]+[Task 3]+[Task 4]
- [Task 6] 依赖 [Task 1]+[Task 5]
- [Task 7] 依赖 [Task 6]
