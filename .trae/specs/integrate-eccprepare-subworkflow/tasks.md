# Tasks

- [x] Task 1: 拷贝 seqkit/concat nf-core 模块到 circdna.nf。

  - [x] SubTask 1.1: 从 `/Users/siyangming/nextflow_nf_core/modules/modules/nf-core/seqkit/concat/` 拷贝到 `circdna.nf/modules/nf-core/seqkit/concat/`（含 main.nf, meta.yml, environment.yml, tests/）。

- [x] Task 2: 创建 3 个 bin 自定义脚本（samplesheet_check 模式：bin/ 纯 Python + process 包装，共用 `quay.io/biocontainers/biopython:1.84`）。

  - [x] SubTask 2.1: 创建 `bin/eccsplorer_prepare_read_length.py`（numpy 替代 scipy 的 `get_best_read_length()`），创建 `modules/local/eccsplorer_prepare_read_length/main.nf`，含 stub。
  - [x] SubTask 2.2: 创建 `bin/eccsplorer_prepare_read_count.py`（`get_max_read_count()`，Bio.SeqIO PE 配对计数），创建 `modules/local/eccsplorer_prepare_read_count/main.nf`，含 stub。
  - [x] SubTask 2.3: 创建 `bin/eccsplorer_prepare_prexing.py`（`prexing_reads()`：PE 配对子采样+截断+前缀 `_#0/1` `_#0/2`），创建 `modules/local/eccsplorer_prepare_prexing/main.nf`，含 stub。

- [x] Task 3: 创建子流程 `subworkflows/local/ecc_preprocessing/main.nf`。

  - [x] SubTask 3.1: include SEQTK_SEQ（现有）、3 个 custom modules、SEQKIT_CONCAT（新拷贝）。
  - [x] SubTask 3.2: 编排数据流：SEQTK_SEQ → ECCSplorer_PREPARE_READ_LENGTH → ECCSplorer_PREPARE_READ_COUNT per sample → ECCSplorer_PREPARE_PREXING per sample → SEQKIT_CONCAT。
  - [x] SubTask 3.3: 定义 take（trimmed_reads）/ emit（prexed_fasta, versions）接口。

- [x] Task 4: 接入 workflows/circdna.nf 并更新配置。

  - [x] SubTask 4.1: `nextflow.config` 新增 `run_eccprepare`（bool, default false）、`eccprepare_fold_cov`（0.1）、`eccprepare_target_read_count`（'auto'）。
  - [x] SubTask 4.2: `conf/modules.config` 为 3 个新 process + SEQKIT_CONCAT 添加资源标签。
  - [x] SubTask 4.3: `workflows/circdna.nf` 新增 include 和条件调用分支（`if (params.run_eccprepare)`）。

- [x] Task 5: 两轮 test_local + resume 验证最小化重跑成本。

  - [x] SubTask 5.0: **冻结其他修改**。项目无 git 仓库，跳过 stash 步骤。
  - [x] SubTask 5.1: 第一轮基线 `nextflow run main.nf -profile test_local,docker --mode eccdna` 51 任务全部成功。
  - [x] SubTask 5.2: 第二轮 `--run_eccprepare -resume` 执行成功。`Pipeline completed successfully`，CPU hours: 4.3 (99% cached)。
  - [x] SubTask 5.3: 检查 trace：已有 FASTQC/TRIMGALORE/BWA*INDEX/BWA_MEM/ECCsplorer/CIRCLE_MAP 等全部 CACHED；新模块 SEQTK_SEQ_R1/R2 + ECCSPLORER_PREPARE*\* + SEQKIT_CONCAT 为 NEW。
  - [x] SubTask 5.4: 无 git 仓库，无需恢复。

- [x] Task 6: 更新 CHANGELOG.md。
  - [x] SubTask 6.1: 记录新增 ecc_preprocessing 子流程及 3 个 bin 脚本（MINOR bump，v4.5.0）。

# Task Dependencies

- Task 2 depends on Task 1
- Task 3 depends on Task 1, Task 2
- Task 4 depends on Task 3
- Task 5 depends on Task 4
- Task 6 depends on Task 5
