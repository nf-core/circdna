# Tasks

- [ ] Task 1: 创建 ECCSPLORER_CLU 模块

  - [ ] SubTask 1.1: 创建 `circdna.nf/modules/local/eccsplorer_clu/main.nf`
  - [ ] SubTask 1.2: 输入: `tuple val(meta), path(treatment_reads), path(control_reads), val(taxon)`
  - [ ] SubTask 1.3: 调用 `ECCsplorer.py <R1> <R2> <C1> <C2> --mode clu --taxon <taxon>`，提取 cluster candidates/cluster table/HTML summary
  - [ ] SubTask 1.4: 创建 `environment.yml`（引用 `yangmingsi::eccsplorer=2022.01.1.1`）+ `meta.yml`
  - [ ] SubTask 1.5: container 指向 `quay.io/bioinfortools/eccsplorer:2022.01.1.1`

- [ ] Task 2: 清理 stub 代码

  - [ ] SubTask 2.1: 删除 `modules/local/eccsplorer_clu_prepare/`
  - [ ] SubTask 2.2: 删除 `modules/local/eccsplorer_clu_core/`
  - [ ] SubTask 2.3: 删除 `modules/local/eccsplorer_clu_candidates_plot/`
  - [ ] SubTask 2.4: 删除 `subworkflows/local/eccsplorer_cluster/`
  - [ ] SubTask 2.5: 删除 `subworkflows/local/eccsplorer_all/`

- [ ] Task 3: 接入 eccdna_mode 子工作流

  - [ ] SubTask 3.1: include `ECCSPLORER_CLU` 模块
  - [ ] SubTask 3.2: 当 `params.eccsplorer_clu` 为 true 且 `run_eccsplorer` 为 true 时，在 map 完成后调用 ECCSPLORER_CLU
  - [ ] SubTask 3.3: ECCSPLORER_CLU 接收 eccDNA 的 FASTQ + 配对的 gDNA FASTQ（按 pair 列）

- [ ] Task 4: 更新配置

  - [ ] SubTask 4.1: `conf/modules.config` — 添加 ECCSPLORER_CLU 资源标签 `process_high` + publishDir
  - [ ] SubTask 4.2: `conf/test_local.config` — 添加 `eccsplorer_database` 指向 PublicDB、`eccsplorer_clu=true`
  - [ ] SubTask 4.3: `nextflow.config` — 确保 `eccsplorer_clu` 参数定义（已有，默认 false）

- [ ] Task 5: 真实数据测试（使用 -resume）
  - [ ] SubTask 5.1: 运行 `nextflow run main.nf -profile test_local,docker --eccsplorer_clu true -resume`（注：完整版 clu 未在本轮运行；slim clu 已通过 tests/eccsplorer_clu/test.nf 与 circle_identifier=eccsplorer_clu_slim 验证，见 test-tidehunter-module spec）
  - [ ] SubTask 5.2: 验证 map 阶段缓存命中（不重跑）
  - [ ] SubTask 5.3: 验证 ECCSPLORER_CLU 成功执行
  - [ ] SubTask 5.4: 验证输出 cluster candidates 文件非空

# Dependencies

- Task 2、Task 4 可与 Task 1 并行
- Task 3 依赖 Task 1、Task 2
- Task 5 依赖 Task 3、Task 4
