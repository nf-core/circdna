# Tasks

- [x] Task 1: 替换 quay.io Python 容器镜像

  - [x] SubTask 1.1: 将 `modules/local/eccsplorer/main.nf` 的 container 从 `python:3.9-slim` 替换为 `quay.io/biocontainers/python:3.9--py39`
  - [x] SubTask 1.2: 将 `modules/local/candidate_merge/main.nf` 的 container 替换为 `quay.io/biocontainers/python:3.9--py39`
  - [x] SubTask 1.3: 将 `modules/local/ecc_score/main.nf` 的 container 替换为 `quay.io/biocontainers/python:3.9--py39`

- [x] Task 2: 为新模块添加 publishDir

  - [x] SubTask 2.1: 在 `modules/nf-core/mosdepth/main.nf` 添加 `publishDir "${params.outdir}/reference_mode/mosdepth", mode:'copy', enabled:true`
  - [x] SubTask 2.2: 在 `modules/local/eccsplorer/main.nf` 添加 `publishDir "${params.outdir}/eccdna_mode/eccsplorer"`
  - [x] SubTask 2.3: 在 `modules/local/candidate_merge/main.nf` 添加 `publishDir "${params.outdir}/eccdna_mode/candidate_merge"`
  - [x] SubTask 2.4: 在 `modules/local/ecc_score/main.nf` 添加 `publishDir "${params.outdir}/integrated_mode/ecc_score"`

- [x] Task 3: 修复 CIRCLE_MAP_PIPELINE 缺失 bed emit

  - [x] SubTask 3.1: 在 `subworkflows/local/circle_map_pipeline/main.nf` 添加 `bed` 输出 channel，合并 CIRCLEMAP_REALIGN 和 CIRCLEMAP_REPEATS 的 bed 输出

- [x] Task 4: 修复 Integrated Mode channel join 逻辑

  - [x] SubTask 4.1: 修改 `subworkflows/local/integrated_mode/main.nf`，将 gDNA mosdepth 作为 value channel 传入（collect），不再依赖 meta.id join

- [x] Task 5: 添加重复序列注释集成

  - [x] SubTask 5.1: 确认 `reference_mode/main.nf` 正确接收 `repeat_gff` 参数并传递
  - [x] SubTask 5.2: 确认 `integrated_mode/main.nf` 正确将 repeat_bed 传递给 ECC_SCORE 模块
  - [x] SubTask 5.3: 确认 `ecc_score/main.nf` 在 repeat_bed 为 null 时跳过 TE penalty

- [x] Task 6: MultiQC 报告集成

  - [x] SubTask 6.1: 在 `workflows/circdna.nf` 的 MultiQC channel 收集逻辑中添加 mosdepth summary 文件
  - [x] SubTask 6.2: 在 `assets/multiqc_config.yml` 添加 mosdepth 和 eccDNA 统计的 section 配置

- [x] Task 7: 生成 gDNA 测试数据

  - [x] SubTask 7.1: 使用 wgsim 从 `testdatasets/reference/genome.fa` 生成 gDNA paired-end reads
  - [x] SubTask 7.2: 压缩为 `gdna_1_R1.fastq.gz` / `gdna_1_R2.fastq.gz` 放入 `testdatasets/testdata/`

- [x] Task 8: 创建测试 samplesheet 和配置

  - [x] SubTask 8.1: 创建 `testdatasets/samplesheet/samplesheet_integrated.csv`，包含 gDNA 和 eccDNA 样本，带 `datatype` 字段
  - [x] SubTask 8.2: 创建 `conf/test_integrated.config` 配置文件
  - [x] SubTask 8.3: 更新 `conf/test_local.config` 添加 `mode` 参数

- [x] Task 9: 端到端测试验证
  - [x] SubTask 9.1: 运行 `nextflow run main.nf -profile test_local,docker --mode reference --outdir <OUTDIR> -stub` 验证 Reference Mode
  - [x] SubTask 9.2: 运行 `nextflow run main.nf -profile test_local,docker --mode eccdna --outdir <OUTDIR> -stub` 验证 eccDNA Mode
  - [x] SubTask 9.3: 运行 `nextflow run main.nf -profile test_integrated,docker --outdir <OUTDIR> -stub` 验证 Integrated Mode

# Task Dependencies

- Task 3 需在 Task 9 之前完成（eccDNA Mode 依赖 bed 输出）
- Task 4 需在 Task 9 之前完成（Integrated Mode 依赖 join 修复）
- Task 7 需在 Task 8 之前完成（samplesheet 引用 gDNA 数据文件）
- Task 1, 2, 5, 6 可并行执行
