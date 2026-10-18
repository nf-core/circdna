# Tasks

## Phase 1 — bio.nf 中完善 ECCsplorer 模块

- [x] Task 1: 检查 kavonrtep/repeatexplorer:2.3.8 镜像结构

  - [x] SubTask 1.1: 拉取镜像，确认 seqclust 安装路径、Perl 依赖目录、R 包位置
  - [x] SubTask 1.2: 记录需复制的文件清单

- [x] Task 2: 重写 ECCsplorer Dockerfile 为多阶段构建

  - [x] SubTask 2.1: 第一阶段 `FROM kavonrtep/repeatexplorer:2.3.8 AS repex`
  - [x] SubTask 2.2: 第二阶段 `condaforge/mambaforge:latest`，COPY --from=repex
  - [x] SubTask 2.3: seqclust → `/opt/conda/envs/eccsplorer/bin/seqclust`
  - [x] SubTask 2.4: Perl/R 依赖 → 对应路径
  - [x] SubTask 2.5: 配置 `lib/config.py` RepeatExplorer 路径
  - [x] SubTask 2.6: 本地构建并验证 `seqclust --help`、`ECCsplorer.py --help`

- [x] Task 3: 推送 Docker 镜像 → `quay.io/bioinfortools/eccsplorer:2022.01.1.1`

- [x] Task 4: 更新 Conda 声明 + nf-test 真实测试
  - [x] SubTask 4.1: 更新 `ECCsplorer/conda-recipe/meta.yaml` 注释
  - [x] SubTask 4.2: 添加非 stub nf-test 用例（testdata: R1/R2/ref.fa, `--mode map`）
  - [x] SubTask 4.3: 运行 `nf-test test` 验证通过

## Phase 2 — circdna.nf 接入 ECCsplorer

- [x] Task 5: 同步 bio.nf ECCsplorer 模块到 circdna.nf

- [x] Task 6: 增加 samplesheet `pair` 列

  - [x] SubTask 6.1: `bin/check_samplesheet.py` OPTIONAL_FIELDS 添加 `pair`
  - [x] SubTask 6.2: `bin/check_samplesheet.py` 读取传递 `pair` 列
  - [x] SubTask 6.3: `subworkflows/local/input_check/main.nf` 填充 `meta.pair`
  - [x] SubTask 6.4: `assets/schema_input.json` 添加 `pair` 字段
  - [x] SubTask 6.5: `subworkflows/local/eccdna_mode/main.nf` 改用 `meta.pair` 配对
  - [x] SubTask 6.6: 更新 `samplesheets/test_local_eccdna.csv` 添加 `pair` 列

- [x] Task 7: ECCsplorer 改由 circle_identifier 控制

  - [x] SubTask 7.1: `nextflow.config` 中 `circle_identifier` 默认值加 `eccsplorer`，移除 `eccsplorer_map_core`
  - [x] SubTask 7.2: `workflows/circdna.nf` 从 circle_identifier 解析 `eccsplorer`→`run_eccsplorer_new`
  - [x] SubTask 7.3: 更新 `conf/test_local.config`、`conf/server.config`、`conf/test.config`

- [x] Task 8: testdata 端到端测试
  - [x] SubTask 8.1: `nextflow run main.nf -profile test_local,docker -resume`
  - [x] SubTask 8.2: 验证 ECCSPLORER 输出 BED/FASTA 非空、流程退出码 0

# Dependencies

Task 2→1, Task 4→3, Task 5→4, Task 6/7→5（可并行）, Task 8→6+7
