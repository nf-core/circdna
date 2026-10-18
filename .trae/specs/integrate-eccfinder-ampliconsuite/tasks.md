# Tasks

- [ ] Task 1: 创建分支并拷贝 ecc_finder_asm_sr 模块

  - [ ] SubTask 1.1: `git checkout -b integrate-eccfinder-ampliconsuite`
  - [ ] SubTask 1.2: 从 `bio.nf/modules/ecc_finder/asm_sr/` 拷贝到 `circdna.nf/modules/local/ecc_finder_asm_sr/`

- [ ] Task 2: 接入 ECCDNA_MODE 子工作流

  - [ ] SubTask 2.1: include `ECC_FINDER_MAP_SR`、`ECC_FINDER_ASM_SR`、`AMPLICONSUITE`
  - [ ] SubTask 2.2: 添加 `run_ecc_finder_map_sr`、`run_ecc_finder_asm_sr`、`run_ampliconsuite` 参数
  - [ ] SubTask 2.3: `run_ecc_finder_map_sr=true` 时调用 MAP_SR（sorted BAM + FASTA）
  - [ ] SubTask 2.4: `run_ecc_finder_asm_sr=true` 时调用 ASM_SR（FASTQ 直接传入）
  - [ ] SubTask 2.5: `run_ampliconsuite=true` 时调用 AMPLICONSUITE（sorted BAM + mosek + aa_data_repo）

- [ ] Task 3: 更新 workflows/circdna.nf

  - [ ] SubTask 3.1: 从 circle_identifier 解析 `ecc_finder_map_sr` / `ecc_finder_asm_sr` / `ampliconsuite`
  - [ ] SubTask 3.2: 传递三个布尔参数到 ECCDNA_MODE 调用

- [ ] Task 4: 更新配置

  - [ ] SubTask 4.1: `conf/modules.config` — 添加三个 publishDir 配置
  - [ ] SubTask 4.2: `conf/test_local.config` — 添加 ampliconsuite 所需参数（aa_data_repo 等）

- [ ] Task 5: test_local 测试（-resume）
  - [ ] SubTask 5.1: `nextflow run main.nf -profile test_local,docker -resume --circle_identifier eccsplorer,ecc_finder_map_sr,ecc_finder_asm_sr`
  - [ ] SubTask 5.2: 验证 MAP_SR 和 ASM_SR 输出非空
  - [ ] SubTask 5.3: 流程退出码 0

# Dependencies

- Task 2 依赖 Task 1；Task 3/4 与 Task 2 可并行；Task 5 依赖 Task 2/3/4
