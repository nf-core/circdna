# Tasks

- [ ] Task 1: 合并 4 个 ecc_finder 子模块到单一 main.nf

  - [ ] SubTask 1.1: 从 bio.nf 拷贝 `map_ont/main.nf`、`asm_ont/main.nf` 内容
  - [ ] SubTask 1.2: 创建 `modules/local/ecc_finder/main.nf`，包含 4 个 process（MAP_SR → ASM_SR → MAP_ONT → ASM_ONT）
  - [ ] SubTask 1.3: 统一修正 tag 和 output meta（MAP_SR 的 `$idx` → `$meta`，output 用 `meta2`）
  - [ ] SubTask 1.4: 创建 `modules/local/ecc_finder/environment.yml`（4 份完全一致，取一份）+ `meta.yml`
  - [ ] SubTask 1.5: 删除 `modules/local/ecc_finder_map_sr/` 和 `modules/local/ecc_finder_asm_sr/`

- [ ] Task 2: 创建 ECC_FINDER_PIPELINE 子工作流

  - [ ] SubTask 2.1: 创建 `subworkflows/local/ecc_finder_pipeline/main.nf`
  - [ ] SubTask 2.2: 输入 reads + bwa_index + fasta_meta + run_map + run_asm + platform('sr')
  - [ ] SubTask 2.3: platform=sr + run_map → ECC_FINDER_MAP_SR（BWA idx+R1+R2+ref）
  - [ ] SubTask 2.4: platform=sr + run_asm → ECC_FINDER_ASM_SR（R1+R2）
  - [ ] SubTask 2.5: platform=ont + run_map → ECC_FINDER_MAP_ONT（预留，需 minimap2 idx）
  - [ ] SubTask 2.6: platform=ont + run_asm → ECC_FINDER_ASM_ONT（预留，单端 query）

- [ ] Task 3: 更新 eccdna_mode + workflows/circdna.nf

  - [ ] SubTask 3.1: eccdna_mode 中用 `include ECC_FINDER_PIPELINE` 替代 2 个独立 include
  - [ ] SubTask 3.2: 传递 reads + bwa_index + fasta_meta + run_map + run_asm 到子工作流
  - [ ] SubTask 3.3: 更新 modules.config 中 withName 引用

- [ ] Task 4: Stub 测试验证
  - [ ] SubTask 4.1: `-profile test_local,docker -stub --circle_identifier ecc_finder_map_sr,ecc_finder_asm_sr` 通过
  - [ ] SubTask 4.2: 输出文件中 `ecc_finder_map_sr/` 和 `ecc_finder_asm_sr/` 目录存在

# Dependencies

- Task 2 依赖 Task 1；Task 3 依赖 Task 2；Task 4 依赖 Task 3
