# Tasks

- [x] Task 1: 迁移 6 个 bio.nf 模块到 local（物理移动）

  - [x] SubTask 1.1: cresil → local/cresil（含 5 子模块 + testdata）
  - [x] SubTask 1.2: flye → local/flye（含 testdata）
  - [x] SubTask 1.3: repeatexplorer2 → local/repeatexplorer2
  - [x] SubTask 1.4: ampliconsuite → local/ampliconsuite_ec（重命名避免覆盖 prepareaa）
  - [x] SubTask 1.5: eccsplorer → local/eccsplorer（黑盒，含 ECCSPLORER + ECCSPLORER_WITH_CONTROL）
  - [x] SubTask 1.6: ecc_finder/{asm_ont,asm_sr,map_ont,map_sr} + testdata → local/ecc_finder/{...}

- [x] Task 2: ecc_finder 拆分为 nf-core 规范子模块

  - [x] SubTask 2.1: 删除旧单文件 local/ecc_finder/{main.nf,environment.yml,meta.yml}
  - [x] SubTask 2.2: 更新 ecc_finder_pipeline/main.nf 的 4 个 include（→ 4 子模块）；修复 versions emit（versions_ecc_finder）与 tag（meta→meta2）

- [x] Task 3: 更新所有跨仓 include 引用

  - [x] SubTask 3.1: eccsplorer_pipeline/main.nf：AMPLICONSUITE（→local/ampliconsuite_ec）、ECCSPLORER(\_WITH_CONTROL)（→local/eccsplorer）
  - [x] SubTask 3.2: eccsplorer_clu_slim/main.nf + tests/eccsplorer_clu/test.nf：REPEATEXPLORER2（→local/repeatexplorer2）
  - [x] SubTask 3.3: 删除 ECCSPLORER_CLU 断链（include + 调用块 + take 参数）
  - [x] SubTask 3.4: `grep bio.nf/modules` 在 \*.nf 无残留（清理 cresil .nf-test 缓存）

- [x] Task 4: 黑盒版 circle_identifier 接入 workflows

  - [x] SubTask 4.1: 恢复黑盒标识符解析（eccsplorer / ecc_finder_map_sr / asm_sr / map_ont / asm_ont）
  - [x] SubTask 4.2: 黑盒标识符触发 ECCSPLORER_PIPELINE（与 slim 并存）
  - [x] SubTask 4.3: stub 编译验证黑盒 + slim 标识符共存（ECC_FINDER_MAP_SR/ECCSPLORER process 创建成功）

- [x] Task 5: scripts 归档清理

  - [x] SubTask 5.1: prepare_eccsplorer_database.sh → testdatasets/eccsplorer_db/
  - [x] SubTask 5.2: update_samplesheets.py → samplesheets/
  - [x] SubTask 5.3: 删除 test_incremental_cache.py
  - [x] SubTask 5.4: 删除 scripts/ 目录（含 convert_sra_to_fastq_parallel.sh/process_fasta_bgzip.py/路径问题.txt/.DS_Store）

- [x] Task 6: 评估并提交未提交更改

  - [x] SubTask 6.1: 保留 PRExing 链（prepare\_\* + seqkit/concat + ecc_preprocessing）
  - [x] SubTask 6.2: 解决冲突文件（nextflow.config/modules.config/workflows.circdna.nf）：master 为底 + PRExing 内容，丢弃旧黑盒 clu 配置
  - [x] SubTask 6.3: 保留 samplesheet testdata→ngs 路径更新 + testdatasets/ngs 重命名
  - [x] SubTask 6.4: 删除 errors.txt 等冗余文件

- [x] Task 7: 验证
  - [x] SubTask 7.1: stub 编译通过（slim 标识符 + 黑盒标识符）
  - [x] SubTask 7.2: 确认无 bio.nf/modules 跨仓引用残留、无冲突标记
  - [x] SubTask 7.3: 提交所有变更（commit 8cd9991）

# Task Dependencies

- Task 2 依赖 Task 1；Task 3 依赖 Task 1；Task 4 依赖 Task 1/3；Task 5/6 独立；Task 7 依赖全部
