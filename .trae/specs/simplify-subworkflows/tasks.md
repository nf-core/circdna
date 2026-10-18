# Tasks

- [ ] Task 1: 重命名 eccdna_mode → eccsplorer_pipeline

  - [ ] SubTask 1.1: `mv subworkflows/local/eccdna_mode subworkflows/local/eccsplorer_pipeline`
  - [ ] SubTask 1.2: 更新 main.nf 内 workflow 名 `ECCDNA_MODE` → `ECCSPLORER_PIPELINE`
  - [ ] SubTask 1.3: 创建 `eccsplorer_pipeline/meta.yml`
  - [ ] SubTask 1.4: 更新 `workflows/circdna.nf` 中所有 `ECCDNA_MODE` 和 `eccdna_mode` 引用

- [ ] Task 2: 合并 reference_mode 到 ECCSPLORER_PIPELINE

  - [ ] SubTask 2.1: 在 `workflows/circdna.nf` 中，mode='reference' 分支改为调用 `ECCSPLORER_PIPELINE`（所有检测工具 off）
  - [ ] SubTask 2.2: 移除 `REFERENCE_MODE` include

- [ ] Task 3: 补全 ecc_finder_pipeline/meta.yml

- [ ] Task 4: 删除冗余子工作流

  - [ ] SubTask 4.1: 删除 `reference_mode/`、`eccsplorer_all/`、`eccsplorer_cluster/`、`eccsplorer_mapping/`、`eccsplorer_prepare/`

- [ ] Task 5: Stub 测试验证（--mode reference 和 --circle_identifier eccsplorer,ecc_finder_map_sr）

# Dependencies

- Task 2 依赖 Task 1；Task 4 依赖 Task 2；Task 5 依赖 Task 3/4
