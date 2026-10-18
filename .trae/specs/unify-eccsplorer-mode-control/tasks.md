# Tasks

- [ ] Task 1: 更新 workflows/circdna.nf

  - [ ] SubTask 1.1: 解析 `eccsplorer_map`（原 `eccsplorer`）→ `run_eccsplorer_new`
  - [ ] SubTask 1.2: 解析 `eccsplorer_clu` → `run_eccsplorer_clu_new`
  - [ ] SubTask 1.3: `eccsplorer_clu` 在 circle_identifier 中时自动 `run_eccsplorer_new = true`
  - [ ] SubTask 1.4: 移除 `params.eccsplorer_clu` 引用

- [ ] Task 2: 更新配置和元数据

  - [ ] SubTask 2.1: `test_local.config` — 改用 `circle_identifier` 统配
  - [ ] SubTask 2.2: `nextflow.config` — 移除 `eccsplorer_clu = false`
  - [ ] SubTask 2.3: `nextflow_schema.json` — 更新参数定义

- [ ] Task 3: Stub 测试
  - [ ] SubTask 3.1: `--circle_identifier eccsplorer_map` → map 运行
  - [ ] SubTask 3.2: `--circle_identifier eccsplorer_clu` → map + clu 运行

# Dependencies

- Task 2 依赖 Task 1；Task 3 依赖 Task 2
