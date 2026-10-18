# Tasks

- [ ] Task 1: 合并 eccsplorer 和 eccsplorer_clu 模块

  - [ ] SubTask 1.1: 将 `ECCSPLORER_CLU` process 代码从 `eccsplorer_clu/main.nf` 复制到 `eccsplorer/main.nf` 末尾
  - [ ] SubTask 1.2: 更新 `subworkflows/local/eccdna_mode/main.nf` 中 include 路径：`eccsplorer_clu/main` → `eccsplorer/main`
  - [ ] SubTask 1.3: 删除 `modules/local/eccsplorer_clu/` 整个目录
  - [ ] SubTask 1.4: 验证语法（`nextflow config . 2>/dev/null` 或 nextflow lint）

- [ ] Task 2: ampliconsuite 迁移到 bio.nf + 官方镜像
  - [ ] SubTask 2.1: 移动 `circdna.nf/modules/local/ampliconsuite/` → `bio.nf/modules/ampliconsuite/`
  - [ ] SubTask 2.2: 更新 `main.nf` 第 6 行 container 为 `quay.io/biocontainers/ampliconsuite:1.6.0--pyh109da93_0`
  - [ ] SubTask 2.3: 更新 `environment.yml`：依赖改为 `bioconda::ampliconsuite=1.6.0`，移除 mosek channel/依赖
  - [ ] SubTask 2.4: 删除 `Dockerfile`
  - [ ] SubTask 2.5: 更新 `subworkflows/local/ampliconarchitect_pipeline/main.nf` 中 include 路径
  - [ ] SubTask 2.6: 删除 `modules/local/ampliconsuite/` 目录（移动后清理）

# Dependencies

- Task 1 和 Task 2 相互独立，可并行
