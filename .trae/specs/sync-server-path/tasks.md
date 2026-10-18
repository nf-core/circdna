# Tasks

- [x] Task 1: 更新 SERVER_RUN_GUIDE.md 中的服务器路径

  - [x] 将所有 `cd /data1/users/siyangming/nextflow_nf_core/circdna.nf/` 替换为 `cd /data1/users/siyangming/PlanteccDNADB/circdna.nf/`（约16处）
  - [x] 将 rsync 源路径中本地 Mac 路径保持不变（`/Users/siyangming/nextflow_nf_core/circdnalr.nf/FASTA/`）
  - [x] 验证所有 `cd` 命令指向新路径，`--outdir` 路径（`/data1/users/siyangming/PlanteccDNADB/eccDNA_results/`）已正确无需修改

- [x] Task 2: 更新 samplesheets/data_issues.txt 中的服务器路径

  - [x] 将第51行 `cd /data1/users/siyangming/nextflow_nf_core/circdna.nf` 替换为新路径

- [x] Task 3: 更新 errors.txt 中的服务器路径

  - [x] 将第126行 `/data1/users/siyangming/nextflow_nf_core/circdna.nf/work/` 替换为新路径

- [x] Task 4: 更新 circdna.nf/AGENTS.md 中的跨项目路径表
  - [x] 更新第518-520行中的项目路径表（circdna.nf、circrna.nf、bio.nf 的服务器路径描述）

# Task Dependencies

- Task 1, 2, 3, 4 互相独立，可并行执行
