# Tasks

- [x] Task 1: 提交 circdna.nf 的 5 个修改文件到 GitHub

  - [x] git add SERVER_RUN_GUIDE.md AGENTS.md conf/server.config errors.txt samplesheets/data_issues.txt
  - [x] git commit -m "fix(docker): add runOptions for user mapping and sync server paths"
  - [x] git push

- [x] Task 2: 生成容器权限审计报告
  - [x] 检查 circdna.nf 的容器权限配置（nextflow.config + conf/server.config）
  - [x] 检查 circrna.nf 的容器权限配置（nextflow.config + conf/server.config）
  - [x] 检查 bio.nf 的容器权限配置
  - [x] 检查 nanoseq.nf 的容器权限配置
  - [x] 检查 isoseq.nf 的容器权限配置
  - [x] 检查 fetchngs.nf 的容器权限配置
  - [x] 检查 riboseq.nf 的容器权限配置
  - [x] 检查 rnaseq 的容器权限配置
  - [x] 汇总生成报告表格

# Task Dependencies

- Task 1 和 Task 2 互相独立，可并行执行
