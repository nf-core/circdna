# Tasks

- [ ] Task 1: 生成 BWA index

  - [ ] SubTask 1.1: `bwa index testdatasets/reference/genome.fa`

- [ ] Task 2: 运行全流程

  - [ ] SubTask 2.1: `nextflow run main.nf -profile test_local,docker --circle_identifier eccsplorer,ecc_finder_map_sr,ecc_finder_asm_sr`
  - [ ] SubTask 2.2: 使用 -resume 加速（如有缓存）

- [ ] Task 3: 验证输出
  - [ ] 检查各检测工具输出文件非空
  - [ ] 确认退出码 0

# Dependencies

- Task 2 依赖 Task 1；Task 3 依赖 Task 2
