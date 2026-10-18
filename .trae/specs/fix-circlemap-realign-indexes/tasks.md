# Tasks

- [x] Task 1: 修正 CIRCLEMAP_REALIGN 的 qname BAI 传递

  - [x] SubTask 1.1: 在 `subworkflows/local/circle_map_pipeline/main.nf` 中将 `ch_qname_sorted_bai` 并入 `ch_cm_realign_in` 通道（在 join 链末尾增加 `.join(ch_qname_sorted_bai)`）
  - [x] SubTask 1.2: 更新 `modules/local/circlemap/realign/main.nf` 的 `input:` tuple 定义，从 6 个路径扩展为 7 个：`tuple val(meta), path(re_bam), path(re_bai), path(qname), path(qname_bai), path(sbam), path(sbai)`

- [x] Task 2: 在 CIRCLEMAP_REALIGN 中创建 FASTA .gzi 索引
  - [x] SubTask 2.1: 在 `modules/local/circlemap/realign/main.nf` 的 `script:` 块中，在 `circle_map.py Realign` 调用前添加 `samtools faidx $fasta`

# Task Dependencies

- [Task 1] 和 [Task 2] 无依赖，可并行执行
