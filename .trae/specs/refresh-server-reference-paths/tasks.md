# Tasks

- [x] Task 1: 核对服务器 profile 的 genome 键与最新参考目录清单，明确哪些物种只需改路径，哪些物种还需要同步文档中的 genome 命名。

  - [x] SubTask 1.1: 对照 `circdna.nf/conf/server.config` 与用户提供的 `/data1/users/siyangming/PublicDB/reference/` 清单，整理每个 genome 键当前指向的文件名与最新目录位置。
  - [x] SubTask 1.2: 标记存在单倍型分支或双组装文件的物种，确认文档最终应该暴露哪些 `--genome` 值。

- [x] Task 2: 更新 `circdna.nf/conf/server.config` 的参考基因组定位，使服务器 profile 只引用最新 PublicDB 路径下实际存在的 `.fa.gz` 文件。

  - [x] SubTask 2.1: 修改 `fasta_base_path` 或各物种 `fasta` 字段的拼接方式，使其适配 `/data1/users/siyangming/PublicDB/reference/<species>/` 目录结构。
  - [x] SubTask 2.2: 逐项检查受支持物种的 FASTA 文件名是否与最新清单一致，避免继续引用旧文件名或旧目录。
  - [x] SubTask 2.3: 如多单倍型物种的 genome 键需要调整，确保配置与后续文档同步收敛到同一套键名。

- [x] Task 3: 更新 `circdna.nf/SERVER_RUN_GUIDE.md` 中与参考基因组位置相关的说明和示例。

  - [x] SubTask 3.1: 将“环境先决条件”“注意事项”“常见错误”中所有旧的 FASTA 目录说明替换为新的 PublicDB 参考目录。
  - [x] SubTask 3.2: 修订受影响物种的运行示例，确保每个 `--genome` 示例都能在 `server.config` 中找到对应键。
  - [x] SubTask 3.3: 保留现有 screen、resume 和大基因组运行说明，仅在与参考文件路径或 genome 键直接相关时做必要调整。

- [x] Task 4: 完成一次配置与文档的一致性验证，确认不存在残留旧路径或失配示例。
  - [x] SubTask 4.1: 检查 `server.config` 中每个受支持物种的 `fasta` 条目都能映射到用户提供的最新文件清单。
  - [x] SubTask 4.2: 搜索并清理 `SERVER_RUN_GUIDE.md` 与 `server.config` 中残留的 `/data1/users/siyangming/FASTA` 旧路径引用。
  - [x] SubTask 4.3: 复核文档中出现的所有 `--genome` 示例值都与 `params.genomes` 的键完全一致。

# Task Dependencies

- Task 2 depends on Task 1
- Task 3 depends on Task 1
- Task 4 depends on Task 2
- Task 4 depends on Task 3
