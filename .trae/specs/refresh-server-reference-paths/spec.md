# 更新服务器参考基因组路径 Spec

## Why

`circdna.nf` 的服务器 profile 仍将参考基因组定位到旧的扁平目录 `/data1/users/siyangming/FASTA/`，而最新参考文件已经整理到 `/data1/users/siyangming/PublicDB/reference/<species>/`。如果不更新配置和运行文档，服务器运行会继续依赖错误路径，并误导使用者传入与现状不一致的 `--genome` 值。

## What Changes

- 更新 `circdna.nf/conf/server.config` 中服务器 profile 的参考基因组文件定位方式，使其指向 `/data1/users/siyangming/PublicDB/reference/<species>/` 下的最新 `.fa.gz` 文件。
- 逐项核对 `params.genomes` 中每个物种键与用户提供的最新参考文件清单，确保配置引用的是实际存在且已建立 `.fai`/`.gzi` 索引的压缩 FASTA。
- 同步更新 `circdna.nf/SERVER_RUN_GUIDE.md` 中关于参考目录、FASTA 存放位置、物种运行命令和注意事项的描述，保证文档与 `server.config` 一致。
- 对多单倍型或双组装物种的文档说明进行收敛，避免指南继续展示与 `server.config` 不一致的 `--genome` 示例。

## Impact

- Affected specs: 服务器 profile 下的本地参考基因组解析、服务器运行文档的物种命令示例
- Affected code: `circdna.nf/conf/server.config`, `circdna.nf/SERVER_RUN_GUIDE.md`

## ADDED Requirements

### Requirement: 服务器 profile 使用 PublicDB 参考目录

系统 SHALL 在 `server` profile 下从 `/data1/users/siyangming/PublicDB/reference/<species>/` 解析每个受支持物种的参考基因组文件，而不是依赖旧的 `/data1/users/siyangming/FASTA/` 扁平目录。

#### Scenario: 单组装物种成功解析新路径

- **WHEN** 用户在服务器上以 `-profile server` 运行任一单组装物种
- **THEN** `server.config` 中对应的 `params.genomes.<species>.fasta` SHALL 指向该物种目录内实际存在的 `.fa.gz` 文件
- **AND** 该路径 SHALL 与用户提供的最新参考文件清单一致

#### Scenario: 多单倍型物种的键和值保持一致

- **WHEN** 用户查看或使用存在 hap1/hap2 参考文件的物种
- **THEN** `server.config` 中暴露给用户的 genome 键 SHALL 与其所映射的具体 FASTA 文件一一对应
- **AND** 文档 SHALL 不再展示与配置不一致的 `--genome` 值

### Requirement: 服务器运行指南与配置保持一致

系统 SHALL 在 `SERVER_RUN_GUIDE.md` 中使用与 `server.config` 完全一致的参考目录说明、物种示例和注意事项。

#### Scenario: 用户按指南准备环境

- **WHEN** 用户阅读“环境先决条件”或“注意事项”
- **THEN** 文档 SHALL 明确指出参考基因组位于 `/data1/users/siyangming/PublicDB/reference/`
- **AND** 不再要求用户将 FASTA 上传到旧目录 `/data1/users/siyangming/FASTA/`

#### Scenario: 用户复制物种运行命令

- **WHEN** 用户从指南复制任一物种的运行命令
- **THEN** 命令中的 `--genome` 示例 SHALL 对应 `server.config` 中存在的 genome 键
- **AND** 与大基因组附加配置相关的说明 SHALL 保持原有适用范围不变，除非本次配置核对结果显示必须调整

## MODIFIED Requirements

### Requirement: server profile 的本地参考基因组映射

`server.config` 中的 `params.genomes` SHALL 以“物种键 -> 最新 PublicDB 参考 FASTA 文件”的方式维护服务器运行所需映射。对于已迁移到物种子目录的参考文件，配置 SHALL 使用物种目录加文件名的完整相对结构，而不是继续拼接旧的根目录加文件名。

### Requirement: 服务器运行文档中的参考目录说明

`SERVER_RUN_GUIDE.md` SHALL 将所有涉及参考基因组位置的说明统一为最新目录结构，并确保“路径约定”“注意事项”“常见错误”中的路径、示例命令和术语与 `server.config` 当前实现一致。
