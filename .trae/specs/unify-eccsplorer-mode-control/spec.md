# 统一 ECCsplorer 模式控制 Spec

## Why

当前 `eccsplorer_clu` 通过独立的 `params.eccsplorer_clu` 布尔值控制，与 `circle_identifier` 字符串体系不统一。改为统一从 `circle_identifier` 解析，移除独立参数。

## What Changes

- 全部模式通过 `circle_identifier` 控制
- 两标识符：`eccsplorer_map` / `eccsplorer_clu`
- 移除 `params.eccsplorer_clu` 独立参数

### 语义（基于 ECCsplorer 原始工具模式）

```
eccsplorer_map  → 仅 map 检测（每样本独立）
eccsplorer_clu  → map + 聚类 + 比较分析（全流程）
```

`clu` 模式在 ECCsplorer 原始工具中包含 seqclust 聚类 + 跨样本比较分析。`ECCSPLORER_CLU` 模块输出 `*_comparative_cluster_table.csv` 即为 compare 产物。`eccsplorer_clu` 自动启用 map。

## Impact

- 修改: `workflows/circdna.nf` — 解析 `eccsplorer_map` / `eccsplorer_clu`
- 修改: `conf/test_local.config` — 改用 `circle_identifier` 统配
- 修改: `nextflow.config` — 移除 `eccsplorer_clu = false`
- 修改: `nextflow_schema.json` — 更新参数定义
- **BREAKING**：`--eccsplorer_clu true` 不再生效，需改为 `--circle_identifier eccsplorer_clu`

## ADDED Requirements

### Requirement: circle_identifier 统一控制

| 值               | 行为              | 产出                                                                                   |
| ---------------- | ----------------- | -------------------------------------------------------------------------------------- |
| `eccsplorer_map` | 仅 map            | `ecc_sequences.fasta`, candidates BED                                                  |
| `eccsplorer_clu` | map + 聚类 + 比较 | 上述 + `cluster_candidates.csv`, `comparative_cluster_table.csv`, `eccCL_summary.html` |
| 无上述值         | ECCsplorer 不运行 | -                                                                                      |

## REMOVED Requirements

### Requirement: params.eccsplorer_clu

**Reason**: 统一到 circle_identifier
**Migration**: `--eccsplorer_clu true` → `--circle_identifier eccsplorer_clu`
