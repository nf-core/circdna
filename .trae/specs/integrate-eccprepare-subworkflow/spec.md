# 将 eccPrepare.py 预处理拆分为 Nextflow 子流程 Spec

## Why

`ECCsplorer/lib/eccPrepare.py` 的预处理逻辑可在不依赖 eccsplorer 大容器的情况下模块化。通过组合 nf-core 模块（seqtk/seq、seqkit/concat）和 3 个轻量自定义脚本（全部共用 `quay.io/biocontainers/biopython:1.84`），以子流程接入并验证最小化重跑成本。

## What Changes

- 复用现有 nf-core 模块：`seqtk/seq`（已有）、`seqkit/concat`（从 modules 仓库拷贝）
- 3 个自定义 `bin/` 脚本，遵循 `samplesheet_check` 模式，**全部使用 `quay.io/biocontainers/biopython:1.84`（106MB，含 numpy）**
- `get_best_read_length()` 中 `scipy.optimize.fmin(xtol=10000)` 精度极低，用 `numpy.median` 或网格搜索替代，消除 scipy 依赖
- 新增 `--run_eccprepare` 开关（默认 false），第二轮 `-resume` 验证已有任务全 CACHED

## 模块镜像选型

| Process                          | Python 依赖       | 镜像                                                             | 大小  |
| -------------------------------- | ----------------- | ---------------------------------------------------------------- | ----- |
| `SEQTK_SEQ`                      | C 工具            | `quay.io/biocontainers/seqtk:1.4--he4a0461_1`（已有）            | ~1MB  |
| `ECCSplorer_PREPARE_READ_LENGTH` | numpy + biopython | `quay.io/biocontainers/biopython:1.84`                           | 106MB |
| `ECCSplorer_PREPARE_READ_COUNT`  | biopython         | `quay.io/biocontainers/biopython:1.84`                           | 106MB |
| `ECCSplorer_PREPARE_PREXING`     | numpy + biopython | `quay.io/biocontainers/biopython:1.84`                           | 106MB |
| `SEQKIT_CONCAT`                  | Go 工具           | `community.wave.seqera.io/library/seqkit:2.13.0`（与源模块一致） | ~20MB |

**全部自定义模块共享同一个 `biopython:1.84` 镜像**，不引入新镜像，不依赖 eccsplorer 大容器。

## Impact

- 拷贝 `modules/nf-core/seqkit/concat/`
- 新建 `bin/eccsplorer_prepare_read_length.py` + `modules/local/eccsplorer_prepare_read_length/main.nf`
  - 新建 `bin/eccsplorer_prepare_read_count.py` + `modules/local/eccsplorer_prepare_read_count/main.nf`
  - 新建 `bin/eccsplorer_prepare_prexing.py` + `modules/local/eccsplorer_prepare_prexing/main.nf`
- 新建 `subworkflows/local/ecc_preprocessing/main.nf`
- 修改 `workflows/circdna.nf`、`conf/modules.config`、`nextflow.config`、`CHANGELOG.md`

## ADDED Requirements

### Requirement: 读长优化模块

`bin/eccsplorer_prepare_read_length.py` 用 numpy 替代 scipy 计算优化读长（原 `ftol=10000, xtol=10000` 精度极低，numpy.median 或网格搜索等价）。输入 `[meta, fasta]`，输出 `[meta, val(optimal_length)]`。镜像：`quay.io/biocontainers/biopython:1.84`。

### Requirement: PE 读段计数模块

`bin/eccsplorer_prepare_read_count.py`（`get_max_read_count()` 逻辑）。输入 `[meta, fasta1, fasta2, val(best_read_length)]`，输出 `[meta, val(count)]`。镜像同上。

### Requirement: 子采样截断前缀模块

`bin/eccsplorer_prepare_prexing.py`（PE 配对子采样+截断+前缀+`_#0/1` `_#0/2` 后缀，种子 12）。输入 `[meta, fasta1, fasta2, val(best_read_length), val(read_count)]`，输出 `[meta, temp_fasta]`。镜像同上。

### Requirement: ECCSplorer 预处理子流程

```
ch_trimmed_reads
  ├──► SEQTK_SEQ ──► PE FASTA
  │       ├──► ECCSplorer_PREPARE_READ_LENGTH ──► 全局 min best_read_length
  │       └──► ECCSplorer_PREPARE_READ_COUNT per sample ──► 全局 min read_count
  │                  └──► ECCSplorer_PREPARE_PREXING per sample ──► temp_fasta
  │                            └──► SEQKIT_CONCAT ──► REPEATEXPLORER_READY.fa
```

全部自定义 process 共用 `biopython:1.84` 镜像。

#### Scenario: 缓存隔离

- **WHEN** 第一轮 baseline，第二轮 `--run_eccprepare -resume`
- **THEN** 已有任务全部 CACHED，新模块 NEW

## 变更隔离原则（冻结其他修改）

本 spec 仅关注 ecc_preprocessing 子流程的新增接入。测试流程：

1. **测试前**：`git stash` 冻结工作区中非本次变更的所有未提交修改，确保工作区仅含 eccprepare 代码，记录被冻结内容清单
2. **测试中**：遇非本次变更引入的问题（其他流程 bug、已有模块失败、配置不兼容等）一律冻结不修，仅记录
3. **测试后**：`git stash pop` 恢复被冻结的修改

## MODIFIED / REMOVED

无。
