# Fix `first` on Value Channel Warning Spec

## Why

运行 `nextflow run main.nf -profile test_local,docker --mode eccdna` 时，Nextflow 输出警告：

```
WARN: The operator `first` is useless when applied to a value channel which returns a single value by definition
```

该警告源于 `subworkflows/local/eccdna_mode/main.nf` 第 50 行对已经是 value channel 的 `fasta_meta.map{...}` 结果再次调用 `.first()`。需要消除该警告，同时不破坏现有的 channel 广播语义。

## What Changes

- 在 `subworkflows/local/eccdna_mode/main.nf` 第 50 行，移除多余的 `.first()` 调用：
  - 修改前：`ch_eccsplorer_fasta = fasta_meta.map { meta, fasta -> fasta }.first()`
  - 修改后：`ch_eccsplorer_fasta = fasta_meta.map { meta, fasta -> fasta }`
- **不修改**其他位置的 `.first()` 调用（`bam_preprocessing/main.nf` 第 39 行、`workflows/circdna.nf` 第 50/81 行），因为它们作用于 queue channel，不会触发警告。
- 更新 `CHANGELOG.md`，记录该修复。

## Impact

- Affected specs: 无（独立修复）
- Affected code:
  - `subworkflows/local/eccdna_mode/main.nf` — 第 50 行 `ch_eccsplorer_fasta` 构造逻辑
  - `CHANGELOG.md` — 新增 PATCH 条目

## Root Cause Analysis

### Channel 类型推导

```
workflows/circdna.nf 第 50 行:
  ch_fasta_meta = ch_fasta.map{...}.first()   ← queue → value channel ✔

subworkflows/local/eccdna_mode/main.nf 第 50 行:
  ch_eccsplorer_fasta = fasta_meta.map { meta, fasta -> fasta }.first()
                         ↑ value channel       ↑ value channel    ↑ 冗余！警告！
```

- `fasta_meta` 由 `workflows/circdna.nf` 传入，已经是 value channel（经过 `.first()` 转换）
- `value.map{...}` 仍然是 value channel
- `value.first()` 是冗余操作 → Nextflow 警告

### 为什么其他 `.first()` 不触发警告

| 位置                           | 代码                                      | 输入类型                         | `.first()` 是否冗余        |
| ------------------------------ | ----------------------------------------- | -------------------------------- | -------------------------- |
| `workflows/circdna.nf:50`      | `ch_fasta.map{...}.first()`               | queue（fromPath）                | 否                         |
| `workflows/circdna.nf:81`      | `channel.fromPath(...).map{...}.first()`  | queue（fromPath）                | 否（且 test_local 不执行） |
| `bam_preprocessing/main.nf:39` | `fasta_meta.join(queue).map{...}.first()` | queue（value.join(queue)=queue） | 否                         |

### 修改安全性分析

`ch_eccsplorer_fasta` 的使用方式：

```groovy
// 第 56 行：BAM 模式
ECCSPLORER(BAM_PREPROCESSING.out.bam_sorted.combine(ch_eccsplorer_fasta))

// 第 61 行：FASTQ 模式
ECCSPLORER(reads.combine(ch_eccsplorer_fasta))
```

- `.combine(queue, value)` 会将 value channel 广播到 queue channel 的每个元素
- 移除 `.first()` 后，`ch_eccsplorer_fasta` 仍然是 value channel（`value.map{...}` = value）
- 广播行为不变，每个样本仍然会与同一个 fasta 配对

## ADDED Requirements

### Requirement: 消除 `first` on value channel 警告

`ECCDNA_MODE` subworkflow 中 `ch_eccsplorer_fasta` 的构造不得对 value channel 调用 `.first()`，以消除 Nextflow 警告。

#### Scenario: test_local eccdna 模式运行无警告

- **WHEN** 用户运行 `nextflow run main.nf -profile test_local,docker --mode eccdna`
- **THEN** Nextflow 日志中不出现 "The operator `first` is useless when applied to a value channel" 警告
- **AND** ECCSPLORER process 正常处理所有样本（每个样本都与参考基因组 fasta 正确配对）

#### Scenario: channel 广播语义保持不变

- **WHEN** `ch_eccsplorer_fasta` 与 `BAM_PREPROCESSING.out.bam_sorted`（queue channel）通过 `.combine()` 组合
- **THEN** 每个样本的 BAM 都与同一个参考基因组 fasta 配对
- **AND** 不发生 one-to-one 匹配导致只处理第一个样本的问题

## MODIFIED Requirements

### Requirement: ch_eccsplorer_fasta 构造方式

`ch_eccsplorer_fasta` 的构造方式从：

```groovy
ch_eccsplorer_fasta = fasta_meta.map { meta, fasta -> fasta }.first()
```

改为：

```groovy
ch_eccsplorer_fasta = fasta_meta.map { meta, fasta -> fasta }
```

- `fasta_meta` 已经是 value channel（由 `workflows/circdna.nf` 第 50 行 `.first()` 转换而来）
- `value.map{...}` 仍然是 value channel，无需再调用 `.first()`
- 移除后，channel 类型和行为不变
