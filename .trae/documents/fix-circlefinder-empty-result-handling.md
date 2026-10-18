# 修复 CIRCLEFINDER 空结果处理 — 计划文档

## 任务摘要

样本 CRR1082977 (Solanum_lycopersicum) 在 Circle-Finder 步骤输出 `circle_finder_exit_log.txt`，内容为误导性的 "ERROR - CIRCLE_FINDER" 信息。经查证：

1. **这不是流程失败** — 任务 `exit: 0` 成功完成（`.nextflow.log` 第 1310 行确认）
2. **这是 Circle-Finder 算法的正常结果** — 该样本没有 reads 满足 "频率==3" 模式（1 concordant + 2 split = 3）
3. **"ERROR" 措辞是 circdna.nf 自定义添加的** — 原版 Circle-Finder 工具根本不做此检查

## 根因分析（基于网络调研与源码对比）

### 三个版本对比

| 版本                                      | 空结果处理方式                                     | "无 eccDNA" 时输出                                |
| ----------------------------------------- | -------------------------------------------------- | ------------------------------------------------- |
| **原版** (`pk7zuva/Circle_finder`)        | 无任何检查，awk/grep 静默产生空文件                | 空 `microDNA-JT.txt`（无报错）                    |
| **改进版** (`suda-huanglab/circlehunter`) | 早期检查 `split_count > 0`，否则 `touch` 空文件    | 空 `microDNA-JT.txt`（无报错）                    |
| **circdna.nf 当前版本**                   | 多处 `file_exists` 检查 → 写 "ERROR" 日志 → `exit` | **仅** `circle_finder_exit_log.txt`（误导性报错） |

### 关键发现

**改进版 circlehunter 的最佳实践**（[源码](https://github.com/suda-huanglab/circlehunter/blob/master/tools/circlefinder-pipeline-bwa-mem-samblaster.sh)）：

```bash
split_count=$(awk 'END{print NR}' $6-$7\.split.txt)
if (($split_count > 0)); then
    # ... 运行完整 Circle-Finder 流程 ...
    awk '...' > $6-$7\.microDNA-JT.txt
else
    touch $6-$7\.microDNA-JT.txt    # ← 创建空文件，无报错
fi
```

**circdna.nf 当前做法**（[circlefinder/main.nf:22-32](file:///Users/siyangming/nextflow_nf_core/circdna.nf/modules/local/circlefinder/main.nf)）：

```bash
file_exists () {
    [[ ! -s \$1 ]] && \
    echo "ERROR - CIRCLE_FINDER - $prefix ... No circular DNA was identified ..." \
        > ${prefix}.circle_finder_exit_log.txt && exit
}
```

**问题**：

1. 把"未检测到 eccDNA"误标为 "ERROR"，实际是有效科学结果
2. `exit` 后不产生 `microDNA-JT.txt`，下游通道收不到该样本的输出文件
3. 共有 7 处 `file_exists` 检查（行 45、53、61、66、72、79、99），任一触发即提前退出
4. 原版工具和改进版都不做此类"报错"——它们只是产生空文件

## 修改方案

**核心思路**：采用 circlehunter 改进版的处理模式 — 把"无 eccDNA"视为有效结果，产生空 `microDNA-JT.txt` 文件，不写误导性 "ERROR" 日志。

### 修改文件 1：[circdna.nf/modules/local/circlefinder/main.nf](file:///Users/siyangming/nextflow_nf_core/circdna.nf/modules/local/circlefinder/main.nf)

**改动要点**：

1. 移除 `file_exists` 函数定义（行 22-32）
2. 在脚本开头检查 `split` 文件是否有内容
3. 若 `split` 为空 → `touch ${prefix}.microDNA-JT.txt` 并正常退出（exit 0）
4. 若 `split` 有内容 → 运行原有 awk 流程
5. 输出声明中的 `circle_finder_exit_log.txt`（行 10）保留但极少触发（仅用于真正的脚本错误，而非"无 eccDNA"）

**改动后的脚本结构**：

```groovy
script:
def prefix = task.ext.prefix ?: "${meta.id}"
"""
#!/usr/bin/env bash
set -uo pipefail

# 检查 split 文件是否有数据（参考 circlehunter 最佳实践）
split_count=\$(awk 'END{print NR}' ${split})
if (( \${split_count:-0} == 0 )); then
    # 无 split reads → 无 eccDNA，产生空输出文件
    touch ${prefix}.microDNA-JT.txt
    echo "INFO - CIRCLE_FINDER - $prefix - No split reads found; no eccDNA detected." > ${prefix}.circle_finder_exit_log.txt
    exit 0
fi

awk '{print \$4}' ${split} | sort -T ./ | uniq -c > ${prefix}.split.id-freq.txt
awk '\$1=="2" {print \$2}' ${prefix}.split.id-freq.txt > ${prefix}.split.id-freq2.txt

awk '{print $4}' ${concordant} | sort -T ./ | uniq -c > ${prefix}.concordant.id-freq.txt
awk '\$1=="3" {print \$2}' ${prefix}.concordant.id-freq.txt > ${prefix}.concordant.id-freq3.txt

# 检查 concordant.id-freq3.txt 是否为空（原 file_exists 检查点）
if [[ ! -s ${prefix}.concordant.id-freq3.txt ]]; then
    touch ${prefix}.microDNA-JT.txt
    echo "INFO - CIRCLE_FINDER - $prefix - No reads with frequency==3 in concordant file; no eccDNA detected." > ${prefix}.circle_finder_exit_log.txt
    exit 0
fi

# ... 后续步骤 7-11 保持不变（删除中间的 file_exists 检查） ...
# Step 7-11 原有 awk/sed/sort 逻辑保留
"""
```

**具体行级改动**：

| 行号  | 当前内容                                                          | 改为                                                 |
| ----- | ----------------------------------------------------------------- | ---------------------------------------------------- |
| 19-32 | `file_exists` 函数定义                                            | 替换为 split_count 早期检查（见上）                  |
| 45    | `file_exists ${prefix}.concordant.id-freq3.txt`                   | 改为 `if [[ ! -s ... ]]; then touch ...; exit 0; fi` |
| 53    | `file_exists ${prefix}.split_freq2.txt`                           | 删除（后续 awk 容忍空输入）                          |
| 61    | `file_exists ${prefix}.concordant_freq3.txt`                      | 删除                                                 |
| 66    | `file_exists ${prefix}.split_freq2.oneline.txt`                   | 删除                                                 |
| 72    | `file_exists ${prefix}.split_freq2.oneline.S-R-S-CHR-S-ST.ID.txt` | 删除                                                 |
| 79    | `file_exists ${prefix}.concordant_freq3.2SPLIT-1M.txt`            | 删除                                                 |
| 99    | `file_exists ${prefix}.concordant_freq3.2SPLIT-1M.inoneline.txt`  | 删除                                                 |
| 10    | `emit: log` 保留                                                  | 保留（INFO 日志仍有价值）                            |

### 修改文件 2：[circdna.nf/CHANGELOG.md](file:///Users/siyangming/nextflow_nf_core/circdna.nf/CHANGELOG.md)

在当前版本（3.2.0）下方新增 3.2.1 段落：

```markdown
## 3.2.1 - [2026-08-02]

### Enhancements & fixes

- **CIRCLEFINDER 空结果处理优化**: 参照 suda-huanglab/circlehunter 最佳实践，
  将"未检测到 eccDNA"从误导性 "ERROR" 改为 "INFO" 提示，并产生空 `microDNA-JT.txt`
  文件以保证下游通道完整性。原版 Circle-Finder 工具本身就静默产生空文件，
  circdna.nf 之前自定义的 file_exists 检查过度报错。
```

### 修改文件 3：[circdna.nf/nextflow.config](file:///Users/siyangming/nextflow_nf_core/circdna.nf/nextflow.config)

```groovy
manifest {
    version = '3.2.1'  // 从 3.2.0 升级
}
```

## 假设与决策

| 项                                         | 决策                      | 理由                                                            |
| ------------------------------------------ | ------------------------- | --------------------------------------------------------------- |
| 是否保留 `circle_finder_exit_log.txt` 输出 | **保留**                  | INFO 级别日志对用户有诊断价值（知道为何无 eccDNA）              |
| 是否取消所有 `file_exists` 检查            | **保留首个、删除中间**    | 首个检查点（concordant.id-freq3.txt）有诊断价值；中间检查点冗余 |
| 是否同时启用 `freqGr3`（频率>3）分支       | **不启用**                | 原版默认不使用 freqGr3，且会改变算法行为，超出本次修复范围      |
| 是否新增 `--circle_finder_min_freq` 参数   | **不新增**                | 保持算法一致性，避免引入假阳性                                  |
| 版本升级类型                               | **PATCH (3.2.0 → 3.2.1)** | 行为优化+错误信息修复，无新功能、无破坏性变更                   |

## 验证步骤

1. **本地测试**：

   ```bash
   cd /Users/siyangming/nextflow_nf_core/circdna.nf
   nextflow run main.nf -profile test_local --circle_identifier circle_finder
   ```

   验证：CIRCLEFINDER 任务正常完成，产出 `microDNA-JT.txt`（可能为空）

2. **检查下游通道完整性**：

   - 确认 CIRCLEFINDER 的 `circdna` emit 总是产出文件（即使是空文件）
   - 下游 `CANDIDATE_MERGE` 或结果汇总不再因缺失文件而报错

3. **验证 schema**：

   ```bash
   conda activate nextflow
   nf-core schema build
   ```

4. **Lint 检查**：

   ```bash
   nf-core lint
   ```

5. **版本一致性检查**：
   - `nextflow.config` 中 `version = '3.2.1'`
   - `CHANGELOG.md` 中 `## 3.2.1 - [2026-08-02]`

## 不做的事项

- ❌ 不修改 Circle-Finder 算法逻辑（频率阈值仍为 `==3`）
- ❌ 不启用被注释的 `freqGr3` 分支
- ❌ 不新增可配置参数
- ❌ 不修改 SAMBLASTER / BEDTOOLS 上游模块
- ❌ 不重新跑服务器上的 CRR1082977 样本（修改后可在下次运行时验证）
