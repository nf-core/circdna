# Checklist

## 阶段一：注释数据库预处理与参数声明

### scripts/prepare_eccsplorer_database.sh

- [ ] 脚本支持 `--species` 参数（fungi/plant/human 等）
- [ ] 脚本解压 Dfam RepeatMasker.lib.gz（若未解压）
- [ ] 脚本根据物种从 RepBase31.07.fasta/ 选择对应 .ref 文件
- [ ] 脚本合并 Dfam + RepBase 物种库为单一 eccsplorer_db.fa
- [ ] 脚本支持 `--out` 参数指定输出路径
- [ ] 脚本默认输出到 `circdna.nf/testdatasets/eccsplorer_db/`

### nextflow.config 参数

- [ ] 新增 `eccsplorer_database = null`
- [ ] **未新增** `eccsplorer_control_fastq_1` 或 `eccsplorer_control_fastq_2` 参数（control 完全由 samplesheet 驱动）

### 数据库生成

- [ ] 执行脚本生成 `testdatasets/eccsplorer_db/eccsplorer_db.fa`
- [ ] 文件非空且为有效 FASTA 格式

## 阶段二：ECCSPLORER 模块修改

### circdna.nf/modules/local/eccsplorer/main.nf

- [ ] input 块新增可选 control_r1/control_r2 输入（或新增 ECCSPLORER_WITH_CONTROL process 变体）
- [ ] script 块（BAM 模式）条件性附加 control 位置参数
- [ ] script 块（FASTQ 模式）条件性附加 control 位置参数
- [ ] stub 块保持不变

### conf/modules.config

- [ ] ECCSPLORER 的 ext.args 动态构造，条件性附加 `-d ${params.eccsplorer_database}`
- [ ] 保持现有 `--mode map` 和 `--trim_reads` 逻辑

### bio.nf 同步

- [ ] bio.nf/modules/eccsplorer/main.nf 与 circdna.nf 一致
- [ ] bio.nf/modules/eccsplorer/meta.yml 与 circdna.nf 一致（新增 control 输入描述）
- [ ] circdna.nf/modules/local/eccsplorer/meta.yml 新增 control 输入描述

## 阶段三：samplesheet 校验脚本与 subworkflow/workflow 修改

### bin/check_samplesheet.py

- [ ] OPTIONAL_FIELDS 同时识别 `data_type` 和 `datatype` 两种列名
- [ ] 值归一化为小写后校验，接受 `eccDNA`/`gDNA`/`eccdna`/`gdna` 任一写法
- [ ] 新增识别可选的 `group` 列（任意非空字符串值）
- [ ] 输出 samplesheet 保留原始列名（输入是 `data_type` 则输出也是 `data_type`；无 `group` 列则输出也不含）
- [ ] meta map 中以 `data_type` 键存储归一化后的小写值，以 `group` 键存储分组值

### subworkflows/local/eccdna_mode/main.nf

- [ ] take 块新增 control_reads 通道（gDNA reads，可为 empty channel）
- [ ] 根据 meta.data_type 过滤出 eccDNA 样本（dataset A）
- [ ] 按 `meta.group` 列值匹配 eccDNA 与 gDNA（同一 group 下配对）
- [ ] 一个 group 下有 1 个 gDNA + N 个 eccDNA → N 个 eccDNA 共用该 gDNA control
- [ ] 一个 group 下无 gDNA → 该 group 的 eccDNA 走无 control 路径
- [ ] 一个 group 下有多个 gDNA → 报错
- [ ] 无 group 列但 gDNA 仅 1 个 → 所有 eccDNA 共用该 control
- [ ] 无 group 列且 gDNA 多于 1 个 → 报错要求用户提供 group 列
- [ ] 匹配的 eccDNA + control combine 后传入 ECCSPLORER（带 control）
- [ ] 无匹配 control 的 eccDNA 走原有 ECCSPLORER 路径（无 control）
- [ ] 两次运行结果合并到统一 ch_eccsplorer_bed 通道

### workflows/circdna.nf

- [ ] 从 samplesheet 解析 data_type 列（若存在）并传递到 meta map（键名 `data_type`）
- [ ] 从 samplesheet 解析 group 列（若存在）并传递到 meta map（键名 `group`）
- [ ] 将 gDNA reads 作为 channel 传入 ECCDNA_MODE（无 gDNA 时为 empty channel）
- [ ] **未新增** `eccsplorer_control_fastq_1/2` 命令行参数

## 阶段四：测试数据与 samplesheet

### gDNA 测试数据

- [ ] gdna_1_R1.fastq.gz、gdna_1_R2.fastq.gz 已存在（确认无需创建）
- [ ] gdna_2_R1.fastq.gz、gdna_2_R2.fastq.gz 已创建（复制自 circdna_2）
- [ ] gdna_3_R1.fastq.gz、gdna_3_R2.fastq.gz 已创建（复制自 circdna_3）

### samplesheets/samplesheet_local.csv（无 gDNA 对比测试）

- [ ] 补充 `data_type` 列
- [ ] 3 个样本（circdna_1/2/3）data_type 均为 `eccDNA`
- [ ] 不包含任何 gDNA 行
- [ ] 不包含 `group` 列（无 gDNA 时无需分组）

### samplesheets/samplesheet_local_with_gdna.csv（有 gDNA control 测试）

- [ ] 3 个 eccDNA 样本（circdna_1/2/3），data_type 为 `eccDNA`
- [ ] 3 个 gDNA 样本（gdna_1/2/3），data_type 为 `gDNA`
- [ ] FASTQ 路径指向 testdatasets/testdata/
- [ ] 新增 `group` 列：circdna_1↔gdna_1 同属 group=A，circdna_2↔gdna_2 同属 group=B，circdna_3↔gdna_3 同属 group=C
- [ ] **不包含** `eccsplorer_control_fastq_1/2` 列（control 由 data_type + group 列驱动）

## 阶段五：测试配置

### conf/test_local.config

- [ ] input 保持 `samplesheets/samplesheet_local.csv`（已补充 data_type 列，无 group 列）
- [ ] 新增 `eccsplorer_database` 默认路径
- [ ] outdir 改为 `results/test_local`
- [ ] **未设置** `eccsplorer_control_fastq_1/2` 参数

### conf/test_local_gdna.config

- [ ] input 设为 `samplesheets/samplesheet_local_with_gdna.csv`
- [ ] eccsplorer_database 同 test_local
- [ ] outdir 设为 `results/test_local_gdna`
- [ ] 其余参数（genome、mode 等）与 test_local 一致
- [ ] **未设置** `eccsplorer_control_fastq_1/2` 参数（control 由 samplesheet data_type + group 列驱动）

### nextflow.config profiles

- [ ] 注册 `test_local_gdna { includeConfig 'conf/test_local_gdna.config' }`

## 阶段六：版本管理与文档

### nextflow.config

- [ ] manifest.version 从 '4.2.3' 改为 '4.3.0'

### CHANGELOG.md

- [ ] 顶部新增 `## v4.3.0 - [2026-08-04]` 版本段
- [ ] Enhancements & fixes 记录：BLAST 注释数据库支持、gDNA control（samplesheet data_type + group 列驱动）、check_samplesheet.py 兼容 data_type 列与新增 group 列、统一输出目录规范
- [ ] Added 部分记录新增参数、脚本、配置、samplesheet、group 列
- [ ] **未记录** `eccsplorer_control_fastq_1/2` 参数（因未新增）

## 阶段七：验证

### stub 模式（无 control）

- [x] `nextflow run main.nf -profile test_local,docker -stub` 执行成功
- [x] `results/test_local/eccsplorer/` 下产出文件存在
- [x] ECCSPLORER 任务无 failed
- [x] 命令中不含 dataset B 参数（无 control）

### stub 模式（有 control）

- [x] `nextflow run main.nf -profile test_local_gdna,docker -stub` 执行成功
- [x] `results/test_local_gdna/eccsplorer/` 下产出文件存在
- [x] ECCSPLORER 任务无 failed
- [x] 命令中包含 gDNA control 位置参数（dataset B）

### 真实模式对比（可选）

- [ ] test_local 真实运行完成，`results/test_local/eccsplorer/*_blast.m6` 非空
- [ ] test_local_gdna 真实运行完成，`results/test_local_gdna/eccsplorer/*_blast.m6` 非空
- [ ] 对比两次运行的 candidates.bed 数量差异（验证 control 降低假阳性效果）
