# Tasks

## 阶段一：注释数据库预处理与参数声明

- [ ] Task 1: 创建注释数据库预处理脚本 `scripts/prepare_eccsplorer_database.sh`

  - [ ] SubTask 1.1: 脚本解压 Dfam RepeatMasker.lib.gz（若未解压），输出 `Dfam-RepeatMasker.lib`
  - [ ] SubTask 1.2: 脚本根据物种参数（--species，如 fungi/plant/human）从 RepBase31.07.fasta/ 选择对应 `.ref` 文件（fngrep.ref / athrep.ref / humrep.ref）
  - [ ] SubTask 1.3: 脚本合并 Dfam + RepBase 对应物种库为单一 `eccsplorer_db.fa`，输出到指定目录
  - [ ] SubTask 1.4: 脚本支持 `--out` 参数指定输出路径，默认输出到 `circdna.nf/testdatasets/eccsplorer_db/`

- [ ] Task 2: 在 `circdna.nf/nextflow.config` 新增 ECCsplorer 注释数据库参数

  - [ ] SubTask 2.1: 新增 `eccsplorer_database = null` 参数（仅此一个命令行参数，**不新增任何 control 参数**）

- [ ] Task 3: 执行 `scripts/prepare_eccsplorer_database.sh --species fungi`，生成 `circdna.nf/testdatasets/eccsplorer_db/eccsplorer_db.fa`（用于酵母测试数据）

## 阶段二：ECCSPLORER 模块修改

- [ ] Task 4: 修改 `circdna.nf/modules/local/eccsplorer/main.nf` 支持可选 control 输入

  - [ ] SubTask 4.1: input 块新增可选 `path(control_r1)` 和 `path(control_r2)`（使用 `StageStrategy` 或通过 when 条件区分两套 input 签名；若 Nextflow 不支持可选 input，则使用两个 process 变体 ECCSPLORER 与 ECCSPLORER_WITH_CONTROL）
  - [ ] SubTask 4.2: script 块（BAM 模式）条件性附加 control 参数：当 control_r1/control_r2 非空时，将 `${control_r1} ${control_r2}` 作为 ECCsplorer 的位置参数 3、4 传入
  - [ ] SubTask 4.3: script 块（FASTQ 模式）同 SubTask 4.2
  - [ ] SubTask 4.4: stub 块无需修改（stub 不执行真实命令）

- [ ] Task 5: 更新 `circdna.nf/conf/modules.config` 的 ECCSPLORER `ext.args`

  - [ ] SubTask 5.1: 将 `ext.args` 改为动态构造，条件性附加 `-d ${params.eccsplorer_database}`（当参数非 null 时）
  - [ ] SubTask 5.2: 保持现有 `--mode map` 和 `--trim_reads` 逻辑

- [ ] Task 6: 同步修改 `bio.nf/modules/eccsplorer/main.nf`（与 circdna.nf 保持一致）

- [ ] Task 7: 更新 `bio.nf/modules/eccsplorer/meta.yml` 与 `circdna.nf/modules/local/eccsplorer/meta.yml`（新增 control 输入描述）

## 阶段三：samplesheet 校验脚本与 subworkflow/workflow 修改

- [ ] Task 8: 修改 `circdna.nf/bin/check_samplesheet.py` 兼容现有 `data_type` 列格式 + 新增 `group` 列

  - [ ] SubTask 8.1: OPTIONAL_FIELDS 同时识别 `data_type` 和 `datatype` 两种列名
  - [ ] SubTask 8.2: 值归一化为小写后校验，接受 `eccDNA`/`gDNA`/`eccdna`/`gdna` 任一写法
  - [ ] SubTask 8.3: 新增识别可选的 `group` 列（任意非空字符串值），输出 samplesheet 保留该列
  - [ ] SubTask 8.4: 输出 samplesheet 保留原始列名（输入是 `data_type` 则输出也是 `data_type`；无 `group` 列则输出也不含）
  - [ ] SubTask 8.5: meta map 中以 `data_type` 键存储归一化后的小写值（`eccdna`/`gdna`），以 `group` 键存储分组值（无 group 列时为空字符串）

- [ ] Task 9: 修改 `circdna.nf/subworkflows/local/eccdna_mode/main.nf`

  - [ ] SubTask 9.1: take 块新增 `control_reads` 通道（gDNA reads，可为 empty channel）
  - [ ] SubTask 9.2: 根据 meta.data_type 过滤出 eccDNA 样本（dataset A），将 gDNA 样本的 reads 作为 control
  - [ ] SubTask 9.3: 按 `meta.group` 列值匹配 eccDNA 与 gDNA：同一 group 下的 eccDNA 样本使用该 group 下的 gDNA 样本作为 control
    - 一个 group 下有 1 个 gDNA + N 个 eccDNA → 这 N 个 eccDNA 共用该 gDNA control
    - 一个 group 下无 gDNA → 该 group 的 eccDNA 走无 control 路径
    - 一个 group 下有多个 gDNA → 报错
  - [ ] SubTask 9.4: 无 `group` 列但 samplesheet 含 gDNA：若 gDNA 仅 1 个，所有 eccDNA 共用该 control；若 gDNA 多于 1 个，报错要求用户提供 `group` 列
  - [ ] SubTask 9.5: 匹配的 eccDNA + control combine 后传入 ECCSPLORER（带 control）
  - [ ] SubTask 9.6: 无匹配 control 的 eccDNA 样本走原有 ECCSPLORER 路径（无 control）
  - [ ] SubTask 9.7: 两次运行的结果合并到统一的 ch_eccsplorer_bed 通道

- [ ] Task 10: 修改 `circdna.nf/workflows/circdna.nf`
  - [ ] SubTask 10.1: 从 samplesheet 解析 data_type 列（若存在）并传递到 meta map（键名 `data_type`）
  - [ ] SubTask 10.2: 从 samplesheet 解析 group 列（若存在）并传递到 meta map（键名 `group`）
  - [ ] SubTask 10.3: 将 gDNA reads 作为 channel 传入 ECCDNA_MODE subworkflow（无 gDNA 时为 empty channel）
  - [ ] SubTask 10.4: 不新增任何 `eccsplorer_control_fastq_1/2` 命令行参数

## 阶段四：测试数据与 samplesheet

- [ ] Task 11: 准备 gDNA 测试数据

  - [ ] SubTask 11.1: 复制 `circdna_2_R1.fastq.gz` → `gdna_2_R1.fastq.gz`，`circdna_2_R2.fastq.gz` → `gdna_2_R2.fastq.gz`
  - [ ] SubTask 11.2: 复制 `circdna_3_R1.fastq.gz` → `gdna_3_R1.fastq.gz`，`circdna_3_R2.fastq.gz` → `gdna_3_R2.fastq.gz`
  - [ ] SubTask 11.3: 确认 `gdna_1_R1/R2.fastq.gz` 已存在（无需复制）

- [ ] Task 12: 修改 `samplesheets/samplesheet_local.csv`（无 gDNA 对比测试）

  - [ ] SubTask 12.1: 补充 `data_type` 列，3 个样本（circdna_1/2/3）均为 `eccDNA`
  - [ ] SubTask 12.2: 不包含任何 gDNA 行
  - [ ] SubTask 12.3: 不包含 `group` 列（无 gDNA 时无需分组）

- [ ] Task 13: 创建 `samplesheets/samplesheet_local_with_gdna.csv`（有 gDNA control 测试）
  - [ ] SubTask 13.1: 3 个 eccDNA 样本（circdna_1/2/3），data_type 为 `eccDNA`
  - [ ] SubTask 13.2: 3 个 gDNA 样本（gdna_1/2/3），data_type 为 `gDNA`，FASTQ 路径指向 testdatasets/testdata/
  - [ ] SubTask 13.3: 新增 `group` 列建立配对：circdna_1↔gdna_1 同属 group=A，circdna_2↔gdna_2 同属 group=B，circdna_3↔gdna_3 同属 group=C
  - [ ] SubTask 13.4: 不包含 `eccsplorer_control_fastq_1/2` 列（control 由 data_type + group 列驱动）

## 阶段五：测试配置

- [ ] Task 14: 修改 `conf/test_local.config`

  - [ ] SubTask 14.1: `input` 保持 `samplesheets/samplesheet_local.csv`（已补充 data_type 列，无 group 列）
  - [ ] SubTask 14.2: 新增 `eccsplorer_database = '/Users/siyangming/nextflow_nf_core/circdna.nf/testdatasets/eccsplorer_db/eccsplorer_db.fa'`
  - [ ] SubTask 14.3: `outdir` 改为 `results/test_local`（相对 circdna.nf/ 目录）

- [ ] Task 15: 新建 `conf/test_local_gdna.config`

  - [ ] SubTask 15.1: `input` 设为 `samplesheets/samplesheet_local_with_gdna.csv`
  - [ ] SubTask 15.2: `eccsplorer_database` 同 test_local
  - [ ] SubTask 15.3: `outdir` 设为 `results/test_local_gdna`
  - [ ] SubTask 15.4: 其余参数（genome、mode 等）与 test_local 一致

- [ ] Task 16: 在 `nextflow.config` 的 profiles 块注册 `test_local_gdna` profile
  - [ ] SubTask 16.1: 新增 `test_local_gdna { includeConfig 'conf/test_local_gdna.config' }`

## 阶段六：版本管理与文档

- [ ] Task 17: 更新 `circdna.nf/nextflow.config` 版本号 4.2.3 → 4.3.0
- [ ] Task 18: 更新 `circdna.nf/CHANGELOG.md`
  - [ ] SubTask 18.1: 新增 `## v4.3.0 - [2026-08-04]` 版本段
  - [ ] SubTask 18.2: Enhancements & fixes 记录：ECCsplorer BLAST 注释数据库支持、gDNA control 支持（由 samplesheet data_type + group 列驱动）、check_samplesheet.py 兼容 data_type 列名与新增 group 列、统一输出目录规范
  - [ ] SubTask 18.3: Added 部分记录：`params.eccsplorer_database`、`scripts/prepare_eccsplorer_database.sh`、`conf/test_local_gdna.config`、`samplesheets/samplesheet_local_with_gdna.csv`、samplesheet 新增可选 `group` 列

## 阶段七：验证

- [x] Task 19: stub 模式验证（无 control）

  - [x] SubTask 19.1: 执行 `nextflow run main.nf -profile test_local,docker -stub` 验证流程通过
  - [x] SubTask 19.2: 验证 `results/test_local/eccsplorer/` 下产出文件存在
  - [x] SubTask 19.3: 验证 ECCSPLORER 任务无 failed，命令中不含 dataset B 参数

- [x] Task 20: stub 模式验证（有 control）

  - [x] SubTask 20.1: 执行 `nextflow run main.nf -profile test_local_gdna,docker -stub` 验证流程通过
  - [x] SubTask 20.2: 验证 `results/test_local_gdna/eccsplorer/` 下产出文件存在
  - [x] SubTask 20.3: 验证 ECCSPLORER 任务无 failed，命令中包含 gDNA control 位置参数

- [ ] Task 21: 真实模式对比测试（可选，耗时长）
  - [ ] SubTask 21.1: 执行 `nextflow run main.nf -profile test_local,docker`（无 control + 注释库）
  - [ ] SubTask 21.2: 执行 `nextflow run main.nf -profile test_local_gdna,docker`（有 control + 注释库）
  - [ ] SubTask 21.3: 对比 `results/test_local/eccsplorer/*_candidates.bed` 与 `results/test_local_gdna/eccsplorer/*_candidates.bed` 的候选数量差异
  - [ ] SubTask 21.4: 验证 `*_blast.m6` 非空（注释库生效）

# Task Dependencies

- Task 3 依赖 Task 1（脚本创建后才能执行生成数据库）
- Task 4-7（模块修改）可与 Task 1-3（数据库准备）并行
- Task 8（check_samplesheet.py）独立，可与 Task 4-7 并行
- Task 9-10（subworkflow/workflow）依赖 Task 4-7（模块接口确定后才能接线）+ Task 8（data_type + group 键名确定）
- Task 11-13（测试数据与 samplesheet）依赖 Task 8（samplesheet 格式确定）
- Task 14-16（测试配置）依赖 Task 3、Task 13
- Task 17-18（版本与文档）依赖 Task 1-16 全部完成
- Task 19-20（stub 验证）依赖 Task 17-18 完成
- Task 21（真实测试）依赖 Task 19-20 通过

# 并行化说明

- 阶段一（Task 1-3）与阶段二（Task 4-7）可并行
- Task 8（check_samplesheet.py）独立，可与阶段一、二并行
- 阶段四（Task 11-13）与阶段五（Task 14-16）可部分并行（samplesheet 与 config 独立）
