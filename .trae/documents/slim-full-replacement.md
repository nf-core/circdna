# slim 完全替代完整版 ECCsplorer 与 ecc_finder 方案

## 一、Summary

目标：让 slim 原子化模块链**完全替代**完整版 ECCsplorer（2022.01.1.1）与 ecc_finder（1.0），达到三方面对齐：

1. **功能覆盖**：map / clu / PRExer / all / comparative 全部模式；ecc_finder 的 map-sr / asm-sr / map-ont / asm-ont 全部子命令
2. **输出契约**：slim 链产出一致命名的原版产物（`eccpipe_results/` 目录树，hiconf/lowconf 双 bed、`_blast.m6`、`ECC-SEQUENCES.fasta`、cluster 三表、三份 summary html 等 20+ 项），下游消费者可无缝切换
3. **参数暴露**：完整版 16 个 CLI 参数 + config 阈值全部可通过 slim params / ext.args 配置

一次性全部实现，共 7 个工作块（见第五节）。完成后完整版黑盒模块可下线，但**保留**（本计划不删除完整版模块，避免破坏已有运行记录）。

## 二、当前状态分析（基于代码级探索）

### 2.1 已实现且等价（无需改动）

- map 比对：nf-core segemehl（`--splits --briefcigar --MEOP --accuracy 95` 与 config.py:189-193 逐项一致）
- split-read 检测：haarz + 长度过滤（35000/100/merge 1000）
- peak detection：scipy find_peaks（`ECCsplorer/bin/peak_detect.py`）
- clu 聚类：`bio.nf/modules/repeatexplorer2`（seqclust，独立镜像 repeatexplorer:2.3.8）+ `clu_candidates` 候选表
- samtools/bedtools/seqtk/minimap2/segemehl 等 nf-core 模块已 vendored 于 `circdna.nf/modules/nf-core/`

### 2.2 ecc_finder 缺口（`ecc_finder/` 完整版 vs `circdna.nf/modules/local/ecc_finder_slim/`）

| 缺口                 | 完整版位置                                                                                                                                                                              | 现状                                                                                                                                    |
| -------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------- |
| **map-ont 完全缺失** | map-ont.py：minimap2 `-x map-ont`→PAF 过滤(`-q`/`-a` 默认 200)→TideHunter `-u` 拆 unit→unit 重比对→Genrich `-yv`→bedtools 过滤(`--min-read`/`--min-bound` 0.8/`--min-cov` 10)→csv/fasta | `TIDEHUNTER` 模块存在但未接入任何子工作流                                                                                               |
| **asm-ont 完全缺失** | asm-ont.py：TideHunter `-c/-e/-p/-P 1000000 -5/-3`→cd-hit-est(`-c/-M/-l`)→去 singleton→fasta                                                                                            | cd-hit 无模块（独立快照 `modules/modules/nf-core/cdhit/` 有）                                                                           |
| map-sr disc 证据链   | map-sr.py run_disc（SAMFLAGS `-G 2`/`-f 83/163`）+ run_intersect **split∩disc inner merge**→6 列无表头 CSV（chr,start,end,num_s,num_d,len）                                             | split_detect 产出 disc.bed 但**从未被消费**；仅 split 单证据；CSV 列/表头完全不同（`supporting_reads/confidence_score` 为 slim 自定义） |
| Genrich 参数         | run_Genrich：`-v -l min_peak(200) -g max_dist(1000) -p max_pvalue(0.05)`                                                                                                                | `conf/modules.config` 中 GENRICH 仅 `ext.args = '-v'`，阈值不可调                                                                       |
| asm-sr 参数          | unicycler `--keep 0 --min_fasta_length`                                                                                                                                                 | nf-core UNICYCLER 用 `--no_rotate`；`asm_filter.py --min-length 500` 硬编码                                                             |
| 分布图               | map-sr.py 行 363-371 `{prefix}.distribution.png`                                                                                                                                        | 缺失                                                                                                                                    |

### 2.3 ECCsplorer 缺口（`ECCsplorer/` 完整版 vs `circdna.nf/modules/local/eccsplorer_slim/`）

| 缺口                       | 完整版位置                                                                                                                                                                          | 现状                                                                                                                                                 |
| -------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------- |
| **BLAST 注释链**           | ECCsplorer.py 605-665：RepeatExplorer 默认库+`-d` 库→makeblastdb→blastdb_aliastool→combinedDB.fas；eccMapper.py 176-180 每候选 blastn `-max_hsps 25` outfmt 14 列→`{cand}_blast.m6` | 无任何 blast 模块（快照 `modules/modules/nf-core/blast/` 有 blastn/makeblastdb）；`blastdb_aliastool` 无 nf-core 模块需自建                          |
| **discordant-read 检测**   | eccMapper.py 418-448：`-G 2`/`-f 83/163`→groupby→genomecov `-bga`→阈值 `max_cov×0.05`→`{pre}_regions-DR.bed`                                                                        | candidate_extract/main.nf:36-38 用 `cp peak_all` 顶替 peak_dr，hiconf/lowconf 语义改变                                                               |
| **control 比较（map 侧）** | eccMapper.py 365-370 六种比对 `{pre}.bed`/`-SR.bed`/`-DR.bed`（TR/CO）；eccDNA_Rcodes.py convgenome/convregion/enrichment（`enrich.all = TR.all/CO.all`，`enrich.alt` 对照兜底）    | `ECCSPLORER_SLIM_PIPELINE` 不接收 control 通道（workflows/circdna.nf:336-342 只传 `ch_eccdna_reads`）                                                |
| **normalize 失真**         | mapped bases 来自 `samtools stats`（eccMapper.py 342-360）                                                                                                                          | `eccsplorer_slim_pipeline/main.nf:131-134` 硬编码 `channel.value(1000000)`；region 模式未启用（normalize/main.nf:27-31 默认 genome）→ 富集分数不可信 |
| **comparative 模式**       | eccComparer.py：簇 contig vs 候选 blastn→3 类图→eccComp_summary.html                                                                                                                | 完全缺失                                                                                                                                             |
| **all 编排**               | ECCsplorer.py 779-797：map+clu+comparative 串联                                                                                                                                     | 缺失（circle_identifier 只有 map_slim/clu_slim 两条独立路径）                                                                                        |
| **PRExer 预处理**          | config.py 176-180 trimmomatic（ILLUMINACLIP:2:30:10 + LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:35）；eccPrepare.py seqtk FASTQ→FASTA                                          | 无裁剪/转换模块；`--read-count` 不支持 `auto`（0.1× 覆盖）；seed 42 vs 原版 12                                                                       |
| **可视化/报告**            | rline_multiplot collage（eccDNA_Rcodes.py 257-383）、clu scatter、eccMap/eccCL/eccComp 三份 summary html                                                                            | visualize/main.nf:29 硬编码 `--mode manhattan`，line collage 未接线；html_report 恒 0/N/A                                                            |
| **输出契约**               | config.py 240-249 目录树 + 20+ 固定命名产物                                                                                                                                         | 各模块独立 publishDir，命名/结构全不同                                                                                                               |
| **参数暴露**               | ECCsplorer.py 100-232 共 16 参数                                                                                                                                                    | 仅 `-ref/-tax/-cpu` 有对应；`-win/-cnt(auto)/-dsa/-dsb/-img/-trm/-d/-rgs` 丢失或硬编码                                                               |

### 2.4 可复用基础

- nf-core 模块已 vendored：samtools(全子工具)、bedtools(全)、seqtk、segemehl、minimap2、bwa、unicycler、trimgalore、picard、mosdepth
- 独立快照可 vendored：`modules/modules/nf-core/{blast,trimmomatic,fastp,cdhit}/`（blast 2.17.0 wave 镜像等）
- `testdatasets/eccsplorer_db/`（Dfam-RepeatMasker.lib + eccsplorer_db.fa）可作 BLAST 库
- `bio.nf/modules/repeatexplorer2` 已接入 clu_slim

## 三、总体设计决策

1. **输出兼容层**：不改动各 slim 模块内部逻辑的前提下，新增 **契约导出（contract export）** 子工作流 `subworkflows/local/slim_contract_export/main.nf`，消费各 slim 模块输出，按原版命名（见附表）重写到 `results/eccpipe_results/{mapping_results,clustering_results,comparative_results,reference_data}` 树。这样 slim 核心链保持原子化，契约层集中维护。
2. **参数暴露**：在 `nextflow.config` 新增 params（见附表 2），经 `conf/modules.config` 的 `ext.args` 透传到模块。
3. **control 通道**：`ECCSPLORER_SLIM_PIPELINE` take 增加 control 通道（对齐 clu_slim 的 pair join 模式，workflows/circdna.nf:347-364 可复用）。
4. **DR 检测**：新建独立模块 `eccsplorer_slim/dr_detect`（不做 genomecov 峰，直接产出 `{pre}_aligned-DR.bed`），candidate_extract 恢复 SR∩all∩DR 三证据逻辑。
5. **验证方式**：同一数据集上 slim 与完整版并行跑，产物逐一 diff；归一化/富集数值对齐以 R 脚本为准。

## 四、缺失数据清单（需要用户提供/确认）

| #   | 数据                                                                           | 用途                                                       | 状态/路径线索                                                     |
| --- | ------------------------------------------------------------------------------ | ---------------------------------------------------------- | ----------------------------------------------------------------- |
| 1   | **ONT 长读测试数据**（fasta/fastq，含串联重复的 read）                         | 验证 map-ont/asm-ont 链（TideHunter 拆 unit、cd-hit 聚类） | 用户已确认有；可从真实数据提取子集（提取规格见 4.1）              |
| 2   | **BLAST 注释数据库**                                                           | 验证 combinedDB 构建 + 每候选 blastn                       | 可复用 `testdatasets/eccsplorer_db/eccsplorer_db.fa` 或用户自备库 |
| 3   | **完整样本集 samplesheet**（含 control 配对：data_type=eccDNA/gDNA + pair 列） | 验证 control 比较、comparative、all 模式                   | 用户已确认有；需要 CSV 路径                                       |
| 4   | **参考基因组**（真核，含 chr 命名）                                            | map 链比对与窗口                                           | 酵母 ref 已有；完整样本集如需新基因组请提供                       |
| 5   | **一致性基准（可选但推荐）**：完整版在同一数据上的输出（`eccpipe_results/`）   | 作为 golden 输出逐文件 diff，验证"完全替代"                | 可现场跑完整版生成                                                |
| 6   | **seed 复现确认**：若需与原版子采样完全一致，PIPELINE_SEED=12 需确认           | clu_prepare 复现性                                         | 原版 config.py:227 = 12                                           |

### 4.1 ONT 测试数据提取规格（本地可运行）

建议从真实 ONT 数据中**随机抽样子集**（保留整条 read，勿截短；提取后按条数控制规模）。基于本地环境（macOS arm64 + Docker 仿真 amd64，性能约降 2-3 倍）与各工具特性（minimap2 快、TideHunter/Genrich 中等、cd-hit 快），给出三档：

| 档位                 | 提取量           | 数据规模                      | 本地预计耗时 | 用途                                                      |
| -------------------- | ---------------- | ----------------------------- | ------------ | --------------------------------------------------------- |
| 冒烟（推荐起步）     | 1,000–2,000 条   | ~10–40 MB（ONT 平均 5–15 kb） | 5–15 分钟    | 验证流水线跑通、模块接线、输出非空                        |
| 常规验证（推荐最终） | 5,000–10,000 条  | ~50–150 MB                    | 30–90 分钟   | 验证能检出串联重复候选（Genrich/TideHunter 有足够 reads） |
| 一致性对比（可选）   | 20,000–50,000 条 | ~200–750 MB                   | 2–6 小时     | slim vs 完整版产物 diff（可过夜）                         |

注意事项：

- 完整版 map-ont 对 read 有最小长度过滤（`-q`/`-a` 默认 200 bp，见 map-ont.py:217-218），真实 ONT read 全长均满足，无需预处理
- 样本量过小（<500 条）可能因串联重复 reads 太少而检不出候选，建议至少 1,000 条起步
- 提取方式建议：`seqtk sample <reads.fastq.gz> <N> > test_ont_N.fastq`（随机抽样保持代表性），或按数据提供方样本拆分逻辑取同一 lane/flowcell 的子集
- 若真实数据中 eccDNA/串联重复富集度低，可结合现有短读测试数据（circdna_1）同源样本抽取，保证对照一致性

### 4.2 testdata 覆盖矩阵（除 ONT 外均可覆盖）

现有 `circdna.nf/testdatasets/` 盘点：

- **短读样本**：`testdata/circdna_{1,2,3}_R1/R2.fastq.gz`（treatment，3 样本）+ `gdna_{1,2,3}_R1/R2.fastq.gz`（control，3 样本）
- **参考基因组**：`reference/genome.fa`（酵母，17 条染色体 I–XVI，约 12 Mb，本地比对/窗口极快）
- **模拟真值**：`testdata/circdna_{1,2,3}_simulated.bed`（酵母坐标候选区间，如 `XI 363849 364088`）
- **BLAST/注释库**：`eccsplorer_db/eccsplorer_db.fa`（85 MB）+ `eccsplorer_db/Dfam-RepeatMasker.lib`（53 MB）+ `annotation/yeast_repeat_anno.fa`（34,865 条 RepeatMasker 注释）
- **真实 samplesheet**：`samplesheet/samplesheet_real_integrated.csv`（6 样本，列 `sample,fastq_1,fastq_2,datatype`）

| 工作块 | 功能                                                 | testdata 能否覆盖 | 数据来源                                                                     |
| ------ | ---------------------------------------------------- | ----------------- | ---------------------------------------------------------------------------- |
| 块 1   | map-sr 证据链修复（disc/Genrich 参数/双证据/分布图） | ✅                | circdna_1/2/3 + gdna_1/2/3（control 配对）；`simulated.bed` 可验证检出准确率 |
| 块 2   | ONT 链（map-ont/asm-ont）                            | ❌ **唯一缺口**   | 需用户按 4.1 规格从真实 ONT 提取（1,000–10,000 条）                          |
| 块 3   | BLAST 注释链                                         | ✅                | `eccsplorer_db/eccsplorer_db.fa` 或 `yeast_repeat_anno.fa` 作 `-d` 库        |
| 块 4   | DR 检测 + control 比较 + normalize 修复              | ✅                | circdna/gdna 配对（六种比对 TR/CO×all/SR/DR 均可生成）                       |
| 块 5   | comparative + all 编排                               | ✅                | map（短读）+ clu（circdna_1/gdna_1 配对，5000 对子采样已验证可行）           |
| 块 6   | PRExer + 输出契约 + 参数暴露                         | ✅                | 短读（trimmomatic/seqtk/read_count auto）                                    |
| 块 7   | 可视化与报告补全                                     | ✅                | 短读 map 产物 + clu 产物                                                     |

**唯一数据缺口：块 2（ONT）**。其余 6 块均可用现有 testdata 覆盖，其中：

- `simulated.bed` 真值可做候选检出率的定量验证（slim vs 真值 overlap 率）
- `samplesheet_real_integrated.csv` 已含 3 对 treatment/control

**注意（需用户确认）**：`samplesheet_real_integrated.csv` 的列是 `sample,fastq_1,fastq_2,datatype`，**无 `pair` 列**；而 clu_slim 的 treatment/control 配对依赖 samplesheet 的 `pair` 字段（workflows/circdna.nf:347-364 按 pair join）。用于验证 control/comparative/all 的 samplesheet 需要补 `pair` 列（如 gdna_1/circdna_1 同 pair=p1），或确认现有 check_samplesheet.py 的配对语义（按样本编号同名匹配）。

## 五、Proposed Changes（7 个工作块）

### 块 1：ecc_finder map-sr 证据链修复

- `circdna.nf/modules/local/ecc_finder_slim/split_detect/main.nf`：确认 disc.bed 正常 emit（已产出）
- `circdna.nf/modules/local/ecc_finder_slim/merge_score/main.nf` + `ecc_finder/bin/merge_score.py`：恢复完整版语义——split∩disc **inner merge**、过滤 `read_x - min_read >= 0`、输出**无表头 6 列** `chr,start,end,num_s,num_d,len`；`--min-read` 默认 5
- `circdna.nf/conf/modules.config` GENRICH：透传 `-l`(200)/`-g`(1000)/`-p`(0.05)（新 params `eccfinder_peak_len/eccfinder_max_dist/eccfinder_max_pvalue`）
- 新建 `ecc_finder_slim/distribution/main.nf`：`{prefix}.distribution.png`（候选长度分布，matplotlib）
- 子工作流 `ecc_finder_slim_pipeline/main.nf`：消费 disc.bed → merge_score

### 块 2：ecc_finder ONT 链（map-ont / asm-ont）

- vendored：`cdhit/cdhitest`（自 `modules/modules/nf-core/cdhit/`）
- `circdna.nf/modules/local/tidehunter/main.nf`：**接入**（当前已被块 1 修复为 reads→cons.fa，map-ont 需 `-u` 拆 unit 变体——新建 `modules/local/tidehunter/unit/main.nf` 或 ext.args 分支）
- 新建子工作流 `circdna.nf/subworkflows/local/ecc_finder_ont_pipeline/main.nf`：
  - map-ont：minimap2(`-x map-ont`)→PAF 过滤(≥200/200)→TideHunter `-u` 拆 unit→unit bam→Genrich `-yv`→bedtools intersect(`-f 0.8`)+coverage(`--min-read`≥3/`--min-cov`≥10)→`{prefix}.csv/fasta`
  - asm-ont：TideHunter(`-c 2 -e 0.25 -p 30 -P 1000000 -5/-3`)→`-l` 输出 consensus→cd-hit-est(`-c 0.8 -M 0 -l 1`)→去 singleton→`{prefix}.fasta`
- `workflows/circdna.nf`：circle_identifier 新增 `ecc_finder_map_ont_slim` / `ecc_finder_asm_ont_slim`

### 块 3：ECCsplorer BLAST 注释链

- vendored：`blast/makeblastdb`、`blast/blastn`（自独立快照，wave blast 2.17.0）
- 新建 `circdna.nf/modules/local/eccsplorer_slim/blast_combineddb/main.nf`：`makeblastdb -in {db} -dbtype nucl`（用户库 + RepeatExplorer 默认库 `dna_database_masked.fasta`）→ `blastdb_aliastool -dblist ... -out combinedDB.fas`
- 新建 `circdna.nf/modules/local/eccsplorer_slim/blast_annotate/main.nf`：每候选 `blastn -query {cand}.fasta -db combinedDB.fas -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" -max_hsps 25` → `{cand}_blast.m6`
- 接入 `eccsplorer_slim_pipeline`（blast_results 通道不再传 `[]`）；`params.eccsplorer_database`（nextflow.config:48 已存在）挂载为 `-d`

### 块 4：ECCsplorer DR 检测 + control 比较 + normalize 修复

- 新建 `circdna.nf/modules/local/eccsplorer_slim/dr_detect/main.nf`：`samtools view -G 2 / -f 83 / -f 163` → groupby count_distinct → `{pre}_aligned-DR.bed`（对照 eccMapper.py 297-335）
- `candidate_extract`：恢复 SR∩all∩DR 三证据 → hiconf(3/3) / lowconf(2/3) 双文件输出（对照 eccMapper.py 462-501）
- `ECCSPLORER_SLIM_PIPELINE` take 增加 control；workflows/circdna.nf 按 pair join 传 control（复用 clu_slim 模式）
- `coverage_profile`：六种比对输入（TR/CO × all/SR/DR）
- `normalize`：`mapped_bases` 改真实值（samtools stats 解析，`_alignment-stats.txt`），region 模式启用，实现 `enrich.all` + `enrich.alt`（对照 eccDNA_Rcodes.py 83-191）
- 新建 `samtools stats` 接入（nf-core 已 vendored）

### 块 5：comparative + all 编排

- 新建 `circdna.nf/subworkflows/local/eccsplorer_comparative/main.nf`：
  - `eccMapper_cand.fa` 构建 + 每簇 contig blastn → `blast_CL{cl}.m6`
  - m6 解析 → map_cand_dict → 3 类图（comp_plot / comp_manhattan / comp_scatter）
  - `eccComp_summary.html`（对照 eccComparer.py 全文）
- 新建 `circdna.nf/subworkflows/local/eccsplorer_all_slim/main.nf`：编排 map_slim + clu_slim + comparative；circle_identifier 新增 `eccsplorer_all_slim`

### 块 6：PRExer + 输出契约 + 参数暴露

- vendored：`trimmomatic`（自独立快照）或复用 fastp；新建 `eccsplorer_slim/trim_reads` 模块（ILLUMINACLIP:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:35，`-trm` 参数 nex/tru2/tru3/tru3-2/custom）
- `seqtk seq -A`（nf-core 已 vendored）FASTQ→FASTA
- `clu_prepare`：`--read-count auto` 支持（`floor(genome_size×0.1/best_len)`）、seed 可配（默认 12）
- 新建 `subworkflows/local/slim_contract_export/main.nf`：消费 slim 全部产物 → 重写为原版命名树（附表 1），publish 到 `results/eccpipe_results/`
- 参数暴露：附表 2 全部 params

### 块 7：可视化与报告补全

- `visualize`：`--mode line` 接线（candidate collage，R 脚本已实现）、clu scatter
- `html_report`：mapped_reads/bases 真实值、blast best hit 链接、collage 图引用 → `eccMap_summary.html`；`eccCL_summary.html`（clu 侧）；`eccComp_summary.html`（块 5）

## 附表 1：输出契约文件清单（slim → 原版命名映射）

| 原版产物                                                                                                                              | slim 来源                                     |
| ------------------------------------------------------------------------------------------------------------------------------------- | --------------------------------------------- |
| `reference_data/{ref}.idx/_chrSize.txt/_win{w}.bed`                                                                                   | segemehl index / faidx / makewindows          |
| `mapping_results/{pre}_alignment-stats.txt`                                                                                           | samtools stats                                |
| `mapping_results/{pre}.bed / _aligned-SR.bed / _aligned-DR.bed`                                                                       | bamtobed / haarz / dr_detect                  |
| `mapping_results/{pre}_regions-SR.bed / _regions-DR.bed`                                                                              | haarz+merge / dr_detect+genomecov             |
| `mapping_results/{pre}_hiconf-ECC-REGIONS.bed / {pre}_lowconf-ECC-regions.bed`                                                        | candidate_extract                             |
| `mapping_results/{pre[-CO]}_summary_coverages.csv(_normalized.csv)`                                                                   | normalize(genome)                             |
| `mapping_results/{pre[-CO]}_summary_region-coverages.csv(_normalized.csv)`                                                            | normalize(region)                             |
| `mapping_results/{pre[-CO]}_ECC-SEQUENCES.fasta`                                                                                      | getfasta 汇总                                 |
| `mapping_results/candidates/{cand}/{cand}.fasta / _coverage-summary_RAW.csv / _NORM.csv / _collage-plot.png / _blast.m6`              | coverage_profile / visualize / blast_annotate |
| `mapping_results/{pre}.trns.txt`                                                                                                      | segemehl junction reads                       |
| `clustering_results/COMPARATIVE_CLUSTER_TABLE_eccCANDIDATES.csv / comparative_cluster_table.csv / comp_cl_tab_eccCandidates_list.csv` | clu_candidates（补列表文件）                  |
| `clustering_results/summary_scatter.png`                                                                                              | visualize(clu)                                |
| `comparative_results/*`（eccMapper_cand.fa、blast_CL{cl}.m6、三图、candidates/）                                                      | comparative 子工作流                          |
| `eccpipe_results/…/eccMap_summary.html / eccCL_summary.html / eccComp_summary.html`                                                   | html_report / clu / comparative               |

## 附表 2：新增 params 清单

```
eccsplorer_trim_reads  (nex|tru2|tru3|tru3-2|custom, null)   → 块6
eccsplorer_window_size (100)                                  → 块6 透传 -w
eccsplorer_read_count  (null|auto|int)                        → 块6 clu_prepare
eccsplorer_read_seed   (12)                                   → 块6
eccsplorer_image_format(png)                                  → 块7
eccsplorer_pre_a/pre_b (TR/CO)                                → 块4/6
eccsplorer_genome_size (null)                                 → 块6
eccsplorer_database    (已存在 nextflow.config:48)            → 块3
eccsplorer_max_cand_cnt(500)                                  → 块4
eccsplorer_ecc_proportion(0.8)                                → clu_candidates
eccfinder_peak_len(200)/eccfinder_max_dist(1000)/eccfinder_max_pvalue(0.05)  → 块1
eccfinder_ont_min_read(3)/eccfinder_ont_min_bound(0.8)/eccfinder_ont_min_cov(10) → 块2
eccfinder_asm_min_length(500)                                 → 块1
```

## 六、假设与决策

- 假设 1：完整版模块**保留**不删除（仅新增 slim 补全），避免破坏既有运行记录
- 决策 1：输出契约用**独立契约导出层**实现（不改造各 slim 模块内部命名），降低回归风险
- 决策 2：BLAST 库默认用 `testdatasets/eccsplorer_db/`，可通过 `-d`/params 覆盖
- 决策 3：比对统一走 segemehl（与完整版一致）；`BAM_PREPROCESSING` 的 BWA 产物仅用于完整版 ecc_finder 路径，不改
- 决策 4：归一化/富集数值对齐以 `ECCsplorer/bin/normalize.R`（已提取）为准，确保公式与 eccDNA_Rcodes.py 一致
- 决策 5：clu 判定 `>=` vs 完整版 `>` 保持 slim 现状（差异可接受），但 seed 改为 12 复现原版

## 七、Verification 步骤

1. **每块独立测试**：块 1/2 用 `tests/ecc_finder_ont/test.nf`、`tests/ecc_finder_sr/test.nf`；块 3/4/5/6/7 用 `tests/eccsplorer_full/test.nf`（slim 全链 + 契约导出）
2. **端到端对比（关键）**：同一数据上 slim（`eccsplorer_all_slim` / `ecc_finder_map_sr_slim,ecc_finder_map_ont_slim,...`）与完整版（`eccsplorer` / `ecc_finder`）并行跑，逐文件 diff：
   - hiconf/lowconf BED 区间数 ≥ 90% 一致
   - `ECC-SEQUENCES.fasta`、`_blast.m6`、cluster 三表、三份 summary html 字段非空且命中一致
   - candidates CSV 6 列（num_s/num_d/len）数值对齐
3. **输出契约清单逐项核对**（附表 1 每行存在且非空）
4. **参数暴露验证**：附表 2 每个 params 变更后结果可复现
5. **主流水线回归**：`-resume` 3 个 slim circle_identifier 仍成功；新增 `eccsplorer_all_slim` 全链运行成功
6. **移除验证**：circle_identifier 只含 slim 标识符时，无任何完整版镜像（quay.io/bioinfortools/eccsplorer:2022.01.1.1、ecc_finder:1.0.0）进程被拉起（grep 运行日志确认）

## 八、风险

- 归一化/富集数值与原版**逐位一致**较难：原版经 pyRserve 执行 R 脚本，slim 为独立 R 脚本；需按公式逐行核对（块 4 重点）
- segemehl 版本（0.3.4 vendored）与完整版镜像内版本可能不同 → 比对结果微小差异
- ONT 链验证依赖用户提供真实 ONT 数据（当前 testdata 无 ONT）
- comparative/all 模式原版仅在特定条件触发（four_read_files + mode all/clu），测试需构造对应场景
