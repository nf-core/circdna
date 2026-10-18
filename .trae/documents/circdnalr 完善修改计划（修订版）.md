# circdnalr 完整修改计划（修订版）

范围：`SiYangming/circdna.nf` 分支 `circdnalr`（v4.6.1）。
对照：代表集 [`metadata.csv`](https://github.com/SiYangming/circdna.nf/blob/circdnalr/samplesheets/metadata.csv) 全部 32 条 `circdna` 样本 + 2 条 `circrna`（拒绝）。
原则：按样本分流；统一 BED 契约；不复活 `integrated`；不平行再造 RCA 管道。

本版相对上一稿补上：长读 WGS 背景、RCA 线性化 / T7、PacBio RS vs HiFi、跨平台 `pair`、双 eccDNA pair、metadata → 流水线字段映射。

---

## 0. 明确不做

- 不把 `params.mode` 加回 `integrated`
- 不把 `calculate_ecc_score.py` / `INTEGRATED_EVAL` / `integrated_mode` 搬回 Nextflow（评分仍在 [eccdna.smk](https://github.com/SiYangming/eccdna.smk)）
- 不新建 `pacbio_ecc.nf` 把 CIDER-seq2 和 TideHunter→CircleSeeker 焊成一条
- 不把 TideHunter consensus 喂给 CircleSeeker
- 不把 CIDER-seq2 monomer FASTA 与植物 eccDNA BED 混流
- 不在 read 层 `samtools fastq -f 4` 丢细胞器 reads
- 不新增第二套 TideHunter 参数（沿用 `eccfinder_ont_copy/divergence/period`）
- 不写厨房水槽 process（minimap2|samtools|bedtools 塞一个 script）
- 不把 `library_strategy` 或 `amplification_method` 直接当路由键
- 不把 Illumina+RCA 送进 TideHunter
- 不让全局 `--protocol` / `--mode` 决定多样本混跑

---

## 1. 元数据模型

### 1.1 流水线路由字段（samplesheet 必有或可推断）

| 字段         | 取值                                                                     | 职责                                                     |
| ------------ | ------------------------------------------------------------------------ | -------------------------------------------------------- |
| `platform`   | `illumina` / `pacbio` / `ont`                                            | 预处理链                                                 |
| `assay`      | `wgs` / `circleseq` / `rca` / `ciderseq` / `enriched`                    | 实验类型 → 引擎集合                                      |
| `datatype`   | `gdna` / `eccdna`                                                        | 背景 vs 检测                                             |
| `pair`       | 字符串，可空                                                             | 同研究分组键，**不**等于 ECCsplorer 的 treatment/control |
| `layout`     | 由 fastq_2 是否存在推断 `single_end`                                     | SE/PE                                                    |
| `entrypoint` | `cleaned_fastq` / `raw_fastq` / `subreads` / `hifi_bam`                  | 仅长读入口                                               |
| `concatemer` | `true` / `false`                                                         | 长读 RCA 是否跑 TideHunter                               |
| `read_type`  | `hifi` / `clr` / `ont` / `pe` / `se`                                     | 比对 preset；`clr` 含 PacBio RS / 非 HiFi Sequel         |
| `enrichment` | `none` / `circleseq` / `mobilome` / `ciderseq` / `t7` / `exov` / `other` | 注释，不分流                                             |

CLI `--protocol short_read|pacbio|ont` 保留为兼容缺省：仅当行上缺 `platform` 时回填。
旧列 `protocol=short_read|long_read`：只做兼容映射，不再当路由键。

`params.mode` 仍只有 `reference|eccdna`，但**不再作为总开关**。混表时以行上 `datatype/assay` 为准；`params.mode` 仅在全表同一角色且未写 `datatype` 时回填。

### 1.2 从代表集 metadata.csv 映射到上述字段

| metadata 列                                                                  | 规则                                                             |
| ---------------------------------------------------------------------------- | ---------------------------------------------------------------- |
| `workflow=circrna` 或 `library_source=TRANSCRIPTOMIC`                        | 拒绝，不入 circdna.nf                                            |
| `platform`                                                                   | `ILLUMINA→illumina`，`PACBIO_SMRT→pacbio`，`OXFORD_NANOPORE→ont` |
| `analysis_role=reference_genome`                                             | `datatype=gdna`                                                  |
| `analysis_role=primary_eccDNA`                                               | `datatype=eccdna`                                                |
| `experimental_type=CIDER-seq`                                                | `assay=ciderseq`，`enrichment=ciderseq`                          |
| 长读 + 富集 + RCA + 非 CIDER                                                 | `assay=rca`                                                      |
| 长读 + 富集 + **无 RCA**（如 T7-only）                                       | `assay=enriched`                                                 |
| 短读 + 富集（eccDNA-seq / mobilome-seq / Circle-Seq / ExoV+Phi29 / circSeq） | `assay=circleseq`                                                |
| 未富集 WGS / genome-skimming / size fractionation 背景                       | `assay=wgs`                                                      |
| `library_strategy`                                                           | **忽略**（`Beta_circSeq` 标了 WGS 但是富集样本）                 |
| `amplification_method=RCA` 且 `platform=illumina`                            | 仍是 `circleseq`，**不**改成 `rca`                               |
| notes 含 linearized / T7 debranching                                         | `concatemer=false`                                               |
| 长读 RCA 且无线性化/debranch 描述                                            | `concatemer=true`                                                |
| `assay=wgs` 或 `assay=enriched`                                              | `concatemer=false`                                               |
| PacBio HiFi / 明确 HiFi WGS                                                  | `read_type=hifi`                                                 |
| PacBio RS / RS II / 非 HiFi Sequel                                           | `read_type=clr`                                                  |
| ONT                                                                          | `read_type=ont`                                                  |
| Illumina PE/SE                                                               | `read_type=pe/se`                                                |
| 同一 `bioproject`（或手工指定）                                              | 填同一 `pair`                                                    |

`enrichment` 别名：`mobilome-seq` / `Mobilome-seq` / `circSeq` → `mobilome` 或 `circleseq`（短读都走 circleseq 引擎）。

### 1.3 合法组合（校验表）

```text
# 短读检测 / 背景
illumina + wgs       + gdna
illumina + circleseq + eccdna     # SE、PE 都允许

# 长读检测
pacbio + ciderseq  + eccdna
pacbio + rca       + eccdna       # concatemer 决定是否 TideHunter
ont    + rca       + eccdna
ont    + enriched  + eccdna       # 无 RCA 富集（向日葵 T7）
pacbio + enriched  + eccdna       # 预留，代表集暂无

# 长读背景（代表集 5 条，上一稿缺失）
pacbio + wgs + gdna
ont    + wgs + gdna

# 拒绝
TRANSCRIPTOMIC / workflow=circrna
ciderseq + illumina
rca + gdna
enriched + gdna
wgs + eccdna 且 enrichment=none     # 未富集却当 eccDNA
circleseq + pacbio/ont
assay=rca + platform=illumina       # Illumina RCA 必须写成 circleseq
```

`concatemer` 只对 `assay=rca` 有意义；其它 assay 校验时若缺省则强制 `false`。

### 1.4 代表集 32 条落点

| 样本                                  | platform | assay        | datatype | concatemer | read_type | pair 建议                                |
| ------------------------------------- | -------- | ------------ | -------- | ---------- | --------- | ---------------------------------------- |
| Oryza_PacBio_WGS                      | pacbio   | wgs          | gdna     | false      | hifi      | PRJEB59090                               |
| Oryza_ILLUMINA_eccDNA                 | illumina | circleseq    | eccdna   | —          | pe        | PRJEB59090                               |
| Oryza_ILLUMINA_WGS                    | illumina | wgs          | gdna     | —          | pe        | PRJEB59090                               |
| Wheat_ONT_eccDNA                      | ont      | rca          | eccdna   | true       | ont       | PRJEB72688                               |
| Wheat_ONT_WGS_size                    | ont      | wgs          | gdna     | false      | ont       | PRJEB72688                               |
| Wheat_ONT_WGS                         | ont      | wgs          | gdna     | false      | ont       | PRJEB72688                               |
| Arabidopsis_ONT_eccDNA                | ont      | rca          | eccdna   | true       | ont       | PRJEB46420                               |
| Arabidopsis_ILLUMINA_eccDNA           | illumina | circleseq    | eccdna   | —          | se        | PRJEB46420                               |
| Arabidopsis_ONT_eccDNA_SRP435029      | ont      | rca          | eccdna   | true       | ont       | SRP435029                                |
| Arabidopsis_ONT_WGS                   | ont      | wgs          | gdna     | false      | ont       | SRP435029                                |
| Arabidopsis_ILLUMINA_eccDNA_SRP435029 | illumina | circleseq    | eccdna   | —          | se        | SRP435029                                |
| Artemisia_ILLUMINA_eccDNA             | illumina | circleseq    | eccdna   | —          | pe        | PRJEB63080                               |
| Artemisia_gDNA                        | illumina | wgs          | gdna     | —          | pe        | PRJEB63080                               |
| Amaranthus_CIDER                      | pacbio   | ciderseq     | eccdna   | —          | clr       | SRP346396（可选再与 HiFi WGS 手工 pair） |
| Amaranthus_HiFi_WGS                   | pacbio   | wgs          | gdna     | false      | hifi      | SRP528384                                |
| Alopecurus_RCA                        | pacbio   | rca          | eccdna   | **false**  | clr       | SRP460731                                |
| Tomato_Za                             | ont      | rca          | eccdna   | **false**  | ont       | SRP490261                                |
| Sunflower_T7                          | ont      | **enriched** | eccdna   | false      | ont       | SRP552546                                |
| Arabidopsis_ddm1                      | ont      | rca          | eccdna   | true       | ont       | SRP439050                                |
| Nicotiana_ONT_eccDNA                  | ont      | rca          | eccdna   | true       | ont       | SRP658070                                |
| Nicotiana_ILLUMINA_WGS                | illumina | wgs          | gdna     | —          | pe        | 手工 pair 到 ONT eccDNA                  |
| Cryptomeria_eccDNA                    | illumina | circleseq    | eccdna   | —          | pe        | DRP012371                                |
| Arabidopsis_WGS1 / WGS2               | illumina | wgs          | gdna     | —          | pe        | 各自 study                               |
| Arabidopsis_eccDNA3                   | illumina | circleseq    | eccdna   | —          | pe        | ERP138467                                |
| Oryza_Japonica_eccDNA                 | illumina | circleseq    | eccdna   | —          | pe        | SRP093854                                |
| Daucus_eccDNA                         | illumina | circleseq    | eccdna   | —          | pe        | SRP310078                                |
| Cynodon_Salt3 / CK1                   | illumina | circleseq    | eccdna   | —          | pe        | **同一 pair**（双 eccDNA，无 gdna）      |
| Lycium_G-goji                         | illumina | circleseq    | eccdna   | —          | pe        | PRJCA020795                              |
| Beta_circSeq                          | illumina | circleseq    | eccdna   | —          | pe        | ERP129643                                |
| Oryza_RNASeq / Wheat_RNASeq           | —        | —            | —        | —          | —         | **拒绝**                                 |

---

## 2. 主流程分流

文件：[`workflows/circdna.nf`](https://github.com/SiYangming/circdna.nf/blob/circdnalr/workflows/circdna.nf)

删掉整段 `if (params.protocol in ["pacbio","ont"])` 包死长读。改为：

1. 始终 `INPUT_CHECK`
2. `ch_reads.branch`：
   - `sr_gdna`：`platform==illumina && datatype==gdna`
   - `sr_ecc`：`platform==illumina && datatype==eccdna`
   - `lr_gdna`：`platform in pacbio,ont && assay==wgs && datatype==gdna`
   - `lr_rca`：`assay==rca`
   - `lr_cider`：`assay==ciderseq`
   - `lr_enriched`：`assay==enriched`
3. 未匹配 → `exit 1`，禁止 `generic: true`
4. 短读 gdna → 现有 `REFERENCE_MODE`（BAM + mosdepth）
5. 短读 eccdna → 现有 `ECCDNA_MODE` / slim / blackbox / legacy identifier
6. 长读 gdna → 新子流程 `LONG_READ_REFERENCE`（§4.4）
7. 长读 rca / enriched / ciderseq → §4
8. `--long_read_identifier` 仍是引擎并跑开关，但受 assay/`read_type`/`concatemer` 约束（§3）
9. 混表允许：同一 run 里 Illumina eccDNA + PacBio WGS + ONT RCA 同时跑
10. ECCsplorer clu 的 join **仅** `illumina×circleseq×eccdna` ⋈ `illumina×wgs×gdna` 且 `pair` 相同。跨平台 pair 不进 clu，只保留 meta.pair 供下游 eccdna.smk

`pair` 语义：

| 场景                    | 例子                         | Nextflow                               |
| ----------------------- | ---------------------------- | -------------------------------------- |
| 短读 eccDNA + 短读 gDNA | 青蒿                         | slim clu / comparative                 |
| 长读 eccDNA + 长读 gDNA | 小麦 ONT                     | 检测 + reference 并行，不 join 引擎    |
| 跨平台 eccDNA + gDNA    | 水稻、烟草                   | 并行，不 join                          |
| 双 eccDNA 无 gDNA       | 狗牙根 Salt/CK、拟南芥 46420 | 各自检测；clu 跳过，不要 `join` 丢样本 |
| 仅 eccDNA               | 看麦娘、番茄、向日葵         | 只跑检测                               |

---

## 3. 引擎约束

`--circle_identifier`：只作用于 `sr_ecc`。
`--long_read_identifier`：只作用于长读 **检测** 样本（rca / ciderseq / enriched），**绝不**作用于 `lr_gdna`。

| 条件                                               | 允许的 identifier                         | 禁止                                    |
| -------------------------------------------------- | ----------------------------------------- | --------------------------------------- |
| `assay=ciderseq`                                   | `ciderseq`（可自动打开）                  | circleseeker / tidehunter / cresil 默认 |
| `assay=rca` + `concatemer=true` + `read_type=hifi` | `circleseeker,cresil,fled,flye,eccfinder` | ciderseq                                |
| `assay=rca` + `concatemer=true` + `read_type=ont`  | `cresil,fled,flye,eccfinder`              | circleseeker 默认关；ciderseq           |
| `assay=rca` + `concatemer=true` + `read_type=clr`  | `cresil,fled,eccfinder`（map 为主）       | circleseeker；map-hifi                  |
| `assay=rca` + `concatemer=false`                   | 同上但 **跳过 TideHunter unit/asm**       | TideHunter 默认                         |
| `assay=enriched`                                   | `cresil,fled,eccfinder`（map）            | TideHunter、ciderseq、circleseeker      |
| `assay=wgs` + 长读                                 | 无检测引擎                                | 任何 long_read_identifier               |

minimap2 preset：

```text
read_type=hifi → -x map-hifi
read_type=clr  → -x map-pb
read_type=ont  → -x map-ont
```

CircleSeeker 仅 `read_type=hifi && assay=rca && concatemer=true`。
CIDER-seq 样本即使 `platform=pacbio` 也不跑 CircleSeeker。

---

## 4. 子流程改动

### 4.1 RCA：接到现有链

复用 [`subworkflows/local/ecc_finder_ont_slim/main.nf`](https://github.com/SiYangming/circdna.nf/blob/circdnalr/subworkflows/local/ecc_finder_ont_slim/main.nf)：

- 提升为 RCA 默认切分/组装路径，不限于 `ecc_finder_map_ont_slim` identifier
- 输入：`LONG_READ_PREPROCESSING` 后的 FASTQ
- preset 跟 `meta.read_type`，不要写死 `map-ont`
- `concatemer=false` 时：不跑 `TIDEHUNTER_UNIT` / `TIDEHUNTER_ASM`；只跑 minimap2 map + 后续 peak/merge（若 identifier 含 eccfinder map）
- TideHunter 参数继续 `eccfinder_ont_*`

CircleSeeker：

- 只吃 **原始** 预处理 FASTQ
- 仅 hifi + rca + concatemer
- 已有 [`bin/circleseeker_to_bed.py`](https://github.com/SiYangming/circdna.nf/blob/circdnalr/bin/circleseeker_to_bed.py) 保留
- 空 CSV 补文件移到 bed 转换 **之前**（[`modules/local/circleseeker/main.nf`](https://github.com/SiYangming/circdna.nf/blob/circdnalr/modules/local/circleseeker/main.nf)）

### 4.2 CIDER-seq 独立

[`subworkflows/local/ciderseq_pipeline/main.nf`](https://github.com/SiYangming/circdna.nf/blob/circdnalr/subworkflows/local/ciderseq_pipeline/main.nf)

- 不与 RCA FASTA `mix`
- 接上已声明未使用的 `params.ciderseq_host_genome`（可选宿主比对/注释）
- `cseccdna` / `cstools` 按需入主流程，默认关
- 产物留在 `long_read/ciderseq/`；无宿主坐标则不进植物 ecc 统一 BED

### 4.3 无 RCA 富集（`assay=enriched`）

不新建平行检测工具。预处理后走 CRESIL/FLED/ecc_finder **map**（由 identifier 选）。
跳过 TideHunter、CircleSeeker、CIDER-seq2。
向日葵 T7 走这条。

### 4.4 长读 WGS 背景（新）

新建 `subworkflows/local/long_read_reference/main.nf`：

```text
LONG_READ_PREPROCESSING → MINIMAP2_ALIGN（preset=read_type）
  → SAMTOOLS_SORT/INDEX → MOSDEPTH
  → 可选 organelle_tag
```

禁止任何 eccDNA identifier。
对应：水稻 PacBio WGS、苋菜 HiFi WGS、小麦两条 ONT WGS、拟南芥 ONT WGS。

短读 gdna 继续 `REFERENCE_MODE`，不要把长读塞进 BWA。

### 4.5 无坐标组装物才 remap

需要：TideHunter consensus、ecc_finder asm FASTA、Flye contig。
不需要：CircleSeeker、CReSIL/FLED BED、CIDER-seq2 病毒 monomer、长读 WGS。

新建 `subworkflows/local/remap_assembled_circles/main.nf`，组合已有模块：

- `MINIMAP2_ALIGN`
- `SAMTOOLS_SORT` / `SAMTOOLS_INDEX`
- `BEDTOOLS_BAMTOBED`
- 薄脚本 `bin/collapse_circle_alignments.py`：按 query 折叠、丢 unmapped/低 MapQ、处理 SA、输出 BED6+`read_count`

preset 跟 `meta.read_type`。
Flye 输出目前没进过滤，remap 后并入 `LONG_READ_FILTERING`。

不要建 `modules/local/minimap2_remap.nf`。

### 4.6 统一植物 eccDNA BED 契约

```text
chr  start  end  name  score  strand  read_count  [te_overlap]  [origin]  [pair]  [engine]
```

| 引擎                  | 动作                                 |
| --------------------- | ------------------------------------ |
| CircleSeeker          | 已有 BED，补可选列                   |
| CReSIL                | 确认 `convert_cresil_to_bed.py` 列名 |
| FLED                  | junctions → 同一 BED                 |
| ecc_finder map        | 用现有坐标                           |
| ecc_finder asm / Flye | §4.5 remap                           |
| CIDER-seq2            | 默认不进；宿主比对出坐标才进         |
| 长读/短读 WGS         | 不进候选 BED，只出 mosdepth          |

全部检测候选 `mix` 进 [`LONG_READ_FILTERING`](https://github.com/SiYangming/circdna.nf/blob/circdnalr/subworkflows/local/long_read_filtering/main.nf)（长读）或短读现有过滤。

`repeats_bed`：**改为 annotate**（`te_overlap`），禁止 `bedtools intersect` 硬删。
`blacklist_bed` 仍可硬过滤。

短读契约保持：`mosdepth_bed` + `eccsplorer_bed` + `circle_map_bed`。
长读契约新增：`long_read_bed` + `mosdepth_bed`（gdna 长读）。
eccdna.smk 消费这些文件做评分，不在 Nextflow 里算 ECC_SCORE。

### 4.7 细胞器：打 origin，默认不丢

新建 `subworkflows/local/organelle_tag/main.nf` + `bin/tag_organelle_origin.py`。

- `params.filter_organelle=false`，`organelle_genome=null`，`drop_organelle_candidates=false`
- 打开且给了 FASTA：短读 BWA、长读 minimap2，比对 organelle（或核+细胞器 concat）
- 在候选 BED 打 `origin=nuclear|pt|mt|ambiguous`
- 仅 `drop_organelle_candidates=true` 时从目录去掉 organelle 候选
- 不做 FastQC 后抽 unmapped 当主输入
- 若以后要做 read 级减速：PE 必须双方 unmapped（`-f 12` + `-1/-2/-s`），organelle BAM 单独 publish

---

## 5. Samplesheet 校验与 channel

### 5.1 [`bin/check_samplesheet.py`](https://github.com/SiYangming/circdna.nf/blob/circdnalr/bin/check_samplesheet.py)

- 增加 `VALID_ASSAYS`、`VALID_READ_TYPES`、`VALID_CONCATENER`
- `VALID_DATATYPES` 保持小写 `gdna|eccdna`
- 删除作为路由键的 `VALID_PROTOCOLS=short_read|long_read`；旧列只映射
- **长读分支必须透传** `datatype,platform,assay,pair,concatemer,read_type,enrichment`
- 实现 §1.3 组合断言
- `library_source=TRANSCRIPTOMIC` 或 `assay` 无法映射时失败
- 长读输出 CSV 表头包含上述列，不要只写 `sample,single_end,fastq_1,input_bam,entrypoint`
- 缺省：
  - 短读无 assay：`datatype=gdna→wgs`，`eccdna→circleseq`
  - 长读无 assay：禁止猜 ciderseq；必须显式写，或仅当全局 identifier **只有** ciderseq 时才回填
  - `concatemer` 缺省：`assay=rca` → `true`，其余 `false`
  - `read_type` 缺省：`ont→ont`，`pacbio→hifi`（**CLR 必须显式写**，避免 RS 被当成 HiFi）

### 5.2 [`assets/schema_input.json`](https://github.com/SiYangming/circdna.nf/blob/circdnalr/assets/schema_input.json)

- 删掉 datatype 文案里的 Integrated Mode
- 新增 `assay`、`concatemer`、`read_type`、`enrichment`、`pair`
- `protocol` 标 deprecated 或移除枚举
- `platform` / `datatype` 与 Python 对齐

### 5.3 [`subworkflows/local/input_check/main.nf`](https://github.com/SiYangming/circdna.nf/blob/circdnalr/subworkflows/local/input_check/main.nf)

- 短读、长读 channel 都写 `meta.platform/assay/datatype/data_type/pair/concatemer/read_type/enrichment`
- 长读仍 emit `[meta, fastq, bam, entrypoint]`
- map 函数再做一遍组合断言
- `meta.platform` 优先行，否则 `params.protocol` 兼容
- 不要未定义的 `files`

### 5.4 [`modules/local/samplesheet_check/main.nf`](https://github.com/SiYangming/circdna.nf/blob/circdnalr/modules/local/samplesheet_check/main.nf)

- 继续传 `params.protocol`；若加 `params.assay` 缺省则一并传入

---

## 6. 参数、配置、文档

### 6.1 [`nextflow.config`](https://github.com/SiYangming/circdna.nf/blob/circdnalr/nextflow.config)

只加：

```groovy
filter_organelle          = false
organelle_genome          = null
drop_organelle_candidates = false
assay                     = null   // 仅当 samplesheet 缺 assay 时回填，默认 null 强制行上填写长读
```

不加：`mode=integrated`、`w_*`、`tidehunter_kmer`、`tidehunter_min_period`。

`mode` 默认仍 `eccdna`。
`long_read_identifier` 默认保持现有；CircleSeeker 不要塞进 ONT 默认列表。

### 6.2 [`nextflow_schema.json`](https://github.com/SiYangming/circdna.nf/blob/circdnalr/nextflow_schema.json)

- 登记新参
- `mode` 枚举 `reference|eccdna`
- `--protocol` help：遗留全局平台开关，优先 samplesheet `platform`

### 6.3 [`conf/modules.config`](https://github.com/SiYangming/circdna.nf/blob/circdnalr/conf/modules.config)

- `MINIMAP2_*` 的 `ext.args` 用 `meta.read_type`（不要只用 `meta.platform`）
- Flye：`hifi→--pacbio-hifi`，`clr→--pacbio-raw`，`ont→--nano-hq`
- TideHunter 的 `withName` 仅绑 `concatemer=true` 的 RCA 任务
- 长读 reference 的 mosdepth publishDir 单独目录，避免和 eccDNA 混

### 6.4 文档

[`AGENTS.md`](https://github.com/SiYangming/circdna.nf/blob/circdnalr/AGENTS.md)

- 写清 `platform × assay × datatype × concatemer × read_type`
- 保持 “integrated 已移除，评分在 eccdna.smk”
- 写 CIDER / RCA / enriched / 长读 WGS 分流
- 写 BED 契约与 `pair` 语义
- 写 Illumina+RCA = circleseq

README / usage

- 长读 samplesheet 示例含 `platform,assay,datatype,pair,concatemer,read_type,entrypoint`
- 说明一张表可混 Illumina Circle-seq + PacBio WGS + ONT RCA
- 删掉 `mode=integrated` 残留
- 给出代表集 4 行示例（水稻三件套、小麦 ONT、苋菜 CIDER、向日葵 T7）

[`CHANGELOG.md`](https://github.com/SiYangming/circdna.nf/blob/circdnalr/CHANGELOG.md)

- breaking：长读新列；`protocol` 列不再路由；repeats 改注释；混表按样本分流
- 明确 CLR 必须写 `read_type=clr`

---

## 7. 代表集 samplesheet 落地

不要把 metadata.csv 整表当 Nextflow 输入。生成/改写运行用 samplesheet：

- 列：`sample,fastq_1,fastq_2,input_bam,entrypoint,platform,assay,datatype,pair,concatemer,read_type,enrichment`
- `sample` 用 `run_id`（与现有 `circdna_ngs_clean.csv` / `circdna_tgs_clean.csv` 一致），`sample_id` 可放 notes 或另列
- 去掉 `/data1/users/siyangming/...` 绝对路径，改文档化根目录或 `${params.input_dir}`
- 新增 `samplesheets/circdna_representative.csv`（32 条 circdna，按 §1.4 填字段）
- `test_ciderseq_lr.csv`：`platform=pacbio,assay=ciderseq,datatype=eccdna,read_type=clr`
- `test_pacbio_lr.csv` / `test_ont_lr.csv`：补列
- `test_integrated.csv`：改名 `test_paired_gdna_eccdna.csv`，避免暗示 integrated 模式
- 物种级 `circdnalr_*_long_read.csv`：补 `platform,assay,datatype,concatemer,read_type`；把现在误放进 long_read 表的 **WGS 背景**（如 `ERR11838731`）标成 `assay=wgs,datatype=gdna`，不要当 RCA
- 2 条 circrna **不要**写入 circdna samplesheet

建议同时提供一个只读映射脚本（可选，非运行时必须）：`bin/metadata_to_samplesheet.py`，把代表集 metadata.csv 转成上述运行表，规则固化为 §1.2，避免手工再漂。

---

## 8. 顺手修的现有缺陷

| 位置                                                 | 改什么                                                                 |
| ---------------------------------------------------- | ---------------------------------------------------------------------- |
| `modules/local/circleseeker/main.nf`                 | 空 CSV 在 bed 转换之前补                                               |
| `workflows/circdna.nf`                               | 接 `ciderseq_host_genome`；按 meta branch                              |
| `subworkflows/local/flye_pipeline/main.nf`           | 用 `read_type`；组装 FASTA 送 remap                                    |
| `subworkflows/local/long_read_preprocessing/main.nf` | 已按 `params.protocol` 切 PacBio/ONT；改为按 `meta.platform`，以便混表 |
| `LONG_READ_PREPROCESSING` 的 LIMA/PBCCS              | 仅 `platform=pacbio` 且 entrypoint 需要时；ONT 行不得进 PBCCS          |
| 大基因组                                             | Cryptomeria / 小麦继续 `-c conf/large_genome.config`（已有，文档点名） |

预处理混表要点：现在 `LONG_READ_PREPROCESSING` 用全局 `params.protocol == "pacbio"`。混 ONT+PacBio 时必须改成 `reads.branch { pacbio: meta.platform=='pacbio'; ont: meta.platform=='ont' }`。

---

## 9. 目标数据流

```text
samplesheet
  → INPUT_CHECK（platform × assay × datatype × concatemer × read_type × pair）
      │
      ├─ illumina + wgs + gdna        → REFERENCE_MODE（BWA + mosdepth）
      ├─ illumina + circleseq + eccdna → 现有短读检测（legacy / slim / blackbox）
      │                                  仅当同 pair 存在 illumina gdna 时才 clu join
      │
      ├─ pacbio/ont + wgs + gdna      → LONG_READ_REFERENCE（minimap2 + mosdepth）
      │
      ├─ ciderseq                     → CIDERSEQ_PIPELINE → 可选 host 注释
      │
      ├─ rca / enriched
      │     → LONG_READ_PREPROCESSING（按 meta.platform 分支）
      │     → rca + concatemer=true  → TideHunter + ecc_finder_ont_slim
      │     → rca + hifi + concatemer → 另：CircleSeeker(原始 reads)
      │     → rca + concatemer=false → 跳过 TideHunter，map 引擎
      │     → enriched               → CRESIL/FLED/ecc_finder map
      │     → 无坐标 FASTA           → remap_assembled_circles
      │     → mix → LONG_READ_FILTERING（support + blacklist + TE 注释）
      │
      └─ 可选 organelle_tag → origin 列
                ↓
         统一 BED + mosdepth  →  publish 给 eccdna.smk
```

---

## 10. 测试

| 测试                                                               | 断言                                                                         |
| ------------------------------------------------------------------ | ---------------------------------------------------------------------------- |
| `test_mixed_assay`                                                 | 一张表：illumina circleseq + pacbio wgs + ont rca + ciderseq，四条支路都启动 |
| 负例 `ciderseq+illumina`                                           | SAMPLESHEET_CHECK 失败                                                       |
| 负例 `rca+gdna`                                                    | 失败                                                                         |
| 负例 circrna 行                                                    | 失败                                                                         |
| `concatemer=false`（看麦娘/番茄 stub）                             | 无 TIDEHUNTER 进程                                                           |
| `assay=enriched`（向日葵 stub）                                    | 无 TideHunter / CircleSeeker / CIDER                                         |
| `assay=wgs` 长读                                                   | 只有 preprocess + minimap2 + mosdepth，无检测引擎                            |
| `read_type=clr`                                                    | minimap2 `-x map-pb`，无 CircleSeeker                                        |
| 青蒿 pair                                                          | slim clu join 成功                                                           |
| 狗牙根双 eccDNA                                                    | 两条都检测，clu 不 join、不丢样本                                            |
| 水稻三件套                                                         | 三条支路并行，clu 只可能 join 两条 Illumina                                  |
| CircleSeeker 空结果                                                | bed 仍产出表头                                                               |
| remap stub                                                         | FASTA → BED6+read_count                                                      |
| organelle_tag stub                                                 | 打 origin，默认不丢                                                          |
| 现有 `test_pacbio_lr` / `test_nanopore_lr` / `test_ciderseq_local` | 补列后仍通                                                                   |
| nf-test                                                            | TideHunter 输出不连接 CircleSeeker 输入                                      |

---

## 11. 文件级清单

| 动作 | 文件                                                   | 做什么                                                               |
| ---- | ------------------------------------------------------ | -------------------------------------------------------------------- |
| 改   | `bin/check_samplesheet.py`                             | assay/concatemer/read_type、组合断言、长读列透传、拒 circrna         |
| 改   | `assets/schema_input.json`                             | 新列；去掉 integrated 文案                                           |
| 改   | `subworkflows/local/input_check/main.nf`               | 短读+长读都写齐 meta                                                 |
| 改   | `modules/local/samplesheet_check/main.nf`              | 传新缺省                                                             |
| 改   | `workflows/circdna.nf`                                 | 按 meta 六路 branch；混表；clu 限制                                  |
| 改   | `subworkflows/local/long_read_preprocessing/main.nf`   | 按 `meta.platform` 分支                                              |
| 改   | `subworkflows/local/ecc_finder_ont_slim/main.nf`       | RCA 默认路径；preset=`read_type`；`concatemer=false` 跳过 TideHunter |
| 改   | `subworkflows/local/ciderseq_pipeline/main.nf`         | 可选宿主锚定；不 mix RCA                                             |
| 改   | `subworkflows/local/long_read_filtering/main.nf`       | repeats → annotate                                                   |
| 改   | `subworkflows/local/flye_pipeline/main.nf`             | `read_type`；FASTA → remap                                           |
| 改   | `modules/local/circleseeker/main.nf`                   | 空文件时序；调用条件在 workflow 层限制                               |
| 改   | `conf/modules.config`                                  | preset / Flye / publishDir                                           |
| 改   | `nextflow.config`                                      | organelle 三开关 + 可选 `assay=null`                                 |
| 改   | `nextflow_schema.json`                                 | 同步                                                                 |
| 改   | `AGENTS.md` / `README.md` / `CHANGELOG.md`             | 路由、契约、breaking                                                 |
| 改   | 全部相关 `samplesheets/*.csv`                          | 新列；WGS 背景标对；去绝对路径                                       |
| 新建 | `subworkflows/local/long_read_reference/main.nf`       | 长读 gdna：map + mosdepth                                            |
| 新建 | `subworkflows/local/remap_assembled_circles/main.nf`   | 组合 nf-core 模块                                                    |
| 新建 | `bin/collapse_circle_alignments.py`                    | 折叠比对 → 统一 BED                                                  |
| 新建 | `subworkflows/local/organelle_tag/main.nf`             | 候选打 origin                                                        |
| 新建 | `bin/tag_organelle_origin.py`                          | 分类                                                                 |
| 新建 | `samplesheets/circdna_representative.csv`              | 32 条已映射运行表                                                    |
| 新建 | `bin/metadata_to_samplesheet.py`                       | metadata.csv → 运行表（可选）                                        |
| 新建 | mixed-assay / 负例 samplesheet + test config + nf-test | §10                                                                  |
| 不建 | `modules/local/minimap2_remap.nf`                      | 用模块组合                                                           |
| 不建 | `modules/local/integrated_eval.nf`                     | 留 eccdna.smk                                                        |
| 不建 | `bin/calculate_ecc_score.py`                           | 同上                                                                 |
| 不建 | `subworkflows/local/pacbio_ecc.nf`                     | 分流进现有子流程                                                     |
| 不建 | read 级 `filter_organelle.nf`                          | 用 tag 子流程                                                        |

---

## 12. 验收标准

1. `metadata.csv` 中 32 条 `circdna` 都能写入运行 samplesheet，且合法组合校验通过。
2. 2 条 `circrna` 被拒绝。
3. 一张混合表一次 `nextflow run` 能同时跑：Illumina Circle-seq、长读 WGS mosdepth、ONT/PacBio RCA、CIDER-seq。
4. `ERR11838731` 等长读 WGS **不**进 CircleSeeker/TideHunter。
5. `SRR16958693` 走 CIDER-seq2，`read_type=clr`，不走 CircleSeeker。
6. `SRR26069818` / `SRR28004411` 不跑 TideHunter。
7. `SRR31773424` 走 `enriched` map 引擎。
8. 青蒿 gDNA/eccDNA 仍能 clu；狗牙根两条 eccDNA 都不丢。
9. 植物 eccDNA 候选都是同一 BED 契约；TE 只注释不硬删。
10. 细胞器默认只打 `origin`。
11. 评分仍在 eccdna.smk；Nextflow 只发布 BED + mosdepth。

按此修订后，代表集的样本组合与引擎选择一一对应，不再出现「长读背景被当 RCA」或「线性化 RCA 被 TideHunter 误切」的缺口。
