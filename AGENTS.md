# AGENTS.md — nf-core/circdna Pipeline Architecture & Change Rules

> **Purpose**: This document defines the complete architecture, data flow, and change management rules for the `nf-core/circdna` Nextflow pipeline. Any AI assistant making changes to this pipeline MUST follow these rules.

> **UNIVERSAL RULES**: Before reading this document, first check `/Users/siyangming/nextflow_nf_core/AGENTS.md` for universal Nextflow pipeline conventions, CHANGELOG rules, and version bumping policies. This document is the pipeline-specific supplement.

---

## 1. Pipeline Overview

**Pipeline**: `nf-core/circdna` v4.2.0
**Main Entry**: `main.nf` → `workflows/circdna.nf`
**Two Modes**:
| Mode | Description |
|------|-------------|
| `reference` | gDNA/WGS depth analysis only (BAM prep → mosdepth) |
| `eccdna` | eccDNA (circular DNA) detection (legacy + slim/blackbox chains) |

**Two Protocols** (via `--protocol`):
| Protocol | Description |
|----------|-------------|
| `short_read` | Illumina NGS short-read detection (default) |
| `pacbio` / `ont` | Long-read detection via `--long_read_identifier` (cresil/fled/flye/eccfinder/circleseeker) |

> **Note**: The `integrated` mode was removed in v4.0.0. ECC_SCORE scoring was moved to `eccdna.smk` (independent post-processing layer, located at `/Users/siyangming/nextflow_nf_core/eccdna.smk`), which consumes circdna.nf detection outputs.

---

## 2. Directory Architecture

```
circdna.nf/
├── main.nf                          # ENTRY POINT: orchestrates workflow lifecycle
├── nextflow.config                  # GLOBAL CONFIG: params, profiles, manifests
├── nextflow_schema.json             # PARAMETER SCHEMA: auto-generated from params
├── modules.json                     # MODULE MANIFEST: all module dependencies
│
├── workflows/
│   └── circdna.nf                   # TOP-LEVEL WORKFLOW: mode selection + data routing
│
├── subworkflows/local/              # CUSTOM SUBWORKFLOWS
│   ├── input_check/main.nf          #   Samplesheet parsing + validation
│   ├── bam_preprocessing/main.nf    #   FASTQC → TrimGalore → BWA → Picard → SAMtools
│   ├── eccdna_mode/main.nf          #   eccDNA mode: BAM prep + mosdepth + Circle-Map
│   ├── reference_mode/main.nf       #   Reference mode: BAM prep + mosdepth only
│   ├── circle_finder_pipeline/main.nf   # Legacy: Circle-Finder (SAMBLASTER → BEDTools → CircleFinder)
│   ├── circle_map_pipeline/main.nf      # Legacy/Active: Circle-Map (ReadExtractor → Repeats → Realign)
│   ├── ampliconarchitect_pipeline/main.nf  # Legacy: AmpliconArchitect + CNVkit
│   ├── unicycler_pipeline/main.nf       # Legacy: Unicycler circular assembly
│   ├── eccsplorer_pipeline/main.nf      # ECCsplorer 检测（黑盒链）
│   ├── eccsplorer_slim_pipeline/main.nf # ECCsplorer slim 链（map）
│   ├── eccsplorer_clu_slim/main.nf      # ECCsplorer slim 链（clu）
│   ├── eccsplorer_all_slim/main.nf      # ECCsplorer slim 链（all）
│   ├── eccsplorer_prexer_slim/main.nf   # ECCsplorer slim 预处理（prexer）
│   ├── ecc_finder_pipeline/main.nf      # ecc_finder 检测（map/asm × SR/ONT 黑盒）
│   ├── ecc_finder_slim_pipeline/main.nf # ecc_finder slim 链
│   ├── ecc_finder_ont_slim/main.nf      # ecc_finder slim 链（ONT）
│   ├── cresil_pipeline/main.nf          # CReSIL 长读检测
│   ├── fled_pipeline/main.nf            # FLED 长读检测
│   ├── flye_pipeline/main.nf            # Flye 长读组装检测
│   ├── circleseeker_pipeline/main.nf    # CircleSeeker 长读检测
│   ├── long_read_preprocessing/main.nf  # 长读预处理（PBCCS/LIMA/Pychopper/Chopper）
│   ├── long_read_mapping/main.nf        # 长读比对（minimap2）
│   ├── long_read_filtering/main.nf      # 长读过滤（按读长/支持数）
│   └── utils_nfcore_circdna_pipeline/main.nf  # Shared utility: methods description
│
├── subworkflows/nf-core/            # NF-CORE SUBWORKFLOWS (vendored)
│   ├── utils_nfcore_pipeline/main.nf      # Version collection, MultiQC summary
│   ├── bam_markduplicates_picard/main.nf  # Picard MarkDuplicates workflow
│   └── bam_stats_samtools/main.nf         # SAMtools stats collection
│
├── modules/local/                   # CUSTOM MODULES
│   ├── eccsplorer/                  #   ECCsplorer 检测（含 tests）
│   ├── eccsplorer_slim/             #   ECCsplorer slim 原子模块（normalize/peak_detect/clu_candidates/html_report/...）
│   ├── ecc_finder/                  #   ecc_finder（map_ont/map_sr/asm_ont/asm_sr）
│   ├── ecc_finder_slim/             #   ecc_finder slim 原子模块（split_detect/paf_filter/merge_score/...）
│   ├── cresil/                      #   CReSIL（trim/identify/annotate/visualize + identify_wgls）
│   ├── circleseeker/                #   CircleSeeker（含 tests）
│   ├── circlefinder/main.nf         #   CircleFinder: circular DNA detection
│   ├── circlemap/                   #     Circle-Map sub-modules
│   │   ├── readextractor/main.nf
│   │   ├── realign/main.nf
│   │   └── repeats/main.nf
│   ├── getcircularreads/main.nf     #   Extract circular reads for Unicycler
│   ├── ampliconsuite/main.nf        #   AmpliconArchitect + AmpliconClassifier
│   ├── ampliconsuite_ec/            #   AmpliconSuite EC
│   ├── bedtools/                    #     BEDTools sub-modules
│   │   ├── sortedbam2bed/main.nf
│   │   └── splitbam2bed/main.nf
│   ├── fled/                        #   FLED 检测
│   ├── flye/                        #   Flye 组装（含 tests）
│   ├── tidehunter/                  #   TideHunter 串联重复检测（长读）
│   ├── segemehl/haarz/              #   Segemehl HaarZ
│   ├── haarz/                       #   HaarZ 独立模块
│   ├── genrich/                     #   Genrich peak calling
│   ├── repeatexplorer2/             #   RepeatExplorer2 预处理
│   ├── filter_eccdna_by_support/    #   按读支持数过滤 eccDNA
│   └── samplesheet_check/main.nf    #   Samplesheet validation
│
├── modules/nf-core/                 # NF-CORE MODULES (vendored from nf-core/modules)
│   ├── bwa/                         #     BWA index + mem
│   ├── fastqc/                      #     FastQC quality check
│   ├── trimgalore/ + trimmomatic/   #     Trimming
│   ├── samtools/                    #     SAMtools (sort, index, flagstat, idxstats, faidx, view, stats)
│   ├── picard/                      #     Picard MarkDuplicates
│   ├── samblaster/                  #     SAMBLASTER split-read detection
│   ├── bedtools/                    #     BEDTools (bamtobed/coverage/genomecov/getfasta/groupby/intersect/makewindows/merge/sort)
│   ├── circexplorer2/parse/         #     CIRCexplorer2 parse
│   ├── cnvkit/                      #     CNVkit batch + segment
│   ├── minimap2/                    #     Minimap2 align + index（长读）
│   ├── pbccs/ + lima/               #     PacBio CCS + 引物拆分（长读）
│   ├── pychopper/ + chopper/        #     ONT 修剪 + 过滤（长读）
│   ├── nanoplot/                    #     NanoPlot 长读 QC
│   ├── cdhit/                       #     cd-hit-est 去冗余（ecc_finder）
│   ├── mosdepth/                    #     mosdepth coverage
│   ├── multiqc/                     #     MultiQC aggregation
│   ├── seqtk/seq/ + cat/fastq/      #     FASTQ 工具 + 合并
│   ├── blast/                       #     BLAST（ECCsplorer 注释）
│   ├── genrich/                     #     Genrich（nf-core 版）
│   ├── flye/ + unicycler/           #     组装工具
│   └── segemehl/                    #     Segemehl（nf-core 版）
│
├── conf/                            # CONFIGURATION FILES
│   ├── base.config                  #   Default resource allocation (cpus/memory/time labels)
│   ├── modules.config               #   Per-module ext.args, ext.prefix, publishDir overrides
│   ├── igenomes.config              #   Reference genome paths (ref genomes)
│   ├── igenomes_ignored.config      #   Skip ref genome paths
│   ├── large_genome.config          #   大基因组 CSI 索引（SAMTOOLS_INDEX -c）
│   ├── server.config                #   Server execution profile
│   ├── test.config                  #   GitHub CI test profile (smallest)
│   ├── test_local.config            #   Local test profile (4 CPU, 8GB)
│   ├── test_local_gdna.config       #   gDNA/reference 模式本地测试
│   ├── test_bam_local.config        #   BAM 输入本地测试
│   ├── test_AA.config               #   AmpliconArchitect CI 测试
│   ├── test_AA_local.config         #   AmpliconArchitect 本地测试
│   ├── test_nanopore_lr.config      #   ONT 长读测试
│   ├── test_pacbio_lr.config        #   PacBio 长读测试
│   └── test_full.config             #   Full test profile
│
├── samplesheets/                    # ALL SAMPLESHEETS (single source of truth)
│   ├── circdna_ngs_clean.csv        #   **MAIN INPUT**: NGS short-read data（`data_type` 列，12 物种）
│   ├── circdna_tgs_clean.csv        #   **MAIN INPUT**: TGS long-read FASTQ
│   ├── circdna_{species}_eccDNA.csv #   Per-species NGS files (auto-generated by script)
│   ├── circdnalr_{species}_long_read.csv # Per-species TGS files (hand-maintained)
│   ├── circrna_*.csv                #   circRNA samplesheets (cross-project)
│   ├── SraRunInfo_eccDNA_all2.csv   #   SRA Run Info metadata (TaxID, ScientificName, download_path)
│   ├── update_samplesheets.py       #   Auto-generate per-species samplesheets
│   ├── test_local_eccdna.csv        #   Local test samplesheet
│   ├── test_*.csv                   #   Additional test samplesheets
│   └── data_issues.txt              #   Known data issues log
│
├── bin/                             # BIN SCRIPTS (available to Nextflow processes)
│   ├── check_samplesheet.py         #   VALIDATOR: samplesheet format checker
│   ├── Coverage.py                  #   Coverage calculation
│   ├── bam2bam.py                   #   BAM format conversion
│   ├── circle_map.py                #   Circle-Map integration
│   ├── extract_circle_SV_reads.py   #   SV read extraction
│   ├── realigner.py                 #   Read realignment
│   ├── repeats.py                   #   Repeat analysis
│   ├── scrape_software_versions.py  #   Version scraping
│   ├── simulations.py               #   Simulation utilities
│   ├── summarise_aa.py              #   AmpliconArchitect summary
│   ├── circleseeker_to_bed.py       #   CircleSeeker 结果转 BED
│   ├── convert_cresil_to_bed.py     #   CReSIL 结果转 BED
│   ├── patch_cresil.py / patch_wgls.py  #   CReSIL/WGLS 兼容补丁
│   ├── filter_by_read_support.py    #   按读支持数过滤
│   └── utils.py                     #   Shared utilities
│
├── lib/                             # GROOVY LIBRARY
│   ├── WorkflowMain.groovy          #   Main workflow helpers (citation, initialise, getGenomeAttribute)
│   ├── WorkflowCircdna.groovy       #   Pipeline-specific helpers
│   ├── NfcoreTemplate.groovy        #   nf-core template utilities (logo, email, summary)
│   └── Utils.groovy                 #   General utilities
│
├── assets/                          # STATIC ASSETS
│   ├── multiqc_config.yml           #   MultiQC custom config
│   ├── schema_input.json            #   Input schema
│   ├── samplesheet.csv              #   Example samplesheet
│   └── methods_description_template.yml  # Methods description template
│
├── testdatasets/                    # TEST DATA (FASTQ files + reference)
│   ├── ngs/ ont/ pacbio/            #   Test FASTQ files (按协议分目录)
│   ├── reference/                   #   Test reference genome
│   ├── annotation/                  #   注释文件（yeast_genes.bed 等）
│   ├── cnvkit/                      #   CNVkit reference
│   ├── eccsplorer_db/               #   ECCsplorer 数据库
│   ├── mosek/                       #   Mosek license
│   └── README.md                    #   Test dataset documentation
│
└── docs/                            # DOCUMENTATION
    ├── README.md                    #   Pipeline documentation
    ├── usage.md                     #   Usage guide
    └── output.md                    #   Output description
```

---

## 3. Data Flow Architecture

```
┌─────────────────────────────────────────────────────────────────────┐
│                        INPUT LAYER                                  │
│                                                                     │
│  samplesheet.csv ──► INPUT_CHECK ──► ch_reads ──► CAT_FASTQ         │
│  (sample,fastq_1, │  (validate +   (meta, reads)  (merge lanes)    │
│   fastq_2,        │  add metadata)                                │
│   data_type)      │                                                │
└────────────────────────────┬────────────────────────────────────────┘
                             │
                             ▼
┌─────────────────────────────────────────────────────────────────────┐
│                     QC & TRIMMING LAYER                              │
│                                                                     │
│  ch_cat_fastq ──► FASTQC ──► TRIMGALORE ──► ch_trimmed_reads        │
│                       │                       │                     │
│                       ▼                       ▼                     │
│                  ch_fastqc_multiqc      ch_trimgalore_multiqc       │
└────────────────────────────┬────────────────────────────────────────┘
                             │
                             ▼
┌─────────────────────────────────────────────────────────────────────┐
│                       ALIGNMENT LAYER                                │
│                                                                     │
│  ch_trimmed_reads ──► BWA_INDEX ──► BWA_MEM ──► BAM_PREPROCESSING    │
│                                               │                     │
│                                               ▼                     │
│                              (Picard MarkDup + SAMtools sort/index) │
│                                               │                     │
│                                               ▼                     │
│                            ch_bam_sorted + ch_bam_sorted_bai         │
└────────────────────────────┬────────────────────────────────────────┘
                             │
         ┌───────────────────┴───────────────────┐
         ▼                                       ▼
┌───────────────┐                      ┌───────────────────────┐
│ REFERENCE_MODE │                      │  ECCDNA_MODE          │
│               │                      │  MOSDEPTH             │
│ MOSDEPTH      │                      │  CIRCLE_MAP_PIPELINE  │
│ (gDNA depth)  │                      │  (short-read 默认链)   │
└───────┬───────┘                      └───────────┬───────────┘
        │                                          │
        ▼                                          ▼
  mosdepth_bed                              merged_bed
  (gDNA)                                   (eccDNA candidates)
                                            (+ legacy: circle_finder /
                                              ampliconarchitect / unicycler /
                                              circexplorer2)
                                            (+ blackbox/slim: eccsplorer /
                                              ecc_finder 链)

┌─────────────────────────────────────────────────────────────────────┐
│                LONG-READ PATH (protocol=pacbio|ont)                 │
│                                                                     │
│  ch_long_reads ──► LONG_READ_PREPROCESSING (PBCCS/LIMA/Pychopper)   │
│         ──► LONG_READ_MAPPING (minimap2)                            │
│         ──► LONG_READ_FILTERING                                     │
│         ──► CRESIL_PIPELINE / FLED_PIPELINE / FLYE_PIPELINE         │
│             / ECC_FINDER_PIPELINE / CIRCLESEEKER_PIPELINE           │
│             （由 --long_read_identifier 选择）                       │
│         ──► FILTER_ECCDNA_BY_SUPPORT                                │
└─────────────────────────────────────────────────────────────────────┘
```

---

## 4. Samplesheet Specification

### 4.1 NGS (FASTQ) Format

```csv
sample,fastq_1,fastq_2,data_type
ERR1830502,/data1/users/siyangming/PlanteccDNADB/eccDNA/Oryza_sativa/ERR1830502_1.fastq.gz,/data1/users/siyangming/PlanteccDNADB/eccDNA/Oryza_sativa/ERR1830502_2.fastq.gz,eccDNA
SRR5051136,/data1/users/siyangming/PlanteccDNADB/eccDNA/Oryza_sativa/SRR5051136_1.fastq.gz,/data1/users/siyangming/PlanteccDNADB/eccDNA/Oryza_sativa/SRR5051136_2.fastq.gz,eccDNA
```

- **Required columns**: `sample`, `fastq_1`, `fastq_2`
- **Optional columns**: `lane`, `datatype`, `platform`, `protocol`
- **`data_type` values**: `eccDNA` or `gDNA` (used by eccdna mode to split channels for ECCsplorer control analysis)
- **`datatype` values**: `eccdna` or `gdna` (used by `check_samplesheet.py` for validation, normalized to lowercase)
- **`platform` values**: `illumina`, `pacbio`, `ont`
- **`protocol` values**: `short_read`, `long_read`
- **File path convention**: `/data1/users/siyangming/PlanteccDNADB/eccDNA/{Species}/{sample}_{1,2}.fastq.gz`

### 4.2 TGS (Long-read) Format — single-end FASTQ with routing columns

```csv
sample,fastq_1,fastq_2,input_bam,entrypoint,platform,assay,datatype,pair,concatemer,read_type,enrichment
ERR11838731,/data1/users/siyangming/PlanteccDNADB/eccDNA/Oryza_sativa/ERR11838731.fastq.gz,,,cleaned_fastq,pacbio,wgs,gdna,PRJEB59090,false,hifi,none
ERR12724336,/data1/users/siyangming/PlanteccDNADB/eccDNA/Triticum_aestivum/ERR12724336.fastq.gz,,,cleaned_fastq,ont,rca,eccdna,PRJEB72688,true,ont,circleseq
```

- 长读样本表必须含 `platform,assay,datatype,concatemer,read_type` 路由列（`assay` 长读必填；`protocol` 列仅做兼容映射，不再是路由键）
- 长读 WGS 背景标 `assay=wgs,datatype=gdna`，绝不当作 RCA 检测
- 与 `circdna_tgs_clean.csv`、`circdnalr_{species}_long_read.csv` 一致

### 4.3 路由模型（v4.7，§1.1）

| 字段         | 取值                                                                     | 职责                                                                                |
| ------------ | ------------------------------------------------------------------------ | ----------------------------------------------------------------------------------- |
| `platform`   | `illumina` / `pacbio` / `ont`                                            | 预处理链（按行，支持混表）                                                          |
| `assay`      | `wgs` / `circleseq` / `rca` / `ciderseq` / `enriched`                    | 实验类型 → 引擎集合                                                                 |
| `datatype`   | `gdna` / `eccdna`                                                        | 背景 vs 检测                                                                        |
| `pair`       | 字符串，可空                                                             | 研究分组键（clu join 仅 `illumina×circleseq×eccdna` ⋈ `illumina×wgs×gdna` 同 pair） |
| `concatemer` | `true` / `false`                                                         | 长读 RCA 是否跑 TideHunter（`false`=线性化/T7，跳过）                               |
| `read_type`  | `hifi` / `clr` / `ont` / `pe` / `se`                                     | minimap2 preset（`clr` 必须显式写）                                                 |
| `enrichment` | `none` / `circleseq` / `mobilome` / `ciderseq` / `t7` / `exov` / `other` | 注释，不分流                                                                        |

- 短读缺省：`datatype=gdna→assay=wgs`，`eccdna→circleseq`；Illumina+RCA 必须写成 `circleseq`
- 长读缺省：`gdna→wgs`；eccdna 必须显式（或 `--assay` 回填），禁止猜 ciderseq
- `mode` 不再作为总开关；仅全表同一角色且行上未写 datatype 时回填（`reference→gdna`）

### 4.4 关键 Samplesheet 文件

| File                                | Purpose              | Notes                                                   |
| ----------------------------------- | -------------------- | ------------------------------------------------------- |
| `circdna_ngs_clean.csv`             | **Master NGS input** | 12 species, `data_type` column                          |
| `circdna_tgs_clean.csv`             | **Master TGS input** | Long-read FASTQ                                         |
| `circdna_{species}_eccDNA.csv`      | Per-species NGS      | Auto-generated by `samplesheets/update_samplesheets.py` |
| `circdnalr_{species}_long_read.csv` | Per-species TGS      | Hand-maintained                                         |
| `SraRunInfo_eccDNA_all2.csv`        | SRA metadata         | TaxID, ScientificName, download_path                    |
| `test_local_eccdna.csv`             | Local test           | Small test dataset                                      |

### 4.4 Data Consistency Rules

1. **Path → Species Consistency**: File paths in the samplesheet MUST match actual filesystem locations

   - **Always verify** with `SraRunInfo_eccDNA_all2.csv` (check `TaxID` and `ScientificName`)
   - Common issue: fastq files downloaded to wrong species directory → fix by moving files, not by changing species in samplesheet

2. **No Duplicates**: `circdna_ngs_clean.csv` must not have duplicate sample names

3. **Master File is Source of Truth**: Species-specific files are auto-generated from master files — never edit them directly

---

## 5. Change Management Rules

### 5.1 When Changing Samplesheet Data

**Files to update (MANDATORY)**:

1. `samplesheets/circdna_ngs_clean.csv` (master NGS)
2. `samplesheets/circdna_tgs_clean.csv` (master TGS, if applicable)
3. Run `samplesheets/update_samplesheets.py` to regenerate all per-species files

**Files to verify**:

- `bin/check_samplesheet.py` — ensure column compatibility (supports `data_type`/`datatype` columns)
- `samplesheets/update_samplesheets.py` — ensure it handles new columns correctly

### 5.2 When Changing Pipeline Parameters

**Files to update (MANDATORY)**:

1. `nextflow.config` (add/modify params)
2. `nextflow_schema.json` (run `nf-core schema build` to regenerate)
3. `assets/schema_input.json` (if input schema changes)
4. `conf/modules.config` (if new process needs resource/publishDir config)
5. `conf/test_local.config` (update test params)
6. 相关 test profile（如 `conf/test_local_gdna.config` / `conf/test_nanopore_lr.config` 等）

### 5.3 When Adding/Modifying Modules

**Files to update (MANDATORY)**:

1. `modules/local/{module}/main.nf` (module logic)
2. `modules/local/{module}/meta.yml` (module metadata)
3. `modules/local/{module}/environment.yml` (conda env, if new dependencies)
4. `modules.json` (module manifest — add entry)
5. `conf/modules.config` (add resource/publishDir config for new processes)
6. `workflows/circdna.nf` (include new module in appropriate subworkflow)

**Files to update (if test profiles changed)**:

- `conf/test_local.config`
- `conf/test.config`

### 5.4 When Adding/Modifying Subworkflows

**Files to update (MANDATORY)**:

1. `subworkflows/local/{subworkflow}/main.nf` (subworkflow logic)
2. `subworkflows/local/{subworkflow}/meta.yml` (metadata)
3. `workflows/circdna.nf` (include + wire into mode selection)
4. `conf/modules.config` (resource overrides for new processes)

### 5.5 When Adding New Pipeline Mode / Identifier

**Files to update (MANDATORY)**:

1. `nextflow.config` (add mode to `mode` param options)
2. `workflows/circdna.nf` (add mode to `valid_modes`, create branch; or add identifier to `circle_identifier`/`long_read_identifier` 解析)
3. `conf/test_local.config` (add test config for new mode)
4. `docs/usage.md` (document new mode)
5. `nextflow_schema.json` (regenerate)

### 5.6 When Changing BAM Processing Logic

**Files to update (MANDATORY)**:

1. `subworkflows/local/bam_preprocessing/main.nf`
2. `subworkflows/local/eccdna_mode/main.nf` (uses BAM_PREPROCESSING)
3. `subworkflows/local/reference_mode/main.nf` (uses BAM_PREPROCESSING)
4. `subworkflows/local/circle_finder_pipeline/main.nf` (legacy path)
5. `subworkflows/local/circle_map_pipeline/main.nf` (legacy path)
6. `conf/modules.config` (update resource/publishDir for affected processes)

### 5.7 When Changing the Check Script

**Files to update (MANDATORY)**:

1. `bin/check_samplesheet.py` — update validation logic
2. `subworkflows/local/input_check/main.nf` — verify compatibility
3. `samplesheets/test_local_eccdna.csv` — ensure test samplesheet passes validation
4. Run: `python bin/check_samplesheet.py samplesheets/test_local_eccdna.csv /tmp/test_out.csv FASTQ`

---

## 6. Mode & Identifier Mapping

### Mode Selection (in `workflows/circdna.nf`)

```groovy
params.mode == 'reference'  → REFERENCE_MODE
params.mode == 'eccdna'     → ECCDNA_MODE (+ optional legacy/blackbox/slim branches)
```

- `valid_modes = ['reference', 'eccdna']`

### Legacy Circle Identifiers (only when `mode='eccdna'` + `circle_identifier` is set)

| Identifier           | Module                        | Description                 |
| -------------------- | ----------------------------- | --------------------------- |
| `circexplorer2`      | CIRCEXPLORER2_PARSE           | CIRCexplorer2 parse         |
| `circle_map_realign` | CIRCLE_MAP_PIPELINE (realign) | Circle-Map realignment      |
| `circle_map_repeats` | CIRCLE_MAP_PIPELINE (repeats) | Circle-Map repeat analysis  |
| `circle_finder`      | CIRCLE_FINDER_PIPELINE        | Circle-Finder detection     |
| `ampliconarchitect`  | AMPLICONARCHITECT_PIPELINE    | AmpliconArchitect + CNVkit  |
| `unicycler`          | UNICYCLER_PIPELINE            | Unicycler circular assembly |

### Blackbox Identifiers (short-read, `circle_identifier`)

| Identifier           | Module                        |
| -------------------- | ----------------------------- |
| `eccsplorer`         | ECCSPLORER_PIPELINE           |
| `ecc_finder_map_sr`  | ECC_FINDER_PIPELINE (map-sr)  |
| `ecc_finder_asm_sr`  | ECC_FINDER_PIPELINE (asm-sr)  |
| `ecc_finder_map_ont` | ECC_FINDER_PIPELINE (map-ont) |
| `ecc_finder_asm_ont` | ECC_FINDER_PIPELINE (asm-ont) |

### Slim Identifiers (原子化链，`circle_identifier` 含 `_slim` 触发)

| Identifier                                            | Subworkflow              |
| ----------------------------------------------------- | ------------------------ |
| `eccsplorer_map_slim`                                 | ECCSPLORER_SLIM_PIPELINE |
| `eccsplorer_clu_slim`                                 | ECCSPLORER_CLU_SLIM      |
| `eccsplorer_all_slim`                                 | ECCSPLORER_ALL_SLIM      |
| `ecc_finder_map_sr_slim` / `ecc_finder_asm_sr_slim`   | ECC_FINDER_SLIM_PIPELINE |
| `ecc_finder_map_ont_slim` / `ecc_finder_asm_ont_slim` | ECC_FINDER_ONT_SLIM      |

### Long-read Identifiers (only when `protocol=pacbio|ont` + `long_read_identifier`，默认 `cresil,fled,flye,eccfinder`)

| Identifier     | Subworkflow                        |
| -------------- | ---------------------------------- |
| `cresil`       | CRESIL_PIPELINE                    |
| `fled`         | FLED_PIPELINE                      |
| `flye`         | FLYE_PIPELINE                      |
| `eccfinder`    | ECC_FINDER_PIPELINE (长读 map/asm) |
| `circleseeker` | CIRCLESEEKER_PIPELINE              |

### Default Mode Behavior (no `circle_identifier`)

| Mode        | Subworkflows Run                          |
| ----------- | ----------------------------------------- |
| `reference` | BAM_PREPROCESSING → MOSDEPTH              |
| `eccdna`    | BAM_PREPROCESSING → MOSDEPTH + CIRCLE_MAP |

---

## 7. Key Conventions

### 7.1 Naming Conventions

- **Process names**: UPPER_CASE (e.g., `BWA_MEM`, `TRIMGALORE`, `CIRCLEFINDER`)
- **Channel names**: `ch_` prefix (e.g., `ch_fastq`, `ch_bam_sorted`)
- **Subworkflow names**: UPPER_CASE (e.g., `BAM_PREPROCESSING`, `ECCDNA_MODE`)
- **Custom modules**: snake_case in path, UPPER_CASE in include
- **Species files**: `circdna_{Species}_eccDNA.csv` (underscore-separated species name)

### 7.2 Resource Labels (from `conf/base.config`)

| Label                  | CPUs     | Memory     | Time     |
| ---------------------- | -------- | ---------- | -------- |
| `process_single`       | 1        | 6 GB       | 4h       |
| `process_low`          | 2        | 12 GB      | 4h       |
| `process_medium`       | 6        | 36 GB      | 8h       |
| `process_high`         | 12       | 72 GB      | 16h      |
| `process_long`         | —        | —          | 20h      |
| `process_maximum_time` | —        | —          | max_time |
| `process_high_memory`  | —        | 200 GB     | —        |
| `process_max`          | max_cpus | max_memory | max_time |

**长读/检测专用资源覆盖**（`conf/base.config` `withName`）：`PBCCS`(12c/48GB)、`CRESIL_IDENTIFY`(12c/64GB/48h)、`FLYE`(12c/100GB/96h)、`ECC_FINDER_ASM_*`(12c/64GB/48h) 等。

### 7.3 Samplesheet Column Names

- **Production files** (`circdna_ngs_clean.csv`): `data_type` (camelCase)
- **`check_samplesheet.py`** validates: accepts both `data_type` and `datatype` column names, values normalized to lowercase (`eccdna`/`gdna`)
- **Both conventions work**: `check_samplesheet.py` reads by column name, and `data_type` is accepted as a valid column name via `OPTIONAL_FIELDS`

### 7.4 Data Type Values

| Column Name | Valid Values     | Usage                                                                           |
| ----------- | ---------------- | ------------------------------------------------------------------------------- |
| `data_type` | `eccDNA`, `gDNA` | Production files, filter in `workflows/circdna.nf`, ECCsplorer control analysis |
| `datatype`  | `eccdna`, `gdna` | Test files, validation in `check_samplesheet.py` (normalized to lowercase)      |

### 7.5 Pipeline Modes

- `reference` → gDNA/WGS only (mosdepth)
- `eccdna` → eccDNA detection (Circle-Map + ECCsplorer + ecc_finder 等可选链)

### 7.6 Circle Identifier

- Only active when `mode='eccdna'` AND `circle_identifier` is explicitly set
- Values（legacy）: `circexplorer2`, `circle_map_realign`, `circle_map_repeats`, `circle_finder`, `ampliconarchitect`, `unicycler`
- Values（blackbox）: `eccsplorer`, `ecc_finder_map_sr`, `ecc_finder_asm_sr`, `ecc_finder_map_ont`, `ecc_finder_asm_ont`
- Values（slim）: `eccsplorer_map_slim`, `eccsplorer_clu_slim`, `eccsplorer_all_slim`, `ecc_finder_map_sr_slim`, `ecc_finder_asm_sr_slim`, `ecc_finder_map_ont_slim`, `ecc_finder_asm_ont_slim`
- Default path (no `circle_identifier`): uses ECCDNA_MODE subworkflow (BAM prep + mosdepth + Circle-Map)

### 7.7 Long-read Identifier

- Only active when `protocol='pacbio'|'ont'`，由 `--long_read_identifier` 选择
- Values: `cresil`, `fled`, `flye`, `eccfinder`, `circleseeker`
- 默认值：`cresil,fled,flye,eccfinder`

### 7.8 ECC_SCORE (moved to eccdna.smk)

- ECC_SCORE scoring was moved to `eccdna.smk` (independent post-processing layer) in v4.0.0
- `eccdna.smk` consumes circdna.nf detection outputs (mosdepth bed, eccsplorer bed, circle_map bed)
- Data flow: `circdna.nf (detection) → eccdna.smk (scoring)`

### 7.9 Platform Values

- `illumina`: Short-read sequencing (NGS)
- `pacbio`: Long-read sequencing (TGS)
- `ont`: Oxford Nanopore

### 7.10 Protocol Values

- `short_read`: Short-read sequencing
- `long_read`: Long-read sequencing

---

## 8. Validation & Testing Commands

### 8.1 Schema Validation

```bash
# Regenerate schema after param changes
conda activate nextflow
nf-core schema build
```

### 8.2 Samplesheet Validation

```bash
# Validate a samplesheet
conda activate nextflow
python bin/check_samplesheet.py samplesheets/test_local_eccdna.csv /tmp/test_out.csv FASTQ

# Validate with data_type column
python bin/check_samplesheet.py samplesheets/circdna_ngs_clean.csv /tmp/test_out.csv FASTQ
```

### 8.3 Regenerate Species Files

```bash
# Auto-generate per-species samplesheets
python samplesheets/update_samplesheets.py
```

### 8.4 Local Test Run

```bash
# Run local test
conda activate nextflow
nextflow run main.nf -profile test_local

# Run gDNA/reference mode test
nextflow run main.nf -profile test_local_gdna

# Run long-read tests
nextflow run main.nf -profile test_nanopore_lr
nextflow run main.nf -profile test_pacbio_lr
```

### 8.5 Full Pipeline Run (Example)

```bash
# eccDNA mode with specific circle identifiers
nextflow run main.nf \
    -profile docker \
    --input samplesheets/circdna_Oryza_sativa_eccDNA.csv \
    --input_format FASTQ \
    --genome Oryza_sativa \
    --mode eccdna \
    --circle_identifier circexplorer2,circle_map_realign,circle_map_repeats \
    --outdir results

# Long-read (ONT) mode with cresil + eccfinder
nextflow run main.nf \
    -profile docker \
    --input samplesheets/circdnalr_Oryza_sativa_long_read.csv \
    --input_format FASTQ \
    --genome Oryza_sativa \
    --protocol ont \
    --long_read_identifier cresil,eccfinder \
    --outdir results
```

---

## 9. Common Issues & Fixes

### 9.1 File Not Found in Species Directory

**Symptom**: `ERROR: /path/to/{species}/{sample}.fastq.gz not found`
**Root Cause**: Fastq files downloaded to wrong species directory
**Fix**:

1. Check `SraRunInfo_eccDNA_all2.csv` for correct `ScientificName` and `TaxID`
2. Move files to correct directory
3. Verify and update samplesheet paths

### 9.2 Duplicate Samples

**Symptom**: `Samplesheet contains duplicate rows!`
**Root Cause**: Same sample listed multiple times in same samplesheet
**Fix**: Check `circdna_ngs_clean.csv` for duplicates, remove them

### 9.3 Column Name Mismatch

**Symptom**: Validation fails on `data_type` vs `datatype`
**Root Cause**: Production files use `data_type`, test files use `datatype`
**Fix**: `check_samplesheet.py` supports both — ensure column name matches expected convention

### 9.4 Module Not Found

**Symptom**: `Include of module failed`
**Root Cause**: Module not properly vendored or path incorrect
**Fix**: Check `modules.json`, verify file paths, ensure `meta.yml` exists

### 9.5 大基因组 SAMtools Index 失败

**Symptom**: `samtools index` 报 `Numerical result out of range`
**Root Cause**: 参考序列超过 BAI 上限（约 512 Mb/染色体）
**Fix**: 附加 `-c conf/large_genome.config`（正则 `.*SAMTOOLS_INDEX.*` 匹配全部 SAMTOOLS_INDEX 实例，启用 CSI 索引）

---

## 10. 关联项目与共享资源

本流程与以下项目/工作流存在数据或资源共享关系：

| 项目           | 路径                                                                                                               | 关联关系                                                  |
| -------------- | ------------------------------------------------------------------------------------------------------------------ | --------------------------------------------------------- |
| **circdna.nf** | 本地：`/Users/siyangming/nextflow_nf_core/circdna.nf` · 服务器：`/data1/users/siyangming/PlanteccDNADB/circdna.nf` | eccDNA 检测（本流程）                                     |
| **circrna.nf** | `/Users/siyangming/nextflow_nf_core/circrna.nf`                                                                    | circRNA 检测；共享 SRA 元数据与部分样本表                 |
| **eccdna.smk** | `/Users/siyangming/nextflow_nf_core/eccdna.smk`                                                                    | 独立后处理层：消费本流程检测产物并产出 ECC_SCORE 评分 BED |

**Shared Resources**:

- `samplesheets/SraRunInfo_eccDNA_all2.csv` — SRA 元数据（TaxID、ScientificName），用于物种目录核对
- `samplesheets/circrna_*.csv` — circRNA 样本表（与 circrna.nf 共享）

---

## 11. Git Branch Strategy

| Branch           | Purpose                            |
| ---------------- | ---------------------------------- |
| `master`         | Production-ready code, stable      |
| `circdnalr`      | Long-read (TGS) development branch |
| `ECCsplorer`     | ECCsplorer 集成开发分支            |
| `circleseeker`   | CircleSeeker 长读检测开发分支      |
| Feature branches | Feature-specific development       |

**Sync Rules**:

- Samplesheet changes: Copy files directly (avoid cherry-pick due to frequent conflicts)
- Only sync `samplesheets/` when explicitly requested
- Always verify with `SraRunInfo_eccDNA_all2.csv` before syncing species classifications

---

_Last updated: 2026-08-23_
_Generated based on actual pipeline structure (circdnalr 分支, v4.2.0) and verified conventions_
