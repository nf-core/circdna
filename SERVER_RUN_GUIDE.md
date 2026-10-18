# circdna.nf 服务器运行指南

## Screen 会话管理

> **首次使用 screen 只需配置一次**，后续通过 `screen -r eccdna` 恢复

```bash
# 创建 screen 会话 (首次)
screen -S eccdna

# 恢复 screen 会话 (每次登录服务器后)
screen -r eccdna

# 列出所有 screen 会话
screen -ls

# 从 screen 会话分离 (保持运行)
# 快捷键: Ctrl+A+D

# 重新连接到已分离的会话
screen -r eccdna
```

## 环境先决条件

> **每次进入 screen 会话后必须执行**，否则会出现 `nextflow: 未找到命令` 错误

```bash
# 1. 激活 conda 环境 (nextflow 命令仅在此环境中可用)
conda activate nextflow

# 2. 进入父目录 (所有操作在此目录下执行)
cd /data1/users/siyangming/PlanteccDNADB/

# 3. 验证环境
which nextflow    # 应输出 conda 环境中的 nextflow 路径
nextflow -version # 应显示 Nextflow 版本

# 4. 确认目录结构
ls
# 应包含: circdna.nf/  eccDNA/  eccDNA_results/
```

### 首次完整设置（仅需一次）

```bash
# Step 1: 创建并进入 screen
screen -S eccdna

# Step 2: 激活环境并进入目录
conda activate nextflow
cd /data1/users/siyangming/PlanteccDNADB/

# Step 3: 保存环境变量到 screen 会话 (可选，避免每次重复)
# conda env config set auto_activate_base false
# echo 'conda activate nextflow' >> ~/.bashrc
# echo 'cd /data1/users/siyangming/PlanteccDNADB/' >> ~/.bashrc

# Step 4: 分离会话
# Ctrl+A+D

# Step 5: 验证恢复
screen -r eccdna
# 应该自动在 nextflow 环境和 PlanteccDNADB 目录中
```

## 路径约定

所有命令在 `/data1/users/siyangming/PlanteccDNADB/` 父目录下执行，使用 `circdna.nf/` 前缀引用 pipeline 文件：

```bash
# 父目录结构
PlanteccDNADB/
├── circdna.nf/          # Pipeline 代码 (nextflow run ./circdna.nf/main.nf)
│   ├── main.nf
│   ├── nextflow.config
│   ├── samplesheets/    # circdna.nf/samplesheets/*.csv
│   ├── conf/            # circdna.nf/conf/server.config
│   └── ...
├── eccDNA/              # 测序数据
└── eccDNA_results/      # 输出结果
```

## 恢复运行

```bash
# 必须指定具体的 run name，否则 -resume 会尝试恢复最近一次运行
nextflow run ./circdna.nf/main.nf \
    --input circdna.nf/samplesheets/circdna_<species>_eccDNA.csv \
    --genome <species> \
    --outdir eccDNA_results/<species> \
    -profile server \
    -resume <run_name>

# 查看历史 run name
nextflow log -q
# 或查看 .nextflow.log 开头
head -5 circdna.nf/.nextflow.log
```

## 物种对照表（中文 ↔ 拉丁文）

- 大穗看麦娘（Alopecurus_myosuroides）
- 长芒苋（Amaranthus_palmeri）
- 拟南芥（Arabidopsis_thaliana）
- 黄花蒿（Artemisia_annua）
- 甜菜（Beta_vulgaris）
- 日本柳杉（Cryptomeria_japonica）
- 狗牙根（Cynodon_dactylon）
- 胡萝卜（Daucus_carota）
- 向日葵（Helianthus_annuus）
- 黑果枸杞（Lycium_ruthenicum）
- 本氏烟草（Nicotiana_benthamiana）
- 水稻（Oryza_sativa；IRGSP-1.0 默认；Oryza_sativa_Huazhan 仅用于 ecc-\* 样本）
- 番茄（Solanum_lycopersicum）
- 婆罗门参（Tragopogon_porrifolius）
- 普通小麦（Triticum_aestivum）

## 按物种运行命令

> 假设已在 `PlanteccDNADB/` 目录内，且 `conda activate nextflow` 已执行

### 拟南芥 (Arabidopsis_thaliana)

```bash
nextflow run ./circdna.nf/main.nf \
    --input circdna.nf/samplesheets/circdna_Arabidopsis_thaliana_eccDNA.csv \
    --genome Arabidopsis_thaliana \
    --outdir eccDNA_results/Arabidopsis_thaliana \
    -profile server
```

### 水稻 (Oryza_sativa)

`Oryza_sativa`（IRGSP-1.0）为水稻默认基因组；`Oryza_sativa_Huazhan` 目前只用于 `ecc-*` 开头的 8 个样本，必须搭配下方 Huazhan 专用表使用，不用于其他水稻样本或长读样本。

```bash
# 常规水稻运行（IRGSP-1.0；若只跑非 ecc-* 样本，先过滤表内 ecc-* 行）
nextflow run ./circdna.nf/main.nf \
    --input circdna.nf/samplesheets/circdna_Oryza_sativa_eccDNA.csv \
    --genome Oryza_sativa \
    --outdir eccDNA_results/Oryza_sativa \
    -profile server

# ecc-* 样本专用（Huazhan）
nextflow run ./circdna.nf/main.nf \
    --input circdna.nf/samplesheets/circdna_Oryza_sativa_Huazhan_eccDNA.csv \
    --genome Oryza_sativa_Huazhan \
    --outdir eccDNA_results/Oryza_sativa_Huazhan \
    -profile server
```

### 番茄 (Solanum_lycopersicum)

```bash
nextflow run ./circdna.nf/main.nf \
    --input circdna.nf/samplesheets/circdna_Solanum_lycopersicum_eccDNA.csv \
    --genome Solanum_lycopersicum \
    --outdir eccDNA_results/Solanum_lycopersicum \
    -profile server
```

### 胡萝卜 (Daucus_carota)

```bash
nextflow run ./circdna.nf/main.nf \
    --input circdna.nf/samplesheets/circdna_Daucus_carota_eccDNA.csv \
    --genome Daucus_carota \
    --outdir eccDNA_results/Daucus_carota \
    -profile server
```

### 向日葵 (Helianthus_annuus)

```bash
nextflow run ./circdna.nf/main.nf \
    --input circdna.nf/samplesheets/circdna_Helianthus_annuus_eccDNA.csv \
    --genome Helianthus_annuus \
    --outdir eccDNA_results/Helianthus_annuus \
    -profile server
```

### 本氏烟草 (Nicotiana_benthamiana)

```bash
nextflow run ./circdna.nf/main.nf \
    --input circdna.nf/samplesheets/circdna_Nicotiana_benthamiana_eccDNA.csv \
    --genome Nicotiana_benthamiana \
    --outdir eccDNA_results/Nicotiana_benthamiana \
    -profile server
```

### 甜菜 (Beta_vulgaris)

```bash
nextflow run ./circdna.nf/main.nf \
    --input circdna.nf/samplesheets/circdna_Beta_vulgaris_eccDNA.csv \
    --genome Beta_vulgaris \
    --outdir eccDNA_results/Beta_vulgaris \
    -profile server
```

### 黑果枸杞 (Lycium_ruthenicum)

```bash
nextflow run ./circdna.nf/main.nf \
    --input circdna.nf/samplesheets/circdna_Lycium_ruthenicum_eccDNA.csv \
    --genome Lycium_ruthenicum \
    --outdir eccDNA_results/Lycium_ruthenicum \
    -profile server
```

### 婆罗门参 (Tragopogon_porrifolius) — hap1/hap2，hap1 为大基因组

> 该物种有两个单倍型：`Tragopogon_porrifolius_hap1` 和 `Tragopogon_porrifolius_hap2`（后者将 `--genome` 换为 `Tragopogon_porrifolius_hap2` 即可）
> hap1 参考序列长度超过 BAI 索引上限（约 512 Mb/染色体），**必须**附加 `-c circdna.nf/conf/large_genome.config`（启用 CSI 索引）；hap2 无需附加

```bash
# hap1 (大基因组，需 CSI 索引)
nextflow run ./circdna.nf/main.nf \
    --input circdna.nf/samplesheets/circdna_Tragopogon_porrifolius_eccDNA.csv \
    --genome Tragopogon_porrifolius_hap1 \
    --outdir eccDNA_results/Tragopogon_porrifolius \
    -profile server \
    -c circdna.nf/conf/large_genome.config

# hap2
nextflow run ./circdna.nf/main.nf \
    --input circdna.nf/samplesheets/circdna_Tragopogon_porrifolius_eccDNA.csv \
    --genome Tragopogon_porrifolius_hap2 \
    --outdir eccDNA_results/Tragopogon_porrifolius \
    -profile server \
    -c circdna.nf/conf/large_genome.config
```

### 狗牙根 (Cynodon_dactylon)

```bash
nextflow run ./circdna.nf/main.nf \
    --input circdna.nf/samplesheets/circdna_Cynodon_dactylon_eccDNA.csv \
    --genome Cynodon_dactylon \
    --outdir eccDNA_results/Cynodon_dactylon \
    -profile server
```

### 青蒿 (Artemisia_annua)

```bash
nextflow run ./circdna.nf/main.nf \
    --input circdna.nf/samplesheets/circdna_Artemisia_annua_eccDNA.csv \
    --genome Artemisia_annua \
    --outdir eccDNA_results/Artemisia_annua \
    -profile server
```

### 苋 (Amaranthus_palmeri) — hap1/hap2

> 该物种有两个单倍型：`Amaranthus_palmeri_hap1` 和 `Amaranthus_palmeri_hap2`（后者将 `--genome` 换为 `Amaranthus_palmeri_hap2` 即可）

```bash
nextflow run ./circdna.nf/main.nf \
    --input circdna.nf/samplesheets/circdna_Amaranthus_palmeri_eccDNA.csv \
    --genome Amaranthus_palmeri_hap1 \
    --outdir eccDNA_results/Amaranthus_palmeri \
    -profile server
```

### 黑麦草 (Alopecurus_myosuroides) — 大基因组

```bash
nextflow run ./circdna.nf/main.nf \
    --input circdna.nf/samplesheets/circdna_Alopecurus_myosuroides_eccDNA.csv \
    --genome Alopecurus_myosuroides \
    --outdir eccDNA_results/Alopecurus_myosuroides \
    -profile server \
    -c circdna.nf/conf/large_genome.config
```

### 日本柳杉 (Cryptomeria_japonica) — 大基因组

```bash
nextflow run ./circdna.nf/main.nf \
    --input circdna.nf/samplesheets/circdna_Cryptomeria_japonica_eccDNA.csv \
    --genome Cryptomeria_japonica \
    --outdir eccDNA_results/Cryptomeria_japonica \
    -profile server \
    -c circdna.nf/conf/large_genome.config
```

### 小麦 (Triticum_aestivum) — 大基因组

```bash
nextflow run ./circdna.nf/main.nf \
    --input circdna.nf/samplesheets/circdna_Triticum_aestivum_eccDNA.csv \
    --genome Triticum_aestivum \
    --outdir eccDNA_results/Triticum_aestivum \
    -profile server \
    -c circdna.nf/conf/large_genome.config
```

## 常用操作

### Screen 会话

```bash
# 查看所有 screen 会话
screen -ls

# 恢复会话
screen -r eccdna

# 新建会话
screen -S <name>

# 从当前会话分离 (保持进程运行)
# Ctrl+A+D

# 终止会话 (在会话内执行)
exit

# 远程 detached 会话 (会话在后台)
screen -d -r eccdna
```

### Pipeline 操作

```bash
# 查看运行日志
tail -f circdna.nf/.nextflow.log

# 查看最近 run name
head -3 circdna.nf/.nextflow.log

# 查看运行进度
nextflow log <run_name> -f

# 清理指定 run 之前的 work 缓存
nextflow clean -f -before <run_name>

# 只查看会被清理的目录，不实际删除
nextflow clean -n -f -before <run_name>

# 快速重启某个物种
nextflow run ./circdna.nf/main.nf \
    --input circdna.nf/samplesheets/circdna_<species>_eccDNA.csv \
    --genome <species> \
    --outdir eccDNA_results/<species> \
    -profile server
```

## 中间产物与 work 清理

`-profile server` 的 `conf/server.config` 已默认开启：

- `save_reference`、`save_trimmed`、`save_merged_fastq`、`save_markduplicates_bam`、`save_sorted_bam`、`save_circle_map_intermediate`、`save_circle_finder_intermediate`、`save_unicycler_intermediate`、`save_long_read_intermediate` 全部为 `true`。
- 没有专门 publish 路径的 process 输出，会复制到 `<outdir>/intermediate/<完整流程路径>/`。
- `cleanup = true`：运行成功后自动删除 `work/` 目录。

因此，正式产出放在 `<outdir>`，`work/` 只作为临时计算目录，不用于长期保存。失败运行不会自动清理；需要手动清理旧缓存时使用 `nextflow clean -f -before <run_name>`，需要删除全部缓存时使用 `nextflow clean -f`。若某个运行仍需 `-resume`，先不要执行清理，或用临时 `cleanup=false` 覆盖 server 配置。

## 大基因组说明

染色体超过 512 Mb 的物种（如小麦、日本柳杉、黑麦草），BAI 索引无法表示超过 2^29 ≈ 536 Mb 的坐标，会导致 `samtools index` 或 `samtools sort --write-index` 报 `Numerical result out of range`。

`large_genome.config` 做两件事：

1. `params.use_csi_index = true` — 子工作流中 `SAMTOOLS_SORT_RE` 等进程的 `samtools sort --write-index` 改用 CSI 格式
2. 精确覆盖全部 4 个 `SAMTOOLS_INDEX` 实例（`SAMTOOLS_INDEX_BAM` / `_FILTERED` / `_RE` / `BAM_MARKDUPLICATES_PICARD:SAMTOOLS_INDEX`）的 `ext.args = '-c'`，使 `samtools index` 输出 CSI

用法：在大基因组物种的命令末尾加 `-c circdna.nf/conf/large_genome.config`，如下方黑麦草/小麦/日本柳杉/婆罗门参 hap1 所示。

## 注意事项

- **必须执行 `conda activate nextflow`**：`nextflow` 命令仅在 `nextflow` conda 环境中可用
- **`-resume` 必须指定 run name**：使用 `-resume`（不带参数）会恢复最近一次运行，可能不是你想要的
- **参考基因组文件**：需已存在于 `/data1/users/siyangming/PublicDB/reference/<species>/` 目录下
- **样本数据**：需存在于 `eccDNA/` 目录
- **大基因组**：染色体超过 512 Mb 的物种需加 `-c circdna.nf/conf/large_genome.config`，详见上方「大基因组说明」
- **`circle_identifier`、`input_format`** 等参数已在 `circdna.nf/conf/server.config` 中配置，无需在命令中指定
- **清理旧 work 目录**：优先使用 `nextflow clean -f -before <run_name>`；如遇 root 权限文件，再退回 Docker 删除

## 常见错误

### `bash: nextflow: 未找到命令`

```bash
# 原因: 未激活 conda 环境
# 解决:
conda activate nextflow
which nextflow  # 验证
```

### `-resume` 恢复了错误的运行

```bash
# 原因: -resume 不带参数时默认恢复最近一次运行
# 解决: 指定 run name
nextflow run ./circdna.nf/main.nf ... -resume <具体的run_name>

# 查看可用的 run names
head -10 circdna.nf/.nextflow.log | grep "Launching\|run name"
```

### `WARN: Access to undefined parameter 'fasta'`

```bash
# 已在 circdna.nf/workflows/circdna.nf#L37 修复:
# 使用 params.containsKey('fasta') 检查而非直接访问
# 如果仍出现警告，确保 --genome 参数正确指定
```

### `Remote resource not found: nextflow-io/circdna.nf`

```bash
# 原因: Nextflow 把 circdna.nf 误解析为 GitHub 远程仓库（org/repo 格式）
# 因为没有路径前缀，被识别为远程仓库名

# 错误命令:
nextflow run ./circdna.nf/main.nf ...

# 解决方案 1: 使用 ./ 前缀（推荐）
nextflow run ./circdna.nf/main.nf \
    --input circdna.nf/samplesheets/circdna_<species>_eccDNA.csv \
    --genome <species> \
    --outdir eccDNA_results/<species> \
    -profile server

# 解决方案 2: 使用绝对路径
nextflow run /data1/users/siyangming/PlanteccDNADB/circdna.nf/main.nf \
    --input circdna.nf/samplesheets/circdna_<species>_eccDNA.csv \
    ...

# 解决方案 3: 进入 circdna.nf 目录内运行
cd circdna.nf/
nextflow run main.nf \
    --input samplesheets/circdna_<species>_eccDNA.csv \
    --genome <species> \
    --outdir ../eccDNA_results/<species> \
    -profile server
```

---

## 三代（TGS）长读数据分析

> 三代数据 = PacBio / ONT 长读。路由完全由 samplesheet 行级字段控制（`platform × assay × datatype × concatemer × read_type`），**无需指定 `--protocol`**；`--genome` 为单基因组参数，因此按物种表（`circdnalr_{species}_long_read.csv`）逐个运行。

### 数据总览（141 个长读样本，`samplesheets/circdna_tgs_clean.csv`）

| 物种                   | 样本数 | 平台   | assay    | datatype | read_type      | 自动路由                    |
| ---------------------- | ------ | ------ | -------- | -------- | -------------- | --------------------------- |
| Oryza_sativa           | 1      | pacbio | wgs      | gdna     | hifi           | LONG_READ_REFERENCE（背景） |
| Triticum_aestivum      | 80     | ont    | rca      | eccdna   | ont            | 检测（默认引擎）            |
| Amaranthus_palmeri     | 9      | pacbio | rca      | eccdna   | clr            | 检测（默认引擎）            |
| Amaranthus_palmeri     | 1      | pacbio | ciderseq | eccdna   | clr            | CIDER-Seq2                  |
| Amaranthus_palmeri     | 6      | pacbio | wgs      | gdna     | hifi           | LONG_READ_REFERENCE（背景） |
| Arabidopsis_thaliana   | 13     | ont    | wgs      | gdna     | ont            | LONG_READ_REFERENCE（背景） |
| Arabidopsis_thaliana   | 15     | ont    | rca      | eccdna   | ont            | 检测（默认引擎）            |
| Alopecurus_myosuroides | 7      | pacbio | rca      | eccdna   | clr(1)/hifi(6) | 检测（大基因组 + CSI）      |
| Solanum_lycopersicum   | 2      | ont    | rca      | eccdna   | ont            | 检测（默认引擎）            |
| Helianthus_annuus      | 1      | ont    | enriched | eccdna   | ont            | 检测（T7 富集）             |
| Helianthus_annuus      | 3      | ont    | wgs      | gdna     | ont            | LONG_READ_REFERENCE（背景） |
| Nicotiana_benthamiana  | 4      | ont    | rca      | eccdna   | ont            | 检测（默认引擎）            |

> 单张表内同时含 `wgs/gdna` 与 `rca/eccdna` 行时（如拟南芥、向日葵、长芒苋），背景行自动走 LONG_READ_REFERENCE、检测行自动走检测引擎，一次 run 并行完成，无需拆表。

### 引擎与参数速查

| 参数                                     | 推荐值                       | 说明                                                       |
| ---------------------------------------- | ---------------------------- | ---------------------------------------------------------- |
| `--long_read_identifier`                 | `cresil,fled,flye,eccfinder` | 默认引擎（PacBio/ONT RCA、PacBio CLR 均适用）              |
| `--long_read_identifier`                 | 默认 + `,circleseeker`       | HiFi + RCA + concatemer=true 样本（黑麦草 6 个 hifi）      |
| `--long_read_identifier`                 | `''`（空）                   | 纯 WGS 背景表，不跑任何检测引擎                            |
| `--long_read_identifier`                 | `ciderseq`                   | CIDER-Seq2（需同时提供 4 个 `ciderseq_*` 参数）            |
| `--eccfinder_mode`                       | `map`（默认）                | `map` / `asm` / `both`；`asm` 会触发组装 + remap，耗时翻倍 |
| `--min_read_support`                     | `2`（默认）                  | 统一 BED 的读支持数过滤阈值                                |
| `-c circdna.nf/conf/large_genome.config` | 小麦 / 黑麦草                | 染色体 >512 Mb，启用 CSI 索引                              |

### 前置检查（每次进入 screen 后）

```bash
conda activate nextflow
cd /data1/users/siyangming/PlanteccDNADB/

# 1. 数据文件（应 ≥ 物种表样本数）
for sp in Oryza_sativa Triticum_aestivum Amaranthus_palmeri Arabidopsis_thaliana \
          Alopecurus_myosuroides Solanum_lycopersicum Helianthus_annuus Nicotiana_benthamiana; do
    echo "$sp: $(ls eccDNA/$sp/*.fastq.gz 2>/dev/null | wc -l) fastq.gz"
done

# 2. 参考基因组（server.config 中 18 个基因组路径已验证全部存在，含 Oryza_sativa_Huazhan）
ls /data1/users/siyangming/PublicDB/reference/*/   # 核对各物种 .fa.gz

# 3. 当前代码版本（应含 v4.7 路由改造）
cd circdna.nf && git log --oneline -1 && cd ..
```

### 运行命令（按物种）

> 冒烟优先：先跑番茄（仅 2 样本）验证链路，再按需全量。所有长读 run 使用 `circdnalr_{species}_long_read.csv`，输出到 `eccDNA_results/{species}_longread/`。

#### ① 冒烟验证 — 番茄（2 样本，最快）

```bash
nextflow run ./circdna.nf/main.nf \
    --input circdna.nf/samplesheets/circdnalr_Solanum_lycopersicum_long_read.csv \
    --genome Solanum_lycopersicum \
    --long_read_identifier cresil,fled,flye,eccfinder \
    --outdir eccDNA_results/Solanum_lycopersicum_longread \
    -profile server
```

#### ② 拟南芥（27 样本：13 WGS 背景 + 15 RCA 检测，双链路并行验证）

```bash
nextflow run ./circdna.nf/main.nf \
    --input circdna.nf/samplesheets/circdnalr_Arabidopsis_thaliana_long_read.csv \
    --genome Arabidopsis_thaliana \
    --long_read_identifier cresil,fled,flye,eccfinder \
    --outdir eccDNA_results/Arabidopsis_thaliana_longread \
    -profile server
```

#### ③ 向日葵（1 T7 富集检测 + 3 WGS 背景）

```bash
nextflow run ./circdna.nf/main.nf \
    --input circdna.nf/samplesheets/circdnalr_Helianthus_annuus_long_read.csv \
    --genome Helianthus_annuus \
    --long_read_identifier cresil,fled,flye,eccfinder \
    --outdir eccDNA_results/Helianthus_annuus_longread \
    -profile server
```

#### ④ 本氏烟草（4 样本 ONT RCA）

```bash
nextflow run ./circdna.nf/main.nf \
    --input circdna.nf/samplesheets/circdnalr_Nicotiana_benthamiana_long_read.csv \
    --genome Nicotiana_benthamiana \
    --long_read_identifier cresil,fled,flye,eccfinder \
    --outdir eccDNA_results/Nicotiana_benthamiana_longread \
    -profile server
```

#### ⑤ 长芒苋（9 RCA CLR + 6 HiFi WGS 背景；含 1 个 CIDER-Seq2 行）

```bash
# RCA 检测 + HiFi WGS 背景（CIDER 行无引擎，仅预处理；如需跑 CIDER 见下方⑨）
nextflow run ./circdna.nf/main.nf \
    --input circdna.nf/samplesheets/circdnalr_Amaranthus_palmeri_long_read.csv \
    --genome Amaranthus_palmeri_hap1 \
    --long_read_identifier cresil,fled,flye,eccfinder \
    --outdir eccDNA_results/Amaranthus_palmeri_longread \
    -profile server
```

#### ⑥ 黑麦草（7 样本，大基因组 + CSI；6 个 HiFi 可加 circleseeker）

```bash
nextflow run ./circdna.nf/main.nf \
    --input circdna.nf/samplesheets/circdnalr_Alopecurus_myosuroides_long_read.csv \
    --genome Alopecurus_myosuroides \
    --long_read_identifier cresil,fled,flye,eccfinder,circleseeker \
    --outdir eccDNA_results/Alopecurus_myosuroides_longread \
    -profile server \
    -c circdna.nf/conf/large_genome.config
```

#### ⑦ 水稻（仅 1 个 HiFi WGS 背景，纯背景无检测引擎）

```bash
nextflow run ./circdna.nf/main.nf \
    --input circdna.nf/samplesheets/circdnalr_Oryza_sativa_long_read.csv \
    --genome Oryza_sativa \
    --long_read_identifier '' \
    --outdir eccDNA_results/Oryza_sativa_longread \
    -profile server
```

> 长读水稻仍使用 `--genome Oryza_sativa`（IRGSP-1.0）。`Oryza_sativa_Huazhan` 目前仅用于 `circdna_Oryza_sativa_Huazhan_eccDNA.csv` 中的短读 `ecc-*` 样本，不用于长读背景或其他样本。

#### ⑧ 小麦（80 样本 ONT RCA，大基因组 + CSI；计算量最大，务必先抽样验证）

```bash
# 第一步：抽取 3 个代表样本验证（含 pair 为空/非空、benchmark 样本各 1）
head -1 circdna.nf/samplesheets/circdnalr_Triticum_aestivum_long_read.csv > /tmp/wheat_smoke.csv
grep -E "^ERR12724336|^ERR12724360|^ERR6326020" circdna.nf/samplesheets/circdnalr_Triticum_aestivum_long_read.csv >> /tmp/wheat_smoke.csv

# 第二步：冒烟跑（验证后删掉 --outdir 或换 run name 全量）
nextflow run ./circdna.nf/main.nf \
    --input /tmp/wheat_smoke.csv \
    --genome Triticum_aestivum \
    --long_read_identifier cresil,fled,flye,eccfinder \
    --outdir eccDNA_results/Triticum_aestivum_longread_smoke \
    -profile server \
    -c circdna.nf/conf/large_genome.config

# 第三步：全量（-resume 必须指定 run name）
nextflow run ./circdna.nf/main.nf \
    --input circdna.nf/samplesheets/circdnalr_Triticum_aestivum_long_read.csv \
    --genome Triticum_aestivum \
    --long_read_identifier cresil,fled,flye,eccfinder \
    --outdir eccDNA_results/Triticum_aestivum_longread \
    -profile server \
    -c circdna.nf/conf/large_genome.config \
    -resume <冒烟_run_name>
```

#### ⑨ CIDER-Seq2（长芒苋 SRR16958693，单独运行）

> 前置准备：CIDER-Seq2 需要 4 个参数指向数据库（服务器尚无生产库，需从官方 CIDER-Seq2 资源准备）：`--ciderseq_config`（JSON）、`--ciderseq_blastdb`、`--ciderseq_align_targets`、`--ciderseq_protein_db`。并准备只含 ciderseq 行的子表：

```bash
head -1 circdna.nf/samplesheets/circdnalr_Amaranthus_palmeri_long_read.csv > circdna.nf/samplesheets/ciderseq_prod.csv
grep "SRR16958693" circdna.nf/samplesheets/circdnalr_Amaranthus_palmeri_long_read.csv >> circdna.nf/samplesheets/ciderseq_prod.csv

nextflow run ./circdna.nf/main.nf \
    --input circdna.nf/samplesheets/ciderseq_prod.csv \
    --genome Amaranthus_palmeri_hap1 \
    --long_read_identifier ciderseq \
    --ciderseq_config /path/to/ciderseq_config.json \
    --ciderseq_blastdb /path/to/ciderseq_blastdb \
    --ciderseq_align_targets /path/to/ciderseq_align_targets \
    --ciderseq_protein_db /path/to/ciderseq_protein_db \
    --outdir eccDNA_results/Amaranthus_palmeri_ciderseq \
    -profile server
```

### 三代运行注意事项

1. **数据放置**：小麦表中 ERR6326020/21 的 FASTQ 位于 `Triticum_aestivum/` 目录，但 metadata 显示其属于拟南芥 ecc_finder benchmark（PRJEB46420）——按"路径决定物种"规则会按小麦处理；如要正确归类，先把文件移到 `Arabidopsis_thaliana/` 并移出小麦表。
2. **引擎资源**：Flye 组装最重（12c/100GB/96h，base.config 已覆盖）；`eccfinder_mode=asm` 会额外触发组装+remap。全量跑 80 样本小麦前先评估队列（executor cpus=96）。
3. **短读参数无影响**：server.config 中 `circle_identifier` 为短读 legacy 值，长读 run 下短读通道为空（数据触发门控已验证），不会执行短读进程。
4. **`-resume` 必须指定 run name**（见上文「常见错误」）。
5. **断点续跑**：长读大表建议拆成多个独立 run（按物种），避免单次排队过久。
