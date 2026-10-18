# slim 链服务器测试运行指南

参考基因组**不下载**：通过 `conf/server.config` 的 `genomes` 块（`fasta_base_path=/data1/users/siyangming/PublicDB/reference`）经 `--genome <species>` 提供。

## 测试文件（相对路径，放仓库根目录即可运行）

| 文件                                    | 内容                                                                        | 用途                        |
| --------------------------------------- | --------------------------------------------------------------------------- | --------------------------- |
| `samplesheets/test_ont.csv`             | ONT 长读 smoke(1500)/regular(7500) reads，`platform=ont,protocol=long_read` | map-ont / asm-ont           |
| `samplesheets/test_pacbio.csv`          | PacBio HiFi smoke(1508)/regular(2512) reads，`platform=pacbio`              | map-ont / asm-ont           |
| `samplesheets/test_eccsplorer_pair.csv` | 短读 6 样本（gdna/circdna × 3），含 `datatype`+`pair` 列                    | map(clu) control 配对 / all |

## 运行命令模板

```bash
cd /path/to/circdna.nf   # 服务器仓库根目录

# 1) ONT 链（map-ont + asm-ont），Arabidopsis 参考走 server.config
nextflow run main.nf -profile server \
  --genome Arabidopsis_thaliana \
  --input samplesheets/test_ont.csv \
  --circle_identifier 'ecc_finder_map_ont_slim,ecc_finder_asm_ont_slim' \
  --outdir results/ont_slim

# 2) PacBio 链（同 ONT 标识符，数据为 HiFi）
nextflow run main.nf -profile server \
  --genome Arabidopsis_thaliana \
  --input samplesheets/test_pacbio.csv \
  --circle_identifier 'ecc_finder_map_ont_slim,ecc_finder_asm_ont_slim' \
  --outdir results/pacbio_slim

# 3) 短读 control 配对（map + clu，按 datatype 分流、pair 配对）
nextflow run main.nf -profile server \
  --genome Oryza_sativa \
  --input samplesheets/test_eccsplorer_pair.csv \
  --circle_identifier 'eccsplorer_map_slim,eccsplorer_clu_slim' \
  --outdir results/pair_slim

# 4) all 模式（map + clu + comparative + 契约导出）
nextflow run main.nf -profile server \
  --genome Oryza_sativa \
  --input samplesheets/test_eccsplorer_pair.csv \
  --circle_identifier 'eccsplorer_all_slim' \
  --outdir results/all_slim

# 5) 短读 ecc_finder map-sr/asm-sr（BAM 由 BAM_PREPROCESSING 生成）
nextflow run main.nf -profile server \
  --genome Oryza_sativa \
  --input samplesheets/test_eccsplorer_pair.csv \
  --circle_identifier 'ecc_finder_map_sr_slim,ecc_finder_asm_sr_slim' \
  --outdir results/sr_slim
```

## 说明

- `test_ont.csv`/`test_pacbio.csv` 的 fastq 路径为相对路径（`testdatasets/ont/...`），仓库置于服务器任意路径即可运行
- 物种选择：ONT/PacBio 数据为 Arabidopsis（SRR24335762）/ rice HiFi，对应 `--genome Arabidopsis_thaliana` / `--genome Oryza_sativa`（server.config 已含 TAIR10 / IRGSP-1.0 链接）
- 数据提取脚本：`testdatasets/extract_test_data.sh`（固定 seed，可复现）
