# Checklist

- [x] master 提交已并入 ECCsplorer 分支（git merge 完成，冲突文件以 master 为准，commit 444ce77）
- [x] merge 后 slim 独有文件保留（stash pop 成功，未丢失 slim 模块链；commit 691fbd6）
- [x] `testdatasets/ont/ont_eccdna_smoke/regular.fastq.gz` 存在且可读（1500/7500 reads）
- [x] `testdatasets/pacbio/pacbio_eccdna_smoke/regular.fastq.gz` 存在且可读（1508/2512 reads）
- [x] `samplesheets/test_ont.csv` 就位，fastq 路径为**相对路径**，含 `platform=ont`/`protocol=long_read`
- [x] `samplesheets/test_pacbio.csv` 就位，fastq 路径为**相对路径**，含 `platform=pacbio`/`protocol=long_read`
- [x] `samplesheets/test_eccsplorer_pair.csv` 已生成（6 样本，含 `datatype`+`pair` 列），check_samplesheet 解析通过（EXIT=0，pair 透传）
- [x] master 版 `conf/server.config` 就位（`fasta_base_path=/data1/users/siyangming/PublicDB/reference` + `genomes` 块 17 物种）
- [x] 服务器运行命令模板已记录到文档（`-profile server` + `--genome <species>` + `--input test_*.csv`）
- [x] ECCsplorer 分支冻结修改保留价值评估报告产出（保留/合并/丢弃清单），已呈现用户确认
- [x] 有价值修改已提交到 ECCsplorer 分支（commit 691fbd6，181 文件）、无价值项已清理
- [x] 提交后主流水线 slim 标识符 stub 编译通过（回归无破坏，编译无错误）
- [x] 交付物汇总（test\_ 文件清单 + 运行命令模板 + 相对路径说明）已写入文档
