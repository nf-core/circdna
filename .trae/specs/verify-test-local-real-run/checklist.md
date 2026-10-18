# Checklist

## 阶段一：eccdna 模式运行

- [ ] `nextflow run main.nf -profile test_local,docker --outdir /tmp/circdna_test_local_eccdna` 执行完成
- [ ] 流程 exit code = 0（或记录失败原因）
- [ ] 记录总任务数和执行时间
- [ ] 无 failed tasks（或记录失败任务详情）

## 阶段二：eccdna 模式产出验证

### FASTQ 预处理

- [ ] `trimgalore/` 目录存在 trimmed FASTQ 文件
- [ ] `fastqc/` 目录存在 QC 报告

### 比对与排序

- [ ] `bwa/` 目录存在 sorted BAM 文件（非空）
- [ ] `bwa/` 目录存在 BAI 索引文件

### 去重

- [ ] `markduplicates/` 目录存在 dedup BAM 文件
- [ ] `markduplicates/` 目录存在 metrics 文件

### 深度计算

- [ ] `mosdepth/` 目录存在区域深度文件（`*.bed.gz` 或 `*.regions.bed`）

### ECCSPLORER 检测

- [ ] `eccsplorer/` 目录存在 `*_candidates.bed` 文件（非 0 字节占位，或确认真实运行后为空）
- [ ] `eccsplorer/` 目录存在 `*_junction_reads.txt` 文件
- [ ] ECCSPLORER 的 `versions.yml` 包含真实版本号（eccsplorer: 2022.01.1.1, samtools: 真实版本）

### Circle-Map

- [ ] `circlemap_realign/` 目录存在 realign 结果文件

### 汇总报告

- [ ] `multiqc/multiqc_report.html` 文件存在
- [ ] MULTIQC 报告包含各模块统计信息

## 阶段三：reference 模式运行

- [ ] `nextflow run main.nf -profile test_local,docker --mode reference --outdir /tmp/circdna_test_local_reference` 执行完成
- [ ] 流程 exit code = 0（或记录失败原因）
- [ ] 记录总任务数和执行时间
- [ ] 无 failed tasks（或记录失败任务详情）

## 阶段四：reference 模式产出验证

### FASTQ 预处理

- [ ] `trimgalore/` 目录存在 trimmed FASTQ 文件
- [ ] `fastqc/` 目录存在 QC 报告

### 比对与排序

- [ ] `bwa/` 目录存在 sorted BAM 文件（非空）
- [ ] `bwa/` 目录存在 BAI 索引文件

### 去重

- [ ] `markduplicates/` 目录存在 dedup BAM 文件
- [ ] `markduplicates/` 目录存在 metrics 文件

### CIRCEXPLORER2 解析

- [ ] `circexplorer2/` 目录存在解析结果文件

### Circle-Finder 检测

- [ ] `circlefinder/` 目录存在检测结果文件（`microDNA-JT.txt` 等）

### AmpliconArchitect

- [ ] `ampliconarchitect/` 目录存在结果文件（或确认因 testdata 限制跳过）

### UNICYCLER 组装

- [ ] `unicycler/` 目录存在组装结果文件（或确认因 testdata 限制跳过）

### 汇总报告

- [ ] `multiqc/multiqc_report.html` 文件存在
- [ ] MULTIQC 报告包含各模块统计信息

## 关键验证点

### eccdna 模式

- [ ] ECCSPLORER 真实执行（非 stub），无 `ambiguous option: -c` 错误
- [ ] ECCSPLORER 调用 `python /opt/eccsplorer/ECCsplorer.py` 成功执行
- [ ] Docker 镜像 `quay.io/bioinfortools/eccsplorer:2022.01.1.1` 在真实模式下工作正常
- [ ] BAM_PREPROCESSING 子工作流正常工作（FASTQ 输入模式）

### reference 模式

- [ ] CIRCEXPLORER2_PARSE 真实执行
- [ ] Circle-Finder 真实执行（产出 `microDNA-JT.txt`，可能为空但非占位）
- [ ] AmpliconArchitect/UNICYCLER 按预期执行或合理跳过

### 通用

- [ ] 两种模式端到端可用，产出真实分析结果
- [ ] MULTIQC 汇总报告正常生成
