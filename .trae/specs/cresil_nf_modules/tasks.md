# CReSIL Nextflow 模块 - The Implementation Plan (Decomposed and Prioritized Task List)

## [x] Task 1: 准备测试数据 (testdata)
- **Priority**: high
- **Depends On**: None
- **Description**:
  - 创建 `modules/cresil/testdata/` 目录
  - 复制 example/exp_reads.fastq 到 testdata
  - 生成小的测试参考基因组（约 50kb 含重复序列的模拟基因组）
  - 用 minimap2 生成参考基因组的 .mmi 索引
  - 用 samtools faidx 生成 .fai 索引
  - 生成简单的 gene.bed、cpg.bed、rmsk.bed 注释文件（用于 annotate 测试）
- **Acceptance Criteria Addressed**: [AC-1, AC-2, AC-3, AC-4, AC-5, AC-8]
- **Test Requirements**:
  - `programmatic` TR-1.1: testdata 目录包含 exp_reads.fastq、ref.fa、ref.fa.fai、ref.mmi
  - `programmatic` TR-1.2: 注释 BED 文件格式正确（BED6 格式）
  - `human-judgement` TR-1.3: 测试数据规模合理（参考基因组大小适中，可在合理时间内完成测试）
- **Notes**: 使用现有的 minimap2 和 samtools 模块或 Docker 容器生成测试数据

## [x] Task 2: 构建 CRESIL_TRIM 模块
- **Priority**: high
- **Depends On**: Task 1
- **Description**:
  - 创建 `modules/cresil/trim/` 目录
  - 创建 main.nf：process CRESIL_TRIM
  - 创建 meta.yml：描述输入输出和工具信息
  - 创建 environment.yml：conda 环境配置
  - 输入：tuple val(meta), path(reads) + tuple val(meta2), path(mmi)
  - 输出：tuple val(meta), path("*.trim.txt") emit: trim + versions.yml
  - 支持 task.ext.args 和 task.ext.prefix
  - 实现 stub 模式
  - 使用 Docker: quay.io/bioinfortools/cresil:1.2.0
  - 使用 Conda: 本地 environment.yml 或 yangmingsi/cresil
- **Acceptance Criteria Addressed**: [AC-1, AC-6, AC-7, AC-9]
- **Test Requirements**:
  - `programmatic` TR-2.1: main.nf 语法正确，可被 Nextflow 解析
  - `programmatic` TR-2.2: 正常运行模式下生成 trim.txt 文件且非空
  - `programmatic` TR-2.3: stub 模式下快速生成占位文件
  - `programmatic` TR-2.4: versions.yml 包含正确版本号
  - `human-judgement` TR-2.5: 代码风格与现有 flair/longshot 模块一致
- **Notes**: 参考 modules/flair/align/main.nf 的写法

## [x] Task 3: 构建 CRESIL_IDENTIFY 模块
- **Priority**: high
- **Depends On**: Task 1, Task 2
- **Description**:
  - 创建 `modules/cresil/identify/` 目录
  - 创建 main.nf：process CRESIL_IDENTIFY
  - 创建 meta.yml：描述输入输出和工具信息
  - 创建 environment.yml
  - 输入：tuple val(meta), path(fasta), path(fai), path(reads), path(trim)（trim 可选）
  - 输出：tuple val(meta), path("*eccDNA_final.txt") emit: identify + versions.yml
  - 支持 -s (skip-variant) 和 -sg (skip-gfa) 参数
  - 支持 -cm medaka 模型参数
  - 实现 stub 模式
- **Acceptance Criteria Addressed**: [AC-2, AC-6, AC-7, AC-9]
- **Test Requirements**:
  - `programmatic` TR-3.1: main.nf 语法正确，可被 Nextflow 解析
  - `programmatic` TR-3.2: 正常运行模式下生成 eccDNA_final.txt
  - `programmatic` TR-3.3: stub 模式下快速生成占位文件
  - `programmatic` TR-3.4: versions.yml 包含正确版本号
  - `human-judgement` TR-3.5: 代码风格与现有模块一致
- **Notes**: 由于 medaka 运行较慢，测试时可使用 -s (skip-variant) 加速

## [x] Task 4: 构建 CRESIL_IDENTIFY_WGLS 模块
- **Priority**: medium
- **Depends On**: Task 1, Task 2
- **Description**:
  - 创建 `modules/cresil/identify_wgls/` 目录
  - 创建 main.nf：process CRESIL_IDENTIFY_WGLS
  - 创建 meta.yml、environment.yml
  - 输入：tuple val(meta), path(mmi), path(fasta), path(fai), path(reads), path(trim)（trim 可选）
  - 输出：tuple val(meta), path("*eccDNA_final.txt") emit: identify_wgls + versions.yml
  - 支持 -m (mode) 参数
  - 实现 stub 模式
- **Acceptance Criteria Addressed**: [AC-3, AC-6, AC-7, AC-9]
- **Test Requirements**:
  - `programmatic` TR-4.1: main.nf 语法正确
  - `programmatic` TR-4.2: 正常运行模式下生成 eccDNA_final.txt
  - `programmatic` TR-4.3: stub 模式下快速生成占位文件
  - `programmatic` TR-4.4: versions.yml 包含正确版本号
  - `human-judgement` TR-4.5: 代码风格一致
- **Notes**: WGLS 模式适用于全基因组测序数据，测试时需注意数据量

## [x] Task 5: 构建 CRESIL_ANNOTATE 模块
- **Priority**: high
- **Depends On**: Task 1
- **Description**:
  - 创建 `modules/cresil/annotate/` 目录
  - 创建 main.nf：process CRESIL_ANNOTATE
  - 创建 meta.yml、environment.yml
  - 输入：tuple val(meta), path(identify_table), path(rmsk_bed) 可选, path(cpg_bed) 可选, path(gene_bed) 可选
  - 输出：tuple val(meta), path("*annotate*") 各注释文件 + versions.yml
  - 实现 stub 模式
- **Acceptance Criteria Addressed**: [AC-4, AC-6, AC-7, AC-9]
- **Test Requirements**:
  - `programmatic` TR-5.1: main.nf 语法正确
  - `programmatic` TR-5.2: 正常运行模式下生成注释结果文件
  - `programmatic` TR-5.3: stub 模式下快速生成占位文件
  - `programmatic` TR-5.4: versions.yml 包含正确版本号
  - `human-judgement` TR-5.5: 代码风格一致
- **Notes**: 三个 BED 输入均为可选，需处理部分输入的情况

## [x] Task 6: 构建 CRESIL_VISUALIZE 模块
- **Priority**: medium
- **Depends On**: Task 1
- **Description**:
  - 创建 `modules/cresil/visualize/` 目录
  - 创建 main.nf：process CRESIL_VISUALIZE
  - 创建 meta.yml、environment.yml
  - 输入：tuple val(meta), path(identify_table), val(eccdna_id)
  - 输出：tuple val(meta), path("for_Circos/") emit: circos_config + versions.yml
  - 支持 -uc (unit circos) 和 -mc (mode circos) 参数
  - 实现 stub 模式
- **Acceptance Criteria Addressed**: [AC-5, AC-6, AC-7, AC-9]
- **Test Requirements**:
  - `programmatic` TR-6.1: main.nf 语法正确
  - `programmatic` TR-6.2: 正常运行模式下生成 Circos 配置目录
  - `programmatic` TR-6.3: stub 模式下快速生成占位文件
  - `programmatic` TR-6.4: versions.yml 包含正确版本号
  - `human-judgement` TR-6.5: 代码风格一致
- **Notes**: 仅生成配置文件，不执行 circos 命令

## [/] Task 7: 编写各模块的 nf-test 测试（真实运行验证）
- **Priority**: high
- **Depends On**: Task 1, Task 2, Task 3, Task 4, Task 5, Task 6
- **Description**:
  - 为 trim 模块创建 tests/main.nf.test 和 snapshot
  - 为 identify 模块创建 tests/main.nf.test 和 snapshot
  - 为 identify_wgls 模块创建 tests/main.nf.test 和 snapshot
  - 为 annotate 模块创建 tests/main.nf.test 和 snapshot
  - 为 visualize 模块创建 tests/main.nf.test 和 snapshot
  - 每个测试包含正常运行测试和 stub 测试
- **Acceptance Criteria Addressed**: [AC-8, AC-9]
- **Test Requirements**:
  - `programmatic` TR-7.1: 所有模块的测试文件语法正确
  - `programmatic` TR-7.2: stub 测试全部通过
  - `programmatic` TR-7.3: 正常运行测试全部通过（使用 testdata）
  - `programmatic` TR-7.4: snapshot 文件正确生成
  - `human-judgement` TR-7.5: 测试用例覆盖关键输入输出场景
- **Notes**: 先跑 stub 测试确保结构正确，再跑完整测试

## [x] Task 8: 端到端验证与文档
- **Priority**: medium
- **Depends On**: Task 7
- **Description**:
  - 验证所有模块的 Docker 环境可正常运行
  - 验证所有模块的 Conda 环境可正常运行（Linux 平台）
  - 检查所有 meta.yml 的完整性和准确性
  - 确保模块间输入输出可串联（trim -> identify -> annotate -> visualize）
  - （可选）在 modules/cresil/ 下添加 README.md 说明使用方法
- **Acceptance Criteria Addressed**: [AC-6, AC-7, AC-8]
- **Test Requirements**:
  - `programmatic` TR-8.1: 模块间数据串联正常（trim 输出可作为 identify 输入等）
  - `programmatic` TR-8.2: 所有 meta.yml 格式正确且信息完整
  - `human-judgement` TR-8.3: 整体架构符合 bio.nf 项目规范
  - `human-judgement` TR-8.4: 代码可读性良好，有适当的空行和缩进
- **Notes**: 串联测试无需单独的 subworkflow，手动验证即可
