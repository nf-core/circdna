# CReSIL Nextflow 模块 - Verification Checklist

## 目录结构检查
- [x] `modules/cresil/trim/` 目录存在，包含 main.nf, meta.yml, environment.yml, tests/
- [x] `modules/cresil/identify/` 目录存在，包含 main.nf, meta.yml, environment.yml, tests/
- [x] `modules/cresil/identify_wgls/` 目录存在，包含 main.nf, meta.yml, environment.yml, tests/
- [x] `modules/cresil/annotate/` 目录存在，包含 main.nf, meta.yml, environment.yml, tests/
- [x] `modules/cresil/visualize/` 目录存在，包含 main.nf, meta.yml, environment.yml, tests/
- [x] `modules/cresil/testdata/` 目录存在，包含测试数据文件

## 测试数据检查
- [x] testdata/exp_reads.fastq 存在且非空
- [x] testdata/ref.fa 存在且为有效 FASTA
- [x] testdata/ref.fa.fai 存在且有效
- [x] testdata/ref.mmi 存在且为有效 minimap2 索引
- [x] testdata/gene.bed 存在且为有效 BED 格式
- [x] testdata/cpg.bed 存在且为有效 BED 格式
- [x] testdata/rmsk.bed 存在且为有效 BED 格式

## TRIM 模块检查
- [x] main.nf 中 process 名为 CRESIL_TRIM（全大写）
- [x] meta.yml 包含完整的 input/output 描述
- [x] environment.yml 包含正确的 conda 配置
- [x] 输入：tuple val(meta), path(reads) + tuple val(meta2), path(mmi)
- [x] 输出：trim.txt + versions
- [x] 支持 task.ext.args 传参
- [x] 支持 task.ext.prefix 自定义前缀
- [x] stub 模式正常工作
- [x] conda 指令指向正确的 environment.yml
- [x] container 指令指向 quay.io/bioinfortools/cresil:1.2.0
- [x] tests/main.nf.test 存在且包含正常测试和 stub 测试
- [x] tests/main.nf.test.snap 存在

## IDENTIFY 模块检查
- [x] main.nf 中 process 名为 CRESIL_IDENTIFY
- [x] meta.yml 包含完整的 input/output 描述
- [x] environment.yml 包含正确的 conda 配置
- [x] 输入：tuple val(meta), path(fasta), path(fai), path(reads), path(trim) 可选
- [x] 输出：eccDNA_final.txt + versions
- [x] 支持 -s (skip-variant) 和 -sg (skip-gfa) 参数（通过 task.ext.args）
- [x] 支持 -cm medaka 模型参数（通过 task.ext.args）
- [x] stub 模式正常工作
- [x] conda 和 container 配置正确
- [x] tests/main.nf.test 存在且包含正常测试和 stub 测试

## IDENTIFY_WGLS 模块检查
- [x] main.nf 中 process 名为 CRESIL_IDENTIFY_WGLS
- [x] meta.yml 包含完整的 input/output 描述
- [x] environment.yml 包含正确的 conda 配置
- [x] 输入：tuple val(meta), path(mmi), path(fasta), path(fai), path(reads), path(trim) 可选
- [x] 输出：eccDNA_final.txt + versions
- [x] 支持 -m (mode) 参数（通过 task.ext.args）
- [x] stub 模式正常工作
- [x] conda 和 container 配置正确（Singularity 后缀已统一）
- [x] tests/main.nf.test 存在

## ANNOTATE 模块检查
- [x] main.nf 中 process 名为 CRESIL_ANNOTATE
- [x] meta.yml 包含完整的 input/output 描述
- [x] environment.yml 包含正确的 conda 配置
- [x] 输入：tuple val(meta), path(identify_table), 三个可选 BED 文件
- [x] 输出：各注释结果文件 + versions
- [x] stub 模式正常工作
- [x] conda 和 container 配置正确
- [x] tests/main.nf.test 存在

## VISUALIZE 模块检查
- [x] main.nf 中 process 名为 CRESIL_VISUALIZE
- [x] meta.yml 包含完整的 input/output 描述（含注释文件输入）
- [x] environment.yml 包含正确的 conda 配置
- [x] 输入：tuple val(meta), path(identify_table), 注释文件 x4, val(eccdna_id)
- [x] 输出：for_Circos/ 目录 + versions
- [x] 支持 -uc 和 -mc 参数（通过 task.ext.args）
- [x] stub 模式正常工作
- [x] conda 和 container 配置正确
- [x] tests/main.nf.test 存在

## 规范一致性检查
- [x] 所有 process 名遵循 CRESIL_<SUBCOMMAND> 全大写命名
- [x] 所有模块的 meta.yml 结构与 flair/longshot 一致
- [x] 所有模块的 versions 输出格式一致
- [x] 所有模块使用相同的 Docker 镜像 tag
- [x] 所有模块使用相同的标签策略（tag "$meta.id"）
- [x] 所有模块都有 label 'process_high'
- [x] 所有模块的 Singularity 镜像后缀统一为 hdfd78af_0

## 测试验证
- [x] TRIM 模块 stub 测试通过
- [x] TRIM 模块正常运行测试通过，输出 trim.txt 非空
- [x] IDENTIFY 模块 stub 测试通过
- [ ] IDENTIFY 模块正常运行测试通过（测试数据无法检测到 eccDNA，需更复杂的数据）
- [x] IDENTIFY_WGLS 模块 stub 测试通过
- [ ] IDENTIFY_WGLS 模块正常运行测试通过（CReSIL 源码 bug: KeyError 'name'）
- [x] ANNOTATE 模块 stub 测试通过
- [x] ANNOTATE 模块正常运行测试通过
- [x] VISUALIZE 模块 stub 测试通过
- [x] VISUALIZE 模块正常运行测试通过

## 端到端串联验证
- [x] trim 输出可作为 identify 输入（文件模式匹配 *.trim.txt）
- [x] identify 输出可作为 annotate 输入（文件模式匹配 *eccDNA_final.txt）
- [x] identify 输出可作为 visualize 输入（文件模式匹配 *eccDNA_final.txt）
- [x] annotate 输出可作为 visualize 输入（4 个注释文件模式匹配）
- [x] 完整流程（trim -> identify -> annotate -> visualize）数据流兼容
