# Flye Nextflow模块 - 验证清单

## 文件结构验证
- [ ] 模块目录 modules/flye/ 存在
- [ ] main.nf 文件存在且语法正确
- [ ] meta.yml 文件存在且格式正确
- [ ] environment.yml 文件存在且包含正确的依赖
- [ ] tests/ 目录存在
- [ ] tests/main.nf.test 文件存在
- [ ] tests/nextflow.config 文件存在
- [ ] flye分支已创建并推送到远程

## Process定义验证
- [ ] process名称为 FLYE（全大写）
- [ ] tag使用 $meta.id
- [ ] label设置为 'process_high'
- [ ] 包含 when 块支持 task.ext.when

## 环境配置验证
- [ ] conda配置指向 ${moduleDir}/environment.yml
- [ ] Docker镜像为 quay.io/biocontainers/flye:2.9.6--py313h7fbb527_1
- [ ] Singularity配置正确（Galaxy镜像或docker镜像转换）
- [ ] environment.yml包含 bioconda::flye=2.9.6
- [ ] tests/nextflow.config 中 docker.runOptions 设置为 '--platform linux/amd64'

## 输入输出验证
- [ ] 输入为 tuple val(meta), path(reads)
- [ ] 输出 assembly 通道（assembly.fasta）带 emit: assembly
- [ ] 输出 gfa 通道（assembly_graph.gfa）带 emit: gfa
- [ ] 输出 info 通道（assembly_info.txt）带 emit: info
- [ ] 输出 versions 通道（versions.yml）带 emit: versions
- [ ] meta.yml中input/output描述与main.nf一致

## 功能验证
- [ ] script块使用 flye 命令
- [ ] 支持 --out-dir 参数设置为当前目录
- [ ] 使用 task.ext.args 传递额外参数
- [ ] 使用 task.ext.prefix 定义输出前缀（默认meta.id）
- [ ] versions.yml文件格式正确（YAML格式，包含版本号）
- [ ] stub块生成模拟的assembly.fasta文件
- [ ] stub块生成模拟的assembly_graph.gfa文件
- [ ] stub块生成模拟的assembly_info.txt文件
- [ ] stub块生成模拟的versions.yml文件

## 测试验证
- [ ] 至少包含一个stub模式测试用例
- [ ] 测试用例包含assertAll断言
- [ ] snapshot测试文件 main.nf.test.snap 存在
- [ ] nf-test stub模式运行成功

## 文档验证
- [ ] meta.yml包含name字段
- [ ] meta.yml包含description字段
- [ ] meta.yml包含keywords字段
- [ ] meta.yml包含tools字段（含homepage, licence等）
- [ ] meta.yml包含完整的input描述
- [ ] meta.yml包含完整的output描述
- [ ] meta.yml包含authors和maintainers

## 代码风格验证
- [ ] 代码风格与现有nf-core模块一致
- [ ] 变量命名规范
- [ ] 缩进正确（4空格）
- [ ] 无多余的注释或调试代码

## 分支与提交验证
- [ ] flye分支已创建
- [ ] 所有模块文件已提交
- [ ] .trae/specs/flye_nf_module/文档已提交
- [ ] 源代码目录（Flye/）未提交到分支
- [ ] 分支已推送到远程仓库
