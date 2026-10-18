# Flye Nextflow模块 - 实现计划（分解并排序的任务列表）

## [ ] Task 1: 创建flye分支和模块目录结构
- **Priority**: high
- **Depends On**: None
- **Description**: 
  - 从主分支创建新分支flye
  - 创建modules/flye/目录结构
  - 创建modules/flye/tests/目录
  - 创建environment.yml文件，配置bioconda::flye=2.9.6
- **Acceptance Criteria Addressed**: AC-1, AC-3
- **Test Requirements**:
  - `programmatic` TR-1.1: 分支flye已创建
  - `programmatic` TR-1.2: 目录modules/flye/存在且包含environment.yml和tests/子目录
  - `programmatic` TR-1.3: environment.yml包含name: flye和bioconda::flye=2.9.6依赖
- **Notes**: 参照现有模块的目录结构

## [ ] Task 2: 准备测试数据（E.coli小型长读长数据集）
- **Priority**: high
- **Depends On**: Task 1
- **Description**: 
  - 从Zenodo或GitHub获取E.coli PacBio/ONT长读长示例数据的子集
  - 下载小型FASTA格式测试数据
  - 保存为modules/flye/testdata/目录下的文件
  - 确保测试数据足够小（stub模式不需要真实数据，正常模式可快速运行）
- **Acceptance Criteria Addressed**: AC-9
- **Test Requirements**:
  - `programmatic` TR-2.1: 测试数据文件存在于modules/flye/testdata/目录
  - `programmatic` TR-2.2: 至少包含一个FASTA格式的长读长测试文件（E.coli数据集子集）
- **Notes**: 数据集参考Flye官方文档中提到的E.coli PacBio或ONT数据集，取一小部分即可

## [ ] Task 3: 创建main.nf模块主文件
- **Priority**: high
- **Depends On**: Task 1
- **Description**: 
  - 定义FLYE process
  - 配置input: tuple val(meta), path(reads)
  - 配置output: assembly(assembly.fasta), gfa(assembly_graph.gfa), info(assembly_info.txt), versions(versions.yml)
  - 配置conda和docker/singularity环境
  - 实现script块，调用flye命令
  - 实现stub块，生成模拟输出文件
  - 支持task.ext.args和task.ext.prefix
  - 支持--out-dir参数指向当前目录
- **Acceptance Criteria Addressed**: AC-2, AC-3, AC-4, AC-5, AC-6, AC-7
- **Test Requirements**:
  - `programmatic` TR-3.1: process名为FLYE
  - `programmatic` TR-3.2: 输入为tuple val(meta), path(reads)
  - `programmatic` TR-3.3: 输出包含assembly、gfa、info、versions四个通道
  - `programmatic` TR-3.4: conda和docker镜像配置正确
  - `programmatic` TR-3.5: script块使用flye命令，支持$args
  - `programmatic` TR-3.6: stub块创建模拟assembly.fasta、assembly_graph.gfa、assembly_info.txt和versions.yml
  - `human-judgement` TR-3.7: 代码风格与现有nf-core模块一致
- **Notes**: 参照minimap2或类似模块的写法

## [ ] Task 4: 创建meta.yml模块文档
- **Priority**: high
- **Depends On**: Task 3
- **Description**: 
  - 编写meta.yml文件，包含完整的模块文档
  - 描述工具信息（name, description, homepage, licence等）
  - 详细描述input和output参数
  - 添加authors和maintainers
- **Acceptance Criteria Addressed**: AC-8
- **Test Requirements**:
  - `human-judgement` TR-4.1: meta.yml包含所有必需字段
  - `human-judgement` TR-4.2: input和output描述准确，与main.nf一致
  - `human-judgement` TR-4.3: tools部分包含flye的主页、文档和许可证信息
- **Notes**: 参照现有模块的meta.yml格式

## [ ] Task 5: 创建测试配置文件和测试用例
- **Priority**: high
- **Depends On**: Task 3, Task 4
- **Description**: 
  - 创建tests/nextflow.config，配置docker.enabled=true和runOptions='--platform linux/amd64'
  - 创建tests/main.nf.test，包含nf-test测试用例
  - 测试用例1: stub模式运行（使用测试数据）
  - 测试用例2: 正常模式（可选，需要较长时间）
- **Acceptance Criteria Addressed**: AC-9, AC-11
- **Test Requirements**:
  - `programmatic` TR-5.1: tests/nextflow.config存在，配置了docker.enabled=true和--platform linux/amd64
  - `programmatic` TR-5.2: tests/main.nf.test存在，包含至少2个测试用例
  - `programmatic` TR-5.3: stub模式测试用例可以成功运行
- **Notes**: 参照现有模块的测试写法

## [ ] Task 6: 运行stub模式测试并生成snapshot
- **Priority**: high
- **Depends On**: Task 5
- **Description**: 
  - 运行nf-test stub模式测试
  - 验证输出文件是否符合预期
  - 生成并保存snapshot文件
  - 修复测试中发现的问题
- **Acceptance Criteria Addressed**: AC-10
- **Test Requirements**:
  - `programmatic` TR-6.1: nf-test --stub模式测试通过
  - `programmatic` TR-6.2: main.nf.test.snap文件生成且包含正确的快照
  - `programmatic` TR-6.3: 所有assertAll断言通过
- **Notes**: 确保stub输出与实际输出格式一致

## [ ] Task 7: Docker环境测试验证（可选）
- **Priority**: medium
- **Depends On**: Task 6
- **Description**: 
  - 使用真实Docker环境运行测试
  - 验证flye可以正常运行（--version）
  - 验证输出文件格式正确
  - 验证--platform linux/amd64配置正常工作
- **Acceptance Criteria Addressed**: AC-3, AC-5, AC-11
- **Test Requirements**:
  - `programmatic` TR-7.1: Docker镜像可以正常拉取和启动（带--platform linux/amd64）
  - `programmatic` TR-7.2: flye --version输出正确版本号（v2.9.6）
  - `programmatic` TR-7.3: 组装输出文件格式正确
- **Notes**: 如果数据量大，可以只验证版本号和基本功能

## [ ] Task 8: 提交并推送flye分支
- **Priority**: high
- **Depends On**: Task 6
- **Description**: 
  - 将所有修改提交到flye分支
  - 推送到远程仓库
  - 不包含源代码（Flye/目录），只包含模块文件和.trae文档
- **Acceptance Criteria Addressed**: AC-1
- **Test Requirements**:
  - `programmatic` TR-8.1: flye分支已推送到远程
  - `programmatic` TR-8.2: 所有模块文件已提交
  - `programmatic` TR-8.3: .trae/specs/flye_nf_module/文档已提交
- **Notes**: 遵循之前cresil/ecc_finder分支的提交规范
