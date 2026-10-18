# Checklist

## 合并与代码完整性

- [ ] circdnalr 分支成功 `git merge master`，无未解决冲突
- [ ] 合并后分支同时包含短读长（master 侧 ECC_SCORE/Snakemake）与长读长（CReSIL/FLED/FLYE/ECCFINDER）逻辑
- [ ] `nextflow config -profile test_local` 可正常解析
- [ ] `bin/check_samplesheet.py` 语法校验通过
- [ ] master 内置 ONT/PacBio 测试数据（testdatasets/ont、testdatasets/pacbio）在合并后分支中存在
- [ ] samplesheet 冲突按约定解决：circdna*\* 系列以 master 为准，circdnalr*\* 长读系列合并两侧改动

## 测试配置

- [ ] `conf/test_nanopore_lr.config` 覆盖四引擎（cresil,fled,flye,eccfinder），指向内置 ont smoke 数据
- [ ] `conf/test_pacbio_lr.config` 覆盖四引擎（cresil,fled,flye,eccfinder），指向内置 pacbio smoke 数据
- [ ] 长读测试样本表路径使用服务器绝对路径且文件存在
- [ ] 本地 `nextflow config -profile test_nanopore_lr` / `test_pacbio_lr` 参数解析通过

## 服务器同步

- [ ] `circdnalr` 已推送到 origin/circdnalr
- [ ] 服务器 `/data1/users/siyangming/PlanteccDNADB/circdna.nf` 已切换到 circdnalr 分支并更新到最新
- [ ] 服务器 nextflow 环境（conda env nextflow + Java）可用

## ONT 测试

- [ ] ONT 长读长流程成功运行（四引擎）
- [ ] CReSIL 输出存在
- [ ] FLED 输出存在
- [ ] FLYE 输出存在
- [ ] ECCFINDER 输出存在
- [ ] 无致命错误，退出状态正常

## PacBio 测试

- [ ] PacBio 长读长流程成功运行（四引擎）
- [ ] CReSIL 输出存在
- [ ] FLED 输出存在
- [ ] FLYE 输出存在
- [ ] ECCFINDER 输出存在
- [ ] 无致命错误，退出状态正常

## 结果记录

- [ ] ONT/PacBio 测试命令、状态、输出路径、产出文件已记录
- [ ] 测试失败项已定位并修复（如有）
