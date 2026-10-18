# ECCsplorer 分支冻结修改保留价值评估（merge master 后）

日期：2026-08-13，分支 ECCsplorer（已并入 master 提交 444ce77）

## 背景

- master 已删除完整版 ECCsplorer 黑盒模块（816918e），走向 slim 原子化链
- ECCsplorer 分支冻结修改经 `stash -u → merge master（冲突以 master 为准）→ stash pop` 后保留如下状态
- 冲突文件处理：`conf/server.config`（master 版）、`conf/modules.config`/`nextflow.config`/`workflows/circdna.nf`（master 版为底 + slim 专属内容合并）

## 保留清单（提交，保留价值核心）

| 类别             | 内容                                                                                                                                                                                                                                                                         | 价值                                                                                                                              |
| ---------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------- | -------------------- |
| slim 模块链      | `modules/local/ecc_finder_slim/`、`eccsplorer_slim/`（peak*detect/candidate_extract/coverage_profile/normalize/visualize/html_report/dr_detect/blast*\_/distribution/ont\_\_/paf*filter/contract_export/comparative*\*）、`genrich/`、`haarz/`、`tidehunter/`、`ecc_finder/` | **核心**：完全替代完整版的原子化链                                                                                                |
| vendored 模块    | `modules/nf-core/blast/`、`cdhit/`、`trimmomatic/`                                                                                                                                                                                                                           | 完整版 BLAST 注释 / ONT 聚类 / PRExer 裁剪                                                                                        |
| 新子工作流       | `ecc_finder_ont_slim/`、`ecc_finder_slim_pipeline/`、`eccsplorer_slim_pipeline/`、`eccsplorer_clu_slim/`、`eccsplorer_all_slim/`、`eccsplorer_prexer_slim/`                                                                                                                  | map/clu/all/PRExer/ONT 全模式编排                                                                                                 |
| workflows 接入   | `workflows/circdna.nf`（master 版为底 + slim include/调用/datatype 分流/pair 配对）                                                                                                                                                                                          | slim 标识符可运行                                                                                                                 |
| 参数暴露         | `nextflow.config`（slim params 段）、`conf/modules.config`（slim 模块配置段）                                                                                                                                                                                                | 附表 2 全部可配                                                                                                                   |
| samplesheet 扩展 | `bin/check_samplesheet.py`（+pair 列透传）、`subworkflows/local/input_check/main.nf`（+meta.pair）                                                                                                                                                                           | clu/control pair 配对必需                                                                                                         |
| 测试数据与文件   | `testdatasets/ont                                                                                                                                                                                                                                                            | pacbio/`（smoke/regular fastq）、`test_ont.csv`/`test_pacbio.csv`（相对路径）、`test_eccsplorer_pair.csv`、`extract_test_data.sh` | 服务器 slim 测试入口 |
| 测试             | `tests/`（ecc_finder_ont/ecc_finder_sr/tidehunter 等）                                                                                                                                                                                                                       | 各链独立验证                                                                                                                      |
| 修复             | `circle_map_pipeline`（+ch_qname_sorted_bai，master 9fdf560 引入的未定义变量 bug）、`unicycler` stub                                                                                                                                                                         | master bug 修复                                                                                                                   |
| 文档             | `.trae/specs/`（新 spec）、`.trae/documents/`（slim 计划/指南/评估）                                                                                                                                                                                                         | 追踪与交付                                                                                                                        |

## 删除清单（已删除，合理清理）

| 内容                                          | 原因                                                          |
| --------------------------------------------- | ------------------------------------------------------------- | --------------------------------------------------------------- | --------------------------------------------------------------------------- |
| `modules/local/eccsplorer_clu_candidates_plot | clu_core                                                      | clu_prepare/`                                                   | 旧黑盒 clu 模块，被 slim 的 clu_prepare/repeatexplorer2/clu_candidates 替代 |
| `subworkflows/local/eccsplorer_all            | eccsplorer_cluster/`                                          | 旧黑盒编排，被 `eccsplorer_all_slim`/`eccsplorer_clu_slim` 替代 |
| `testdatasets/samplesheet/`                   | samplesheet 已统一至 `samplesheets/`（test\_ 前缀，相对路径） |

## 丢弃清单

无（删除项均为合理清理，无冗余残留）。

## 结论

**全部保留**：slim 模块链 + vendored + 新子工作流 + workflows 接入 + params + samplesheet 扩展 + 测试数据为"完全替代"核心资产，全部提交到 ECCsplorer 分支。冲突文件已按用户决策以 master 为准（master 删除的黑盒模块配置/调用不恢复）。
