# Checklist

- [x] TIDEHUNTER 模块按官方 CLI 调用（输入 reads，输出 consensus fasta）
- [x] clu_prepare.py 提取成功并加入 slim 镜像
- [x] REPEATEXPLORER2 模块存在于 bio.nf，使用独立 repeatexplorer:2.3.8 镜像并可调用 seqclust
- [x] clu_candidates.py 提取成功并加入 slim 镜像
- [x] eccsplorer_clu_slim 子工作流可编译（stub），circle_identifier 支持 eccsplorer_clu_slim
- [x] tests/tidehunter/test.nf 运行成功，consensus fasta 非空含重复单元
- [x] tests/eccsplorer_clu/test.nf 运行成功，cluster_candidates.csv 非空
- [x] results/slim_run/tidehunter_test/ 与 results/slim_run/eccsplorer_clu_test/ 有产物
- [x] 主流水线 -resume 回归成功（3 个 slim circle_identifier）
- [x] modularize-eccsplorer-modes Task 14/15 已勾选（all 骨架 14.3 保留待后续；完整版 implement Task 5 未单独运行，已在 tasks.md 注释说明由 slim clu 验证覆盖）
