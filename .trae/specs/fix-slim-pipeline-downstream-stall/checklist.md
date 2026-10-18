# Checklist

- [x] SEGEMEHL_ALIGN 已配置 `--splits --briefcigar --MEOP --accuracy 95`
- [x] SEGEMEHL_ALIGN publishDir 排除 .sam
- [x] HAARZ 包含 awk start>end 行修复
- [x] Docker 全流程 -resume 运行成功（3 个 slim 模式，不清理缓存）
- [x] SEGEMEHL_ALIGN 因 ext.args 变更自然重跑，3 样本重新 publish
- [x] HAARZ + 下游 5 个进程（CANDIDATE_EXTRACT/COVERAGE_PROFILE/NORMALIZE/VISUALIZE/HTML_REPORT）全部执行
- [x] bam 目录仅含 .bam，3 个样本齐全，无 .sam
