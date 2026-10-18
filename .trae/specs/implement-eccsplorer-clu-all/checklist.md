# Checklist

## ECCSPLORER_CLU 模块

- [ ] `main.nf` 进程定义正确，container 指向 `quay.io/bioinfortools/eccsplorer:2022.01.1.1`
- [ ] 输入签名: `tuple val(meta), path(treatment_reads), path(control_reads), val(taxon)`
- [ ] 调用 `ECCsplorer.py <R1> <R2> <C1> <C2> --mode clu --taxon <taxon>`
- [ ] 正确提取 `CLUSTER_TABLE.csv`、`COMPARATIVE_CLUSTER_TABLE.csv`、`eccCL_summary.html`
- [ ] `environment.yml` + `meta.yml` 存在
- [ ] Stub 模式正常

## 清理

- [ ] `eccsplorer_clu_prepare/` 已删除
- [ ] `eccsplorer_clu_core/` 已删除
- [ ] `eccsplorer_clu_candidates_plot/` 已删除
- [ ] `subworkflows/local/eccsplorer_cluster/` 已删除
- [ ] `subworkflows/local/eccsplorer_all/` 已删除

## 流程接入

- [ ] `eccdna_mode/main.nf` 中 include ECCSPLORER_CLU
- [ ] 仅当 `run_eccsplorer=true` 且 `eccsplorer_clu=true` 时执行
- [ ] ECCSPLORER_CLU 与现有 map process 正确编排（map 后运行）

## 配置

- [ ] `conf/modules.config` 中 ECCSPLORER_CLU 配置正确
- [ ] `eccsplorer_clu` 参数默认 false
- [ ] PublicDB 路径可访问

## 测试

- [ ] `-profile test_local,docker --eccsplorer_clu true -resume` map 阶段缓存命中
- [ ] ECCSPLORER_CLU 执行成功
- [ ] Cluster candidates 输出非空
- [ ] 流程退出码 0
