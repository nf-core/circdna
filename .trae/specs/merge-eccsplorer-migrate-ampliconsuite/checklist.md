# Checklist

## ECCsplorer 合并

- [ ] `eccsplorer/main.nf` 包含三个 process：ECCSPLORER、ECCSPLORER_WITH_CONTROL、ECCSPLORER_CLU
- [ ] `eccdna_mode/main.nf` include 路径已更新
- [ ] `modules/local/eccsplorer_clu/` 已删除
- [ ] 语法检查通过

## ampliconsuite 迁移

- [ ] `bio.nf/modules/ampliconsuite/` 存在且包含 main.nf、meta.yml、environment.yml
- [ ] main.nf container 为 `quay.io/biocontainers/ampliconsuite:1.6.0--pyh109da93_0`
- [ ] environment.yml 使用 `bioconda::ampliconsuite=1.6.0`
- [ ] Dockerfile 已删除
- [ ] `ampliconarchitect_pipeline/main.nf` include 路径已更新
- [ ] `circdna.nf/modules/local/ampliconsuite/` 已删除
