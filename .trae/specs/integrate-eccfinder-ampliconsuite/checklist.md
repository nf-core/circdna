# Checklist

## 分支

- [ ] `integrate-eccfinder-ampliconsuite` 分支已创建

## ecc_finder_asm_sr 模块

- [ ] 已从 bio.nf 拷贝到 `modules/local/ecc_finder_asm_sr/`

## ECCDNA_MODE 接入

- [ ] include ECC_FINDER_MAP_SR、ECC_FINDER_ASM_SR、AMPLICONSUITE
- [ ] `run_ecc_finder_map_sr` / `run_ecc_finder_asm_sr` / `run_ampliconsuite` 参数正确
- [ ] MAP_SR 输入 BAM+FASTA，输出 CSV+FASTA
- [ ] ASM_SR 输入 FASTQ，输出 FASTA
- [ ] AMPLICONSUITE 输入 BAM+mosek+aa_data_repo

## circle_identifier 独立控制

- [ ] `ecc_finder_map_sr` → 仅 MAP_SR 运行
- [ ] `ecc_finder_asm_sr` → 仅 ASM_SR 运行
- [ ] `ampliconsuite` → 仅 AMPLICONSUITE 运行
- [ ] 互不影响

## 配置

- [ ] modules.config 中三个 publishDir 正确
- [ ] test_local.config 参数完整

## 测试

- [ ] -profile test_local,docker -resume 成功
- [ ] MAP_SR 输出非空
- [ ] ASM_SR 输出非空
- [ ] 流程退出码 0
