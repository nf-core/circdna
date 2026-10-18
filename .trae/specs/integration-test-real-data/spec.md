# Spec

## Why

验证最近所有更新（ECCsplorer + ecc_finder + 子工作流简并）在真实数据上能否正常产出结果。

## What Changes

无代码变更。仅执行以下测试命令并观察输出。

**测试命令**：

```bash
nextflow run main.nf -profile test_local,docker \
  --circle_identifier eccsplorer,ecc_finder_map_sr,ecc_finder_asm_sr
```

**检测工具**：ECCsplorer (map) + ECC_FINDER_MAP_SR + ECC_FINDER_ASM_SR（三者并行）
**排除项**：eccsplorer CLU（聚类耗时长）、ampliconsuite（缺 mosek 许可证实际数据）

## Expected Results

- 所有 3 个检测工具对 3 个样本各产生输出
- ECCsplorer: `results/test_local/eccsplorer/circdna_N_ecc_sequences.fasta` 非空
- ecc_finder MAP_SR: `results/test_local/ecc_finder/map_sr/circdna_N.fasta` 非空
- ecc_finder ASM_SR: `results/test_local/ecc_finder/asm_sr/circdna_N.fasta` 非空
- 退出码 0

## Notes

- ECCSPLORER_CLU 会因 `eccsplorer_clu=false` 跳过
- AMPLICONSUITE 会因 `aa_data_repo` 空目录跳过
