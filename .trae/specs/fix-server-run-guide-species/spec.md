# 修正 SERVER_RUN_GUIDE 物种列表 Spec

## Why

SERVER*RUN_GUIDE.md 中包含了 3 个不存在 NGS 短读长 samplesheet 的物种，其命令引用了不存在的 `circdna*{species}\_eccDNA.csv`（实际仅存在于 long_read 目录）。命令无法执行，需根据实际 samplesheets 重新归类。

## What Changes

- 从"按物种运行命令" NGS 短读长章节中，将 3 个仅有 long_read 数据的物种移至新增的 circdnalr 三代长读长章节（占位）：
  - 向日葵 (Helianthus_annuus) → `circdnalr_Helianthus_annuus_long_read.csv`
  - 苋 (Amaranthus_palmeri) → `circdnalr_Amaranthus_palmeri_long_read.csv`
  - 黑麦草 (Alopecurus_myosuroides) → `circdnalr_Alopecurus_myosuroides_long_read.csv`，标注大基因组
- 新增「三代长读长物种 (circdnalr)」章节，以占位命令列出，标记为待实现
- 三代章节中黑麦草保留大基因组标注，需添加 `-c circdna.nf/conf/large_genome.config`
- 更新"注意事项"中"大基因组"列表：NGS 短读长仅保留"小麦、日本柳杉"；另在三代章节中单独标注黑麦草
- NGS 短读长章节保留 samplesheets 中实际存在的 12 个物种

## Impact

- Affected specs: 无
- Affected code: `circdna.nf/SERVER_RUN_GUIDE.md`
