# Checklist

- [x] 6 个 bio.nf 模块已迁入 local（cresil/ecc_finder/eccsplorer/flye/repeatexplorer2/ampliconsuite_ec）
- [x] ecc_finder 拆分为 4 个子模块（每 process 一目录），旧单文件 main.nf 已删除
- [x] ecc_finder_pipeline 已改引 4 个独立子模块（含 versions_ecc_finder 与 tag meta2 修复）
- [x] ampliconsuite 已重命名为 local/ampliconsuite_ec，未覆盖 local/ampliconsuite（prepareaa）
- [x] eccsplorer 黑盒已就位，eccsplorer_pipeline 引用不再断链（ECCSPLORER_CLU 断链已删除）
- [x] `grep bio.nf/modules` 在 \*.nf 无残留
- [x] 黑盒版与 slim 版 circle_identifier 在 workflows 中共存并按选择运行（stub 编译通过）
- [x] prepare_eccsplorer_database.sh 已归档到 testdatasets/eccsplorer_db/
- [x] update_samplesheets.py 已归档到 samplesheets/
- [x] test_incremental_cache.py 已删除，scripts/ 目录已删除
- [x] PRExing 预处理链（prepare\_\* + seqkit + ecc_preprocessing）已保留
- [x] 冲突文件已解决（master 为底 + PRExing 内容，黑盒旧 clu 配置丢弃）
- [x] samplesheet testdata→ngs 路径更新 + testdatasets/ngs 重命名已保留
- [x] 冗余文件（errors.txt 等）已删除
- [x] stub 编译通过（slim + 黑盒标识符），所有变更已提交（commit 8cd9991）
