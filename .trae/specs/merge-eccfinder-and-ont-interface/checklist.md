# Checklist

## 模块合并

- [ ] `modules/local/ecc_finder/main.nf` 包含 4 个 process
- [ ] `environment.yml` 和 `meta.yml` 存在
- [ ] `ecc_finder_map_sr/` 和 `ecc_finder_asm_sr/` 已删除
- [ ] 所有 process 共享同一 container/conda

## 子工作流

- [ ] `subworkflows/local/ecc_finder_pipeline/main.nf` 存在
- [ ] SR map → ECC_FINDER_MAP_SR 路由正确
- [ ] SR asm → ECC_FINDER_ASM_SR 路由正确
- [ ] ONT map → ECC_FINDER_MAP_ONT 代码就绪（默认关闭）
- [ ] ONT asm → ECC_FINDER_ASM_ONT 代码就绪（默认关闭）

## 流程接入

- [ ] eccdna_mode 使用 ECC_FINDER_PIPELINE 替代独立 include
- [ ] modules.config 更新 withName 路径
- [ ] circle_identifier 控制不变（ecc_finder_map_sr / ecc_finder_asm_sr）

## Stub 测试

- [ ] `-stub` 测试 55 tasks 成功
- [ ] ecc_finder 输出目录存在
