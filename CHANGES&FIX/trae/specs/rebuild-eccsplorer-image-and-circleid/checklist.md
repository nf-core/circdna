# Checklist

## Phase 1 — bio.nf

### Docker 镜像

- [x] 多阶段 Dockerfile 语法正确
- [x] `seqclust --help` 正常输出
- [x] `ECCsplorer.py --help` 正常输出
- [x] 镜像推送至 quay.io

### Conda 声明

- [x] meta.yaml 注释准确

### nf-test

- [x] 非 stub 测试用例存在
- [x] `nf-test test` 全部通过
- [x] 输出 BED/FASTA 非空

## Phase 2 — circdna.nf

### 模块同步

- [x] bio.nf ↔ circdna.nf main.nf 一致

### pair 列

- [x] check_samplesheet.py 支持 pair
- [x] input_check 填充 meta.pair
- [x] schema_input.json 含 pair
- [x] eccdna_mode 用 meta.pair 配对（不用 group）
- [x] 同 pair 值的 eccDNA+gDNA → ECCSPLORER_WITH_CONTROL
- [x] 无 pair 值 → ECCSPLORER（无 control）
- [x] group 列保留用于通用分组，不影响 ECCsplorer 配对

### circle_identifier

- [x] circle_identifier 默认值含 eccsplorer
- [x] eccsplorer_map_core 已移除
- [x] 含 eccsplorer → 运行，不含 → 不运行
- [x] 仅 eccDNA → ECCSPLORER
- [x] eccDNA+gDNA + 同 pair → ECCSPLORER_WITH_CONTROL
- [x] 各配置文件中 circle_identifier 已更新
- [x] 无额外 mode/control 参数

### 端到端测试

- [x] `-profile test_local,docker -resume` 成功
- [x] ECCSPLORER 退出码 0
- [x] BED/FASTA 非空
- [x] 流程退出码 0
- [x] 文件归属正确
