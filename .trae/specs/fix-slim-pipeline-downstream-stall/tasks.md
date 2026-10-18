# Tasks

- [x] Task 1: SEGEMEHL_ALIGN 模块补充 segemehl 参数

  - [x] 1.1 modules.config 中为 SEGEMEHL_ALIGN 添加 `ext.args = '--splits --briefcigar --MEOP --accuracy 95'`
  - [x] 1.2 SEGEMEHL_ALIGN publishDir 的 saveAs 排除 `.sam` 文件
  - [x] 1.3 验证：stub 编译通过

- [x] Task 2: HAARZ 模块增加 awk 行修复（复刻 eccMapper.py start>end 纠错）

  - [x] 2.1 在 haarz.x 输出后添加 awk 修复并重命名
  - [x] 2.2 验证：stub 编译通过

- [x] Task 3: 修改流程后 -resume 全流程验证
  - [x] 3.1 直接以 -resume 运行 Docker 完整测试（test_3.csv + 3 个 slim circle_identifier）
  - [x] 3.2 验证 SEGEMEHL_ALIGN 因 ext.args 变更自然重跑，3 个样本全部重新 publish
  - [x] 3.3 验证 HAARZ/CANDIDATE_EXTRACT/COVERAGE_PROFILE/NORMALIZE/VISUALIZE/HTML_REPORT 全部运行
  - [x] 3.4 验证 bam 目录仅含 .bam，3 个样本齐全（不手动清理缓存）

# Task Dependencies

- [Task 2] depends on [Task 1]
- [Task 3] depends on [Task 1] + [Task 2]
