# 为所有流程添加 Docker 用户映射权限配置 - The Implementation Plan (Decomposed and Prioritized Task List)

## [ ] Task 1: bio.nf 11 个模块测试配置补充 `-u $(id -u):$(id -g)`
- **Priority**: high
- **Depends On**: None
- **Description**: 
  - bio.nf 下 11 个模块测试配置文件（`modules/*/tests/nextflow.config`）目前只有 `runOptions = '--platform linux/amd64'`
  - 统一在 `--platform` 之前插入 `-u $(id -u):$(id -g)`，保持 `--platform linux/amd64` 不变
  - 涉及的 11 个文件：
    1. `modules/flye/tests/nextflow.config`
    2. `modules/ecc_finder/map_ont/tests/nextflow.config`
    3. `modules/ecc_finder/asm_ont/tests/nextflow.config`
    4. `modules/ecc_finder/map_sr/tests/nextflow.config`
    5. `modules/ecc_finder/asm_sr/tests/nextflow.config`
    6. `modules/cresil/identify_wgls/tests/nextflow.config`
    7. `modules/cresil/trim/tests/nextflow.config`
    8. `modules/cresil/visualize/tests/nextflow.config`
    9. `modules/cresil/annotate/tests/nextflow.config`
    10. `modules/cresil/identify/tests/nextflow.config`
    11. `modules/fastqdl/download/tests/nextflow.config`
  - **实现方式**：批量搜索替换，将 `runOptions = '--platform linux/amd64'` 统一替换为 `runOptions = '-u $(id -u):$(id -g) --platform linux/amd64'`
- **Acceptance Criteria Addressed**: AC-1, AC-3, AC-4
- **Test Requirements**:
  - `programmatic` TR-1.1: `grep -n 'runOptions' bio.nf/modules/*/tests/nextflow.config` 返回的所有 11 行内容均包含 `-u $(id -u):$(id -g) --platform linux/amd64`
  - `programmatic` TR-1.2: git diff 显示仅修改了上述 11 个文件，且每处改动仅在 runOptions 行中插入了 `-u` 参数
  - `human-judgement` TR-1.3: 抽查 3 个文件（例如 flye、cresil/identify、ecc_finder/asm_sr），确认 runOptions 值与预期完全一致，且其他行未被误改
- **Notes**: 11 处改动完全相同，建议用 sed 批量处理或逐个文件 Edit 替换以保证精确

## [ ] Task 2: nanoseq.nf docker profile 添加显式 `docker.runOptions`
- **Priority**: high
- **Depends On**: None
- **Description**: 
  - `nanoseq.nf/nextflow.config` 的 docker profile（L146-153）当前只有 `docker.userEmulation = true`，缺少显式 `docker.runOptions`
  - 在 docker profile 中添加一行：`docker.runOptions = '-u $(id -u):$(id -g)'`
  - **不删除** 现有的 `docker.userEmulation = true`，两者同时存在以提供双重保险（userEmulation 由 Nextflow 自动注入，runOptions 是显式设置，两者不冲突）
  - 改动位置：在 docker profile 中，紧跟在 `docker.userEmulation = true` 之后添加 runOptions 行
- **Acceptance Criteria Addressed**: AC-2, AC-3
- **Test Requirements**:
  - `programmatic` TR-2.1: `grep -n 'docker.runOptions\|docker.userEmulation' nanoseq.nf/nextflow.config` 显示 docker profile 区域同时包含 `docker.runOptions = '-u $(id -u):$(id -g)'` 和 `docker.userEmulation = true`
  - `programmatic` TR-2.2: git diff 显示仅 nanoseq.nf/nextflow.config 的 docker profile 块中增加了 1 行 runOptions
  - `human-judgement` TR-2.3: 检查该行缩进是否与 docker profile 中其他行一致（4 个空格缩进）
- **Notes**: 不修改 arm profile 或 base.config（已在之前对话中修改过 base.config 的 GPU 条件块）

## [ ] Task 3: 批量一致性验证
- **Priority**: medium
- **Depends On**: Task 1, Task 2
- **Description**: 
  - 运行一次性 grep 命令，统计所有流程配置文件中 `runOptions` 包含 `-u` 的覆盖率
  - 生成一份"修改后"的快速验证表格，确认：
    - bio.nf 11 个模块测试：全部有 `-u`
    - nanoseq.nf docker profile：有 `-u`
    - circdna.nf、circrna.nf、isoseq.nf、fetchngs.nf、riboseq.nf、rnaseq 的 docker profile：保持未变（有 `-u`）
- **Acceptance Criteria Addressed**: AC-3, AC-4
- **Test Requirements**:
  - `programmatic` TR-3.1: 输出 8 个流程 + bio.nf 模块测试的配置覆盖率汇总表格
  - `human-judgement` TR-3.2: 人工确认汇总表中所有标记为"应有 `-u`"的项均为 ✅
