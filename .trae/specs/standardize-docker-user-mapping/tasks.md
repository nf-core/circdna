# 所有流程统一 Docker 用户映射 A+B+C 方案 - The Implementation Plan

## [ ] Task 1: 6 个标准流程 nextflow.config docker profile 补齐 B+C
- **Priority**: high
- **Depends On**: None
- **Description**: 
  - 为 circdna.nf、circrna.nf、isoseq.nf、fetchngs.nf、riboseq.nf、rnaseq 的 `nextflow.config → profiles → docker { }` 块添加 B（`docker.userEmulation = true`）和 C（`docker.fixOwnership = true`）
  - 这 6 个流程已有 A（`docker.runOptions = '-u $(id -u):$(id -g)'`），只需在 runOptions 行之后插入 B 和 C 两行
  - 标准插入位置：在 `docker.runOptions = '...'` 行之后、其他引擎 `xxx.enabled = false` 行之前
  - 每个文件的 docker profile 块改动后应类似：
    ```groovy
    docker {
        docker.enabled          = true
        conda.enabled           = false
        singularity.enabled     = false
        podman.enabled          = false
        shifter.enabled         = false
        charliecloud.enabled    = false
        apptainer.enabled       = false
        docker.runOptions       = '-u $(id -u):$(id -g)'
        docker.userEmulation    = true                   // 新增 B
        docker.fixOwnership     = true                   // 新增 C
    }
    ```
  - 6 个文件：
    1. `circdna.nf/nextflow.config`（docker profile L154-163，runOptions 在 L162）
    2. `circrna.nf/nextflow.config`（docker profile L164-173，runOptions 在 L172）
    3. `isoseq.nf/nextflow.config`（docker profile L128-137，runOptions 在 L136）
    4. `fetchngs.nf/nextflow.config`（docker profile L94-103，runOptions 在 L102）
    5. `riboseq.nf/nextflow.config`（docker profile L171-180，runOptions 在 L179）
    6. `rnaseq/nextflow.config`（docker profile L184-193，runOptions 在 L192）
- **Acceptance Criteria Addressed**: AC-1
- **Test Requirements**:
  - `programmatic` TR-1.1: 对 6 个文件逐一 Grep，docker profile 块内同时匹配 `docker.runOptions`、`docker.userEmulation`、`docker.fixOwnership` 三个关键词
  - `programmatic` TR-1.2: `docker.userEmulation` 的值为 `true`（不是 false 或缺失）
  - `programmatic` TR-1.3: `docker.fixOwnership` 的值为 `true`
  - `human-judgement` TR-1.4: 抽查 2 个文件确认缩进与周围行一致（8 空格，对齐其他 `docker.*` 行）
- **Notes**: 不修改 emulate_amd64 profile 中的 runOptions（那些已以 `-u` 开头）；不修改其他引擎 enabled 行

## [ ] Task 2: nanoseq.nf nextflow.config docker profile 补齐 A+C
- **Priority**: high
- **Depends On**: None
- **Description**: 
  - 文件：`nanoseq.nf/nextflow.config` 的 `profiles { docker { } }` 块（L146-153）
  - 当前状态：已有 B（`docker.userEmulation = true`，L148），缺少 A 和 C
  - 在 `docker.userEmulation = true` 行之前插入 A（`docker.runOptions = '-u $(id -u):$(id -g)'`）
  - 在 `docker.userEmulation = true` 行之后插入 C（`docker.fixOwnership = true`）
  - 改动后 docker profile 块应类似：
    ```groovy
    docker {
        docker.enabled         = true
        docker.runOptions      = '-u $(id -u):$(id -g)'   // 新增 A
        docker.userEmulation   = true                      // 已有 B
        docker.fixOwnership    = true                      // 新增 C
        singularity.enabled    = false
        podman.enabled         = false
        shifter.enabled        = false
        charliecloud.enabled   = false
    }
    ```
  - 不修改 arm profile（L154-156，已有 `-u` + `--platform`）
- **Acceptance Criteria Addressed**: AC-1
- **Test Requirements**:
  - `programmatic` TR-2.1: docker profile 块内同时匹配 `docker.runOptions`、`docker.userEmulation`、`docker.fixOwnership`
  - `programmatic` TR-2.2: runOptions 值严格为 `'-u $(id -u):$(id -g)'`（不含 `--platform`）
  - `human-judgement` TR-2.3: 缩进与周围行一致
- **Notes**: nanoseq 没有 apptainer.enabled 行，不要额外添加

## [ ] Task 3: bio.nf 顶层 nextflow.config 补齐 B+C
- **Priority**: high
- **Depends On**: None
- **Description**: 
  - 文件：`bio.nf/nextflow.config`（L5-8，顶层 docker 块，无 profiles）
  - 当前状态：已有 A（`runOptions = '-u $(id -u):$(id -g) --platform linux/amd64'`，L7）
  - 在 `enabled = true` 行之后、`runOptions` 行之前插入 B（`userEmulation = true`）
  - 在 `runOptions` 行之后插入 C（`fixOwnership = true`）
  - 改动后：
    ```groovy
    docker {
        enabled = true
        userEmulation = true                               // 新增 B
        runOptions = '-u $(id -u):$(id -g) --platform linux/amd64'
        fixOwnership = true                                // 新增 C
    }
    ```
- **Acceptance Criteria Addressed**: AC-2
- **Test Requirements**:
  - `programmatic` TR-3.1: docker 块内同时匹配 `runOptions`、`userEmulation`、`fixOwnership`
  - `programmatic` TR-3.2: runOptions 值保持 `'-u $(id -u):$(id -g) --platform linux/amd64'` 不变
  - `human-judgement` TR-3.3: 缩进与周围行一致（4 空格）
- **Notes**: bio.nf 无 profiles 块，docker 配置在顶层

## [ ] Task 4: bio.nf 11 个模块测试配置补齐 A+B+C
- **Priority**: high
- **Depends On**: None
- **Description**: 
  - bio.nf 下 11 个模块测试配置文件，路径均为 `bio.nf/modules/{module}/{sub}/tests/nextflow.config`
  - 当前状态：仅有 `runOptions = '--platform linux/amd64'`（缺少 A 的 `-u`，缺少 B 和 C）
  - 改动：将 docker 块从：
    ```groovy
    docker {
        enabled = true
        runOptions = '--platform linux/amd64'
    }
    ```
    改为：
    ```groovy
    docker {
        enabled = true
        userEmulation = true
        runOptions = '-u $(id -u):$(id -g) --platform linux/amd64'
        fixOwnership = true
    }
    ```
  - 11 个文件：
    1. `bio.nf/modules/flye/tests/nextflow.config`
    2. `bio.nf/modules/ecc_finder/map_ont/tests/nextflow.config`
    3. `bio.nf/modules/ecc_finder/asm_ont/tests/nextflow.config`
    4. `bio.nf/modules/ecc_finder/map_sr/tests/nextflow.config`
    5. `bio.nf/modules/ecc_finder/asm_sr/tests/nextflow.config`
    6. `bio.nf/modules/cresil/identify_wgls/tests/nextflow.config`
    7. `bio.nf/modules/cresil/trim/tests/nextflow.config`
    8. `bio.nf/modules/cresil/visualize/tests/nextflow.config`
    9. `bio.nf/modules/cresil/annotate/tests/nextflow.config`
    10. `bio.nf/modules/cresil/identify/tests/nextflow.config`
    11. `bio.nf/modules/fastqdl/download/tests/nextflow.config`
- **Acceptance Criteria Addressed**: AC-2
- **Test Requirements**:
  - `programmatic` TR-4.1: 11 个文件逐一 Grep，docker 块内同时匹配 `runOptions`、`userEmulation`、`fixOwnership`
  - `programmatic` TR-4.2: runOptions 值为 `'-u $(id -u):$(id -g) --platform linux/amd64'`（`-u` 在前，`--platform` 在后）
  - `human-judgement` TR-4.3: 抽查 3 个文件确认格式一致
- **Notes**: 11 处改动完全相同

## [ ] Task 5: 3 个 server.config 删除重复 A/B/C 设置
- **Priority**: medium
- **Depends On**: None
- **Description**: 
  - circdna.nf、circrna.nf、isoseq.nf 的 `conf/server.config` 中顶层 docker 块已有 A+C（与 nextflow.config docker profile 重复）
  - 删除 `fixOwnership = true` 行和 `runOptions = '-u $(id -u):$(id -g)'` 行
  - **保留** `enabled = true`（确保 `-profile server` 仍能启用 Docker）
  - 改动后 docker 块：
    ```groovy
    docker {
        enabled = true
    }
    ```
  - 3 个文件：
    1. `circdna.nf/conf/server.config`（docker 块 L178-182，删除 L180-181）
    2. `circrna.nf/conf/server.config`（docker 块 L244-248，删除 L246-247）
    3. `isoseq.nf/conf/server.config`（docker 块 L214-218，删除 L216-217）
- **Acceptance Criteria Addressed**: AC-3
- **Test Requirements**:
  - `programmatic` TR-5.1: 3 个文件的 docker 块内**不**匹配 `runOptions`、`userEmulation`、`fixOwnership`（这三个关键词不应出现）
  - `programmatic` TR-5.2: 3 个文件的 docker 块仍匹配 `enabled = true`
  - `human-judgement` TR-5.3: 确认只删除了 A/C 行，其他配置（singularity、conda 等）未受影响
- **Notes**: 删除后 `-profile server` 单独使用时 Docker 仍启用（enabled=true），但 A+B+C 权限设置需 `-profile docker,server` 组合；server.config 只负责服务器资源/执行器配置

## [ ] Task 6: AGENTS.md 新增"Docker 用户映射 A+B+C 标准设置规范"章节
- **Priority**: high
- **Depends On**: None
- **Description**: 
  - 文件：`/Users/siyangming/nextflow_nf_core/AGENTS.md`（项目根目录）
  - 在适当位置新增一节，内容结构如下：
    ```
    ## Docker 用户映射 A+B+C 标准设置规范

    为彻底避免容器以 root 身份运行导致输出文件无法删除的问题，所有流程必须在**同一标准位置**同时设置三层 Docker 权限方案：

    ### 三层方案

    | 层 | 参数 | 作用 | 时机 |
    |----|------|------|------|
    | A（预防） | `docker.runOptions = '-u $(id -u):$(id -g)'` | 容器启动时切换为宿主机用户 | 启动前 |
    | B（预防·补充） | `docker.userEmulation = true` | 挂载 /etc/passwd，提供用户名 | 启动前 |
    | C（兜底） | `docker.fixOwnership = true` | 任务结束后 chown 修复所有权 | 结束后 |

    三层共存不冲突：A+B 设置同一 UID（Docker 取后者，值相同等价），C 的 chown 在 A+B 生效时为 no-op。

    ### 标准设置位置

    对有 `profiles {}` 块的标准 nf-core 流程：
    - `nextflow.config` → `profiles { docker { <HERE> } }` 块内

    对无 `profiles` 块的（bio.nf 顶层、模块测试）：
    - 配置文件顶层 `docker { <HERE> }` 块内

    ### 标准配置块

    ```groovy
    docker {
        docker.enabled         = true
        docker.runOptions      = '-u $(id -u):$(id -g)'   // A
        docker.userEmulation   = true                      // B
        docker.fixOwnership    = true                      // C
        // ... 其他引擎设置
    }
    ```

    ### 扩展写法（追加其他参数时 `-u` 在最前）

    - emulate_amd64 profile：`docker.runOptions = '-u $(id -u):$(id -g) --platform=linux/amd64'`
    - arm profile：`docker.runOptions = '-u $(id -u):$(id -g) --platform=linux/amd64'`
    - gpu profile：`docker.runOptions = '-u $(id -u):$(id -g) --gpus all'`
    - bio.nf 模块测试：`runOptions = '-u $(id -u):$(id -g) --platform linux/amd64'`

    ### server.config 规则

    有 `conf/server.config` 的流程，server.config 的 docker 块**只保留** `enabled = true`，不设置 A/B/C（避免与 docker profile 重复）。A+B+C 权限设置统一由 `nextflow.config → profiles → docker { }` 提供。服务器使用方式：`-profile docker,server`。

    ### 不采用的方案

    - `userns-remap`（Docker daemon 级别）：需改 daemon 配置，HPC 服务器不允许
    ```
- **Acceptance Criteria Addressed**: AC-4
- **Test Requirements**:
  - `programmatic` TR-6.1: `grep -c "Docker 用户映射 A+B+C 标准设置规范" AGENTS.md` 返回 1
  - `human-judgement` TR-6.2: 人工检查该节覆盖：三层方案表格、标准位置（两类）、标准配置块、扩展写法、server.config 补充、不采用方案
- **Notes**: 写在根目录 AGENTS.md（全局通用），不是 circdna.nf/AGENTS.md

## [ ] Task 7: 跨流程批量一致性验证
- **Priority**: medium
- **Depends On**: Task 1, Task 2, Task 3, Task 4, Task 5, Task 6
- **Description**: 
  - 一次性批量检查所有流程在标准位置的 A+B+C 三参数是否齐全
  - 用 git diff 验证 emulate_amd64/arm/gpu profile 未被修改
- **Acceptance Criteria Addressed**: AC-5, NFR-1
- **Test Requirements**:
  - `programmatic` TR-7.1: 生成汇总表：
    ```
    circdna.nf  [docker profile]: A+B+C ✅  [server.config]: enabled-only ✅
    circrna.nf  [docker profile]: A+B+C ✅  [server.config]: enabled-only ✅
    isoseq.nf   [docker profile]: A+B+C ✅  [server.config]: enabled-only ✅
    fetchngs.nf [docker profile]: A+B+C ✅  [server.config]: N/A
    riboseq.nf  [docker profile]: A+B+C ✅  [server.config]: N/A
    nanoseq.nf  [docker profile]: A+B+C ✅  [server.config]: N/A
    rnaseq      [docker profile]: A+B+C ✅  [server.config]: N/A
    bio.nf 顶层 [docker block]:  A+B+C ✅
    bio.nf 11测试:               11/11 A+B+C ✅
    ```
  - `programmatic` TR-7.2: `git diff` 中不包含 emulate_amd64 / arm / gpu profile 的 runOptions 行改动
  - `human-judgement` TR-7.3: 人工审阅 AGENTS.md 新增章节内容完整

## [ ] Task 8: 仅同步 circdna.nf 修改到 GitHub
- **Priority**: medium
- **Depends On**: Task 7
- **Description**: 
  - 所有 Task 1-7 完成后，只将 circdna.nf 目录下的修改同步到 GitHub
  - 不同步其他流程（circrna.nf、isoseq.nf、fetchngs.nf、riboseq.nf、nanoseq.nf、rnaseq、bio.nf）的修改
  - 不同步 AGENTS.md 的修改（除非用户另行要求）
  - 操作步骤：
    1. `git add circdna.nf/`（仅暂存 circdna.nf 目录）
    2. 确认 `git status` 中 circdna.nf 的改动文件：`circdna.nf/nextflow.config`（加 B+C）、`circdna.nf/conf/server.config`（删 A+C）
    3. 提交：`git commit -m "feat(docker): add A+B+C user mapping to docker profile, remove duplicate from server.config"`
    4. 推送：`git push origin master`
- **Acceptance Criteria Addressed**: Goal 4
- **Test Requirements**:
  - `programmatic` TR-8.1: `git log --oneline -1` 显示最新 commit 只包含 circdna.nf 路径的文件
  - `programmatic` TR-8.2: `git diff --name-only HEAD~1 HEAD` 输出中只有 circdna.nf/ 开头的文件
  - `human-judgement` TR-8.3: 确认其他流程的修改仍在本地工作区（`git status` 显示未暂存）
- **Notes**: 只同步 circdna.nf 是用户明确要求；其他流程的修改保留在本地，后续按需单独同步

# Task Dependencies
- Task 1-6 互相独立，可并行执行
- Task 7 必须在 Task 1-6 全部完成后执行
- Task 8 必须在 Task 7 完成后执行
