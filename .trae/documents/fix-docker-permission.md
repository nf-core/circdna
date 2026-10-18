# 修复 Docker 容器文件权限问题计划

## 问题分析

### 现象
服务器上 `rm -rf work/` 失败，提示"权限不够"，因为 SPAdes 等工具在 Docker 容器内以 root 身份运行，生成的文件归属于 root。

### 根本原因
`circdna.nf/conf/server.config` 的 docker 配置块缺少 `docker.runOptions = '-u $(id -u):$(id -g)'`，导致 Docker 容器以 root 身份运行。`fixOwnership = true` 只在 Nextflow 层面修复 work 目录所有权，但无法修复容器内进程直接写入的中间文件（如 SPAdes 的 `.bin_reads/` 目录）。

### 全流程排查结果

| 流程 | 配置文件 | docker.runOptions 状态 |
|------|---------|----------------------|
| **circdna.nf** | nextflow.config (docker profile) | ✅ 已有 |
| **circdna.nf** | **conf/server.config** | ❌ **缺失** |
| **circrna.nf** | nextflow.config (docker profile) | ✅ 已有 |
| **circrna.nf** | **conf/server.config** | ❌ **缺失** |
| **bio.nf** | **nextflow.config** | ⚠️ 只有 `--platform linux/amd64`，缺少 `-u` |
| **nanoseq.nf** | **conf/base.config (GPU 条件块)** | ⚠️ 只有 `--gpus all`，缺少 `-u` |
| isoseq.nf | conf/server.config | ✅ 已有 |
| fetchngs.nf | nextflow.config (docker profile) | ✅ 已有 |
| riboseq.nf | nextflow.config (docker profile) | ✅ 已有 |
| rnaseq | nextflow.config (docker profile) | ✅ 已有 |

## 当前缓存删除方法

服务器上已有 root 权限的 work 目录无法用普通 `rm` 删除，有以下三种方法：

### 方法 1: sudo rm（需要 sudo 权限）
```bash
cd /data1/users/siyangming/PlanteccDNADB/circdna.nf/
sudo rm -rf work/
```

### 方法 2: Docker 容器删除（无需 sudo，推荐）
```bash
cd /data1/users/siyangming/PlanteccDNADB/circdna.nf/
docker run --rm -v $(pwd):/work -w /work alpine rm -rf work/
```

### 方法 3: nextflow clean（Nextflow 内置命令）
```bash
cd /data1/users/siyangming/PlanteccDNADB/circdna.nf/
nextflow clean -f
```
注意：`nextflow clean -f` 可能仍无法删除 root 权限文件，建议优先使用方法 1 或 2。

## 修改方案

### 1. circdna.nf/conf/server.config (L178-181)

**文件**: `/Users/siyangming/nextflow_nf_core/circdna.nf/conf/server.config`

```groovy
// 修改前
docker {
    enabled      = true
    fixOwnership = true
}

// 修改后
docker {
    enabled      = true
    fixOwnership = true
    runOptions   = '-u $(id -u):$(id -g)'
}
```

### 2. circrna.nf/conf/server.config (L244-247)

**文件**: `/Users/siyangming/nextflow_nf_core/circrna.nf/conf/server.config`

```groovy
// 修改前
docker {
    enabled      = true
    fixOwnership = true
}

// 修改后
docker {
    enabled      = true
    fixOwnership = true
    runOptions   = '-u $(id -u):$(id -g)'
}
```

### 3. bio.nf/nextflow.config (L5-8)

**文件**: `/Users/siyangming/nextflow_nf_core/bio.nf/nextflow.config`

```groovy
// 修改前
docker {
    enabled = true
    runOptions = '--platform linux/amd64'
}

// 修改后
docker {
    enabled = true
    runOptions = '-u $(id -u):$(id -g) --platform linux/amd64'
}
```

### 4. nanoseq.nf/conf/base.config (L11-14)

**文件**: `/Users/siyangming/nextflow_nf_core/nanoseq.nf/conf/base.config`

```groovy
// 修改前
if (params.deepvariant_gpu) {
    docker.runOptions      = '--gpus all'
    singularity.runOptions = '--nv'
}

// 修改后
if (params.deepvariant_gpu) {
    docker.runOptions      = '-u $(id -u):$(id -g) --gpus all'
    singularity.runOptions = '--nv'
}
```

### 5. circdna.nf/SERVER_RUN_GUIDE.md — 添加清理说明

**文件**: `/Users/siyangming/nextflow_nf_core/circdna.nf/SERVER_RUN_GUIDE.md`

在"注意事项"部分添加清理已有 root 权限文件的说明：

```markdown
- **清理旧 work 目录权限问题**: 如果之前运行未配置 `-u` 参数，work 目录中可能存在 root 权限文件无法删除。使用以下命令清理：
  ```bash
  sudo rm -rf work/
  # 或无需 sudo（通过 Docker 容器删除）：
  docker run --rm -v $(pwd):/work -w /work alpine rm -rf work/
  ```
```

## 不修改的文件（已正确配置）

- circdna.nf/nextflow.config (docker/emulate_amd64/gpu profiles) — 已有 `-u`
- circrna.nf/nextflow.config (docker profile) — 已有 `-u`
- isoseq.nf/conf/server.config — 已有 `-u`
- fetchngs.nf/nextflow.config — 已有 `-u`
- riboseq.nf/nextflow.config — 已有 `-u`
- rnaseq/nextflow.config — 已有 `-u`

## 验证步骤

1. 检查 circdna.nf/conf/server.config 是否包含 `runOptions = '-u $(id -u):$(id -g)'`
2. 检查 circrna.nf/conf/server.config 是否包含 `runOptions = '-u $(id -u):$(id -g)'`
3. 检查 bio.nf/nextflow.config 的 runOptions 是否同时包含 `-u $(id -u):$(id -g)` 和 `--platform linux/amd64`
4. 检查 nanoseq.nf/conf/base.config 的 GPU 条件块是否包含 `-u $(id -u):$(id -g)`
5. 检查 SERVER_RUN_GUIDE.md 是否添加了清理说明
