# 同步 circdna.nf 修改到 GitHub + 容器权限审计报告 Spec

## Why

本次对话中对 circdna.nf 做了路径同步和 Docker 权限修复，需要将这些修改提交到 GitHub。同时需要审计当前目录下其他流程的容器权限解决方案，形成报告供参考。

## What Changes

- 仅提交本次对话中 circdna.nf 的 5 个修改文件到 GitHub
- 不提交其他流程（circrna.nf、bio.nf、nanoseq.nf）的修改
- 生成一份容器权限解决方案审计报告

## Impact

- Affected code: circdna.nf 仓库（SERVER_RUN_GUIDE.md, AGENTS.md, conf/server.config, errors.txt, samplesheets/data_issues.txt）
- 其他流程修改保留在工作区，不提交

## ADDED Requirements

### Requirement: 选择性提交 circdna.nf 修改

系统 SHALL 仅提交本次对话中的 5 个文件到 GitHub，不包含其他修改。

#### Scenario: 提交指定文件

- **WHEN** 执行 git 提交
- **THEN** 仅包含 SERVER_RUN_GUIDE.md, AGENTS.md, conf/server.config, errors.txt, samplesheets/data_issues.txt
- **AND** 不包含其他流程的任何修改

### Requirement: 容器权限审计报告

系统 SHALL 生成一份报告，记录当前目录下所有 Nextflow 流程的 Docker 容器权限配置方案。

#### Scenario: 报告包含所有流程

- **WHEN** 用户查看报告
- **THEN** 报告包含 circdna.nf, circrna.nf, bio.nf, nanoseq.nf, isoseq.nf, fetchngs.nf, riboseq.nf, rnaseq 的容器权限配置
- **AND** 标注每个流程是否配置了 `-u $(id -u):$(id -g)` 参数
