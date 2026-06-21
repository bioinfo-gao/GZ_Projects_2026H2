# 关键操作记录 — nf-core/sarek 安装与测试（2026-06-21）

## 环境探查（确认已有条件，避免重复安装）

```bash
which nextflow nf-core java docker singularity apptainer conda mamba
nextflow -version     # 25.10.4，来自 /Work_bio/gao/configs/.conda/envs/regular_bioinfo
docker ps              # 确认免 sudo 可用
df -h /home /Work_bio   # /home 707G 可用，/Work_bio 1.2T 可用
nproc; free -h          # 64 核，125G 内存
```

结论：nextflow / docker / singularity / apptainer 均已就位，只缺 `nf-core` CLI 工具本身。

## 1. 安装 nf-core/tools

```bash
/home/gao/.conda/envs/regular_bioinfo/bin/pip install nf-core
# -> 成功安装 nf-core-4.0.2（及依赖）
nf-core --version                       # nf-core, version 4.0.2
nf-core pipelines list | grep sarek     # sarek 3.8.1（最新稳定版，4 个月前发布）
```

## 2. 目录与配置搭建

```
/home/gao/projects_2026H2/4_wgs_human_immu/
├── scripts/1_run_sarek_test_profile.sh   # 启动脚本
├── configs/local_resources.config        # 本机资源上限配置
├── test_run/{outdir,work}                 # 测试运行的输出与工作目录
├── docs/README_sarek_tutorial.md          # 完整教程文档
└── logs/{sarek_test_run.log, key_operations_log.md}
```

`configs/local_resources.config`：local executor 限制 48 CPU / 96GB（总资源 64C/125G，预留给系统）。

## 3. 测试运行（官方 test profile，全自动测试数据）

```bash
bash scripts/1_run_sarek_test_profile.sh
# 等价于在 tmux 会话 sarek_test 中执行：
nextflow run nf-core/sarek -r 3.8.1 -profile test,singularity \
  -c configs/local_resources.config \
  --outdir test_run/outdir -work-dir test_run/work -resume
```

**结果：Pipeline completed successfully**，耗时 7m37s，23/23 进程成功（chr22 测试数据集，
germline 单样本：FastQC → BWA-MEM 比对 → MarkDuplicates → BQSR(BaseRecalibrator+ApplyBQSR) →
Strelka germline 变异检测 → bcftools/vcftools 统计 → MultiQC）。

输出验证：
- `test_run/outdir/variant_calling/strelka/test/test.strelka.variants.vcf.gz` 存在且非空
- `test_run/outdir/multiqc/multiqc_report.html` 生成成功
- `test_run/outdir/pipeline_info/execution_trace_*.txt` 显示全部 23 个任务 status=COMPLETED, exit=0

## 4. 调试过程中发现并修复的问题

**问题**：第一次运行时设置了 `export NXF_SINGULARITY_CACHEDIR=...`，但写在调用 `tmux new` 的外层
脚本进程里，而不是 tmux session 内部命令字符串里。因为 tmux server 是常驻后台进程，新建的 session
继承的是 **server 启动时**的环境，不是当前 shell 临时 export 的变量 —— 所以这个 cache 路径没生效，
日志反复出现：

```
WARN: Singularity cache directory has not been defined -- Remote image will be stored in the path:
.../test_run/work/singularity
```

**后果**：本次测试拉取的 ~4.8GB 容器镜像被缓存进了项目内的 `test_run/work/singularity/`，而不是
全局共享缓存目录 `/Work_bio/gao/configs/.singularity/`。

**处理**：
1. 流程仍然跑通（不影响本次测试结果），等流程结束后手动把这批 `.img` 文件 `mv` 到
   `/Work_bio/gao/configs/.singularity/`，避免以后真实样本运行时重新下载这些常用容器
   (fastqc / samtools / strelka / mosdepth / vcftools / gawk 等)。
2. 修正了 `scripts/1_run_sarek_test_profile.sh`：把 `export NXF_OPTS` 和
   `export NXF_SINGULARITY_CACHEDIR` 挪到 tmux 命令字符串内部执行，以后调用此脚本不会再触发这个问题。

这是一个值得记住的通用教训：**用 `tmux new -d -s <name> "<command>"` 启动后台任务时，
所有需要生效的环境变量必须写在 `<command>` 字符串内部，而不能在外层脚本里提前 export**
（除非确定这是 tmux server 第一次启动）。

## 5. 当前磁盘占用

```bash
du -sh /home/gao/projects_2026H2/4_wgs_human_immu          # 222M（容器已移出）
du -sh /Work_bio/gao/configs/.singularity                   # 18G（含本次新增的 sarek 测试容器）
```

## 6. 尚待你确认后才能进行的下一步

见 `docs/README_sarek_tutorial.md` 第 6.4 节 —— "human_immu" 项目具体的免疫相关分析目标
（HLA 分型 / 免疫基因 panel 注释 / TCR-BCR repertoire），决定了真实样本要不要在 sarek 之外
再接其他专用工具。另外需要你提供：真实样本的 fastq 路径、是否有 tumor-normal 配对（somatic）
还是纯 germline、以及预期跑多少个样本（用于评估磁盘与计算时间）。
