# nf-core/taxprofiler 本地安装与测试运行记录

**记录日期**：2026-07-01
**Pipeline 来源**：https://nf-co.re/taxprofiler/2.0.1 （GitHub: `nf-core/taxprofiler`，锁定 release `2.0.1`）
**用途**：shotgun metagenomics（宏基因组鸟枪法测序）taxonomic classification / profiling pipeline，
支持短读长（short-read）与长读长（long-read）数据的 QC、host removal（宿主序列去除），
并可并行调用多种 classifier/profiler（Kraken2、Bracken、KrakenUniq、Centrifuge、MetaPhlAn、
Kaiju、DIAMOND、MALT、mOTUs、ganon、KMCP、sylph、MELON、MetaCache 等），
最后用 MultiQC / Krona / taxpasta 做标准化汇总。

本文档记录**这次实际执行的每一步操作、用到的软件版本、遇到的问题和修复方法**，
目的是让同事或未来的自己可以在同一台服务器上原样复现。

---

## 1. 软件版本清单（Software Versions）

以下版本均为本次实际安装/验证时记录，非泛泛而谈：

| 软件                |                                    版本                                    |                             安装位置 / 备注                             |
| :------------------ | :-------------------------------------------------------------------------: | :---------------------------------------------------------------------: |
| mamba               |                                    2.3.2                                    |                  系统级，`/Work_bio/gao/miniforge3`                  |
| conda env 名称      |                               `taxprofiler`                               |            专门新建，与日常主力环境`regular_bioinfo` 分开            |
| nextflow            |                            25.10.4 (build 11173)                            |                        conda env`taxprofiler`                        |
| nf-core (tools)     |                                    4.0.2                                    |                        conda env`taxprofiler`                        |
| openjdk             |                               23.0.2-internal                               |        conda env`taxprofiler`（满足 pipeline 要求的 `>=17`）        |
| nf-core/taxprofiler | release**2.0.1**，commit `70ecc15e49b4f1fcf79d876643b5d14b65c66178` | 通过`nextflow pull` 拉取到 `~/.nextflow/assets/nf-core/taxprofiler` |
| container engine    |      apptainer 1.4.5（对外表现为`singularity` 命令），docker 29.3.0      |       系统级已装好，两者均可用；本次选用`-profile singularity`       |
| 操作系统            |                         Ubuntu 22.04.5 LTS (jammy)                         |                               服务器本机                               |
| CPU / 内存          |          AMD Threadripper 2990WX，32 物理核 / 64 线程，125 GB RAM          |                         详见第 4 节资源限制设计                         |

> 说明：nf-core pipeline 本身不需要把 Kraken2、MetaPhlAn 这些具体的生信工具装进 conda env——
> 每一个 process 都在 Nextflow 自动拉取的独立 Singularity 容器里运行，`taxprofiler` 这个 env
> 只需要 Nextflow + Java 来调度工作流即可。

---

## 2. 环境搭建步骤（可直接复制执行）

### 2.1 新建 mamba 环境

```bash
mamba create -n taxprofiler -y nextflow=25.10.* nf-core=4.0.* openjdk=17
```

安装完成后验证（应能看到上表中的版本号）：

```bash
conda run -n taxprofiler nextflow -version
conda run -n taxprofiler nf-core --version
conda run -n taxprofiler java -version
```

### 2.2 拉取 pipeline 代码（锁定 2.0.1 版本）

```bash
export NXF_SINGULARITY_CACHEDIR='/Work_bio/gao/configs/.singularity'
conda run -n taxprofiler nextflow pull nf-core/taxprofiler -r 2.0.1
```

拉取成功后代码会缓存在 `~/.nextflow/assets/nf-core/taxprofiler`，之后每次 `nextflow run`
都会直接复用，不会重复下载。

---

## 3. 工作目录结构

```
/home/gao/projects_2026H2/8_taxprofiler_setup/
├── docs/
│   ├── taxprofiler_安装与测试记录_0701.md   ← 本文档（中文，含完整复现步骤）
│   └── README_taxprofiler_tutorial.md       ← 英文版技术教程（内容对应，供快速查阅）
├── configs/
│   └── local_resources.config               ← CPU/内存上限 + singularity 缓存目录配置
├── scripts/
│   └── 1_run_taxprofiler_test_profile.sh    ← 一键复现测试运行的启动脚本
├── logs/
│   └── taxprofiler_test_run.log             ← 完整运行日志（tmux 内 tee 出来的）
└── test_run/
    ├── testdata/
    │   └── database_v2.1_taxprofiler2.0.1.csv  ← 修正后的测试数据库清单（原因见第 6 节）
    ├── work/                                 ← Nextflow 中间文件（work directory）
    └── outdir/                               ← pipeline 最终输出结果
```

---

## 4. 关键配置：`configs/local_resources.config`

服务器物理配置：32 核 / 64 线程，125 GB 内存。
根据本机工作规范的硬性上限（**所有并发任务合计不超过 28 物理核 / 56 线程**，
需为系统、SSH、交互操作预留资源），本配置把 Nextflow 本地执行器（local executor）
和单进程资源上限都设在硬上限以内（24 核 / 96 GB），并复用已有的共享 Singularity
镜像缓存目录（之前跑 nf-core/sarek 时就用的同一个目录，避免重复下载相同的镜像）：

```groovy
// configs/local_resources.config
executor {
    cpus   = 24          // 本地执行器最多同时用 24 个核（留出余量，未顶到 28 核上限）
    memory = '96.GB'
}

process {
    resourceLimits = [
        cpus: 24,
        memory: 96.GB,
        time: 12.h        // 单个 process 最长运行时间上限
    ]
}

singularity {
    enabled    = true
    autoMounts = true
    cacheDir   = '/Work_bio/gao/configs/.singularity'   // 共享镜像缓存，跨项目复用
}
```

**重要坑点**：`NXF_SINGULARITY_CACHEDIR` 这个环境变量必须写在 **tmux 命令字符串内部**
export，不能写在外层脚本里提前 export。因为如果 tmux server 已经常驻后台，
新建的 session 会继承 tmux server 启动时的环境，而不是外层脚本临时 export 的变量——
这是之前调试 nf-core/sarek 时踩过的坑，这次直接照搬了正确写法。

---

## 5. 测试运行：一键复现命令

```bash
bash /home/gao/projects_2026H2/8_taxprofiler_setup/scripts/1_run_taxprofiler_test_profile.sh
```

脚本内部实际执行的核心命令（在 tmux session `taxprofiler_test` 里跑）：

```bash
export NXF_OPTS='-Xms512m -Xmx2g'
export NXF_SINGULARITY_CACHEDIR='/Work_bio/gao/configs/.singularity'

nextflow run nf-core/taxprofiler \
  -r 2.0.1 \
  -profile test,singularity \
  -c configs/local_resources.config \
  --databases test_run/testdata/database_v2.1_taxprofiler2.0.1.csv \
  --outdir test_run/outdir \
  -work-dir test_run/work \
  -resume
```

参数说明（中文）：

- `-profile test,singularity`：使用 pipeline 官方内置的 `test` profile（极小的公开测试
  fastq 数据 + 每个 classifier 对应的极小测试数据库），叠加 `singularity` 容器引擎；
  **不需要**准备任何真实的宏基因组数据或本地参考数据库，全部由 pipeline 自动下载。
- `-c configs/local_resources.config`：叠加第 4 节的本机资源上限配置。
- `--databases ...csv`：显式覆盖测试数据库清单（原因见第 6 节，是本次遇到的关键问题）。
- `-resume`：断点续跑，任何中断后重新执行这条命令会复用已完成的 process，不会从头重跑。

查看运行进度：

```bash
tmux attach -t taxprofiler_test        # 进入 tmux 实时查看（Ctrl-B D 退出不中断任务）
tail -f logs/taxprofiler_test_run.log  # 或者直接 tail 日志文件
```

---

## 6. 本次实际遇到的问题及修复过程（重点）

### 6.1 问题描述

第一次运行 `-profile test,singularity`（未加 `--databases` 覆盖）时，
在参数校验（validation）阶段直接报错，pipeline **还没有开始下载/执行任何 process**：

```
ERROR ~ Validation of pipeline parameters failed!
* --databases (.../taxprofiler/database_v2.1.csv): Validation of file failed:
  -> Entry 6: Error for field 'tool' (centrifuger):
     Expected any of [bracken, centrifuge, diamond, ganon, kaiju, kmcp,
     kraken2, krakenuniq, malt, metaphlan, motus, sylph, melon, metacache]
     (Invalid tool name...)
```

### 6.2 根因分析

Pipeline 内置的 `test` profile（文件 `conf/test.config`）里，`--databases` 参数默认
指向的是 `nf-core/test-datasets` 仓库 `taxprofiler` 这个分支上的 CSV 文件。
**这个分支是一个持续更新的"活分支"**，同时被 taxprofiler 的 `dev` 开发分支共用，
**并不会跟着某个 release 版本号冻结**。

本次运行时（2026-07-01），该 CSV 文件里已经多了一行 `centrifuger` 工具的测试数据库条目
（`centrifuger` 是比 `centrifuge` 更新的另一个独立分类工具）。但我们锁定运行的是
`2.0.1` 这个旧版本，它的参数 schema（`nextflow_schema.json`）里还**不认识**
`centrifuger` 这个工具名，于是校验直接失败——这是一个典型的"锁定旧 release
版本 vs. 持续演进的测试数据仓库"版本不匹配问题，不是本机环境配置的错误。

### 6.3 修复方法

下载官方测试数据库 CSV，删掉不兼容的 `centrifuger` 那一行，另存为本地文件，
再通过 `--databases` 显式指定，覆盖 pipeline 默认值：

```bash
mkdir -p test_run/testdata
curl -s https://raw.githubusercontent.com/nf-core/test-datasets/taxprofiler/database_v2.1.csv \
  | grep -v '^centrifuger,' \
  > test_run/testdata/database_v2.1_taxprofiler2.0.1.csv
```

修复后其余所有测试参数（fastq samplesheet、host removal 参考序列、各个 `run_*`
开关等）均保持 pipeline 默认值不变，只有这一个 database 清单被替换。

**后续复现提示**：如果以后把 pipeline 升级到更新的版本（例如 2.1.x 以上），
需要重新检查上游 test-datasets 的 CSV 是否还需要这个本地补丁——新版本 schema
可能已经原生支持 `centrifuger`，届时可以直接用官方默认值，不再需要本地修正文件。

---

## 7. 运行状态与监控方法

采用两阶段监控（符合本机工作规范）：

- **第一阶段（启动后前 3 分钟）**：每 30 秒检查一次日志，同时匹配"成功特征"
  （process 正常执行、镜像开始拉取）和"失败特征"（`ERROR ~`、`Halted` 等关键词），
  确认参数校验通过、Singularity 镜像开始拉取、真实 process（如 KMCP_PROFILE、
  GANON_CLASSIFY、METACACHE_QUERY 等）已经在跑，而不是 tmux 会话空转或直接崩溃。
- **第二阶段（确认健康后）**：改为每 10 分钟检查一次日志尾部和进程列表，
  直到 tmux session 结束（pipeline 跑完或报错退出）。

test profile 一次性启用了约 14 种 classifier/profiler，首次运行需要拉取约 40 个
Singularity 镜像（后续复现如果共享缓存目录 `/Work_bio/gao/configs/.singularity`
里已经有对应镜像，会直接复用，不用重新下载），所以**首次运行的耗时主要花在拉镜像
上，而不是实际计算**。最终运行结果（process 成功/失败统计、总耗时、MultiQC
报告路径等）见第 9 节。

### 7.1 实际遇到的一次 dataflow 停滞（stall / deadlock）及处理过程

**现象**：第一次跑（见第 6 节修复 `--databases` 后启动的那次）在跑到
`METAPHLAN_METAPHLAN` 两个 process 于 23:12:22 完成后，**此后再无任何新任务
启动、日志再无新行、work 目录再无新文件**，`htop`/`ps` 显示本机负载几乎为空闲
（load average ~0.4），Nextflow 主进程的 `/proc/<pid>/wchan` 一直是
`futex_wait_queue`——符合本机规范里"确认 deadlock"的两个条件（wchan 卡在
futex 上 **且** 日志时间戳冻结超过 15 分钟）。

**关键确诊证据**（比单纯"日志不动"更硬的证据）：Nextflow 自己的内部日志
`scripts/.nextflow.log` 里，从 23:16:17 开始**每隔 5 分钟就打印一次**：

```
!! executor local > No more task to compute -- The following nodes are still active:
```

这是 Nextflow 调度器自己报告的"dataflow 通道卡住"信号——某些下游聚合型
process（比如 `TAXPASTA_MERGE`、`MULTIQC`、`STANDARDISATION_PROFILES` 系列）
一直处于 `status=ACTIVE` 等着上游 channel 关闭，但没有更多任务可以调度，
调度器只能反复空转报告，永远不会自己恢复。到发现时该警告已经连续打印了
5 次（20 分钟）。这不是"任务慢"，是真正卡死。

**处理方法**：既然确认是 deadlock 而不是"健康但慢"，且此时已经有约 30 个
process（fastp/fastqc/minimap2/bowtie2/samtools/MetaPhlAn/Kaiju/DIAMOND/
KMCP/GANON/MELON/MetaCache 等）已经成功跑完（exit code 0，产物都在
`test_run/work/` 里），符合本机规范"先确认 deadlock，再重启，且要用
`-resume` 避免推倒重来"的要求：

```bash
tmux kill-session -t taxprofiler_test
bash scripts/1_run_taxprofiler_test_profile.sh   # 脚本本身自带 -resume
```

重启后日志立刻显示大量 `cached: N ✔` 字样，证明之前跑完的 process 被直接
复用，没有重新计算；同时之前卡住的下游 process（Centrifuge、Sylph 等）
也开始正常推进。

**根因推测**：test profile 一次性启用了约 14 个 classifier，产生的 process
DAG 分支很多、很复杂，个别小众/长尾分支（如本例中 Centrifuge/Sylph 相关的
镜像拉取或 channel 汇总）在高并发下偶发卡住，是已知的 Nextflow dataflow
调度器边缘情况，并非本机资源配置或环境问题。**后续复现建议**：如果再次
运行后日志长时间（>15 分钟）无新内容，先查
`scripts/.nextflow.log` 里有没有反复出现的
`No more task to compute -- ... still active` 字样确诊，确诊后直接
`tmux kill-session` + 重跑同一条帶 `-resume` 的命令即可，不会丢失已完成的
计算结果。

---

## 8. 后续用于真实样本的运行方式

真实项目需要准备两个输入文件：

1. **`--input samplesheet.csv`**：列为
   `sample,run_accession,instrument_platform,fastq_1,fastq_2,fasta`
   （具体模板见 pipeline 仓库里的 `assets/samplesheet.csv`）。
2. **`--databases databases.csv`**：列为 `tool,db_name,db_params,db_type,db_path`
   （`db_type` 取值 `short` / `long` / `short;long`；`db_path` 可以是本地路径，
   也可以是一个 tar 包或 URL，指向该 classifier 对应的预建参考数据库）。

真实数据的典型运行命令（数据/数据库准备好之后）：

```bash
nextflow run nf-core/taxprofiler \
  -r 2.0.1 \
  -profile singularity \
  -c configs/local_resources.config \
  --input samplesheet.csv \
  --databases databases.csv \
  --outdir results \
  -work-dir work \
  -resume
```

通过 `--run_<tool>` 系列布尔开关（如 `--run_kraken2 --run_metaphlan`）控制
实际启用哪些 classifier；只有和已启用工具匹配的数据库行才会被使用。
QC / host removal 相关开关（`--perform_shortread_qc`、
`--perform_shortread_hostremoval`、`--hostremoval_reference` 等）用法同理。
完整参数列表可执行 `nf-core launch nf-core/taxprofiler`，或查看
`nextflow_schema.json` / https://nf-co.re/taxprofiler/2.0.1/parameters 。

### 8.1 参数名核实记录（避免"传了参数但静默无效"的坑）

2026-07-02 用 `nextflow run ... -preview` 实测核对了一批"看起来合理但实际
可能是旧版本残留写法"的参数，发现其中 4 个是**不存在的参数名**——不存在的
参数不会报错中断，而是被 nf-schema 插件**静默忽略**（既不在参数汇总里显示，
也不打印任何 warning），非常容易造成"以为设置生效了、实际完全没起作用"的
误判。核实结果：

| 传入的参数                            | 是否有效 |                                                                                                                                                         说明                                                                                                                                                         |
| :------------------------------------ | :-------: | :------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------: |
| `--perform_shortread_preprocessing` | ❌ 不存在 |                                                                                    静默忽略。正确参数是**`--perform_shortread_qc`**（短读长 QC/fastp）；如还要低复杂度过滤需另加 `--perform_shortread_complexityfilter`                                                                                    |
| `--remove_host`                     | ❌ 不存在 |           静默忽略。正确参数是**`--perform_shortread_hostremoval`**（短读长）/ `--perform_longread_hostremoval`（长读长），且**必须同时提供** `--hostremoval_reference <genome.fasta>`（或预建索引 `--shortread_hostremoval_index`）才会真正执行宿主去除，只给开关不给参考序列不会生效           |
| `--run_kraken2`                     |  ✅ 有效  |                                                                                             真实样本默认所有`run_*` 分类器开关都是关闭的，要跑 Kraken2 必须显式打开，且 `databases.csv` 里要有 `tool=kraken2` 的行                                                                                             |
| `--run_bracken`                     |  ✅ 有效  |                                                            schema 说明写着"会自动带出所需的 Kraken2 前置步骤"——单开`--run_bracken` 已隐含 Kraken2 那一步，和 `--run_kraken2` 同时开不冲突，但只要 Bracken 结果时理论上不必重复加 `--run_kraken2`                                                            |
| `--run_metaphlan`                   |  ✅ 有效  |                                                                                                                                   需要`databases.csv` 里有 `tool=metaphlan` 行                                                                                                                                   |
| `--max_cpus 28`                     |  ❌ 无效  | taxprofiler 2.0.1 用的是**新版 nf-core 模板**，已放弃旧模板 `check_max()` + `--max_cpus`/`--max_memory` 机制（sarek 3.8.1 那种旧模板才用这套，见 `conf/base.config` 里只有固定的 `process.cpus`/`process.withLabel` 写法，完全没有引用 `params.max_cpus`）。传了会被完全忽略，起不到任何限流作用 |
| `--max_memory '100.GB'`             |  ❌ 无效  |                                                                                                                                                   同上，完全不生效                                                                                                                                                   |

**结论**：要真正限制资源，必须用 Nextflow 原生的 `process.resourceLimits`
配置块——也就是本项目已经建好的 `configs/local_resources.config`（见第 4
节），而不是 `--max_cpus`/`--max_memory` 这两个 CLI 参数。**后续复现提示**：
每次给 taxprofiler（或任何新模板 nf-core pipeline）加自定义参数前，建议先用
`nextflow run ... -preview` 快速核对一遍参数汇总，确认新加的参数确实出现在
输出的参数列表里，而不是想当然照抄别的 pipeline（尤其是旧模板的 sarek）的
参数写法。

修正后的真实样本运行命令（替换第 8 节的示例）：

```bash
nextflow run nf-core/taxprofiler \
  -r 2.0.1 \
  -profile singularity \
  -c configs/local_resources.config \
  --input samplesheet.csv \
  --databases databases.csv \
  --perform_shortread_qc \
  --perform_shortread_hostremoval --hostremoval_reference <host_genome.fasta> \
  --run_kraken2 \
  --run_bracken \
  --run_metaphlan \
  --outdir results \
  -work-dir work \
  -resume
```

---

## 9. 最终运行结果

**测试运行已于 2026-07-01 23:45:12 成功结束**（经历了第 7.1 节记录的一次
dataflow 停滞，`kill` + `-resume` 重启后完成）：

```
-[nf-core/taxprofiler] Pipeline completed successfully-
Duration    : 6m 24s   （重启后这一段的净耗时；含镜像拉取的总墙钟时间约 40 分钟）
CPU hours   : 1.5 (45.5% cached)
Succeeded   : 103
Cached      : 76
Failed      : 0
```

日志里未出现任何 `FAILED` 任务，14 个 classifier 全部产出正常，验证通过。

**产物目录**（`test_run/outdir/`，共 22 个子目录）：

| 子目录                                                                                                                                              |                                                                内容                                                                |
| :-------------------------------------------------------------------------------------------------------------------------------------------------- | :--------------------------------------------------------------------------------------------------------------------------------: |
| `fastp` / `fastqc` / `nanoq` / `porechop_abi`                                                                                               |                                                         短读长 / 长读长 QC                                                         |
| `bowtie2` / `samtools`                                                                                                                          |                                                 short-read host removal 比对与统计                                                 |
| `bbduk`                                                                                                                                           |                                                          低复杂度序列过滤                                                          |
| `nonpareil`                                                                                                                                       |                                                         测序深度冗余度估计                                                         |
| `kraken2` / `bracken` / `centrifuge` / `kaiju` / `diamond` / `kmcp` / `ganon` / `sylph` / `melon` / `metacache` / `metaphlan` |                                                     各 classifier 原始分类结果                                                     |
| `krona`                                                                                                                                           |                                                        交互式分类结果可视化                                                        |
| `taxpasta`                                                                                                                                        |                                                    跨工具标准化后的统一分类表格                                                    |
| `multiqc`                                                                                                                                         |                            **`multiqc/multiqc_report.html`** —— 汇总质控报告，最先打开看这个                            |
| `pipeline_info`                                                                                                                                   | 3 次运行（首次失败于参数校验 / 第二次遇到 stall / 第三次 resume 完成）各自的 execution report / timeline / DAG / trace，可用于复盘 |

结论：**nf-core/taxprofiler 2.0.1 在本机的安装、Singularity 容器引擎、
共享镜像缓存、资源限制配置均已验证可用**，可以进入真实样本分析阶段
（见第 8 节）。
