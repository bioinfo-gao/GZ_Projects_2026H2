# 关键操作记录 — 项目13/14 资源并跑（2026-07-09）

> 📌 **本文件横跨两个项目**：内容同时涉及项目13（13_Ellen_knockin_wgs）和项目14
> （14_geneedit_lats12_wgs）的资源共享决策。**权威原件在本路径**
> （`13_Ellen_knockin_wgs/logs/key_operations_log.md`，git 跟踪）；
> `14_geneedit_lats12_wgs/logs/sLk_of_key_operations_log_in_project13.md` 是指向本文件的
> **软链接**（symlink，非独立副本；文件名前缀 `sLk_of_..._in_project13` 就是"这是软链接、
> 原件在项目13"的自解释命名，IDE文件浏览器看不出symlink箭头时也能一眼认出），修改请在
> 本文件（权威原件）进行——软链接会自动同步。

## 背景：发现项目14 Study A 提速空间

项目14（Jinpeng Ruan）Study A sarek somatic 跑了约10小时，比对(BWAMEM2_MEM)只完成34/72子任务。
排查发现两处配置问题叠加：

1. `local_resources.config` 的 `executor.queueSize=2` 是**全局并发上限**（不是"2个样本并行"），
   导致72个区间子任务大部分时间只有1-2个真正在跑。
2. `BWAMEM2_MEM` 配置内存申报 50GB/任务，但**实测真实用量仅约18GB**——nextflow 按申报内存核算
   并发席位，严重高估导致本可并发的任务被挤到串行。

**处理**：

```bash
# 1) 精调配置：BWAMEM2_MEM 内存 50GB→24GB（贴近实测），queueSize 2→3
#    （3×16=48线程，仍在 CLAUDE.md 56线程上限内）
# 2) 优雅停止（非强杀）：
tmux send-keys -t p14_sarek_A C-c
# 确认 nextflow 干净退出（"Execution complete -- Goodbye"），无残留 bwa-mem2/java 进程
# 3) 明确带 -resume 重启（不走脚本默认"先跑一次不带resume"逻辑，避免丢弃已完成子任务缓存）：
nextflow run nf-core/sarek ... -resume
```

**效果验证**：缓存命中38/72子任务（比中断前的34还多，关机瞬间又完成几个，零工作丢失）；
重启后确认3个bwa-mem2真实并发运行，各~18.8GB。提交：`3274090`。

**严谨的提速评估**（★ 用两次运行各自的 `execution_trace_*.txt` 真实提交时间戳算出，不是
对话中的即时印象）：两次运行各有独立的 trace 文件（nextflow 每次启动生成一份），可精确
对比 BWAMEM2_MEM 子任务的**提交时间戳**（`submit` 列）：

| | 完成子任务数 | 提交时间跨度 | 速率 |
| :--- | :---: | :---: | :---: |
| **修复前**（`execution_trace_2026-07-08_10-37-12.txt`，全部 COMPLETED） | 38 | 8.06小时 | 4.72 个/小时 |
| **修复后**（`execution_trace_2026-07-08_21-38-32.txt`，仅统计真正新跑的 COMPLETED，排除38个CACHED） | 34 | 6.32小时 | 5.38 个/小时 |
| **实际提速** | | | **仅 1.14倍** |

（对比窗口经核实完全在项目13启动〈2026-07-09 10:29〉之前，不受项目13并跑影响，是干净对比。）

**⚠️ 如实记录：提速远小于预期**。实测并发确实从1-2个提升到3个真实并行的bwa-mem2进程
（`ps aux` 直接验证），但吞吐只提升了14%，不是配置调整前口头讨论时暗示的"明显加速"。
可能原因（未逐一验证，留待以后有需要再查）：
1. BWAMEM2_MEM 按染色体区间分片，各分片大小/耗时差异大（如chr1/chr2这类大染色体远慢于
   小contig）；修复后剩余队列里恰好轮到较慢的区间，部分抵消了并发数增益。
2. 真正瓶颈可能是磁盘I/O或内存带宽等硬件层面，而非queueSize这个逻辑并发槽位数——3个
   bwa-mem2同时读写大文件时会互相争抢物理I/O，"3倍并发槽位"不等于"3倍实际吞吐"。

**结论**：这次精调方向没错（消除了一个不合理的人为节流），但不是万能药；如果以后再遇到
类似"进度慢"的情况，应优先怀疑硬件I/O瓶颈，而不是默认调大queueSize就能线性提速。

重启操作的具体命令已存成脚本 `../../14_geneedit_lats12_wgs/scripts/study_A/A2b_restart_with_resume.sh`
（原先只是临时Bash命令，未第一时间存成可复查脚本，事后补记）。

## 决定：趁项目14收尾、系统利用率低，并跑项目13

项目14 Study A 比对阶段（最吃线程的部分）已跑完，进入去重(MarkDuplicates,单线程为主)+TIDDIT的
收尾阶段。实测：整机负载 ~5（64核），内存仅用17GB/125GB，判断有空间同时跑项目13的正式6样。

**关键认识**：两个 nextflow 实例各自独立调度、互不知晓对方——CLAUDE.md 的"全服务器≤56线程"
预算必须由人工在两个项目间显式切分，不能让两边都按各自的满配（各48线程）跑，否则理论峰值叠加
会超限。

**处理**：

1. 建混合参考（GRCm39+RAGH+MTTH+CD1A，64 contigs）——CPU/内存需求低，与项目14无冲突，直接执行。
2. **新建降配专用配置** `scripts/local_resources_concurrent_with_proj14.config`：
   - 24线程/80GB上限（原满配版48线程/108GB保留不动，供项目14结束后切回）
   - BWAMEM2_MEM 12cpu/24GB（queueSize=2，2×12=24线程=上限）
   - 内存值同样按项目14实测的~18GB经验校准，不用原先偏高的申报值
   - 文件内注明：项目14跑完后应切回 `local_resources.config` 满配版，不要长期用降配版
3. `2_run_sarek.sh` 增加 `SAREK_CONFIG` 环境变量覆盖支持，**显式写入 tmux 命令字符串内部**
   （env var 不会自动传入新 tmux session，是已记录过的坑）。
4. 生成6样 samplesheet，启动：
   ```bash
   SAREK_CONFIG="scripts/local_resources_concurrent_with_proj14.config" \
     bash scripts/2_run_sarek.sh scripts/samplesheet_full.csv
   ```
5. 3分钟早期失败检测通过，无报错。启动后核实整机合计负载 7.46、内存31GB/125GB，安全。
   提交：`2e34190`。

## 并跑到底加快了没有?（2026-07-10 实测评估,用两条 trace 的 realtime×%cpu 算核·小时）

**结论:加快了,而且是"有效"的加快——但完全靠错峰调度,不是并跑本身。**

核心证据(核占用 = Σ realtime×%cpu/100 ÷ 墙钟跨度):

| 窗口 | 平均占用核数 | 含义 |
| :--- | :---: | :--- |
| 项目14 **独跑满配**(比对期,07-08 21:39→07-09 10:29,12.8h) | **28.2 核** | 一个作业在比对期就独自吃满 28 核预算,毫无余量 |
| 项目14 **并跑期**(变异调用尾巴,07-09 10:29→今,32.9h) | **9.9 核** | Mutect2/去重 CPU 轻,独跑会白白空出~18 核 |
| 项目13 **并跑期**(比对) | **~22 核** | 正好填进项目14 空出的~18 核 |
| 两作业**合计**(重叠期) | **~29 核** | ≈28 核预算天花板,机器没闲着 |

**为什么有效**:项目13 恰在项目14 从"比对(CPU重、独自占满28核)"切到"变异调用(CPU轻、只~10核)"
的**拐点**启动。项目13 至今~630 核·小时的比对进度,**几乎全是在项目14 尾巴的空闲产能里"免费"跑出来的**。
串行的话这630核·小时只能等项目14 全跑完再排队。

**量级**:若项目13 能在项目14 剩余尾巴(还有~52个 Mutect2 区块)跑完前结束,则近乎省下项目13 整段墙钟
→ 这一对接近 **1.5–2× 吞吐**。若项目14 先完、项目13 还剩大段 CPU 重活要独跑,收益缩水。

**前提与边界(否则毫无意义)**:收益 100% 来自**两作业处于互补资源画像**(一个CPU重比对 vs 一个CPU轻变异调用)。
- 反例:若在项目14 **比对期**(独自已占满28核)启动项目13,两者抢同一批核 → 纯时间切片、零净收益,各慢一倍。
- 当前两作业都已进入低CPU尾段(项目14 Mutect2 ~10核 + 项目13 转入 MarkDup),CPU 已非瓶颈,内存成为要盯的点。

详细每步资源画像见 `../scripts/sarek_wgs_perstep_timing_and_resources.md`。

## 进展快照(2026-07-11 12:00)

两条 tmux(`p14_sarek_A`/`ellen_sarek`)均 UP,`.nextflow.log` 2 分钟内刷新,健康。

| 项目 | 状态 | 剩余 | ETA |
| :--- | :---: | :---: | :---: |
| **14** somatic(`--tools mutect2,tiddit`) | Mutect2 **72/90(80%)**,去重/SV/深度全完成 | 18 Mutect2 + 合并/污染估计/FilterMutectCalls/MultiQC | **~14h → 07-12 凌晨** |
| **13** germline(`--tools tiddit` 仅结构变异,无 SNV 调用) | 去重6/6、mosdepth6/6、samtools5/6、**TIDDIT 4/6** | 2 TIDDIT + SVDB合并/stats/MultiQC | **~3–5h → 07-11 下午–傍晚** |

- **关键澄清**:P13 只跑 TIDDIT(检测敲入整合位点的结构变异),不做 germline SNV 调用——故其剩余量远小于 07-10 的过度估计,会**先于 P14 完成**。
- P14 Mutect2 近端速率回升到 **1.9 个/h**(07-10 为 1.5,P13 CPU 压力退去后加快);受 `queueSize=3` 限,提不了更高。
- P13 收尾后会生成本轮真正的 6 样本 MultiQC(07-10 归档掉的是旧的)。

## 输出目录清理(2026-07-10)

`output_results/` 里混有 07-06/07 一次旧运行(以 RAGH_153 单样本为主)的残留(旧 `multiqc_report.html`、
旧 RAGH_153 的 markdup/stats/mosdepth/CRAM、manta/tiddit 变异结果、汇总 CSV,共 192 文件 9.6G),违反
"标准输出目录只放本轮 canonical 结果"。已按 mtime<本轮启动(07-09 10:29)判定并整体搬到
`OLD/superseded_run_0706_output/`(含 `README_archive.md` 说明),`OLD/` 已加入 `.gitignore`。
当前 6 样本正式运行未受影响(管线读 work-dir,搬完两 tmux 仍 UP、日志正常),跑到相应阶段会重新生成 RAGH_153 各输出。

## 后续跟踪要点

- **项目14跑完（Study A + B）后，项目13务必切回 `local_resources.config` 满配版**，
  否则会一直背降配速度跑，不必要地慢。
- 项目13的 Cas9/iHPV 混合参考构建体序列（核准中）不影响当前进度——当前跑的是模式B主流程
  （比对/去重/TIDDIT/拷贝数/整合位点），不依赖那批序列。
- 两个项目的 tmux 会话：项目14=`p14_sarek_A`，项目13=`ellen_sarek`。

---

*Zhen Gao, PhD, Principal Bioinformatics Scientist, Athenomics*
