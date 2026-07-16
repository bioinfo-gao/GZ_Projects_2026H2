# 事故复盘教材：一次"提速"如何烧掉 15 小时 —— 提速反噬 / `-resume` 失效 / 监控缺位

- **文档日期（创建，immutable）**: 2026-07-16
- **项目**: `16_Wenliang_Pan_wgs`（Homo sapiens GRCh38，germline WGS，2 samples）
- **作者**: Zhen Gao, PhD — Athenomics
- **性质**: 内部教材 / 事故复盘（client 不可见）
- **一句话**: 为省 3.5 h 的 markdup，加了一个 flag，净亏 ~15 h —— 而其中 **6.4 h 是在无人监控下白烧的**。

> **本文档的定位**：项目内 `analysis_plan_0715.md` §7 已有三条「反面教训」的分散记录；本文把三个问题
> **串成一条因果链**并抽出可复用判据，供教学与复盘。**三者不是并列的三个 bug，而是一个 bug 被两层
> 防线失效逐级放大的过程。**

---

## 0. 三个问题的因果链（先看这个）

```
   [问题 1] 加 --use_gatk_spark markduplicates 想提速
        │        ↓ 触发单线程 ELC 空转（本该 ~1h 内被发现并回退）
        │
   [问题 3] 监控缺位 ← 第一层防线失效
        │        ↓ 11h 无人巡检 → ELC 空转 6.4h 无人察觉（放大 ~6×）
        │
   [问题 2] -resume 救不回来 ← 第二层防线失效
        │        ↓ 以为"回退就好"，实际 24 个比对只命中 1 个 cache，再赔 4.4h
        ↓
     净亏 ~15 h，提速收益 = 0
```

**教学要点**：问题 1 是**技术判断错误**（可原谅，谁都会踩新坑）；问题 3 是**纪律错误**（不可原谅，规则早已写死在
CLAUDE.md 里）；问题 2 是**认知错误**（对工具的错误心智模型，以为 `-resume` 是"后悔药"）。
**真正致命的是 3，不是 1。** 一个被及时发现的坑只值 1 小时；一个夜里无人看管的坑值 6 小时。

---

## 1. 问题一：为什么极其慢 —— `--use_gatk_spark markduplicates` 的双重反噬

### 1.1 动机与操作

- **动机**（合理）：常规 `GATK4_MarkDuplicates` 近单线程，单 WGS 样 ~3.5 h，是流程明显瓶颈。
  `MarkDuplicatesSpark` 号称多线程，看起来是免费的提速。
- **操作**：2026-07-15 **23:50:53** 给 `2_run_sarek.sh` 加 `--use_gatk_spark markduplicates`，带 `-resume` 重启。
  （首跑 run1 于 **22:09:03** 启动，当时**运行完全正常**。）

### 1.2 反噬 A：Spark markdup 不产 metrics → sarek 自动补一个单线程 ELC

**这是最反直觉的一点，也是本次的主杀手。**

- `MarkDuplicatesSpark` **不输出 duplicate metrics**。sarek 为了 MultiQC 的 QC 报告完整，会**自动追加**一个
  独立进程 `GATK4_ESTIMATELIBRARYCOMPLEXITY`（ELC）来补算这个指标。
- ELC 是**单线程**的，且在高深度 WGS 上呈**内存 thrashing + 持续减速**（实测 ~40 groups/min 且越跑越慢）。
- **实测证据**（`.nextflow.log.1`，本项目实际日志）：

  | 事件 | 时间戳 | 证据 |
  | :--- | :---: | :--- |
  | ELC task 129 提交 | `Jul-16 04:30:46` | `Submitted process > ...GATK4_ESTIMATELIBRARYCOMPLEXITY` |
  | ELC task 130 提交 | `Jul-16 04:30:57` | 同上（每样本一个） |
  | 两个 ELC 被终止 | `Jul-16 10:52:30` | `Task completed > TaskHandler[id: 129; ... exit: 129]` |

  **⚠ 关键细节：`exit: 129` = `128 + 1` = SIGHUP —— 它不是"跑完了"，是 10:52 我回退时被信号杀掉的。**
  nextflow 日志里 `Task completed` 只表示"任务结束了"，**不表示成功**。看到 `Task completed` 就以为
  正常完成，是读 nextflow 日志的常见误读 —— **必须看 `exit:` 字段。**
  **ELC 实际运行 04:30:46 → 10:52:30 = 6 h 21 m（6.36 h），且到死都没算完。**

### 1.3 反噬 B：`queueSize=2` 让 ELC 从"慢"升级为"全线阻塞"

- 本项目 `queueSize=2`（受 bwa-mem2 每 task ~40 GB index 的 RAM 限制，见 plan §5）。
- 两个 ELC 各占一个 slot → **两个 slot 全被单线程任务霸占** → HaplotypeCaller 等一切下游**完全无法启动**。
- 结果：**load 从 ~32 掉到个位数，~53 个线程全程空转 6.4 小时**，机器看起来"活着"，实际什么也没产出。

### 1.4 判据（可复用）

| 判据 | 内容 |
| :--- | :--- |
| **默认选择** | germline WGS 一律用**常规 `GATK4_MARKDUPLICATES`**（proj4/proj13 已验证）：单线程但**有界、可预期**、metrics inline、无 ELC。 |
| **何时才考虑 Spark** | 真有 Spark/GPU 集群，**且已在小样本上验证过端到端耗时**。本机单节点场景：不要。 |
| **通用原则** | **"多线程"不等于"更快"**——要看它是否连带引入了别的单线程步骤。评估一个提速 flag，必须看**整条流程的墙钟**，而不是被优化的那一步。 |
| **⛔ 铁律** | **绝不对健康运行中的 pipeline 做未经验证的"提速"改动。** run1 从 22:09 起跑得好好的；这次"优化"是纯粹的自伤。 |

---

## 2. 问题二：为什么 `-resume` 救不回大部分损失

### 2.1 错误的心智模型

回退时的想法是：「去掉 flag，`-resume` 重跑，比对结果还在 cache 里，只重做 markdup 就行。」
**这个想法基于一个错误假设：`-resume` 是"后悔药"，能让我免费撤销一个决定。**

### 2.2 实际发生了什么

**实测证据**（run3 = 去 Spark 后的重启，`.nextflow.log`）：`Cached process` 共 **31** 条，构成如下 ——

| 被复用的 process | 数量 |
| :--- | :---: |
| `TABIX_BGZIPTABIX_INTERVAL_SPLIT` | 21 |
| `FASTQC` | 2 |
| `FASTP` | 2 |
| `TABIX_BGZIPTABIX_INTERVAL_COMBINED` / `ENSEMBLVEP_DOWNLOAD` / `CREATE_INTERVALS_BED` / `CNVKIT_REFERENCE` / `CNVKIT_ANTITARGET` | 各 1 |
| **`BWAMEM2_MEM`** | **1** ← **全部 24 个比对，只命中 1 个** |

**31 这个数字极具欺骗性**：乍看"复用了 31 个任务，`-resume` 生效了"，但其中 21 个是几秒钟的 TABIX 小任务，
**真正值钱的 24 个 BWAMEM2 比对只复用了 1 个**（且那 1 个还是 run1 的遗留产物）。
**教训：看 cache 命中要看"命中了哪些"，不能只看数量。**

### 2.3 根因（已用 `.command.sh` diff 实证）

`--use_gatk_spark` **不只影响 markdup 进程，它会改写 `BWAMEM2_MEM` 自己的命令行**：

- Spark markdup 需要 **queryname-sorted** 输入 → sarek 给 BWAMEM2 发出 `samtools sort -n -@ 16`
- 常规 markdup 需要 coordinate-sorted → `samtools sort -@ 16`（**无 `-n`**）

`diff` 新旧 work dir 的 `.command.sh` 实证：run2(Spark) = `sort -n`，run1/run3(常规) = `sort`。

> **nextflow 的 cache 键是 task 的 script 文本 hash。**
> script 文本变 → hash 变 → **cache 全废**。而且**开、关各废一次**：
> 加 Spark 时废掉 run1 已完成的比对，去 Spark 时又废掉 run2 刚重做的比对。

### 2.4 真实代价核算

| 阶段 | 耗时 | 说明 |
| :--- | :---: | :--- |
| 加 Spark 后重比对 24 chunk | ~4.7 h | run1 的比对成果作废 |
| ELC 空转 | **6.4 h** | 04:30:46 → 10:52:30，且未算完 |
| 去 Spark 后再重比对 23 chunk | ~4.4 h | run2 的比对成果又作废 |
| **合计** | **≈ 15.4 h** | **换来的提速：0** |

### 2.5 判据（可复用）

| 判据 | 内容 |
| :--- | :--- |
| **核心规则** | **`-resume` 只对「不改 script 文本」的改动便宜。** |
| **安全的改动** | `cpus` / `memory` / `queueSize` 等 **directive**（不进 `.command.sh`）→ cache 保留。 |
| **昂贵的改动** | 任何会进入 `.command.sh` 的 flag：**Spark 开关、`ext.args`、`--aligner`、工具参数** → 该 process **及其下游**全部重跑。 |
| **⚠ 隐蔽处** | **一个 flag 可能改写"看似无关"的上游进程的命令行**（本例：markdup 的 flag 改了 alignment 的 `sort -n`）。**不要假设影响面只限于你改的那个工具。** |
| **操作纪律** | **改参数重启前，先 `diff` 新旧 work dir 的 `.command.sh`**，确认影响面；别假设"上游能复用"。 |
| **心智模型修正** | **`-resume` 不是后悔药，是断点续跑。** 它只在"环境/参数没变、只是中断了"时救你；用它来"撤销一个决定"通常是全价重跑。 |

---

## 3. 问题三：监控缺位 —— 让上面两个问题从 1 小时放大成 15 小时

**这一条是元教训，也是本次最该被记住的一条。**

### 3.1 事实（`logs/phase1_monitor.log` 实证）

```
[1 x30s] cached=0 submitted=0
[2 x30s] cached=4 submitted=1
[3 x30s] cached=7 submitted=23
[4 x30s] cached=7 submitted=23
[5 x30s] cached=7 submitted=23
monitor done            ← mtime 07-15 23:53，跑了 3 分钟就自行退出
```

| 时刻 | 事件 |
| :--- | :--- |
| `07-15 23:50:53` | 带 Spark 重启 |
| `07-15 23:53` | Phase-1 监控打出 `monitor done`，**退出** |
| `07-16 04:30:46` | ELC 开始空转（load 从 ~32 掉到个位数） |
| `07-16 10:52:30` | 人工发现并回退 |

- **监控空窗：`23:53 → 10:52` = 11.0 小时。**
- **ELC 无人察觉：`04:30 → 10:52` = 6.4 小时。**
- 按 CLAUDE.md 要求的 Phase-2（10–20 min 巡检），**最晚 05:00 就该发现** → **多烧约 6 小时**。

### 3.2 违反的两条 MANDATORY 规则

1. **「Phase 1 通过 ≠ 监控完成」**
   CLAUDE.md 明写：Phase 1 之后要 *"switch to longer-interval background checks (every 10–20 min)
   **until completion**"*。**Phase 2 从未启动。** 把"启动没炸"当成了"监控已完成"。
2. **「长任务监控必须在 tmux」**
   Phase-1 监控是用 Bash `run_in_background` 跑的 —— **与 agent 会话绑定，会话结束即消失**。
   `scripts/3_monitor.sh` 只是**一次性** `tail` 查询脚本，不是常驻循环。项目里**从头到尾没有任何常驻监控进程**。

### 3.3 更难看的是：工具早就有，只是没用

通用看门狗 `facilities/Server/nextflow_watchdog.sh` **早已存在**，`/corun` skill 也写明了用法，
**本项目压根没调用**。**不是没工具，是没纪律。**

### 3.4 但"记得跑看门狗"也不够 —— 老看门狗抓不到这个事故

**这不是替自己开脱，而是决定了修复方案。** 通用看门狗只有三条触发规则：

| 规则 | ELC 事故时的状态 | 是否命中 |
| :--- | :---: | :---: |
| `SESSION_END` | tmux 会话好好活着 | ❌ |
| `LOW_MEM` | 内存充足（ELC 不吃内存） | ❌ |
| `HEARTBEAT` | 90 min 时 `exit 12` 等**上层复检** | ⚠ 夜里没有上层 = 等于没监控 |

**ELC 空转属于「会话活着 + 内存充足 + 单线程磨洋工」的「健康但病态慢」状态 —— 三条规则一条都不命中。**
老看门狗只会在日志行里默默记下 `nflog+300m`，却**不告警**。

> **教学要点：监控规则必须覆盖"没死但也没干活"这一类状态。**
> 绝大多数监控只查"死没死"（进程在否、内存够否），而**真实世界里最贵的故障往往是"活着但空转"** ——
> 因为它不触发任何告警，可以安静地烧一整夜。

### 3.5 修复（2026-07-16 14:01 已生效并实测）

新增 `scripts/10_watchdog.sh`，**常驻 tmux 会话 `pan_watch`**：

| 特性 | 说明 |
| :--- | :--- |
| **永不退出** | 自我循环，不靠上层复检（老看门狗 heartbeat 就 `exit` 的设计在无人夜里失效） |
| `POLL=600s` | 10 min 巡检，符合 CLAUDE.md Phase-2 要求 |
| **规则 A `STALLED`** | `.nextflow.log` 停滞 > 45 min → 告警（流程没在推进） |
| **规则 B `LOW_CPU_UTILIZATION`** | load 连续 30 min < 10 → 告警。**这条就是照着本次 ELC 事故的特征写的** |
| 保留 | `SESSION_END` / `LOW_MEM` |
| **只告警不 kill** | 杀不杀是需要判断的运维决策，留给人（CLAUDE.md 规则） |
| 告警分离落盘 | `logs/watchdog_ALERTS.log`（异常）与 `logs/watchdog.log`（流水）分开，一眼看异常 |

**⚠ 已实测两条新规则确实触发**（造 3 小时假停滞 nflog + 压低 load 阈值，两条如期报出，且 STALLED 去重不刷屏）。
**"看门狗启动了" ≠ "看门狗会告警"** —— 这正是本次事故的同型错误（把"启动"当"生效"），故必须实测。

### 3.6 判据（可复用）

| 判据 | 内容 |
| :--- | :--- |
| **纪律** | **`tmux ls` 里每个跑长任务的会话，都必须有一个配对的常驻监控会话。** 本项目：`pan_wgs`(sarek) + `pan_watch`(看门狗) + `pan_down`(下游编排器)。 |
| **Phase 1 的定位** | Phase-1（3 min）只回答"启动没炸"，**它退出的那一刻必须无缝接上 Phase-2**，否则等于没监控。 |
| **监控进程的宿主** | **必须 tmux**。Bash `run_in_background` 与 agent 会话同生共死 —— agent 一断，监控就没了，而长任务还在跑。 |
| **规则设计** | 必须有**停滞检测**：(A) 日志停滞；(B) **低 CPU 占用 = 活着但空转**。只查"死没死"的监控抓不到最贵的故障。 |
| **验证** | 监控本身也要**实测告警会触发**，不能只看它启动了。 |

---

## 4. 三条可迁移的通用教训（脱离本项目也成立）

1. **「优化」的默认答案是"不做"。**
   对一个**正在健康运行**的长任务做未验证的提速改动，期望收益（省 3.5 h）远小于风险敞口（赔 15 h）。
   run1 在 22:09 跑得好好的 —— 不动它，本项目今天上午就该出结果了。

2. **工具的"便宜"是有前提的，前提没吃透就会当成免费。**
   `-resume` 便宜的前提是 script 文本不变；`MarkDuplicatesSpark` 快的前提是不引入 ELC。
   **把前提当成无条件保证，是这两个坑的共同结构。**

3. **监控不是"运维琐事"，它是损失的乘数。**
   同一个 bug，有监控值 1 小时，没监控值 6 小时。**投入产出比最高的工程动作，往往不是让流程更快，
   而是让故障更早被看见。** 而监控规则必须覆盖"活着但空转"——最贵的故障不触发任何传统告警。

---

## 5. 本文档的教训已落地到哪里

| 落地位置 | 内容 |
| :--- | :--- |
| `analysis_plan_0715.md` §7「反面教训 1/2/3」 | 项目内的分散记录（本文是其串联版） |
| `scripts/10_watchdog.sh` + tmux `pan_watch` | 问题三的**代码修复**（已跑，已实测） |
| `analysis_plan_0715.md` §6 运行顺序 | 起 sarek 后**必须**起 `pan_watch`，写进标准流程 |
| `/wgs` skill 实测教训段 | 问题一、二（Spark 陷阱 + `-resume` 判据）→ 影响所有 WGS 项目 |
| memory `reference_sarek_spark_markdup_elc_trap` | 问题一、二 |
| memory `feedback_phase2_monitoring_tmux_mandatory` | 问题三 |
| `facilities/Server/nextflow_pipeline/` | 跨项目共享版（见该目录同名主题文档） |

---

*Zhen Gao, PhD, Principal Bioinformatics Scientist, Athenomics — 内部教材，非 client 交付物。*
