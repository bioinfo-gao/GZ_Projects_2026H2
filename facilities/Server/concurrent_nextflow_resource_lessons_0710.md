# 双 nf-core/sarek 管线同时运行——资源实测与运维教训（2026-07-10）

> 场景：系统资源有限，但待跑工作很多。本文记录一次**两个 sarek WGS 管线在同一台服务器上并行运行**
> 的全过程实测（14.5 小时逐点监控）、资源使用曲线、以及由此固化的运维教训。
> 对应可复用 skill：`/corun`（`~/.claude/skills/corun/SKILL.md`）。
> 服务器背景见同目录 `hardware.ipynb`、`1_add_SSD_swap.sh`、`1A_add_SSD_swap_explain.ipynb`。

---

## 1. 硬件与两个作业

| 项 | 值 |
| :--- | :---: |
| 服务器 | AMD Threadripper 2990WX，32 物理核 / 64 线程，125 GiB RAM |
| Swap | 2 GB 系统 swapfile `/swapfile`（在 **OS NVMe 系统盘** Samsung 970 EVO 1TB 上，即 `/`，prio -2 后用）+ **64 GB SSD swapfile**（`/mnt/ex_8T_SSD`，在 **7.3T SATA 数据盘** Desk SSD 上，prio 10 先用），`swappiness=10`。两块 swap 均在 SSD、无一落在机械盘；优先级刻意"先用大容量 SATA SSD、把 NVMe 带宽留给 OS" |
| 策略上限 | **稳态 ≤ 28 物理核 / 56 线程**；繁忙期可短时冲 **30 核 / 60 线程**(不可长期维持,须快回落),永留 ≥2核/4线程给 OS/SSH/交互 |
| 作业 A（proj14） | `p14_sarek_A`，2026-07-08 21:38 启动；Study A 体细胞 WGS（GRCm39，6 样本，Mutect2 + TIDDIT） |
| 作业 B（proj13） | `ellen_sarek`，2026-07-09 10:29 启动；6 样本敲入 WGS（混合参考，bwa-mem2 + TIDDIT） |
| 重叠窗口 | 两作业自 07-09 10:29 起同时运行；本文监控覆盖其中 07-09 21:15 → 07-10 11:51 约 14.5 h |

**关键：两个作业各自用一份"降配" `resourceLimits` 配置**，把单作业的 executor/process 核数与内存压到
全机上限的一半档位，使两者之和仍落在 28 核 / 56 线程、且内存峰值可控。proj13 明确使用
`scripts/local_resources_concurrent_with_proj14.config`（与 proj14 并跑专用），proj14 结束后再换回独立满配版。

---

## 2. 资源实测（304 个采样点，每 180 s 一次）

| 指标 | 最小 | 最大 | 稳态区间 |
| :--- | :---: | :---: | :---: |
| MemAvailable | **32.6 GiB** | 86.0 GiB | 大部分时间 52–72 GiB |
| Swap 已用 | 7.5 GiB | **14.6 GiB** | 基线 ~7.5 GiB，稳态 ~9 GiB |

**峰值时刻（唯一一次真正吃紧）：07-10 04:28–04:40**
- 04:28 available 跌到 32.6 GiB，04:31 swap 冲到 **14.6 GiB**（比基线多换出 ~7 GiB）。
- 无人干预，**约 15 分钟内自行恢复**：04:52 available 回到 68 GiB、swap 回落至 11 GiB，之后继续下降到 ~9 GiB 稳态。
- 成因研判：两作业**同时撞上内存重阶段**——proj13 的比对/去重（bwa-mem2 ~19 GB RSS + samtools sort 缓冲）叠加 proj14 的 Mutect2（多个 JVM）。
- **即便峰值，125 GiB 里仍剩 32 GiB available，离 OOM 很远**；SSD swap 作为缓冲吸收了这波瞬时压力。

**每步实测耗时 + CPU核/峰值RSS 的量化表**见 `13_Ellen_knockin_wgs/scripts/sarek_wgs_perstep_timing_and_resources.md`
（本节只讲整机内存曲线与错峰逻辑，逐步数字不在此重复）。

**两作业的资源画像不同，正是能共存的原因：**
- **proj13 = 比对阶段**（bwa-mem2）：CPU + 内存双重，单进程 RSS ~19 GB，重。全程 14.5 h 都卡在 BWA-MEM2（72 分块只推进了 ~40 个），因为和 proj14 抢 CPU + 人源区（HTT / CD1 基因簇）本身难比对。
- **proj14 = 变异调用阶段**（Mutect2）：许多小 CPU 任务（90 区块 × 肿瘤-正常对），单任务轻，但总量大、耗时长（>15 h 都在变异调用）。
- 一个"重而少"、一个"轻而多"，错峰使用 CPU，内存峰值只在两者偶然同时进入重阶段时短暂出现。

---

## 3. 监控方法（事件驱动看门狗，非被动等成功）

用一个后台 bash 看门狗（`facilities/Server/nextflow_watchdog.sh`，本次实例见 skill），每 180 s 巡检一次：
- 两个 tmux session 是否存活（`tmux has-session`）；
- 各自 `.nextflow.log` 距今多久刷新（判断是否卡住）；
- `free -m` 的 available 与 swap。

**只在发生实际事件时退出并回报**（退出即触发 agent 复检）：
1. 某个 session 消失 → 立刻读该作业日志尾部，判定**正常完成**（`Pipeline completed successfully` 等标记）还是**失败**；
2. available < 2.5 GB 连续两次 → OOM 前兆，提前告警；
3. 每 90 分钟健康心跳 → 即使一切正常也周期性确认，防止监控自己悄悄死掉而不自知。

失败/完成检测是**实时**的，不受心跳间隔影响；心跳只决定"没事时多久确认一次"。稳定数小时后可把心跳从
45 min 放宽到 90 min，减少无谓唤醒。

---

## 4. 运维教训（这才是长期价值）

1. **看"available"不是"free"。** `free` 纯空闲常年只有几百 MB～1 GB（Linux 把内存拿去做 cache 了），
   真正指标是 **available**（含可回收 cache）。本次 free 常显示 <1.5 GB 但 available 稳在 50 GB+，毫无压力。

2. **健康但慢 ≠ 故障，绝不因慢而杀。** proj13 比对 14.5 h 只走了 40/72 分块——极慢，但日志每 0–4 min 都在刷新、0 失败。
   慢是"共享 28 核预算"的**正常代价**，不是问题。杀掉健康但慢的作业会毁掉数小时不可逆的中间结果。

3. **瞬时内存峰值会自行恢复，别恐慌、别抢先 drop_caches / swapoff。** 04:30 那次 available→32 GB、swap→14.6 GB，
   15 min 自愈。此时若手贱 `echo 3 > drop_caches` 会逼两个作业重读参考、更慢；`swapoff` 更可能把冷页硬拉回内存触发 OOM。

4. **SSD swap 是并跑的安全垫，不是浪费。** 64 GB SSD swap（`swappiness=10`）在峰值吸收了 ~7 GB 瞬时溢出，
   让作业平稳越过重阶段而非被 OOM 打死。`swappiness=10` 保证平时不乱换出、只在真紧时才用。

5. **每个并跑作业用"降配" `resourceLimits`，让核数之和 ≤ 28 物理核。** 不要让单作业按满机配置跑两份——
   核数会翻倍超限、内存峰值不可预测。nf-core 现代模板用 `process.resourceLimits`（不再是 `--max_cpus/--max_memory`）。

6. **错峰调度胜过硬挤。** 让先启动的作业领先进入下游阶段（proj14 已在变异调用），后启动的（proj13）还在比对——
   两者资源画像互补时冲突最小。**先完成的那个会释放 CPU，让另一个自动提速**，无需人工干预。

7. **失败恢复必须 `-resume`，且要显式。** nextflow 重启若不带 `-resume` 会丢弃 work-dir 里已缓存的完成任务、白跑数小时。
   proj14 备了显式 `-resume` 重启脚本（`A2b_restart_with_resume.sh`）；proj13 的运行脚本首次调用**不带** resume，
   故整段 session 若被 OOM 打没、需手动重启时**不能直接重跑该脚本**，要用带 `-resume` 的显式命令。

8. **判定"真卡死"要看 `.nextflow.log` 的 `No more task to compute ... still active`，而非仅"日志看着不动"。**
   本次 proj14 画面长时间停在同一屏（Mutect2 是最长的一环、下游在等它），但 trace 的 COMPLETED 持续增长、无该死锁信号 →
   真在推进。只有"进程 wchan 卡在 futex + 日志 mtime 冻结 >15 min + 出现该信号"三者齐备才算确认死锁，届时 kill + `-resume` 才安全。

9. **验证进度用 trace 的 COMPLETED 计数，不看单屏。** nextflow 进度屏只显示部分行，长作业的关键阶段可能滚出可视区。
   `grep -c -w COMPLETED execution_trace_*.txt` 是最可靠的"到底跑了多少"。

---

## 5. 并跑两个 sarek 的经验配额（供后续排程参考）

- **两个 sarek（一个在比对、一个在变异调用）在 125 GiB / 32C 机上可稳定共存**，峰值内存 ~93 GiB（含 cache），
  留 32 GiB available 余量，配 64 GB SSD swap 作缓冲。
- **代价是各自变慢约一倍**（CPU 对半分）。若赶时间且内存允许，宁可**排队串行**跑满配、也不要为"看起来在并行"而两个都拖慢。
- **不要叠第三个重活。** 两个 sarek 已用满 28 核预算；第三个并发会把三者都拖入 CPU 饥饿，且内存峰值可能顶穿。
- 重 work-dir 放 `/mnt/ex_8T_SSD`（SSD、可弃、~3.9 TB），不要放 `/Work_bio`（HDD、慢、空间紧）。

---

*记录人：Zhen Gao / Athenomics。数据来源：2026-07-09-10 双 sarek 并跑实测监控日志（304 采样点）。*
