# 关键操作记录 — 项目13/14 资源并跑（2026-07-09）

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

## 后续跟踪要点

- **项目14跑完（Study A + B）后，项目13务必切回 `local_resources.config` 满配版**，
  否则会一直背降配速度跑，不必要地慢。
- 项目13的 Cas9/iHPV 混合参考构建体序列（核准中）不影响当前进度——当前跑的是模式B主流程
  （比对/去重/TIDDIT/拷贝数/整合位点），不依赖那批序列。
- 两个项目的 tmux 会话：项目14=`p14_sarek_A`，项目13=`ellen_sarek`。

---
*Zhen Gao, PhD, Principal Bioinformatics Scientist, Athenomics*
