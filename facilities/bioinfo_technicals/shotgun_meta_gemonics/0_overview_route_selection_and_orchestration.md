# Shotgun 宏基因组：三条路线的选型与编排总览

> 本目录（`shotgun_meta_gemonics/`）四篇的**入口 / 索引**。回答一个问题：**一个 shotgun 宏基因组项目该跑哪几条路线、什么顺序、谁 gated 在谁后面。**
> 组合不是新分析类型——"全套 shotgun = 组成 + 功能 + MAG"本就是标准全服务交付，只是三条已有路线**按序编排**。

---

## 三条正交路线（一张表看清）

| 路线 | 回答 | 工具 / 管线 | 组装? | 速度 | 对应 skill | 原理篇 |
| :--- | :--- | :--- | :---: | :---: | :--- | :--- |
| **组成 taxonomy** | 谁在这里、各占多少 | nf-core/taxprofiler（Kraken2+Bracken / MetaPhlAn） | 否 | 快（分钟/样本） | `/taxnom` | `1_introduction_and_assembly_or_not` |
| **功能 function** | 它们能干什么（pathway/gene family 丰度、KO/EC） | HUMAnN 3（mag_biobakery 原生） | 否 | 慢（~8h/样本，diamond 主导） | `/tax-functional-humann` | `3_taxonomy_and_functional_profiling_HUMAnN` |
| **基因组 MAG** | 能否拼出基因组草图、它们是谁、编码哪些基因 | nf-core/mag | 是 | 重（GTDB r226 库 ~271GB、组装/分箱耗内存） | `/tax-assembly-mag` | `2_MAG` |

**三者正交、互补，不是二选一**：组成相近可功能不同（株差异），物种不同可功能冗余；MAG 给基因组级证据、抓库里没有的新菌。

---

## 决策：这个项目跑哪几条？

```
送样单 "std analysis"（标准 shotgun）
   └─ 业界惯例 = 组成 + 多样性 + 功能（assembly-free）→ 跑 /taxnom + /tax-functional-humann
                                                          （这就完整覆盖 std 下单内容）
   └─ MAG 是 advanced/增值（单独报价）→ 仅在以下情况加做：
        · 客户/你要基因组级证据、发现新/未培养菌、pangenome
        · 想吃满机器（有空闲算力）
      且数据够深：per-sample ≥10 Gbp 才稳出高质量 MAG；
      深度薄（如 ~4 Gbp/样本）→ 必须 group co-assembly（同臂合并）提升低丰度菌信号，
      MAG 天然是"二期增值"定位，不是回答主问题的入口。
```

- **组成 + 功能** 覆盖绝大多数干预对照类问题（IF vs AL、疾病 vs 健康…）的主诉求。
- **MAG** gated 在"值不值得 + 深度够不够 + 有没有算力"三问之后，且通常等前两条出完主交付、机器空出来再放行（组装/GTDB-Tk 比 HUMAnN 更吃 CPU/RAM，别叠加）。

---

## 推荐编排顺序（含 gating）

```
Phase 1  ── /taxnom (taxprofiler)  → QC/去宿主 + Kraken2/Bracken + MetaPhlAn + taxpasta 合表
   │            └─ 去宿主 unmapped reads 直接复用给下一步（别重复去宿主）
   │            └─ 下游 R：α/β 多样性 + 差异丰度 → 出主交付
   │
Phase 1b ── /tax-functional-humann (HUMAnN)  ← 复用 Phase 1 的 clean reads
   │            └─ diamond ~8h/样本、3 并行；功能表 join/renorm/regroup → 补进同一交付
   │            └─ 【与 Phase 1 一起构成 std analysis 完整交付】
   │
   ▼  （用户 review Phase 1 主交付 + 决定是否放行 MAG）
   │
Phase 2  ── /tax-assembly-mag (nf-core/mag)  ← 仅获批 + 深度够 + 有算力才跑
                └─ 先备好 GTDB r226（~271GB）；group co-assembly；GTDB 分类 + bin 丰度
```

**资源不叠加原则**：HUMAnN（54 线程满载）跑着时不要同时启 MAG 组装；等 functional 收尾、机器空出再上 MAG。多管线并跑参考 `/corun` playbook。

---

## 关键前提与共用资源

- **去宿主一次、多路复用**：Phase 1 taxprofiler 的 `*.unmapped_{1,2}.fastq.gz` 同时喂 HUMAnN 和 MAG，别各自重复去宿主。
- **数据量以 fastp 实测为准**：gzip 大小反推会低估近一倍（宏基因组高复杂度、压缩率低），直接影响 diamond 耗时和 MAG 深度判断。
- **环境**：taxprofiler/mag 走 nextflow（`mag_biobakery` env 调度 + 容器）；HUMAnN 走 `mag_biobakery` 原生二进制。R 下游包装进 `regular_bioinfo`。
- **worked examples**：`8_taxprofiler_setup`（taxprofiler）、`15_mag_setup`（mag 5.4.2）、`17_Daniel_Mendes_gut_metagenomics`（三条全跑的实战：组成+功能已交付，MAG 二期）。

---

*索引 · 本目录其余三篇是各路线的原理/操作详解；三条路线的可执行 setup 分别是 `/taxnom`、`/tax-functional-humann`、`/tax-assembly-mag`。*
