# 项目 14 — 基因编辑肿瘤细胞 + Lats1/2 Hippo 小鼠 WGS 分析方案

- **Project**: 14_geneedit_lats12_wgs
- **Plan Date（创建）**: 2026-07-07　·　文件名固定为创建日，永不更名（研究时间线追踪）
- **文档更新记录**:
  - 2026-07-12 — §1 样本表补入每样本 R1/R2 实测文件大小 + 数据集总量（462 GiB / 463 G）
  - 2026-07-13 — Study B germline 出现 load 峰值152(CNNScoreVariants超订)，实测+根因+配置修正另记于
    `docs/resource_pressure_incident_0713.md`；新增负载修正终版配置
    `scripts/0_common/local_resources_studyB_load_adjusted_FINAL.config`(executor.cpus=48 + CNN cpus=9，
    目标 load ≤56)，取代 `full_machine.config`(后者标记 SUPERSEDED，保留作历史证物)
  - 2026-07-15 — §1 新增「B1TP/B2TP 命名歧义排除」小节：显式关闭 "B2TP = Brca1+Pten 第二个重复" 的读法
    （客户注释列 + 数据双向否定），并如实记录 B2TP 的 Brca2 敲除尚未正面证实（缺 Brca2 sgRNA）
  - 2026-07-16 — **客户补齐 Brca2×3 sgRNA（邮件附图）→ B2TP 的 Brca2 敲除已正面证实，§1 的 caveat 关闭**：
    三条 guide 均唯一命中 Brca2 exon 3 CDS + 完美 NGG PAM（切点 chr5:150452957/961/989）；B2TP 在切点
    **0/14 条 WT read**（唯一无野生型等位的样本）+ 31bp 双切除（150452958-150452988，frameshift @codon~31/3329）
    + g9 切点 7 条 clip 断点堆叠 → **双等位 KO**。tumor3 仍保留 WT Brca2 → **确定不来自 B2TP**
    （B2TP 无 WT 等位可传），tumor3「未编辑亚克隆」结论由 Moderate 升 High。
    新增 `scripts/study_A/A3b_brca2_allele_quant.py`；交付 `custom_research_report_20260716/`（Addendum 1）
  - 2026-07-15 — §3 补入「Study B 品系背景 = 实测发现」重要更正：WGS 实测每样本 ~5–6M 变异 vs GRCm39，
    ~90% 命中 MGP 品系目录 → 组织**并非同源纯 C57BL/6J**，携带全基因组范围非-6J（最可能 129 型）背景。
    这**推翻了 §3 line 125 早前"背景=C57BL/6≈参考≈其正常"的假设**（原假设仅据客户口述）。量级论证排除
    "培养/世代/编辑漂变"解释（少 2–3 个数量级），educated guess = 129 型 ES 细胞背景；精确株系待 strain-assignment 确认
- **Prepared by**: Zhen Gao, PhD, Principal Bioinformatics Scientist, Athenomics
- **样本**: 12 个 WGS 样本，两组各 6 个（详见 §1）
- **数据来源**: 待客户/Jing 提供 fastq 路径

---

## 1. 样本说明与研究目的（据客户 2026-07-07 邮件）

**Study A — 前 6 个样本：基因编辑细胞 → 体内成瘤**

。研究目的：

- A1：基因编辑后，细胞的基因是否**真实发生了预期改变**（编辑验证）。
- A2：后来长出的 **tumor cell 的基因组变化**情况（肿瘤基因组表征）。

**Study B — 后 6 个样本：Lats1/2 flox 小鼠组织**
两个品系：**Lats1/2 flox/flox（L1L2）** 和 **Lats1/2 flox/flox-iHPV（L1L2H）**。在这两个品系里发现一些**不寻常的东西**，想知道其基因组相比正常小鼠**有无异常**、是否某些基因异常导致了该现象。

**客户与单位（MGH 实验室页面核实 2026-07-08）：**

- **Jinpeng Ruan**（金鹏），Research Assistant / 厦门大学在读博士。官方拼写 "Jinpeng"（非 Jingpeng）。
- **Wang Lab**，PI **Cheng Wang, PhD** · Obstetrics & Gynecology · **Vincent Center for Reproductive Biology** · **Massachusetts General Hospital (MGH)**。
- 页面：https://www.massgeneral.org/obgyn/vcrb/research/wang-lab
- 实验室方向：**高级别浆液性卵巢癌 (HGSOC)** 与 **多囊卵巢综合征 (PCOS)**。→ HGSOC 起源于输卵管上皮、由 TP53+BRCA1/2 缺失(HRD)驱动，与本项目 Study A(Trp53;Brca1/2;Pten)、Study B(输卵管异常) 高度吻合——**整个项目是卵巢癌 HGSOC 模型**。

**样本表（客户附图 `docs/client_materials/` + 数据 `/home/gao/Dropbox/JinPeng/` 核实）：** 12 样 PE150 NovaSeq X Plus WGS，gzip fastq 合计 **462 GiB**（`du -shc` 报 463 G）。文件大小 = `/home/gao/Dropbox/JinPeng/` 实测（2026-07-12）。

| # | 生物学标签          |  fastq 前缀 (S#)  | R1/R2 (GiB) |  Type  | 组 | 说明                                                          |
| :-: | :------------------ | :----------------: | :---------: | :----: | :-: | :------------------------------------------------------------ |
| 1 | **RO_origin** | 3852R0origin (S65) | 15.5 / 15.4 |  Cell  | A | **Trp53⁺/⁻;Cas9 未编辑亲本 = Study A matched normal** |
| 2 | RO_B1TP             |  3852R0B1TP (S63)  | 19.9 / 19.3 |  Cell  | A | 电穿孔 sgRNA**Brca1+Pten** KO                           |
| 3 | RO_B2TP             |  3852R0B2TP (S64)  | 19.7 / 19.4 |  Cell  | A | 电穿孔 sgRNA**Brca2+Pten** KO                           |
| 4 | 1st tumor           |     3868 (S66)     | 20.9 / 20.9 |  Cell  | A | B1TP/B2TP 成瘤消化                                            |
| 5 | 2nd tumor           |     7352 (S67)     | 15.7 / 15.8 |  Cell  | A | 同上                                                          |
| 6 | 3rd tumor           |     8599 (S70)     | 26.3 / 26.4 |  Cell  | A | 同上（最大样，~53 GiB）                                       |
| 7 | L1L2_3M             |    A2371 (S72)    | 16.4 / 16.3 | Tissue | B |                                                               |
| 8 | L1L2H_3M            |    A2353 (S71)    | 15.9 / 16.0 | Tissue | B | +iHPV                                                         |
| 9 | L1L2_12M            |     8464 (S68)     | 24.7 / 24.9 | Tissue | B |                                                               |
| 10 | L1L2H_12M           |     8465 (S69)     | 21.4 / 21.2 | Tissue | B | +iHPV                                                         |
| 11 | L1L2_18M            |     3685 (S61)     | 19.1 / 19.0 | Tissue | B |                                                               |
| 12 | L1L2H_18M           |     3689 (S62)     | 16.0 / 15.8 | Tissue | B | +iHPV                                                         |

> **数据来源**：分组以客户附图（`docs/client_materials/5c617d…png`）+ 实际 fastq 文件名为准；两者一致，附图信息正确、无需修正。数据文件夹内的 `Sample_Sheet_for_interal_use.xlsx`（实验员填写）客户已声明完全错误，不采用。
> （更正记录：本方案 0707 早前一版曾误称附图把 18M L1L2 标为 "3852" 并"修正"为 3685——那是**我方转录笔误**（把附图的 `3685` 看成 `3852`），附图本身一直是正确的 `3685`。已撤销该"修正"说法。请你对照附图再确认一下 18M 两行的 `3685`/`3689` 与 L1L2/L1L2H 的对应，因为我已暴露过看错数字。）

- **Study A = Trp53/Brca/Pten CRISPR 成瘤模型（客户 2026-07-07 补充确认）**：
  - **亲本 `RO_origin`**：来自 **Trp53⁺/⁻ ; Cas9** 转基因小鼠的原代细胞（自带一个 Trp53 缺失等位 + 组成型表达 Cas9），**未体外编辑**。
  - **RO_origin 即配对正常**，且它**已含 Trp53⁺/⁻ + Cas9 这两个工程特征** → 拿它作 normal 能把这些背景正确扣掉，只留 sgRNA 编辑(Brca/Pten)与肿瘤获得性改变。**Study A 不缺对照。**
  - **B1TP / B2TP**：以亲本为基础，**电穿孔导入 sgRNA**（Cas9 已在细胞内，只递送 sgRNA）——B1TP 靶 **Brca1+Pten**、B2TP 靶 **Brca2+Pten**（敲除）。
  - **4/5/6 tumor**：把 B1TP/B2TP 打入小鼠 → 长出实体瘤 → 消化**不同**肿瘤得到的细胞株（哪个肿瘤来自 B1TP 还是 B2TP 未指明，**可由 WGS 反推**：Brca1-KO→B1TP、Brca2-KO→B2TP）。

**⚠ B1TP / B2TP 命名歧义排除（2026-07-15 补记）**：标签 `RO_B2TP` 存在两种可能读法——(a) **Brca2+Pten**，(b) **Brca1+Pten 的第二个重复**（把 `B1/B2` 当作重复编号）。经核实，**采用 (a)，排除 (b)**，双向证据如下：

- **客户注释列（原始来源，非我方推断）**：客户附图样本表 `docs/client_materials/5c617d666e25fdc38f6264130f6f77ee.png` 的描述列白纸黑字分别写着 `RO_B1TP → sgRNA-Brca1-Pten`、`RO_B2TP → sgRNA-Brca2-Pten`。即 `B` 后的 `1/2` 对应 **Brca1/Brca2**（各自显式命名靶基因），**不是** replicate 1/2。读法 (b) 与客户自己的注释列直接冲突。
- **WGS 数据独立否定 (b)**：编辑验证（切点基因型）中 **B2TP 在 Brca1 切点为 0/0 野生型**、Pten 切点为纯合 indel。若 B2TP 真是"Brca1+Pten 的第二个重复"，它应像 B1TP 一样携带 Brca1 indel——但它没动 Brca1。故数据本身即排除 (b)。
- ~~**须如实保留的 caveat**：本批客户只提供 Brca1+Pten 的 sgRNA、**未给 Brca2 的 3 条 guide**，因此 B2TP 目前只能确认"**Pten 敲除 + Brca1 野生型**"——此结果**与** Brca2+Pten 一致，但 **Brca2 敲除本身尚未被正面证实**（未见到 Brca2 indel 的直接证据）。现有证据链是"排除 Brca1、符合 Brca2 设计"，而非"看到了 Brca2 的编辑"。~~
  → **✅ 2026-07-16 caveat 已关闭**：客户补齐 Brca2×3 sgRNA，切点定向验证**直接看到了 Brca2 的编辑**——
  B2TP 在 chr5:150452957/961/989 三切点 **0/14 WT read**（六样本中唯一无野生型等位者）、带 31bp 双切除
  （150452958-150452988，两端正好落在 guide 切点上，frameshift @ codon~31/3329 → null）+ g9 切点 7 条 clip
  断点堆叠（第二个被破坏的等位）。证据链已从"排除 Brca1、符合 Brca2 设计"升级为"**直接观察到 Brca2 双等位敲除**"。
  读法 (a) `B2TP = Brca2+Pten` 至此由**直接观察**确证，不再依赖推断。详见 `custom_research_report_20260716/`。
- **⚠ 方法学教训（2026-07-16，务必记住）**：Brca2 切点上 `bcftools call` **一条 indel 都没报 = 假阴性**
  （B2TP 25 条 read 里只有 1 条完美 match，编辑等位全被 soft-clip/31bp-del 漏掉）。若当时信了 caller 的
  "无 indel"，结论会变成"Brca2 没编辑"——与事实完全相反。**复杂编辑位点必须落到 read/CIGAR 层面看**，
  不能只信 variant caller。另：本位点 clip 背景高达 ~30%（连未编辑的 origin 都是），故 **clip 比例本身不特异**，
  真信号靠**断点堆叠在同一切点碱基**；初版判据过松曾把 origin 判成 50% edited（已加硬断言防回归）。
- **Study B = 基因型 × 年龄时间序列（客户 2026-07-08 详答）**：L1L2 vs L1L2H × 3M/12M/18M，各 n=1。**关键澄清**：
  - **背景 = C57BL/6**（≈ GRCm39 参考本身 → 背景扣除几乎免费，见 §3）。
  - **送样小鼠只有 flox、无 Cre → Lats1/2 未删除**（只带 loxP）；**iHPV 是 Cre 依赖的 lox-stop-lox**（CAG–loxP–EGFP–pA–loxP–E6/E7–IRES–Luc），**无 Cre 时 E6/E7 不表达**（表达 EGFP）。→ **两个工程元件都"上膛未击发"。**
  - **iHPV 源** = Addgene #13712（pB-actin E6E7，Munger）；模型见 **PMC4662542**。
  - **表型 = 输卵管异常，年龄依赖（>10 月龄出现）**。这解释了 3/12/18M 年龄序列（3M 在表型前）。

==== 07-15：手工添加

核心逻辑拆解与 WGS 预期表型
根据描述，这只小鼠携带着两颗“定时炸弹”，但由于没有 Cre 酶，它们都处于休眠/伪装状态。

1. 炸弹一：Lats1/2 flox/flox (L1L2)
   生物学状态： 只有 flox（即 loxP 位点），无 Cre。Lats1 和 Lats2 基因功能完全正常，没有被删除。
   WGS 数据中该看到的 (What you should see):
   完整覆盖度： Lats1 和 Lats2 的靶向外显子（Target Exons）的 Read 覆盖度 (Coverage) 应该与上下游区域一致，绝对不能有 Deletion。
   loxP 插入验证： 你能找到比对到 小鼠内含子 -- loxP (34bp) -- 小鼠内含子 的 Junction Reads，证明 loxP 位点确实存在且位置正确。
2. 炸弹二：iHPV (<mark>CAG–loxP–EGFP–pA</mark>–loxP–E6/E7–IRES–Luc)
   这是一个经典的 LSL (lox-Stop-lox) 表达盒，由强大的普遍性启动子 CAG 驱动。
   生物学状态 (无 Cre 时)： CAG 启动子正在疯狂工作，但它只读到了 EGFP 和 pA（提前终止信号）。因此，小鼠目前浑身发绿光 (表达 EGFP)，但致癌的 E6/E7 基因被前面的 lox-Stop-lox 阻挡，处于静默状态（未表达）。
   生物学状态 (假设有 Cre 时)： Cre 会切掉两个 loxP 之间的片段（即切掉 EGFP-pA），此时 CAG 启动子直接连上 E6/E7，致癌基因爆发，小鼠同时发出荧光素酶（Luc）的光。
   WGS 数据中该看到的 (What you should see):
   完整的构建体序列： 你在测序数据里，必须能拼出或比对出完整的 CAG–loxP–EGFP–pA–loxP–E6/E7–IRES–Luc 序列。
   无切除证据： EGFP-pA 区域必须有满额的 Read 覆盖。绝对不能出现直接连接 CAG 和 E6/E7 的跨 loxP 拼接 Reads。

==== 以上 07-15：手工添加

**仍需/可取的材料（Study A/B 基本齐）：**

- ~~Study A：Brca2 的 sgRNA（Pten/Brca1 已到）；三个 tumor 谱系（可数据反推）。~~
  → **✅ 2026-07-16 全部了结**：Brca2 sgRNA 已补齐并完成验证；三个 tumor 谱系已反推完毕
  （tumor1/2 = B1TP；tumor3 = 未编辑亚克隆，非 B1TP/B2TP）。**Study A 唯一剩余开放项 = tumor3 身份**，
  且已不是编辑位点问题而是**样本 identity 问题** → 待办：tumor3 vs RO_origin 的 SNP-fingerprint 比对，
  区分"亲本的未编辑逃逸亚克隆"(指纹一致) vs "样本混淆/调包"(指纹不一致)。现有数据即可跑，无需新测序。
- Study B：**iHPV 全构建体序列**（Addgene #13712 E6/E7 可公开下 + PMC4662542 的 CAG-LSL-EGFP-Luc 载体图）；**Lats1/2 loxP 打靶位点**（判断 floxed 外显子位置，查是否有体细胞重组）。均可自行从公开来源+数据补足。

---

## 2. 生物学背景（公开文献，指导分析与判读）

- **Study A（编辑激活）**：Trp53⁺/⁻ 背景 + Brca1/2 + Pten CRISPR KO——经典强驱动组合：Brca1/2 缺失→同源重组缺陷(HRD)、大尺度重排/拷贝数不稳定；Pten→PI3K；Trp53 常经 LOH 丢第二等位。肿瘤基因组预期显著改变，用 origin 配对可干净鉴定。
- **⚠️ Study B（工程元件"上膛未击发"，需重构假设）**：**关键更正**——送样小鼠 Lats1/2 **只 flox 未删**（无 Cre）、iHPV 是 **lox-stop-lox，E6/E7 未激活**（无 Cre，表达 EGFP）。因此：

  - **"Lats1/2 缺失→CIN/非整倍体"这条机制按设计并不成立**（基因未删）；它只在**存在体细胞/渗漏的 loxP 重组**把 floxed 外显子删掉时才相关。
  - E6/E7（E6→p53、E7→Rb）**按设计也未表达**；同理只在渗漏重组切除 stop 盒后才激活。
  - 参考背景：Lats1/2⁻/⁻ 确实致中心体扩增/非整倍体，E6/E7 致 p53/Rb 失活——这些是"若被激活会怎样"的知识，用于判读，不是这批小鼠的既定状态。
- **重构后的 Study B 核心假设（解释输卵管 >10M 异常）**，按优先级：

  1. **iHPV 构建体插入突变**：CAG-LSL-EGFP-E6E7-Luc 整合在何处？若打断/影响输卵管相关基因 = 表型之因（整合位点分析，B3）。
  2. **渗漏/体细胞 loxP 重组**：随龄(>10M)在输卵管细胞发生 Cre 非依赖重组 → 局部删 Lats1/2 或激活 E6/E7 → 病灶（查 floxed 外显子/stop 盒的体细胞丢失、局部拷贝数/嵌合）。
  3. **de novo 基因组病灶**（CNV/SV/突变）与工程无关的自发改变。

  - L1L2（仅 loxP）也报"异常"这点本身可疑 → 需比 L1L2 vs L1L2H、比年龄，看异常是否与 iHPV 或年龄相关。

---

## 3. 核心问题：没有对照小鼠，能不能做？——能，且大部分高价值问题可做

> 这是本项目最重要的方法学判断。**客户信息到齐后结论更明确：**
>
> - **Study A 自带配对正常**（`RO_origin` 亲本细胞）→ 编辑验证、肿瘤体细胞突变/CNV 全部可**干净配对调用**，无对照问题。
> - ~~**Study B 背景 = C57BL/6，而 GRCm39 参考本身就是 C57BL/6J** → **参考基因组 ≈ 它们的正常**，偏离参考的 CNV/SV/de novo 变异可直接判读；残余的只有 C57BL/6J vs 6N 亚系差异，用 Sanger MGP/dbSNP 目录扣掉即可。**这比背景发散的品系好办得多，专门野生型对照的必要性进一步下降。**~~ **⚠ 此假设（仅据客户口述）已被 WGS 实测推翻——见下方「Study B 品系背景 = 实测发现」。** 组织并非同源纯 6J，携带全基因组非-6J 背景（残余远不止 6N 亚系差异）；不过用 MGP 扣背景这一做法本身仍成立，只是扣掉的量级大得多。
> - 再加 **基因型(L1L2 vs L1L2H)+ 年龄(3/12/18M)** 内部对照轴。
>   总体：**核心问题都能做**（CNV/SV 参考-free、编辑验证配对、基因型×年龄内轴均不受品系背景影响）；仅"de novo 变异直接读参考"这一条需先扣掉大额品系背景（MGP）。

---

### 3bis. Study B 品系背景 = 实测发现（2026-07-15，据 germline WGS）

> **一句话结论**：Study B 六个组织**并非同源纯 C57BL/6J**，而是携带**全基因组范围的非-6J 经典近交系背景**（最可能 129 型，educated guess）。这更正了上面 line 125 早前仅据客户口述的"背景=C57BL/6≈参考"假设。

**证据（每样本 germline vs GRCm39，见客户报告 §6.5 表）**：

| 样本 | 总变异 vs GRCm39 | 扣 MGP 后 private |
| :--- | :---: | :---: |
| L1L2_3M | 5,101,717 | 573,954 |
| L1L2H_3M | 5,334,890 | 562,308 |
| L1L2_12M | 6,110,749 | 658,142 |
| L1L2H_12M | 5,612,998 | 589,237 |
| L1L2_18M | 5,899,863 | 625,093 |
| L1L2H_18M | 5,884,909 | 615,276 |

**为什么是"品系背景"而非"培养/世代/编辑漂变"（关键判读，防审稿人质疑）**：

1. **量级差 2–3 个数量级**：小鼠种系突变率 ~5×10⁻⁹/bp/世代 → 每世代仅 ~10–15 个新变异；即便数百世代也就几千个。同株系长期漂变的经验上限 = 6J vs 6N（分家 ~70 年、数百世代）也**仅数万**。培养体细胞突变数千级，且 Study B 是**组织非长培养细胞系**；CRISPR/floxing 只动少数位点 + 个位到百位脱靶。三者合计 ≤ ~10⁴，比实测 ~5.5×10⁶ **少 2–3 个数量级**。只有**独立近交系谱系**（祖先多态在两谱系间分别固定）才能产生百万级差异。
2. **~90% 命中 MGP（决定性）**：~5–6M 中约 4.5–5.4M（≈90%）能被 Sanger MGP 品系目录扣掉。变异**能进 MGP ⇔ 是与已知品系共享的遗传性多态**；培养/世代/编辑产生的 de novo 突变是**样本私有、不会**被 MGP 扣掉。故数据本身证明主体是遗传性品系背景。
3. **private 余量全基因组均匀**（~0.57–0.66M，∝ 染色体长度）→ 品系背景随基因组分离，而非富集在编辑位点附近。

**株系推断（educated guess）**：~5–6M ≈ 129 vs 6J 的**整段**距离（~4–6M SNP/indel），且 90% 被 MGP 解释 → 背景是**大范围（近乎全基因组）129 型**，**不是**"小残段/未backcross干净"的少量残余，而是该系**大体就在非-6J 背景上**。129 是首选解释因为基因打靶/floxed 等位常在 **129 来源 ES 细胞**里做；本系 ES 细胞来源未见文档，故为 leading explanation 非定论。**精确株系（纯 129 vs 129×B6 混）需专门 strain-assignment 分析**（逐株系 MGP 目录比对每样本 private 变异），可按需运行。

**对分析的影响**：CNV/非整倍体（参考-free）、编辑验证（配对）、基因型×年龄内轴均**不受影响**；仅 Study B de novo 候选挖掘须先扣大额品系背景（MGP，已做）+ 高影响功能限制 + L1L2 vs L1L2H / 随龄复现过滤（见 §B4、客户报告 §6.5 待办）。

**逐问题可行性矩阵：**

| 分析问题                                   |  是否必须配对正常  |           无对照可行性           | 怎么做                                                                                                                                              |
| :----------------------------------------- | :----------------: | :-------------------------------: | :-------------------------------------------------------------------------------------------------------------------------------------------------- |
| 拷贝数 / 非整倍体 / 倍性                   |         否         |          **可以**          | 基因组内 coverage 比值即得拷贝数谱，天然免对照（Control-FREEC tumor-only / CNVkit flat-ref / mosdepth 分箱）。**这是 Study B 最可能的答案。** |
| 大结构变异 SV                              |       弱依赖       |        **基本可以**        | delly/tiddit 调用；用公共正常小鼠(MGP/DGV)过滤常见 SV                                                                                               |
| 基因编辑验证 (A1)                          | 否（但需编辑意图） |          **可以**          | 直接看靶位点 vs GRCm39 参考 + 预期改变；已知 gRNA/donor 即可判 on-target/indel/HDR                                                                  |
| 工程等位验证 (Lats1/2 loxP、iHPV 整合位点) |         否         |          **可以**          | 混合参考 + 结合部检测（复用 /wgs 模式B、13 号项目流程）；需 loxP/HPV/构建体序列                                                                     |
| 肿瘤体细胞**点突变**干净鉴定 (A2)    |      理想需要      | **Study A 可（有 origin）** | Study A 用`RO_origin` 配对 Mutect2；Study B 若要体细胞级则用 MGP + 年龄内基线将就                                                                 |
| 区分"致病变异" vs "品系背景 SNP"           |      理想需要      |         **部分可**         | 用 Sanger MGP 品系目录作替代扣除；聚焦高影响 + 多样本复现变异                                                                                       |

**无对照的四条缓解策略（本方案据此设计）：**

1. **批内互为对照**：Study A **有 `RO_origin` 亲本** → 作编辑细胞与 3 个 tumor 的 matched normal（配对体细胞调用，干净）；Study B **L1L2 vs L1L2H 互比**（HPV 效应）+ **年龄进程 3M→12M→18M**（不稳定是否累积）+（若 L1L2 未删 Lats1/2）以 L1L2 作 L1L2H 的相对基线。
2. **公共替代对照（关键）**：**Sanger 小鼠基因组计划(MGP)** 36 品系相对 C57BL/6J 的 SNP/indel/SV VCF + 小鼠 dbSNP → 扣除背景 germline，代替"正常小鼠"。这是无对照时的核心替身。
3. **倍性/CNV 分析天然免对照**：以基因组自身 coverage 比值判拷贝数，不需要 matched normal——直击 Study B 最可能的"不寻常"。
4. **复现过滤**：多个工程样本共有、而参考/公共正常无的病灶 = 强候选。

**建议但非必需**：若客户能补 **1 个配对正常**（Study A 的未编辑亲本细胞；Study B 的野生型同窝/背景品系一只 WGS），体细胞点突变鉴定会显著变干净。**方案在无对照下即可产出核心结论，有对照则锦上添花**——会把这条写清楚提给客户，但不阻塞开工。

---

## 4. Study A 分析设计（编辑细胞 → 肿瘤）

- **比对/QC/去重**：nf-core/sarek（GRCm39，bwa-mem2，skip BQSR）——沿用 13 号项目验证过的配置。
- **A1 编辑验证（靶点+sgRNA 序列已知）**：**9 条 sgRNA 全部到手**（`docs/client_materials/sgRNA_guides.md`：
  Pten×3、Brca1×3、**Brca2×3（2026-07-16 补齐）**）。流程：在 GRCm39 定位每条 spacer（+反向互补，配 NGG PAM）→ 预测切点(PAM 上游 3bp) → 在切点比 B1TP·B2TP·三个 tumor **vs RO_origin** 查 indel/frameshift（KO 是否发生、等位比例/合子型），IGV 目视。**用 tumor 里 Brca1 还是 Brca2 被 KO 反推其谱系(B1TP/B2TP)**。
  **✅ 已完成**：B1TP=Brca1+Pten KO；B2TP=**Brca2 双等位 KO**+Pten KO；tumor1/2=B1TP 谱系；
  tumor3 三个靶基因全 WT → 未编辑亚克隆（High）。Brca2 切点定量见 `A3b_brca2_allele_quant.py`。
- **A2 肿瘤基因组（RO_origin 作 matched normal，配对干净）**：
  - **点突变/indel**：Mutect2 **tumor-vs-origin 配对**（真正体细胞调用；origin 已含 Trp53⁺/⁻+Cas9，背景被正确扣掉）。
  - 
  - C / CNVkit（以 origin 为对照，更准）。
  - SV：delly + tiddit（origin 过滤 germline）。
  - **HRD 相关**：Brca1/2 缺失 → 同源重组缺陷 → 关注大片段重排、拷贝数不稳定、SBS3/HRD 样特征；Pten 缺失 → PI3K 轴。聚焦癌基因。
  - 三个 tumor 互比 → 共有 vs 私有变异、克隆演化。
  - **注意工程特征**：所有 A 样本含 **Cas9 转基因**（可顺带定位其整合位点）和 Trp53 缺失等位；因都在 origin 里，配对分析不会误报为"获得性"。

## 5. Study B 分析设计（找输卵管异常的基因组根源）—— C57BL/6，参考≈正常

**目标表型**：输卵管异常，>10 月龄出现。**工程元件按设计未激活**（Lats1/2 未删、E6/E7 未表达）→ 分析按 §2 三假设排优先级。比较轴：同龄 **L1L2 vs L1L2H**、同基因型 **3M→12M→18M**。

- **比对/QC**：GRCm39（=C57BL/6J，即它们的近等基因正常）。
- **B1 iHPV 构建体整合位点（首要假设：插入突变）**：建**混合参考**（GRCm39 + CAG-LSL-EGFP-E6E7-Luc 构建体，序列取自 Addgene #13712 + PMC4662542）→ 结合部检测（复用 /wgs 模式B、13 号项目脚本）→ **整合落在哪个基因**、是否打断输卵管相关基因。仅 L1L2H 有此构建体。
- **B2 工程等位状态 + 体细胞重组（第二假设）**：
  - Lats1/2 **loxP 是否存在**、floxed 外显子覆盖是否完整（**查是否有渗漏/体细胞重组把外显子删掉** → 局部覆盖下降/嵌合）；
  - iHPV 的 **EGFP-pA "stop" 盒是否有体细胞丢失**（= E6/E7 被渗漏激活的证据）；
  - 老龄样本(18M)尤其关注这类随龄累积的嵌合事件。
- **B3 拷贝数 / 非整倍体 / SV（第三假设 + 通用扫描）**：Control-FREEC(tumor-only)/CNVkit/mosdepth 分箱 + delly/tiddit → 全基因组拷贝数谱、倍性、大 SV；MGP/DGV 过滤 C57BL/6 背景。
- **B4 de novo 变异找候选**：样本 vs GRCm39 → 用 **Sanger MGP(6J/6N) + dbSNP** 扣 C57BL/6 亚系背景 → 剩余高影响 + 多样本/随龄复现变异 → 候选；**L1L2 vs L1L2H 差异**单列。
- **B5 关联表型**：候选是否落在输卵管发育/纤毛/上皮、Hippo/YAP、DNA 修复、基因组稳定性相关基因；结合 L1L2 也报异常这一疑点判读。

---

## 5.5 Pipeline 架构、计算资源与执行顺序

### 架构（一个合并混合参考 + sarek 双模式 + 定制层）

- **合并混合参考**（一次建、12 样共用，复用 13 号项目做法）：`GRCm39 + SpCas9 + iHPV构建体`。Study A 样本在 Cas9 contig 有覆盖（可定位 Cas9 整合位点）、Study B L1L2H 在 iHPV contig 有覆盖（整合位点）；互为阴性对照。
- **sarek 双模式**（`--aligner bwa-mem2 --skip_tools baserecalibrator --tools tiddit`，Manta 已弃用见 13 号教训）：
  - **Study A = somatic**：samplesheet 一个 patient=`RO`，`RO_origin` 为 normal(status 0)，其余 5 个（B1TP/B2TP/3868/7352/8599）为 tumor(status 1) → sarek 对每个 tumor 跑 **Mutect2 配对**（vs origin）+ SV。干净体细胞调用。
  - **Study B = germline**：6 样各自 germline（HaplotypeCaller/Strelka）→ 变异 → MGP/dbSNP 扣背景。
- **定制层（sarek 之上）**：
  - **CNV/倍性**：Control-FREEC（A 用配对、B 用 tumor-only）/ CNVkit + mosdepth 分箱。
  - **整合位点**：脚本复用 13 号项目 4/5（嵌合读段 + 人源/构建体区，MAPQ≥20，artifact 黑名单）——定 Cas9(A) 与 iHPV(B) 整合位点。
  - **编辑验证**：在 Brca1/Brca2/Pten sgRNA 切点、Lats1/2 loxP/floxed 外显子、iHPV stop 盒处，定向查 indel / 体细胞重组 / 覆盖丢失（配 IGV）。
  - **断点/候选注释**：GENCODE vM35。

### 计算资源与时间（12 样 WGS，~463GB 输入）

- 每样 ~30–40×；据 13 号试跑，单样比对+去重+SV ~7h（不含 Manta）。**queueSize=2** → 12 样 6 波 × ~7h ≈ **~42h 比对层**；变异/CNV 另计。**周五挂、下周一收**。
- 资源沿用 13 号 `local_resources.config`（≤48 线程/108GB，bwa-mem2 16 线程/50GB×2）。tmux 常驻 + 失败自动 resume + 前 3 分钟早期失败检测（CLAUDE.md）。
- 磁盘：输入 463GB + 12 CRAM(~10GB×12) + work 中间文件；`/home/gao/projects_2026H2` 现 4.2T 空余，够但需盯着（work 可分批清）。

### 执行顺序（分阶段，审阅通过后跑）

```
0  取序列建混合参考：SpCas9 + iHPV(Addgene13712/PMC4662542) → GRCm39_plus_Cas9_iHPV.fa (faidx/dict)
   下载 Sanger MGP + dbSNP；补装 control-freec/cnvkit
1  samplesheet：Study A（somatic, patient=RO, origin=normal）；Study B（germline, 6 样）
2  sarek 跑 A（somatic）+ B（germline）  ← 主比对/变异，~2 天，queueSize=2
3  CNV/倍性（Control-FREEC/CNVkit）：A 配对、B tumor-only
4  整合位点（脚本4/5 复用）：Cas9(A) / iHPV(B)
5  编辑验证 + 工程等位/体细胞重组（Brca/Pten/Lats1/2/stop 盒定向查）
6  MGP/dbSNP 背景扣除 → Study B de novo 候选；L1L2 vs L1L2H、年龄进程
7  注释 + 英文报告
```

---

## 6. 参考与工具

| 组件                                                                                                              | 状态/来源                                                       |
| :---------------------------------------------------------------------------------------------------------------- | :-------------------------------------------------------------- |
| GRCm39 (GENCODE M35) + bwa-mem2 索引                                                                              | ✅ 本机`/Work_bio/references/Mus_musculus/GRCm39/...`         |
| nf-core/sarek 3.8.1（somatic + germline）                                                                         | ✅ 已验证（4、13 号项目）                                       |
| delly 2.3.0 / mosdepth 0.3.14                                                                                     | ✅ regular_bioinfo                                              |
| Mutect2 (GATK) 配对/tumor-only                                                                                    | sarek 内置                                                      |
| Control-FREEC / CNVkit / ASCAT（CNV/倍性）                                                                        | ❌ 需补装（`mamba install -c bioconda control-freec cnvkit`） |
| **Sanger MGP 小鼠 VCF（6J/6N SNP/indel/SV，替代对照）**                                                     | ❌ 需下载（公开无限制）                                         |
| 小鼠 dbSNP                                                                                                        | 需确认本机有无                                                  |
| 混合参考构建体序列：**SpCas9**（公开）、**iHPV** CAG-LSL-EGFP-E6E7-Luc（Addgene #13712 + PMC4662542） | ❌ 我方自取                                                     |

---

## 7. 客户信息状态（已基本齐，可开工）

- ✅ **Study A 齐**：Trp53⁺/⁻;Cas9 亲本、电穿孔 sgRNA、B1TP=Brca1+Pten/B2TP=Brca2+Pten、4/5/6 肿瘤细胞株、RO_origin 配对正常；**9 条 sgRNA 全到手（Pten×3 + Brca1×3 + Brca2×3，Brca2 于 2026-07-16 补齐）→ 客户端已无待补项**。
- ✅ **Study B 齐**（客户 2026-07-08）：C57BL/6；flox-only 无 Cre；iHPV=Cre 依赖 LSL（Addgene #13712 + PMC4662542）；表型=输卵管异常 >10M。
- ✅ **数据**：`/home/gao/Dropbox/JinPeng/`（12 样，见 §1；Excel 忽略）。
- **我方自取（不卡客户）**：Brca2 sgRNA（无则扫全基因）；iHPV 构建体序列（Addgene 13712 + PMC4662542）；Cas9(SpCas9) 序列（公开）；Lats1/2 loxP 坐标（文献/数据自查）；Sanger MGP + 小鼠 dbSNP 下载。
- 可选：请客户补 1 只野生型同窝作 Study B 绝对正常（非必需，C57BL/6 已≈参考）。

---

## 8. 交付物（拟）

`custom_research_report_YYYYMMDD/`，两个子报告或合并：

- **Study A**：`edit_verification/`（Brca1/2/Pten KO 确认 + 谱系反推）、`somatic_variants/`（tumor-vs-origin，含 Trp53 LOH）、`cnv/`、`sv/`、`cas9_integration/`。
- **Study B**：`ihpv_integration/`（整合位点，首要）、`engineered_alleles/`（loxP/floxed/stop 盒状态 + 体细胞重组）、`cnv_ploidy/`、`sv/`、`candidates/`（MGP 扣背景后 de novo，L1L2 vs L1L2H / 年龄）。
- 英文报告，结构按 CLAUDE.md（Objectives / Key Findings / Sample Info / Rationale / Methods / Results / Conclusions / Deliverables）；CNV/倍性图必配文字解读。

---

## 参考文献（公开）

- **iHPV 模型（本项目所用）**：PMC4662542（Cre 依赖 CAG-LSL-EGFP-E6E7-Luc）；E6/E7 源 Addgene #13712（pB-actin E6E7, Munger et al, J Virol 1989）。
- 背景知识（"若激活会怎样"）：Lats1/2⁻/⁻ 致中心体扩增/非整倍体（Lats2/Kpm, EMBO J；Lats1 centrosome, Sci Rep 2016；p53–Lats2, Genes Dev 2006）；E6→p53、E7→Rb。
- Brca1/2 缺失 → HRD/基因组不稳定；Pten → PI3K。
- 替代对照：Sanger Mouse Genomes Project（vs C57BL/6J，SNP/indel/SV VCF，公开）。

---

*客户信息与数据已齐（Study A/B 均可开工）。待你审阅本方案；通过后按 §5.5 执行顺序，先建混合参考+补装工具，再 sarek 双模式跑（~2 天），随后定制层。核心结论（编辑验证、肿瘤配对突变、CNV/整合位点）均可产出，Study B 因 C57BL/6 而尤其干净。*
