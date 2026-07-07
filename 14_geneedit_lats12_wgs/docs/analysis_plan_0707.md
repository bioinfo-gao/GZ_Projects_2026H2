# 项目 14 — 基因编辑肿瘤细胞 + Lats1/2 Hippo 小鼠 WGS 分析方案

- **Project**: 14_geneedit_lats12_wgs
- **Plan Date**: 2026-07-07
- **Prepared by**: Zhen Gao, PhD, Principal Bioinformatics Scientist, Athenomics
- **样本**: 12 个 WGS 样本，两组各 6 个（详见 §1）
- **数据来源**: 待客户/Jing 提供 fastq 路径

---

## 1. 样本说明与研究目的（据客户 2026-07-07 邮件）

**Study A — 前 6 个样本：基因编辑细胞 → 体内成瘤**
消化小鼠组织 → 细胞培养 → **基因编辑** → 打入小鼠体内长出 tumor。研究目的：
- A1：基因编辑后，细胞的基因是否**真实发生了预期改变**（编辑验证）。
- A2：后来长出的 **tumor cell 的基因组变化**情况（肿瘤基因组表征）。

**Study B — 后 6 个样本：Lats1/2 flox 小鼠组织**
两个品系：**Lats1/2 flox/flox（L1L2）** 和 **Lats1/2 flox/flox-iHPV（L1L2H）**。在这两个品系里发现一些**不寻常的东西**，想知道其基因组相比正常小鼠**有无异常**、是否某些基因异常导致了该现象。

**样本表（客户附图 `docs/client_materials/`，2026-07-07）：**

| # | Sample Name | Type | 组 | 说明 |
| :---: | :--- | :---: | :---: | :--- |
| 1 | **3852 RO_origin** | Cell | A | **未编辑亲本 primary cell（= Study A 的 matched normal！）** |
| 2 | 3852 RO_B1TP | Cell | A | 基因编辑 **sgRNA-Brca1 + Pten**（CRISPR KO） |
| 3 | 3852 RO_B2TP | Cell | A | 基因编辑 **sgRNA-Brca2 + Pten**（CRISPR KO） |
| 4 | 3868-1st_RO tumor | Cell | A | 1st tumor |
| 5 | 7352-2nd_RO tumor | Cell | A | 2nd tumor |
| 6 | 8599-3rd_RO tumor | Cell | A | 3rd tumor |
| 7 | A2371_L1L2_3M | Tissue | B | L1L2，3 月龄 |
| 8 | A2353_L1L2H_3M | Tissue | B | L1L2H（+iHPV），3 月龄 |
| 9 | 8464_L1L2_12M | Tissue | B | L1L2，12 月龄 |
| 10 | 8465_L1L2H_12M | Tissue | B | L1L2H，12 月龄 |
| 11 | 3852_L1L2_18M | Tissue | B | L1L2，18 月龄 |
| 12 | 3689_L1L2H_18M | Tissue | B | L1L2H，18 月龄 |

- **Study A = Trp53/Brca/Pten CRISPR 成瘤模型（客户 2026-07-07 补充确认）**：
  - **亲本 `RO_origin`**：来自 **Trp53⁺/⁻ ; Cas9** 转基因小鼠的原代细胞（自带一个 Trp53 缺失等位 + 组成型表达 Cas9），**未体外编辑**。
  - **B1TP / B2TP**：以亲本为基础，**电穿孔导入 sgRNA**（Cas9 已在细胞内，只递送 sgRNA）——B1TP 靶 **Brca1+Pten**、B2TP 靶 **Brca2+Pten**（敲除）。
  - **4/5/6 tumor**：把 B1TP/B2TP 打入小鼠 → 长出实体瘤 → 消化**不同**肿瘤得到的细胞株（哪个肿瘤来自 B1TP 还是 B2TP 未指明，**可由 WGS 反推**：Brca1-KO→B1TP、Brca2-KO→B2TP）。
  - **RO_origin 即配对正常**，且它**已含 Trp53⁺/⁻ + Cas9 这两个工程特征** → 拿它作 normal 能把这些背景正确扣掉，只留 sgRNA 编辑(Brca/Pten)与肿瘤获得性改变。**Study A 不缺对照。**
- **Study B = 基因型 × 年龄时间序列**：L1L2 vs L1L2H（iHPV 效应）× 3M/12M/18M（进程）。各基因型/年龄 n=1。

**仍需向客户确认：**
- Study A：三个 tumor 是否都源自 B1TP/B2TP 谱系（哪个肿瘤来自哪种编辑）？sgRNA 具体序列/靶位点坐标（便于精确查编辑）。RO=何种组织来源的 primary cell。
- Study B：背景品系；**Cre 是否诱导**（L1L2 若未删 Lats1/2 则近野生型，可作 L1L2H 的相对基线）；**iHPV(E6/E7) 构建体 + loxP 打靶序列**；观察到的"不寻常"具体表型。
- 是否愿补 1 只**野生型同窝**作 Study B 的绝对正常（非必需，见 §3）。

---

## 2. 生物学背景（公开文献，指导分析与判读）

- **Lats1/Lats2 = Hippo 通路核心激酶、抑癌基因**：缺失 → YAP/TAZ 激活 → 促瘤（如基底样乳腺癌、肾发育异常等）。
- **Lats1/2 缺失是公认的基因组不稳定驱动因素**：Lats1−/− / Lats2−/− 细胞出现**中心体过度复制、多极纺锤体、染色体错排、微核、胞质分裂失败、非整倍体/四倍体**；Lats2 通过 p53 正反馈防止四倍体化。→ **客户所说"基因组不寻常"最可能是非整倍体 / 拷贝数异常 / 染色体不稳定(CIN)**，而这正是 WGS 拷贝数/倍性分析的强项。
- **iHPV = 可诱导 HPV16 E6/E7**（常 K14-CreER 驱动）：E6 降解 **p53**、E7 灭活 **Rb** → 解除检查点，与 Lats1/2 缺失协同放大基因组不稳定。L1L2H 相比 L1L2 多了这层"检查点解除"。
- **推论**：Study B 的核心发现很可能是**大尺度拷贝数改变/非整倍体**；L1L2H 应比 L1L2 更严重（HPV 解除 p53/Rb）。这为"无对照下也能做"提供了关键支点（见 §3）。

---

## 3. 核心问题：没有对照小鼠，能不能做？——能，且大部分高价值问题可做

> 这是本项目最重要的方法学判断。**看到样本表后结论进一步明确：**
> - **Study A 其实自带配对正常**（`3852 RO_origin` 亲本细胞）→ 编辑验证、肿瘤体细胞突变/CNV 全部可**干净配对调用**，无对照问题。
> - **Study B 无专门野生型**，但有**基因型(L1L2 vs L1L2H)+ 年龄(3/12/18M)** 的内部对照结构；且最可能的答案（非整倍体/CNV）**天然免对照**。
> 总体：**最有价值的问题都能做**；真正被削弱的仅"Study B 绝对意义上的品系背景扣除"，用 Sanger MGP 替身 + 基因型/年龄互比缓解。

**逐问题可行性矩阵：**

| 分析问题 | 是否必须配对正常 | 无对照可行性 | 怎么做 |
| :--- | :---: | :---: | :--- |
| 拷贝数 / 非整倍体 / 倍性 | 否 | **可以** | 基因组内 coverage 比值即得拷贝数谱，天然免对照（Control-FREEC tumor-only / CNVkit flat-ref / mosdepth 分箱）。**这是 Study B 最可能的答案。** |
| 大结构变异 SV | 弱依赖 | **基本可以** | delly/tiddit 调用；用公共正常小鼠(MGP/DGV)过滤常见 SV |
| 基因编辑验证 (A1) | 否（但需编辑意图） | **可以** | 直接看靶位点 vs GRCm39 参考 + 预期改变；已知 gRNA/donor 即可判 on-target/indel/HDR |
| 工程等位验证 (Lats1/2 loxP、iHPV 整合位点) | 否 | **可以** | 混合参考 + 结合部检测（复用 /wgs 模式B、13 号项目流程）；需 loxP/HPV/构建体序列 |
| 肿瘤体细胞**点突变**干净鉴定 (A2) | 理想需要 | **Study A 可（有 origin）** | Study A 用 `RO_origin` 配对 Mutect2；Study B 若要体细胞级则用 MGP + 年龄内基线将就 |
| 区分"致病变异" vs "品系背景 SNP" | 理想需要 | **部分可** | 用 Sanger MGP 品系目录作替代扣除；聚焦高影响 + 多样本复现变异 |

**无对照的四条缓解策略（本方案据此设计）：**
1. **批内互为对照**：Study A **有 `RO_origin` 亲本** → 作编辑细胞与 3 个 tumor 的 matched normal（配对体细胞调用，干净）；Study B **L1L2 vs L1L2H 互比**（HPV 效应）+ **年龄进程 3M→12M→18M**（不稳定是否累积）+（若 L1L2 未删 Lats1/2）以 L1L2 作 L1L2H 的相对基线。
2. **公共替代对照（关键）**：**Sanger 小鼠基因组计划(MGP)** 36 品系相对 C57BL/6J 的 SNP/indel/SV VCF + 小鼠 dbSNP → 扣除背景 germline，代替"正常小鼠"。这是无对照时的核心替身。
3. **倍性/CNV 分析天然免对照**：以基因组自身 coverage 比值判拷贝数，不需要 matched normal——直击 Study B 最可能的"不寻常"。
4. **复现过滤**：多个工程样本共有、而参考/公共正常无的病灶 = 强候选。

**建议但非必需**：若客户能补 **1 个配对正常**（Study A 的未编辑亲本细胞；Study B 的野生型同窝/背景品系一只 WGS），体细胞点突变鉴定会显著变干净。**方案在无对照下即可产出核心结论，有对照则锦上添花**——会把这条写清楚提给客户，但不阻塞开工。

---

## 4. Study A 分析设计（编辑细胞 → 肿瘤）

- **比对/QC/去重**：nf-core/sarek（GRCm39，bwa-mem2，skip BQSR）——沿用 13 号项目验证过的配置。
- **A1 编辑验证（靶点已知：Brca1/Brca2/Pten）**：在 **Brca1**(B1TP)/**Brca2**(B2TP)/**Pten** 位点，比 B1TP·B2TP·三个 tumor **vs RO_origin**，查 CRISPR sgRNA 诱导的 indel/frameshift（on-target KO 是否发生、等位比例/合子型）。**用 tumor 里哪个 Brca 被 KO 反推其来自 B1TP 还是 B2TP**（客户未指明谱系）。sgRNA 精确坐标到位可查切点；无坐标先按基因扫 indel。
- **A2 肿瘤基因组（RO_origin 作 matched normal，配对干净）**：
  - **点突变/indel**：Mutect2 **tumor-vs-origin 配对**（真正体细胞调用；origin 已含 Trp53⁺/⁻+Cas9，背景被正确扣掉）。
  - **Trp53 第二次打击（重点）**：亲本是 Trp53⁺/⁻，**肿瘤常丢失剩余野生型等位(LOH)→ Trp53⁻/⁻**；查 Trp53 位点 LOH/缺失，是经典的成瘤关键事件。
  - 拷贝数/非整倍体：Control-FREEC / CNVkit（以 origin 为对照，更准）。
  - SV：delly + tiddit（origin 过滤 germline）。
  - **HRD 相关**：Brca1/2 缺失 → 同源重组缺陷 → 关注大片段重排、拷贝数不稳定、SBS3/HRD 样特征；Pten 缺失 → PI3K 轴。聚焦癌基因。
  - 三个 tumor 互比 → 共有 vs 私有变异、克隆演化。
  - **注意工程特征**：所有 A 样本含 **Cas9 转基因**（可顺带定位其整合位点）和 Trp53 缺失等位；因都在 origin 里，配对分析不会误报为"获得性"。

## 5. Study B 分析设计（Lats1/2 品系找异常）—— 基因型×年龄，无绝对对照

**设计**：L1L2 vs L1L2H（iHPV 效应）× 3M/12M/18M（年龄进程），各 n=1。核心比较轴：
(a) **同龄 L1L2 vs L1L2H** → HPV(E6/E7=p53/Rb 失活) 带来的额外基因组不稳定；
(b) **同基因型跨年龄 3M→12M→18M** → 不稳定/病灶是否随龄累积（Lats1/2 缺失驱动 CIN 的时间演化）。

- **比对/QC**：GRCm39。
- **B1 拷贝数 / 非整倍体 / 倍性（首要，免对照）**：Control-FREEC / CNVkit / mosdepth 分箱 → 全基因组拷贝数谱、染色体臂级增删、倍性估计。**预期 L1L2H > L1L2 的不稳定程度。** 这是"不寻常基因组"最可能的落点。
- **B2 SV**：delly/tiddit → 大缺失/重复/易位/复杂重排；MGP/DGV 过滤背景。
- **B3 工程等位验证**：Lats1/2 位点的 **loxP 是否存在/是否已被 Cre 删除**（判断 Lats1/2 实际状态）；**iHPV(E6/E7) 整合位点**（混合参考 + 结合部，需 HPV/构建体序列）。
- **B4 背景扣除找致病候选**：品系 vs GRCm39 变异 → 用 **Sanger MGP + dbSNP** 扣掉背景 germline → 剩余高影响 + 多样本复现变异 → 候选。**L1L2 vs L1L2H 差异集**单列（HPV 相关）。
- **B5 关联生物学**：候选是否落在 Hippo/YAP、p53/Rb、有丝分裂/中心体、基因组稳定性相关基因上（呼应 Lats1/2 生物学）。

---

## 6. 参考与工具

| 组件 | 状态/来源 |
| :--- | :--- |
| GRCm39 (GENCODE M35) + bwa-mem2 索引 | ✅ 本机 `/Work_bio/references/Mus_musculus/GRCm39/...` |
| nf-core/sarek 3.8.1 | ✅ 已验证（4、13 号项目） |
| delly 2.3.0 / mosdepth 0.3.14 | ✅ regular_bioinfo |
| Control-FREEC / CNVkit（tumor-only CNV） | 需确认/补装 |
| Mutect2 (GATK) tumor-only + germline-resource | GATK 在 sarek 内 |
| **Sanger MGP 小鼠品系 VCF（SNP/indel/SV，替代对照）** | 需下载（公开无限制） |
| 小鼠 dbSNP | 需确认本机有无 |
| HPV16 / loxP / 构建体序列（模式B 整合分析用） | **需客户提供** |

---

## 7. 待客户确认清单（Study A 已基本答全；Study B 未答）

- ✅ **Study A 已澄清**（客户 2026-07-07 补充）：亲本 Trp53⁺/⁻;Cas9、电穿孔 sgRNA、B1TP=Brca1+Pten / B2TP=Brca2+Pten、4/5/6 为肿瘤消化细胞株。tumor 谱系可由数据反推。仅 sgRNA 精确坐标是 nice-to-have（无也能按基因扫 indel）。
- ⏳ **Study B 仍待客户答**：背景品系；**Cre 是否诱导过**（决定 Lats1/2 是否真删；也可由 WGS 自查 floxed 外显子覆盖）；**iHPV(E6/E7) 构建体 + Lats1/2 loxP 打靶序列**（仅卡 B3 整合/等位验证）；"不寻常"具体表型。
- **数据 fastq 路径**（操作性硬需求）。
- 是否愿补 1 只野生型同窝作 Study B 绝对正常（非必需）。

---

## 8. 交付物（拟）

`custom_research_report_YYYYMMDD/`：`qc/`、`copy_number/`（CNV/倍性谱，核心）、`sv/`、`variants/`（含 MGP 过滤后候选）、`edit_verification/`、`integration/`（iHPV/loxP）。英文报告，结构按 CLAUDE.md（Objectives / Key Findings / … / Conclusions / Deliverables），CNV/倍性图必配文字解读。

---

## 参考文献（公开）
- Lats1/2 与 CIN：Lats2/Kpm required for genomic integrity（EMBO J）；Lats1 suppresses centrosome overduplication（Sci Rep 2016）；p53–Lats2 feedback prevents tetraploidization（Genes Dev 2006）。
- Lats1/2 促瘤：Inactivation of LATS1/2 drives basal-like mammary carcinomas（PMC9705439）。
- iHPV：K14-CreER × HPV16 E6/E7 可诱导模型。
- 替代对照：Sanger Mouse Genomes Project（36 品系 vs C57BL/6J，SNP/indel/SV VCF，公开）。

---

*方案待客户确认样本结构与序列后细化并开工。核心结论（CNV/非整倍体/工程等位）在无对照下即可产出。*
