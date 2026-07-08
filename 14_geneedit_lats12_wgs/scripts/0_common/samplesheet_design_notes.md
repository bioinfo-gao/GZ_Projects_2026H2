# 为什么 sarek samplesheet 只能这么写（一眼看懂版）

对象：`A_somatic.csv`、`B_germline.csv`（由 `1_make_samplesheets.py` 从 `../../docs/sample_info.tsv` 生成）。

> **一句话**：这两个 CSV 不是给人看的样本表，而是 **nf-core/sarek 的机器输入**，列名和 `patient` 列的含义都被 sarek 写死，不能按客户原表改。客户一致的命名放在 `../../docs/sample_info.tsv`（含逐字照抄的 `Sample_Name`/`Sample_Type`）。

## 1. 列名是固定 schema，改了 sarek 直接报错
sarek 要求且只认这 6 列：`patient, sample, status, lane, fastq_1, fastq_2`。
- 改成 `Sample_Name`/`Sample_Type` 之类 → schema 校验失败、流程拒跑。
- 所以列名不是我的风格选择，是 sarek 的硬性要求。

## 2. `patient` 是"配对分组键"，不是"病人"，而且**必须共享**
sarek somatic 靠 `patient` 把 **tumor 和它的 normal 归到一组**来配对：
- 同一个 `patient` 下，`status=0`(normal) 与 `status=1`(tumor) 会被 sarek 自动配对，跑 tumor-vs-normal 体细胞调用（Mutect2）。
- **Study A 全部 `patient=RO`**，是**故意**的：RO 谱系里 `RO_origin`=normal(0)，`B1TP/B2TP/3个tumor`=tumor(1)，sarek 就会拿每个 tumor 去和 `RO_origin` 配对。
- 若把 `patient` 改成每样唯一（3852_RO_origin、3868…）→ 每个样本各自成组、**组里没有 normal 可配 → 体细胞调用直接失效**。
- 补充：`RO` 本就是**客户自己的谱系标签**（客户名里 `RO_origin`、`1st_RO tumor` 都带 RO），所以 `patient=RO` 既满足配对、也和客户一致。

## 3. 行的区分靠 `sample` 列，不是 `patient`
每行由 **`sample`**（RO_origin/RO_B1TP/…，各不相同）区分并命名输出文件夹。`patient` 相同不代表"无法区分"——区分是 `sample` 的职责，`patient` 的职责是"谁和谁配对"。

## 4. Study B 为什么每样一个 patient
Study B 是 **germline**（无配对正常），不需要配对分组，所以每样 `patient=自身`、`status=0`。

## 5. 客户一致的命名放哪
`../../docs/sample_info.tsv` 是权威人读表，已含客户原表逐字两列 **`Sample_Name`**（如 `3852_RO_origin`）、**`Sample_Type`**（Cell/Tissue）。报告与对外沟通用这些；sarek CSV 只管喂机器。

---
*简言之：sarek CSV = 机器输入（schema 固定、patient=配对组）；客户一致命名 = sample_info.tsv。两者各司其职。*
