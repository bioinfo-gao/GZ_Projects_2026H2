# CRISPR sgRNA 指导序列（客户 Jinpeng 提供，邮件附图转录）

来源图片：`8c593d5666fc5c36cf80e9c25bb56170.png`(Pten)、`55e5ea53d89200c141e5a53ca52b057e.png`(Brca1)、
`328911310515d5285fd5af28c9cf8b73.png`(**Brca2**，邮件 2026-07-16 补齐)。
产品：CRISPRevolution sgRNA EZ Kit (Modified)，每条 1.5 nmol。
序列图上为 RNA（含 U）；下表转成 DNA 靶序列（20-nt spacer，用于在 GRCm39 定位切点，需配 NGG PAM，正反链均搜）。

⚠ **本表是 `scripts/study_A/A3_edit_verification.sh` 的输入**（脚本 grep 本表行、抽 `[ACGT]{20}` spacer）。
改列序/改格式会静默破坏编辑验证流程——改动请同步脚本。

| 靶基因 | 名称/ID | spacer (DNA, 5'→3') | 用于 |
| :--- | :--- | :--- | :---: |
| Pten | Pten-32799878 | `GGTGGGTTATGGTCTTCAAA` | B1TP + B2TP |
| Pten | Pten-32799895 | `TGATAAGTTCTAGCTGTGGT` | B1TP + B2TP |
| Pten | Pten-32799899 | `GGTTTGATAAGTTCTAGCTG` | B1TP + B2TP |
| Brca1 | (guide 1) | `GGTTCCGGTAGCCCACGCTC` | B1TP |
| Brca1 | (guide 2) | `GGCGTCGATCATCCAGAGCG` | B1TP |
| Brca1 | (guide 3) | `TTCTTGTGAGCGTTTGAATG` | B1TP |
| Brca2 | Brca2+150529497 | `GATAAGCCTCAATTGGTTTG` | B2TP |
| Brca2 | Brca2-150529492 | `AAAGCTCCTCAAACCAATTG` | B2TP |
| Brca2 | Brca2-150529524 | `AGGTTCAGAATTGTATGGGG` | B2TP |

**用法**：编辑验证时，在 GRCm39 上 BLAST/grep 每条 spacer（+反向互补）定位靶点，预测切点 = PAM 上游 3 bp；在该位点比 B1TP/B2TP/肿瘤 vs RO_origin 查 indel/frameshift（CRISPR KO 是否发生、等位比例）。

---

## Brca2 三条 guide 的落位核实（2026-07-16，我方实测，非客户口述）

三条 spacer 在 GRCm39 上**各命中且仅命中一次**，链向与客户命名的 `+`/`-` 完全一致，且**三条都带完美 NGG PAM** →
转录无误、guide 真实有效：

| guide | GRCm39 命中 (1-based) | 链 | protospacer + PAM | 预测切点 (PAM 上游 3bp) |
| :--- | :---: | :---: | :---: | :---: |
| Brca2+150529497 | chr5:150452945-150452964 | + | `GATAAGCCTCAATTGGTTTG`+`AGG` | chr5:150452961 |
| Brca2-150529492 | chr5:150452954-150452973 | − | `AAAGCTCCTCAAACCAATTG`+`AGG` | chr5:150452957 |
| Brca2-150529524 | chr5:150452986-150453005 | − | `AGGTTCAGAATTGTATGGGG`+`GGG` | chr5:150452989 |

**设计解读**：g7(+) 与 g8(−) 相互重叠、反向靶同一位点（切点仅差 4 bp），g9(−) 在其下游 ~32 bp。
三个切点挤在 **~33 bp 窗口**内（150452957–150452989）→ 典型 **multi-guide 同位点轰击**设计，
期望产出 indel 或小片段切除（g7/g8 与 g9 之间 ~30 bp 缺失），而非大片段删除。

**⚠ 客户命名坐标 = GRCm38/mm10，不是 GRCm39**：客户 ID 里的 `150529497` 等比 GRCm39 实测位置
（`150452945`）高 **~76.5 kb**，与 Brca2 位点 GRCm38→GRCm39 的坐标平移一致（mm10 Brca2 ≈ chr5:150.53 Mb，
GRCm39 ≈ chr5:150.45 Mb）。**本项目一律按序列定位、不按客户坐标**，故此差异不影响分析——
但若日后有人直接拿 `150529497` 去 GRCm39 取序列，会取到错误位点，特此记录。
