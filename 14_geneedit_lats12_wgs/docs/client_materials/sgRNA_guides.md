# CRISPR sgRNA 指导序列（客户 Jinpeng 提供，邮件 2026-07-07 附图转录）

来源图片：`8c593d5666fc5c36cf80e9c25bb56170.png`(Pten)、`55e5ea53d89200c141e5a53ca52b057e.png`(Brca1)。
产品：CRISPRevolution sgRNA EZ Kit (Modified)，每条 1.5 nmol。
序列图上为 RNA（含 U）；下表转成 DNA 靶序列（20-nt spacer，用于在 GRCm39 定位切点，需配 NGG PAM，正反链均搜）。

| 靶基因 | 名称/ID | spacer (DNA, 5'→3') | 用于 |
| :--- | :--- | :--- | :---: |
| Pten | Pten-32799878 | `GGTGGGTTATGGTCTTCAAA` | B1TP + B2TP |
| Pten | Pten-32799895 | `TGATAAGTTCTAGCTGTGGT` | B1TP + B2TP |
| Pten | Pten-32799899 | `GGTTTGATAAGTTCTAGCTG` | B1TP + B2TP |
| Brca1 | (guide 1) | `GGTTCCGGTAGCCCACGCTC` | B1TP |
| Brca1 | (guide 2) | `GGCGTCGATCATCCAGAGCG` | B1TP |
| Brca1 | (guide 3) | `TTCTTGTGAGCGTTTGAATG` | B1TP |
| **Brca2** | — | **未提供**（此封只给 Brca1+Pten） | B2TP |

**用法**：编辑验证时，在 GRCm39 上 BLAST/grep 每条 spacer（+反向互补）定位靶点，预测切点 = PAM 上游 3 bp；在该位点比 B1TP/B2TP/肿瘤 vs RO_origin 查 indel/frameshift（CRISPR KO 是否发生、等位比例）。
**待补**：Brca2 的 3 条 sgRNA（B2TP 的 Brca2 精确切点需要；无则先扫 Brca2 全基因 indel）。
