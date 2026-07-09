你自己跑的 demo 流程还在运行中（还没生成 outs/ 最终结果目录），目前的输出路径结构是：

总目录：/Work_bio/gao_test/cellranger_demo_rerun/
QC结果：/Work_bio/gao_test/cellranger_demo_rerun/01.QC/N2-SI_TT_D1_23575YLT4/（已完成，含 html/json报告和过滤后的fastq）
Cell Ranger VDJ结果（运行中，完成后会生成）：/Work_bio/gao_test/cellranger_demo_rerun/N2-SI_TT_D1_23575YLT4/outs/
即最终 web_summary.html 会在：/Work_bio/gao_test/cellranger_demo_rerun/N2-SI_TT_D1_23575YLT4/outs/web_summary.html
实时日志：/Work_bio/gao_test/cellranger_demo_rerun/run.log