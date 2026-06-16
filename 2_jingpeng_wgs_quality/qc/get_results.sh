for j in /home/gao/projects_2026H2/2_jingpeng_wgs_quality/qc_out/*_fastp.json; do
  sample=$(basename "$j" _fastp.json)
  echo "=== $sample ==="
  python3 -c "
import json
d = json.load(open('$j'))
s = d['summary']
bf = d['read1_before_filtering']
print(f\"  Total reads (sampled): {s['before_filtering']['total_reads']}\")
print(f\"  Q20 rate: {s['before_filtering']['q20_rate']}\")
print(f\"  Q30 rate: {s['before_filtering']['q30_rate']}\")
print(f\"  GC content: {s['before_filtering']['gc_content']}\")
print(f\"  Duplication rate: {d.get('duplication',{}).get('rate','N/A')}\")
print(f\"  Adapter trimmed: {d.get('adapter_cutting',{}).get('adapter_trimmed_reads','N/A')}\")
print(f\"  Insert size peak: {d.get('insert_size',{}).get('peak','N/A')}\")
print(f\"  After filter reads: {s['after_filtering']['total_reads']}\")
print(f\"  After Q30 rate: {s['after_filtering']['q30_rate']}\")
"
  echo ""
done