#!/bin/bash
# ============================================================================
# Study B HaplotypeCaller 资源监控(tmux 运行,解耦连接)
#   每 5 min 记录:整机负载/内存/swap + HC 实测并发数 + 已完成 HC 任务的 %cpu/峰值RSS。
#   目的:验证满机配置里 HAPLOTYPECALLER cpus=4/memory=7GB/≈14并发 的假设是否成立;
#         若实测 RSS 明显 >7GB 或并发被内存卡住,据此调配置。
#   自动停止:B 的 session 消失(跑完/退出)时退出。
# ============================================================================
PROJ="/home/gao/projects_2026H2/14_geneedit_lats12_wgs"
MLOG="$PROJ/logs/hc_resource_monitor_0712.log"
trace(){ ls -t "$PROJ"/output_B/pipeline_info/execution_trace_*.txt 2>/dev/null | head -1; }
echo "[$(date '+%m-%d %H:%M')] HC 资源监控启动" | tee -a "$MLOG"
while tmux has-session -t p14_sarek_B 2>/dev/null; do
  ld=$(cut -d' ' -f1-3 /proc/loadavg)
  mem=$(free -g | awk '/Mem/{print "avail="$7"G used="$3"G"} /Swap/{printf " swap="$3"G"}')
  # 活跃 HC:GATK HaplotypeCaller 的 java 进程(命令行含 HaplotypeCaller)
  hc_live=$(pgrep -fc "HaplotypeCaller" 2>/dev/null || echo 0)
  # 已完成 HC 任务的资源实测(从 trace)
  t=$(trace)
  hcstat=""
  if [ -n "$t" ]; then
    hcstat=$(awk -F'\t' 'NR>1 && $4 ~ /HAPLOTYPECALLER/ && $5=="COMPLETED"{
        cpu=$10; sub(/%/,"",cpu); r=$11; gsub(/[^0-9.]/,"",r); u=$11; gsub(/[0-9. ]/,"",u);
        gb=(u=="GB")?r:(u=="MB")?r/1024:r/1048576;
        n++; sc+=cpu; if(gb>mx)mx=gb }
      END{ if(n>0) printf "HC完成=%d 均%%cpu=%.0f(≈%.1f核) 峰值RSS_max=%.1fG", n, sc/n, sc/n/100, mx; else print "HC未完成任务" }' "$t")
  fi
  echo "[$(date '+%m-%d %H:%M')] load=$ld | $mem | HC活跃=$hc_live | $hcstat" | tee -a "$MLOG"
  sleep 300
done
echo "[$(date '+%m-%d %H:%M')] B session 结束,监控退出" | tee -a "$MLOG"
