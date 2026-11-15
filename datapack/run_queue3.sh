#!/bin/bash
# continuation queue (replaces run_queue2.sh): waits for the running ME step, then Table 1, then MI WI MN (FAST_P5=1),
# then the fig2 combine and SOURCES.md. One worker, nice 10. Skips any step whose output exists.
source ~/perseus_me/run/env.sh
cd ~/perseus_me/tracks/ncasi_ms/tools
while pgrep -f 's2_components.py ME' >/dev/null; do sleep 30; done
[ -f ../data/_work/comp_ME.csv ] || nice -n 10 python s2_components.py ME > logs/s2_ME.log 2>&1
[ -f ../data/table1_data_overview.csv ] || nice -n 10 python s4_table1.py > logs/s4.log 2>&1
for st in MI WI MN; do
  [ -f ../data/_work/comp_${st}.csv ] || FAST_P5=1 nice -n 10 python s2_components.py $st > logs/s2_${st}.log 2>&1
done
n=$(ls ../data/_work/comp_??.csv 2>/dev/null | wc -l)
if [ "$n" = 6 ]; then nice -n 10 python s6_fig2_combine.py > logs/s6.log 2>&1; fi
nice -n 10 python s7_sources.py > logs/s7.log 2>&1
echo QUEUE3_DONE $(date) >> logs/queue3.log
