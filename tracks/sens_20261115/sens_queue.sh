#!/bin/bash
# sequential lane: one driver process at a time, nice 10, load gate 30; then the post step (no v3 write)
S=~/perseus_me/tracks/ncasi_ms/sens_20261115; source ~/perseus_me/run/env.sh >/dev/null 2>&1
export OMP_NUM_THREADS=1 TMPDIR=$S/restricted/tmp; unset RC_REGISTER; cd ~/perseus_me/tracks/multistate
while pgrep -f "sens_driver.py" >/dev/null; do sleep 20; done
for st in "$@"; do
  [ -s $S/scores/$st/checks.json ] && continue
  while awk '{exit !($1>30)}' /proc/loadavg; do sleep 60; done
  nice -n 10 python $S/sens_driver.py $st > $S/restricted/logs/sens_$st.log 2>&1; rc=$?
  echo "$(date -Is) driver $st rc=$rc" >> $S/STATUS.md
  ~/perseus_me/venv_fig/bin/python $S/sens_post.py > $S/logs/post_after_$st.log 2>&1
done
echo "$(date -Is) queue done ($*)" >> $S/STATUS.md
