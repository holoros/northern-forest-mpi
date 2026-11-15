#!/bin/bash
# trkC_gain: rerun the CONDPROP stress lane (fixed SIT pipeline, D-C10) after the LANDIS GA h1 lane, with the CEM lane paused.
T=~/perseus_me/tracks/trkC_gain; CEMPG=772744
while pgrep -f "tools/c_landis_h1_lane.sh" > /dev/null; do sleep 60; done
sleep 5; kill -STOP -- -$CEMPG 2>/dev/null; echo "CEM paused $(date -Iseconds)"
bash $T/tools/c_stress_cp.sh > $T/logs/stress_lane.log 2>&1; echo "stress rc=$? $(date -Iseconds)"; echo "- stress lane v2 ended $(date -Iseconds) (logs/stress_lane.log)" >> $T/release/STATUS.md
kill -CONT -- -$CEMPG 2>/dev/null; echo "CEM resumed $(date -Iseconds)"
