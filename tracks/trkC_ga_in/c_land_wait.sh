#!/bin/bash
# trkC_gain: score GA with the LANDIS h1 slot (P2) after the stress lane, CEM paused.
T=~/perseus_me/tracks/trkC_gain; CEMPG=772744; source ~/perseus_me/run/env.sh >/dev/null 2>&1
export OMP_NUM_THREADS=1 NCASI_RC_STATES=$HOME/perseus_me/tracks/multistate/config/multistate/states_ms_inputs.yaml
while pgrep -f "tools/c_stress_wait.sh" > /dev/null; do sleep 60; done
kill -STOP -- -$CEMPG 2>/dev/null; echo "CEM paused $(date -Iseconds)"
cd ~/perseus_me/tracks/multistate; nice -n 10 env HOME=$T/restricted/home python $T/tools/c_score_x2.py GA $T/restricted/score/score_landis_GA $T/restricted/slots > $T/logs/score_landis_GA.log 2>&1; echo "score landis GA rc=$? $(date -Iseconds)"
echo "- GA scored with LANDIS h1 slot $(date -Iseconds) (restricted/score/score_landis_GA)" >> $T/release/STATUS.md
kill -CONT -- -$CEMPG 2>/dev/null; echo "CEM resumed $(date -Iseconds)"
