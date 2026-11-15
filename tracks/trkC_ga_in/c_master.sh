#!/bin/bash
# trkC_gain track queue (one worker): after the CBM lane, keep the CEM lane paused and run the CONDPROP stress variants then LANDIS GA h1
# variants; then resume CEM; when the CEM lane ends, build the CEM slots (treated rows of this track; C2 from the multistate lane
# if its extraction exists) and score GA/IN with every row.
set -uo pipefail; T=~/perseus_me/tracks/trkC_gain; CEMPG=772744
while pgrep -f c_cbm_lane.sh > /dev/null; do sleep 60; done
kill -STOP -- -$CEMPG 2>/dev/null && echo "CEM paused $(date -Iseconds)"
bash $T/tools/c_stress_cp.sh > $T/logs/stress_lane.log 2>&1; echo "stress lane rc=$? $(date -Iseconds)"
bash $T/tools/c_landis_h1_lane.sh > $T/logs/landis_lane.log 2>&1; echo "landis lane rc=$? $(date -Iseconds)"
kill -CONT -- -$CEMPG 2>/dev/null; echo "CEM resumed $(date -Iseconds)"
while pgrep -f c_cem_lane.sh > /dev/null; do sleep 300; done
source ~/perseus_me/run/env.sh >/dev/null 2>&1; export OMP_NUM_THREADS=1 NCASI_RC_STATES=$HOME/perseus_me/tracks/multistate/config/multistate/states_ms_inputs.yaml
R=$T/restricted/cem; S=$T/restricted/slots/cem
for st in GA IN; do [ -s $R/$st/treated/extract/pair_rows.csv ] || continue; mkdir -p $S/$st/restricted
  [ -e $S/$st/restricted/treated ] || ln -s $R/$st/treated $S/$st/restricted/treated
  L2=~/perseus_me/multistate/runs/cem/$st/restricted/pairs_all; [ -s $L2/extract/pair_rows.csv ] && { [ -e $S/$st/restricted/pairs_all ] || ln -s $L2 $S/$st/restricted/pairs_all; }
done
cd ~/perseus_me/tracks/multistate
for st in IN GA; do nice -n 10 env HOME=$T/restricted/home python $T/tools/c_score_x2.py $st $T/restricted/score/score_all_$st $T/restricted/slots > $T/logs/score_all_$st.log 2>&1; echo "score_all $st rc=$? $(date -Iseconds)"; done
echo MASTER_DONE $(date -Iseconds)
