#!/bin/bash
# trkC_gain track queue v2 (one worker; CEM lane paused while each step runs): G2 narrow shuffle rescoring (IN, GA) -> CONDPROP stress
# variants -> LANDIS GA h1 -> resume CEM; when the CEM lane ends, CEM slots and final scoring of every row.
set -uo pipefail; T=~/perseus_me/tracks/trkC_gain; CEMPG=772744; source ~/perseus_me/run/env.sh >/dev/null 2>&1
export OMP_NUM_THREADS=1 NCASI_RC_STATES=$HOME/perseus_me/tracks/multistate/config/multistate/states_ms_inputs.yaml
trap 'kill -CONT -- -$CEMPG 2>/dev/null' EXIT
kill -STOP -- -$CEMPG 2>/dev/null && echo "CEM paused $(date -Iseconds)"
cd ~/perseus_me/tracks/multistate
for st in IN GA; do [ -f $T/work/done_g2_$st ] && continue; nice -n 10 env HOME=$T/restricted/home python $T/tools/c_score_x2.py $st $T/restricted/score/shuffle_$st $T/restricted/slots --shuffle 20261114 > $T/logs/g2_$st.log 2>&1; echo "g2 $st rc=$? $(date -Iseconds)"
  ~/perseus_me/venv/bin/python $T/tools/c_checks_files.py >> $T/logs/g2_$st.log 2>&1 && touch $T/work/done_g2_$st; echo "- g2 $st done $(date -Iseconds) (work/g2_$st.json)" >> $T/release/STATUS.md; done
bash $T/tools/c_stress_cp.sh > $T/logs/stress_lane.log 2>&1; echo "stress lane rc=$? $(date -Iseconds)"; echo "- stress lane ended $(date -Iseconds) (logs/stress_lane.log)" >> $T/release/STATUS.md
bash $T/tools/c_landis_h1_lane.sh > $T/logs/landis_lane.log 2>&1; echo "landis lane rc=$? $(date -Iseconds)"; echo "- LANDIS GA h1 lane ended $(date -Iseconds)" >> $T/release/STATUS.md
kill -CONT -- -$CEMPG 2>/dev/null; echo "CEM resumed $(date -Iseconds)"
while kill -0 $CEMPG 2>/dev/null; do sleep 300; done
R=$T/restricted/cem; S=$T/restricted/slots/cem
for st in GA IN; do [ -s $R/$st/treated/extract/pair_rows.csv ] || continue; mkdir -p $S/$st/restricted
  [ -e $S/$st/restricted/treated ] || ln -s $R/$st/treated $S/$st/restricted/treated
  L2=~/perseus_me/multistate/runs/cem/$st/restricted/pairs_all; [ -s $L2/extract/pair_rows.csv ] && { [ -e $S/$st/restricted/pairs_all ] || ln -s $L2 $S/$st/restricted/pairs_all; }
done
for st in IN GA; do nice -n 10 env HOME=$T/restricted/home python $T/tools/c_score_x2.py $st $T/restricted/score/score_all_$st $T/restricted/slots > $T/logs/score_all_$st.log 2>&1; echo "score_all $st rc=$? $(date -Iseconds)"; done
~/perseus_me/venv/bin/python $T/tools/c_release.py > /dev/null 2>&1; ~/perseus_me/venv/bin/python $T/tools/c_release_md.py > /dev/null 2>&1; echo "- final scoring with CEM/LANDIS ended $(date -Iseconds); TABLES.md regenerated" >> $T/release/STATUS.md
echo MASTER2_DONE $(date -Iseconds)
