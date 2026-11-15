#!/bin/bash
# trkC_gain: CBM gated GA/IN finishing lane on firebreather, one worker for the track: the CEM lane (pgid in $1) is paused (SIGSTOP)
# while this runs and resumed on exit. Steps: ESI at public coordinates, statewide C1/P5, X1 fits (+ permuted refit, G2 for fits),
# scoring GA and IN (c_score_x2.py), G2 row-shuffle rescoring.
set -uo pipefail; umask 077; source ~/perseus_me/run/env.sh >/dev/null 2>&1
CEMPG=${1:-0}; T=~/perseus_me/tracks/trkC_gain; S=$T/restricted/slots; L=$T/logs
export OMP_NUM_THREADS=1 NCASI_RC_STATES=$HOME/perseus_me/tracks/multistate/config/multistate/states_ms_inputs.yaml
resume() { [ "$CEMPG" != 0 ] && kill -CONT -- -$CEMPG 2>/dev/null; echo "CEM resumed $(date -Iseconds)"; }
trap resume EXIT
[ "$CEMPG" != 0 ] && kill -STOP -- -$CEMPG && echo "CEM paused $(date -Iseconds)"
step() { n=$1; shift; [ -f $T/work/done_$n ] && { echo "skip $n"; return 0; }; echo "start $n $(date -Iseconds)"; nice -n 10 "$@" > $L/cbm_$n.log 2>&1; rc=$?; echo "end $n rc=$rc $(date -Iseconds)"; [ $rc = 0 ] && touch $T/work/done_$n; return $rc; }
step esi env R_LD_LIBRARY_PATH=$HOME/opt/geos/lib:/usr/lib/R/lib:/usr/lib/x86_64-linux-gnu Rscript $T/tools/c_esi.R || exit 1
for st in GA IN; do step sw_$st python $T/tools/c_statewide_curves.py $st $S/statewide/$st || exit 1; done
for st in GA IN; do step x1_$st python $T/tools/c_x1.py $st $S || exit 1; step x1perm_$st python $T/tools/c_x1.py $st $S --permute 20261114; done
cd ~/perseus_me/tracks/multistate
for st in IN GA; do step score_$st env HOME=$T/restricted/home python $T/tools/c_score_x2.py $st $T/restricted/score/score_$st $S; done
for st in IN GA; do step shuffle_$st env HOME=$T/restricted/home python $T/tools/c_score_x2.py $st $T/restricted/score/shuffle_$st $S --shuffle 20261114; done
echo CBM_LANE_DONE $(date -Iseconds)
