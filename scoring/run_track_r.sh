#!/bin/bash
# Track R scoring driver (do not start before the runs finish; the driver also refuses states with live runs unless --force).
#   bash tools/ncasi_rc/score/run_track_r.sh --check-only          # readiness per state, no scoring
#   setsid nohup bash tools/ncasi_rc/score/run_track_r.sh --states NH,VT > ~/perseus_me/ncasi_rc/scores/drv_$(date +%Y%m%d_%H%M).log 2>&1 < /dev/null &
set -uo pipefail; umask 077
source ~/perseus_me/run/env.sh >/dev/null 2>&1
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
mkdir -p ~/perseus_me/ncasi_rc/scores; chmod 700 ~/perseus_me/ncasi_rc/scores
cd ~/perseus_me/code && python tools/ncasi_rc/score/run_track_r.py "$@"
