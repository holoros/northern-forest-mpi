#!/bin/bash
cd ~/perseus_me/tracks/multistate; source ~/perseus_me/run/env.sh >/dev/null 2>&1; MS=~/perseus_me/multistate; SR=$MS/state_roots/GA
export OMP_NUM_THREADS=1 TMPDIR=$MS/restricted/tmp NCASI_RC_STATES=$PWD/config/multistate/states_ms_inputs.yaml PYTHONPATH=~/perseus_me/tracks/trkC_gain/tools
nice -n 10 env PME_STATE=GA python tools/multistate/ms_run.py c_curvepts_wp -- --inv $SR/restricted/inventory --aoi $SR/aoi/aoi.json --host $MS/work/host_ms.yaml --state GA --ncyc 10 --workers 1 --limit 3 --out ~/perseus_me/tracks/trkC_gain/restricted/smoke_wp
echo rc=$?
