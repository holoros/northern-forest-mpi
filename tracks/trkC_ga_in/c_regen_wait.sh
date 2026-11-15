#!/bin/bash
# trkC_gain: regenerate release tables, file checks, second-path verification and the claim ledger after the stress and LANDIS steps.
T=~/perseus_me/tracks/trkC_gain; PY=~/perseus_me/venv/bin/python
while pgrep -f "tools/c_land_wait.sh" > /dev/null || pgrep -f "tools/c_stress_wait.sh" > /dev/null; do sleep 60; done
cd $T; nice -n 10 $PY tools/c_checks_files.py > logs/regen.log 2>&1; nice -n 10 $PY tools/c_verify.py GA IN >> logs/regen.log 2>&1
nice -n 10 $PY tools/c_release.py >> logs/regen.log 2>&1; nice -n 10 $PY tools/c_release_md.py >> logs/regen.log 2>&1; nice -n 10 $PY tools/c_claims.py >> logs/regen.log 2>&1
echo "- release tables, checks, verifier and claims regenerated $(date -Iseconds)" >> release/STATUS.md
