#!/bin/bash
# G1: reproduce the released multistate IN rows (work/final/score_IN/ms_scores.csv) with the frozen ms_score_x.py under the overlay HOME.
set -uo pipefail; umask 077
T=/home/aaron/perseus_me/tracks/trkC_gain; source /home/aaron/perseus_me/run/env.sh >/dev/null 2>&1
export HOME=$T/restricted/home OMP_NUM_THREADS=1 NCASI_RC_STATES=/home/aaron/perseus_me/tracks/multistate/config/multistate/states_ms_inputs.yaml
cd /home/aaron/perseus_me/tracks/multistate
md5sum tools/multistate/ms_score_x.py tools/multistate/ms_score.py
nice -n 10 python tools/multistate/ms_score_x.py IN $T/restricted/g1_IN; echo "rc=$?"
python - <<'PY'
import pandas as pd, numpy as np, json
a=pd.read_csv("/home/aaron/perseus_me/multistate/work/final/score_IN/ms_scores.csv"); b=pd.read_csv("/home/aaron/perseus_me/tracks/trkC_gain/restricted/g1_IN/ms_scores.csv")
k=["rule_set","basis","view","row"]; J=a.merge(b,on=k,how="outer",suffixes=("","_r"),indicator=True)
num=[c for c in a.columns if a[c].dtype.kind=="f"]; d=0.0
for c in num:
    x,y=J[c].values,J[c+"_r"].values; m=np.isfinite(x)&np.isfinite(y); d=max(d,float(np.abs(x[m]-y[m]).max()) if m.any() else 0); d=max(d, np.inf if (np.isfinite(x)!=np.isfinite(y)).any() else 0)
r=dict(rows_released=len(a),rows_rerun=len(b),only_one_side=int((J._merge!="both").sum()),max_abs_diff=d,passed=bool(d<=1e-9 and (J._merge=="both").all()))
open("/home/aaron/perseus_me/tracks/trkC_gain/work/g1_IN.json","w").write(json.dumps(r,indent=1)); print(r)
PY
