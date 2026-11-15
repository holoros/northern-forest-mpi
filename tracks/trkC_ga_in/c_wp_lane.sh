#!/bin/bash
# trkC_gain Amendment 1 (md5 bca7a84d): whole-plot CBM gated rebuild and scoring, chained after the track queue (c_master2.sh, which
# ends after the CEM lane and its scoring). One worker, nice 10, firebreather only. Each step appends to release/STATUS.md.
set -uo pipefail; umask 077
T=~/perseus_me/tracks/trkC_gain; R=$T/restricted; W=$R/wp; L=$T/logs; PY=$HOME/perseus_me/ncasi_rc/engines/pycbm; VPY=~/perseus_me/venv/bin/python
MS=~/perseus_me/multistate; WT=~/perseus_me/tracks/multistate
while pgrep -f "tools/c_master2.sh" > /dev/null; do sleep 300; done
source ~/perseus_me/run/env.sh >/dev/null 2>&1
export OMP_NUM_THREADS=1 TMPDIR=$MS/restricted/tmp NCASI_RC_STATES=$WT/config/multistate/states_ms_inputs.yaml PYTHONPATH=$T/tools
[ "$(md5sum < $T/PREREG_AMENDMENT_1.md | cut -c1-32)" = bca7a84d1746f6f95c6239df6d187138 ] || { echo "amendment md5 changed: stop"; exit 9; }
mkdir -p $W; chmod 700 $W
st_() { echo "- $(date -Iseconds) wp: $*" >> $T/release/STATUS.md; }
step() { n=$1; shift; [ -f $T/work/done_wp_$n ] && { echo "skip $n"; return 0; }; echo "start $n $(date -Iseconds)"; nice -n 10 "$@" > $L/wp_$n.log 2>&1; rc=$?; echo "end $n rc=$rc $(date -Iseconds)"; [ $rc = 0 ] && touch $T/work/done_wp_$n && st_ "$n done"; [ $rc = 0 ] || st_ "$n FAILED rc=$rc (logs/wp_$n.log)"; return $rc; }
st_ "lane started (after the track queue)"
cd $WT
for s in GA IN; do SR=$MS/state_roots/$s
  step pts_$s bash -c "$( [ $s = IN ] && echo "source $MS/work/ms_env_cs.sh;" ) env PME_STATE=$s python $WT/tools/multistate/ms_run.py c_curvepts_wp -- --inv $SR/restricted/inventory --aoi $SR/aoi/aoi.json --host $MS/work/host_ms.yaml --state $s --ncyc 10 --workers 1 --out $W/pts_$s" || exit 1
done
head -1 $W/pts_GA/fvs_growth_points_GA.csv > $W/pts_GAIN.csv; tail -n +2 -q $W/pts_GA/fvs_growth_points_GA.csv $W/pts_IN/fvs_growth_points_IN.csv >> $W/pts_GAIN.csv
mkdir -p $W/sit; step fit Rscript --vanilla $T/tools/c_fit_sit.R $T/tools/build_sit_yield_from_fvs_DRAFT_2026-09-17.R $W/pts_GAIN.csv $W/sit || exit 1
for s in GA IN; do f=$W/sit/sit_yield_${s}_default.csv; awk -F, 'BEGIN{OFS=","} NR==1{print;next} {if($2=="")$2="NA"; print}' $f > $f.tmp && mv $f.tmp $f
  step sit_$s bash $T/tools/c_sit_2x2_fb.sh $s $f $W/sit2x2_$s $R/cbm_work_wp || exit 1
  step curves_$s bash -c "cd ~/perseus_me/code && $PY tools/complete/cbm_curves.py cbm_gated $W/sit2x2_$s/adapted_FVS_FIAW_NEW $W/sit2x2_$s/cbm_defaults_$s.db $R/cbm_curves_wp/gated_$s 300" || exit 1
done
for f in c_curves_pairs.py c_curves_hindcast.py c_statewide_curves.py; do sed -e "s#trkC_gain/restricted/cbm_curves/#trkC_gain/restricted/cbm_curves_wp/#" -e "s#trkC_gain/restricted/curve_roots#trkC_gain/restricted/curve_roots_wp#" $T/tools/$f > $T/tools/${f%.py}_wp.py; done
S=$R/slots_wp; mkdir -p $S; [ -e $S/landis ] || ln -s $R/slots/landis $S/landis; [ -e $R/slots/cem ] && { [ -e $S/cem ] || ln -s $R/slots/cem $S/cem; }
for s in GA IN; do OO=$S/cbm_gated/$s/restricted; X=$W/tmp_$s; mkdir -p $OO $X
  phi=$(awk -F, -v t=$s '$1==t{print $5}' $MS/work/common/phi_states.csv)
  step units_$s bash -c "python $T/tools/c_curves_pairs_wp.py --state $s --set pairs_all --engines cbm_gated --out $X/pairs_level && python $T/tools/c_curves_pairs_wp.py --state $s --set treated --engines cbm_gated --out $X/treated_level && python $T/tools/c_curves_pairs_wp.py --state $s --set treated --engines cbm_gated --rule ageq --phi $phi --out $X/treated_ageq && for d in pairs_level treated_level treated_ageq; do for f in $X/\$d/cbm_gated__*pair_rows.csv; do cp \$f $OO/\${d}__\$(basename \$f); done; cp $X/\$d/run_meta.json $OO/\${d}__run_meta.json; done && python $T/tools/c_curves_hindcast_wp.py $s $OO/h1 && python $T/tools/c_statewide_curves_wp.py $s $S/statewide/$s" || exit 1
  step x1_$s python $T/tools/c_x1_wp.py $s $S || exit 1
  step x1perm_$s python $T/tools/c_x1_wp.py $s $S --permute 20261114
done
for s in IN GA; do step score_$s env HOME=$R/home python $T/tools/c_score_x2.py $s $R/score/score_wp_$s $S
  step shuffle_$s env HOME=$R/home python $T/tools/c_score_x2.py $s $R/score/shufflewp_$s $S --shuffle 20261114; done
$VPY - <<'P' >> $L/wp_checks.log 2>&1
import json, numpy as np, pandas as pd
from pathlib import Path
T = Path.home() / "perseus_me/tracks/trkC_gain"; K = ["rule_set", "basis", "view", "row"]; out = {}
for st in ("GA", "IN"):
    a, b = T / f"restricted/score/score_wp_{st}/ms_scores.csv", T / f"restricted/score/shufflewp_{st}/ms_scores.csv"
    if a.exists() and b.exists():
        x, y = pd.read_csv(a), pd.read_csv(b); J = x.merge(y, on=K, suffixes=("", "_s"), validate="one_to_one"); d = 0.0
        for c in [c for c in x.columns if x[c].dtype.kind == "f" and not c.startswith("P_rank1")]:
            u, v = J[c].values, J[c + "_s"].values; m = np.isfinite(u) & np.isfinite(v)
            d = max(d, float(np.abs(u[m] - v[m]).max()) if m.any() else 0.0, np.inf if (np.isfinite(u) != np.isfinite(v)).any() else 0.0)
        out[st] = dict(rows=len(J), max_abs_diff=d, passed=bool(d <= 1e-9 and len(J) == len(x) == len(y)))
(T / "work/g2_wp.json").write_text(json.dumps(out, indent=1)); print(out)
P
st_ "G2 whole-plot: $(cat $T/work/g2_wp.json | tr -d '\n ')"
cd $T; TRKC_TAG=wp nice -n 10 $VPY tools/c_verify.py GA IN > $L/wp_verify.log 2>&1; st_ "second path whole-plot rc=$? (work/verifier_wp_ST.json)"
nice -n 10 $VPY tools/c_checks_files.py >> $L/wp_checks.log 2>&1; nice -n 10 $VPY tools/c_release.py > $L/wp_release.log 2>&1; nice -n 10 $VPY tools/c_release_md.py >> $L/wp_release.log 2>&1
nice -n 10 $VPY tools/c_claims.py >> $L/wp_release.log 2>&1; st_ "release tables (TABLES.md, trkC_before_after.csv, trkC_c5_three_ratios.csv) and claims.csv regenerated with whole-plot rows"
echo WP_LANE_DONE $(date -Iseconds)
