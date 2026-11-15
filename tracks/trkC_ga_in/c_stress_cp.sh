#!/bin/bash
# trkC_gain stress variants of the D-C2 growth-point rule (coordinator instruction 14 Nov): cp100 = CONDPROP_UNADJ = 1 only;
# cpw = all conditions, nlsLM weights = CONDPROP_UNADJ. First a cross-host check: the firebreather SIT pipeline on the D-C2 sit_yield
# must reproduce the Cardinal adapted growth_and_yield.csv and inventory.csv. Then per variant: fit, SIT, curve table, unit rows,
# statewide, scoring (c_score_x2.py with only that variant's cbm_gated). One process at a time, nice 10.
set -uo pipefail; umask 077; source ~/perseus_me/run/env.sh >/dev/null 2>&1
T=~/perseus_me/tracks/trkC_gain; R=$T/restricted; CP=$R/cardinal_pull/out; L=$T/logs; PY=$HOME/perseus_me/ncasi_rc/engines/pycbm
export OMP_NUM_THREADS=1 NCASI_RC_STATES=$HOME/perseus_me/tracks/multistate/config/multistate/states_ms_inputs.yaml
step() { n=$1; shift; [ -f $T/work/done_$n ] && { echo "skip $n"; return 0; }; echo "start $n $(date -Iseconds)"; nice -n 10 "$@" > $L/stress_$n.log 2>&1; rc=$?; echo "end $n rc=$rc $(date -Iseconds)"; [ $rc = 0 ] && touch $T/work/done_$n; return $rc; }
# cross-host check (D-C2 curves)
for st in GA IN; do step xhost_$st bash $T/tools/c_sit_2x2_fb.sh $st $CP/sit/sit_yield_${st}_default.csv $R/stress/xhost_$st $R/cbm_work_xhost; done
python - <<'P' > $T/work/xhost_sit.json
import hashlib, json, pathlib
R = pathlib.Path.home() / "perseus_me/tracks/trkC_gain/restricted"; out = {}
for st in ("GA", "IN"):
    for f in ("growth_and_yield.csv", "inventory.csv", "classifiers.csv"):
        a = R / f"cardinal_pull/out/sit2x2_{st}/adapted_FVS_FIAW_NEW/{f}"; b = R / f"stress/xhost_{st}/adapted_FVS_FIAW_NEW/{f}"
        out[f"{st}/{f}"] = (hashlib.md5(a.read_bytes()).hexdigest() == hashlib.md5(b.read_bytes()).hexdigest()) if b.exists() else "missing"
print(json.dumps(out, indent=1))
P
cat $T/work/xhost_sit.json
python -c "import json,sys; d=json.load(open('$T/work/xhost_sit.json')); sys.exit(0 if all(v is True for v in d.values()) else 1)" || { echo "XHOST CHECK FAILED: stop"; exit 7; }
# variant points
python - <<'P'
import pandas as pd, json, pathlib
R = pathlib.Path.home() / "perseus_me/tracks/trkC_gain/restricted"; M = pathlib.Path.home() / "perseus_me/multistate/state_roots"; rep = {}
for V in ("cp100", "cpw"):
    out = []
    for st in ("GA", "IN"):
        P = pd.read_csv(R / f"cardinal_pull/out/pts_{st}/fvs_growth_points_{st}.csv"); K = pd.read_csv(R / f"cardinal_pull/out/pts_{st}/stand_keys_{st}.csv", dtype=str)
        C = pd.read_parquet(M / st / "restricted/inventory/conds.parquet", columns=["PLT_CN", "CONDID", "CONDPROP_UNADJ"])
        C["cond_key"] = C.PLT_CN.astype(str) + "_" + C.CONDID.astype(str); C = C.drop_duplicates("cond_key")
        K = K.merge(C[["cond_key", "CONDPROP_UNADJ"]], on="cond_key", how="left", validate="one_to_one"); assert K.CONDPROP_UNADJ.notna().all()
        P = P.merge(K[["stand_id", "CONDPROP_UNADJ"]], on="stand_id", how="left", validate="many_to_one")
        Q = P[P.CONDPROP_UNADJ >= 1.0 - 1e-9] if V == "cp100" else P.assign(w=P.CONDPROP_UNADJ)
        out.append(Q.drop(columns=["CONDPROP_UNADJ"])); rep[f"{V}/{st}"] = dict(stands=int(Q.stand_id.nunique()), points=len(Q))
    pd.concat(out).to_csv(R / f"stress/pts_{V}.csv", index=False, float_format="%.10g")
print(json.dumps(rep)); (R.parent / "work/stress_points.json").write_text(json.dumps(rep, indent=1))
P
for V in cp100 cpw; do O=$R/stress/$V; mkdir -p $O/sit
  step fit_$V Rscript --vanilla $T/tools/c_fit_sit.R $T/tools/build_sit_yield_from_fvs_DRAFT_2026-09-17.R $R/stress/pts_$V.csv $O/sit $( [ $V = cpw ] && echo weighted )
  for st in GA IN; do f=$O/sit/sit_yield_${st}_default.csv; [ -f $f ] && awk -F, 'BEGIN{OFS=","} NR==1{print;next} {if($2=="")$2="NA"; print}' $f > $f.tmp && mv $f.tmp $f
    step sit_${V}_$st bash $T/tools/c_sit_2x2_fb.sh $st $f $O/sit2x2_$st $R/cbm_work_$V
    step curves_${V}_$st bash -c "cd ~/perseus_me/code && $PY tools/complete/cbm_curves.py cbm_gated $O/sit2x2_$st/adapted_FVS_FIAW_NEW $O/sit2x2_$st/cbm_defaults_$st.db $R/cbm_curves_$V/gated_$st 300"
  done
  for f in c_curves_pairs.py c_curves_hindcast.py c_statewide_curves.py; do sed -e "s#trkC_gain/restricted/cbm_curves/#trkC_gain/restricted/cbm_curves_$V/#" -e "s#trkC_gain/restricted/curve_roots#trkC_gain/restricted/curve_roots_$V#" $T/tools/$f > $T/tools/${f%.py}_$V.py; done
  S=$R/slots_$V
  for st in GA IN; do OO=$S/cbm_gated/$st/restricted; X=$R/stress/tmp_${V}_$st; mkdir -p $OO $X
    phi=$(awk -F, -v s=$st '$1==s{print $5}' ~/perseus_me/multistate/work/common/phi_states.csv)
    step units_${V}_$st bash -c "python $T/tools/c_curves_pairs_$V.py --state $st --set pairs_all --engines cbm_gated --out $X/pairs_level && python $T/tools/c_curves_pairs_$V.py --state $st --set treated --engines cbm_gated --out $X/treated_level && python $T/tools/c_curves_pairs_$V.py --state $st --set treated --engines cbm_gated --rule ageq --phi $phi --out $X/treated_ageq && for d in pairs_level treated_level treated_ageq; do for f in $X/\$d/cbm_gated__*pair_rows.csv; do cp \$f $OO/\${d}__\$(basename \$f); done; cp $X/\$d/run_meta.json $OO/\${d}__run_meta.json; done && python $T/tools/c_curves_hindcast_$V.py $st $OO/h1 && python $T/tools/c_statewide_curves_$V.py $st $S/statewide/$st"
  done
  cd ~/perseus_me/tracks/multistate
  for st in IN GA; do step score_${V}_$st env HOME=$T/restricted/home python $T/tools/c_score_x2.py $st $R/score/score_${V}_$st $S; done
done
echo STRESS_DONE $(date -Iseconds)
