#!/bin/bash
# NCASI report card (task #5 setup): CEM (shared NCASI build) one-cycle runs on the t1 plots of FIA pairs, any state. No scoring.
# Generalized from tools/v24_H/cem_h_run.sh and cem_h_drive.sh (PERSEUS-ME v2.3 track J / v2.4 track H).
# Usage: cem_pairs_state.sh ST SET OUT [LIMIT] [STRICT] [YEARS]
#   SET pairs_hold | pairs_all | treated; STRICT 1 = cem_v3 strict untreated donor filter (hook 30, default 1); YEARS = optional
#   comma list of t1 years to run (smoke). Per t1 year: subject list (t1 PLT_CN) and same-location donor exclusion (every
#   remeasured pair at a subject's plot location, PLOT.PREV_PLT_CN chain), one 1-cycle run of baseline_window 1, 20 sims.
# Engine: engines/cem/ncasi_rc_patches/assemble_shared.sh (engine_c F1/F2 base + seed block, module list, F4, hooks 30/31,
# converter v3), assembled once into engines/restricted/cem_engine_shared. Flags: NCASI registered s41 gate flags for non-ME
# states (no Maine economic, county harvest or v4 productivity flags; climate off: --no_harvest_all one cycle as track J);
# Maine keeps its PERSEUS-ME extra flags. Treated set adds hook 32 (treated donors for C5) when CEM_H_DIR is given.
set -uo pipefail; umask 077
ST=$1; SET=$2; OUT=$3; LIMIT=${4:-0}; STRICT=${5:-1}; YEARS=${6:-}
export LC_ALL=C.UTF-8 LC_COLLATE=C.UTF-8 R_LIBS=$HOME/R/ncasi_s31_lib
CODE=$HOME/perseus_me/code; ENGR=$HOME/perseus_me/ncasi_rc/engines/restricted; SHARED=$ENGR/cem_engine_shared
IN=$HOME/perseus_me/multistate/inputs/$ST/restricted; FIA=$HOME/restricted/ncasi-modeleval/fia_data_rc/$ST
CFG=$HOME/perseus_me/me_state/v23_E/restricted/inputs/config; RDSDIR=$HOME/perseus_me/tracks/trkC_gain/restricted/cem/rds
[ "$ST" = ME ] && RDSDIR=$HOME/perseus_me/me_state/v23_E/restricted/inputs/fia
case "$OUT" in *restricted*) ;; *) echo "ABORT OUT outside a restricted folder"; exit 1;; esac
[ -s $RDSDIR/fia_db_$ST.rds ] || { echo "ABORT no CEM FIA database fia_db_$ST.rds"; exit 3; }
[ -s $SHARED/MANIFEST.md5 ] || bash $CODE/engines/cem/ncasi_rc_patches/assemble_shared.sh $SHARED
mkdir -p $OUT/lists; chmod 700 $OUT $OUT/lists
python3 - "$IN" "$SET" "$FIA/${ST}_PLOT.csv" "$OUT/lists" "$LIMIT" <<'PY' || exit 1
import sys, os, pandas as pd
IN, SET, PLOT, O, LIM = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], int(sys.argv[5])
sys.path.insert(0, os.path.expanduser("~/perseus_me/code/engines/fvs")); import pairs as PG
P = pd.read_csv(f"{IN}/treated/treated_pairs.csv" if SET == "treated" else f"{IN}/pairs_std/pairs_built.csv", dtype={"pair_id": str, "t1_id": str})
if SET == "pairs_hold": P = P[[PG.half_of(p) == "hold" for p in P.pair_id]]
P = P.sort_values("pair_id");  P = P.head(LIM) if LIM else P
P["t1"] = P.t1_id.str.strip()
pl = pd.read_csv(PLOT, usecols=["CN", "PREV_PLT_CN", "STATECD", "UNITCD", "COUNTYCD", "PLOT"], dtype={"CN": str, "PREV_PLT_CN": str})
pl["CN"] = pl.CN.str.strip(); pl["PREV_PLT_CN"] = pl.PREV_PLT_CN.astype(str).str.replace(r"\.0+$", "", regex=True).str.strip()
pl["loc"] = pl[["STATECD", "UNITCD", "COUNTYCD", "PLOT"]].astype(str).agg("_".join, axis=1); loc = dict(zip(pl.CN, pl["loc"]))
rem = pl[pl.PREV_PLT_CN.ne("nan") & pl.PREV_PLT_CN.ne("")][["PREV_PLT_CN", "CN", "loc"]]
P["cycle"] = [PG.cycle_of(x) for x in P.REMPER]
for y, g in P.groupby("INVYR1"):
    pd.DataFrame(dict(PLT_CN=g.t1.unique())).to_csv(f"{O}/subjects_{y}.csv", index=False)
    L = {loc.get(t) for t in g.t1}; ex = rem[rem["loc"].isin(L)]
    pd.DataFrame(dict(T1_PLT_CN=ex.PREV_PLT_CN, T2_PLT_CN=ex.CN)).to_csv(f"{O}/exclude_{y}.csv", index=False)
cd = pd.read_csv(PLOT.replace("_PLOT.csv", "_COND.csv"), usecols=["PLT_CN", "COND_STATUS_CD"], dtype={"PLT_CN": str}); cd["PLT_CN"] = cd.PLT_CN.str.strip()
P["n_fcond"] = P.t1.map(cd[cd.COND_STATUS_CD == 1].groupby("PLT_CN").size()).fillna(0).astype(int)
P[["pair_id", "t1", "INVYR1", "cycle", "n_fcond"]].to_csv(f"{O}/pairs_in.csv", index=False)
print(f"pairs {len(P)} t1 years {P.INVYR1.nunique()}")
PY
EXTRA=""; [ "$ST" = ME ] && EXTRA="--use_decoupled_climate --use_maine_econ --use_county_harvest --use_v4_prod_mult --v4_prod_mult_strength 1.0"
GEN="--untreated_donors --include_remeasured --use_brms_sdimax --use_disturbance --use_potter_vcc --save_per_plot --skip_supply --skip_report_figs --no_econ --use_owner_stratification --use_owner_balanced --no_harvest_all"
YL=$(ls $OUT/lists/subjects_*.csv | sed 's/.*subjects_\([0-9]*\)\.csv/\1/' | sort -n); [ -n "$YEARS" ] && YL=$(echo $YEARS | tr ',' ' ')
t0=$(date +%s)
for y in $YL; do RUN=$OUT/run_$y; rm -rf $RUN; mkdir -p $RUN; cd $RUN
  cp -r $SHARED/R $SHARED/run_projection.R .
  for f in 13_two_round_figure.R 14_final_figure.R 15_cross_validation.R; do echo "# stubbed (figure script, no projection code)" > R/$f; done
  ln -s $CFG config; mkdir -p data
  A="^  data_list <- prepare_fia_data(CONFIG)$"
  H32=""; [ "$SET" = treated ] && H32=$SHARED/hooks/32_treated_c5_donors.R
  for h in $H32 $CODE/engines/cem/patches/subject_plots.R $SHARED/hooks/31_donor_exclusion.R $( [ "$STRICT" = 1 ] && echo $SHARED/hooks/30_untreated_strict_donor_filter.R ); do
    [ "$(grep -c "$A" run_projection.R)" = 1 ] || { echo "ABORT anchor"; exit 2; }; sed -i "/$A/r $h" run_projection.R; done
  export CEM_SUBJECT_PLOTS=$OUT/lists/subjects_$y.csv CEM_EXCLUDE_PAIRS=$OUT/lists/exclude_$y.csv CEM_RNG_SEED=20261029 FIA_DATA_DIR=$RDSDIR
  [ "$STRICT" = 1 ] && export CEM_V3_STRICT_DONORS=1 || unset CEM_V3_STRICT_DONORS
  if [ "$SET" = treated ]; then REG=NE; case $ST in MI|WI|MN) REG=LS;; GA) REG=SN;; IN) REG=CS;; esac
    export CEM_H_TREATED=1 CEM_H_MIN_TREATED=20 CEM_H_DIR=$HOME/perseus_me/tracks/trkC_gain/restricted/cem/pool/$REG
    for f in cut_t1_trees.csv treated_donors.csv subject_rsc.csv; do [ -s $CEM_H_DIR/$f ] || { echo "ABORT $CEM_H_DIR/$f"; exit 3; }; done
  else unset CEM_H_TREATED CEM_H_DIR; fi
  md5sum R/*.R run_projection.R > ENGINE_MD5.txt
  while [ "$(cut -d' ' -f1 /proc/loadavg | cut -d. -f1)" -ge ${CEM_LOADMAX:-7} ]; do sleep 60; done
  echo "CEM $ST pairs t1 $y set=$SET strict=$STRICT seed=20261029 shared=$(md5sum < $SHARED/MANIFEST.md5 | cut -c1-8) $(date)" > RUN_SPEC.txt
  Rscript run_projection.R --state $ST --n_sims 20 --cycles 1 --baseline_year $y --baseline_window 1 --cores ${CEM_CORES:-2} --scenario_set bau --tag pair_$y $GEN $EXTRA 2>&1 \
    | grep -v -E "LAT|LON|PLT_CN" > run.log
  echo "rc=${PIPESTATUS[0]} $(date)" >> RUN_SPEC.txt; chmod -R go-rwx $RUN; echo "$ST $y $(tail -1 RUN_SPEC.txt)"
done
echo "{\"state\":\"$ST\",\"set\":\"$SET\",\"years\":\"$(echo $YL)\",\"wall_s\":$(( $(date +%s)-t0 )),\"scoring\":\"none (setup runner)\"}" > $OUT/run_meta.json; cat $OUT/run_meta.json
