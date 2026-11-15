#!/bin/bash
# trkC_gain: c_sit_2x2.sh (the registered 2x2 driver's FVS x FIAW x NEW cell) on firebreather, against a copy of the cbm_states files
# pulled from Cardinal (restricted/cbm_pull/cbm_root, md5 manifest cbm_root.md5). python3 = venv_cbm through engines/pycbm. The only
# code edit: the hard-coded Cardinal sys.path line of build_libcbm_state_bundle.py points at the copy's shared/tools.
# usage: c_sit_2x2_fb.sh ST SIT_YIELD_CSV OUTDIR WORKROOT
set -uo pipefail
ST=$1; SY_IN=$2; OUT=$3; ROOT=$4; T=~/perseus_me/tracks/trkC_gain; R0=$T/restricted/cbm_pull/cbm_root; PY=$HOME/perseus_me/ncasi_rc/engines/pycbm
[ -d $ROOT ] || { cp -a $R0 $ROOT; sed -i "s#_sys.path.insert(0, '/users/PUOM0008/crsfaaron/cbm_states/shared/tools')#_sys.path.insert(0, '$ROOT/shared/tools')#" $ROOT/cross_state/libcbm/build_libcbm_state_bundle.py; }
grep -c "$ROOT/shared/tools" $ROOT/cross_state/libcbm/build_libcbm_state_bundle.py
mkdir -p $OUT; export CBM_STATES_ROOT=$ROOT CBM_MAINE_ROOT=$T/restricted/cbm_pull/cbm_maine FIA_REF_FOREST_TYPE=$T/restricted/cbm_pull/REF_FOREST_TYPE.csv; cd $ROOT   # GCBM template db pulled from Cardinal (md5 fb37f046)
PROC=$ROOT/states/$ST/data/processed; PIX=$(ls -t $PROC/pixel_attributes_fia_*.csv | head -1)
AREA=$(awk -F, 'NR>1{s+=$2}END{printf "%.0f", s}' "$PROC/fia_expns_grp_areas_${ST}.csv"); echo "### $ST: area ${AREA} ha"
cp -p "$SY_IN" "$PROC/sit_yield.csv" || exit 3; rm -f $PROC/gcbm_input_${ST}.db; rm -rf $ROOT/cross_state/libcbm/$ST $OUT/adapted_FVS_FIAW_NEW   # no stale build can be reused
$PY -u shared/tools/build_gcbm_input_state.py "$ST" > $OUT/gcbm_${ST}_FVS_NEW.log 2>&1; rc=$?; echo "gcbm rc=$rc"; [ $rc = 0 ] && [ -s $PROC/gcbm_input_${ST}.db ] || { echo "GCBM BUILD FAILED"; exit 4; }
grep -E "specific|state_mean|FT-class|state-pool|falling back" $OUT/gcbm_${ST}_FVS_NEW.log | sed 's/^/     /'
BDIR=$OUT/bundle_FVS_FIAW_NEW; mkdir -p $BDIR
$PY -u cross_state/libcbm/build_libcbm_state_bundle.py "$ST" --total-area-ha "$AREA" --age-distribution auto --pixel-attrs-csv "$PIX" \
  --fia-expns-csv "$PROC/fia_expns_grp_areas_${ST}.csv" > $OUT/bundle_FVS_FIAW_NEW.log 2>&1; rc=$?; echo "bundle rc=$rc"; [ $rc = 0 ] || exit 5
for f in sit_classifiers.csv sit_yield.csv sit_events.csv sit_inventory.csv; do [ -f $ROOT/cross_state/libcbm/$ST/$f ] && cp -p $ROOT/cross_state/libcbm/$ST/$f $BDIR/$f; done
ADIR=$OUT/adapted_FVS_FIAW_NEW
$PY -u cross_state/libcbm/adapt_state_to_libcbm.py "$ST" "$ADIR" --bundle-dir "$BDIR" --total-area-ha "$AREA" > $OUT/adapt_FVS_FIAW_NEW.log 2>&1; rc=$?; echo "adapt rc=$rc"; [ $rc = 0 ] && [ -s $ADIR/growth_and_yield.csv ] || exit 6
cp -p $PROC/cbm_defaults_${ST}.db $OUT/; md5sum $ADIR/* > $OUT/MANIFEST.md5; echo "SIT_DONE $ST"
