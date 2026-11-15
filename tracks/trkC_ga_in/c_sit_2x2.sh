#!/bin/bash
# trkC_gain PREREG 2.3: the registered 2x2 driver (cbm_states/shared/tools/run_2x2_curves_inventory_2026-09-22.sh) reduced to its
# FVS x FIAW x NEW cell, run against a scratch copy of cbm_states (CBM_STATES_ROOT) so no shared states/ST file is written.
# The commands of that cell are copied verbatim (build_gcbm_input_state.py, build_libcbm_state_bundle.py --age-distribution auto
# --pixel-attrs-csv --fia-expns-csv, adapt_state_to_libcbm.py); NEW binning = the current files (GA and IN have no _PREV copies).
# usage: c_sit_2x2.sh ST SIT_YIELD_CSV OUTDIR
set -uo pipefail
ST=$1; SY_IN=$2; OUT=$3; SRC=/users/PUOM0008/crsfaaron/cbm_states; X=/fs/scratch/PUOM0008/crsfaaron/trkC_gain; ROOT=$X/cbm_root
mkdir -p $ROOT/states/$ST/data $ROOT/cross_state/libcbm $OUT
[ -d $ROOT/shared ] || cp -rp $SRC/shared $ROOT/shared
cp -p $SRC/cross_state/libcbm/*.py $ROOT/cross_state/libcbm/
rm -rf $ROOT/states/$ST/config $ROOT/states/$ST/data/processed; cp -rp $SRC/states/$ST/config $ROOT/states/$ST/config; mkdir -p $ROOT/states/$ST/data/processed
for f in $SRC/states/$ST/data/processed/*; do b=$(basename $f); if [ -d $f ]; then ln -s $f $ROOT/states/$ST/data/processed/$b; else cp -p $f $ROOT/states/$ST/data/processed/$b; fi; done
for d in raw external; do [ -e $ROOT/states/$ST/data/$d ] || ln -s $SRC/states/$ST/data/$d $ROOT/states/$ST/data/$d; done
for d in $(ls $SRC/states/$ST); do [ -e $ROOT/states/$ST/$d ] || ln -s $SRC/states/$ST/$d $ROOT/states/$ST/$d; done
export CBM_STATES_ROOT=$ROOT; cd $ROOT
PROC=$ROOT/states/$ST/data/processed; PIX=$(ls -t $PROC/pixel_attributes_fia_*.csv | head -1)
AREA=$(awk -F, 'NR>1{s+=$2}END{printf "%.0f", s}' "$PROC/fia_expns_grp_areas_${ST}.csv"); echo "### $ST: corrected FIA forest area ${AREA} ha; pix $(basename $PIX)"
cp -p "$SY_IN" "$PROC/sit_yield.csv"; md5sum "$SY_IN" "$PROC/sit_yield.csv" $PROC/cbm_defaults_${ST}.db $PIX
python3 -u shared/tools/build_gcbm_input_state.py "$ST" > $OUT/gcbm_${ST}_FVS_NEW.log 2>&1; echo "gcbm rc=$?"
grep -E "specific|state_mean|FT-class|state-pool|falling back" $OUT/gcbm_${ST}_FVS_NEW.log | sed 's/^/     /'
BDIR=$OUT/bundle_FVS_FIAW_NEW; mkdir -p $BDIR
python3 -u cross_state/libcbm/build_libcbm_state_bundle.py "$ST" --total-area-ha "$AREA" --age-distribution auto --pixel-attrs-csv "$PIX" \
  --fia-expns-csv "$PROC/fia_expns_grp_areas_${ST}.csv" > $OUT/bundle_FVS_FIAW_NEW.log 2>&1; echo "bundle rc=$?"
grep -o "Loaded age distribution for [0-9]* FT groups" $OUT/bundle_FVS_FIAW_NEW.log
for f in sit_classifiers.csv sit_yield.csv sit_events.csv sit_inventory.csv; do [ -f $ROOT/cross_state/libcbm/$ST/$f ] && cp -p $ROOT/cross_state/libcbm/$ST/$f $BDIR/$f; done
ADIR=$OUT/adapted_FVS_FIAW_NEW
python3 -u cross_state/libcbm/adapt_state_to_libcbm.py "$ST" "$ADIR" --bundle-dir "$BDIR" --total-area-ha "$AREA" > $OUT/adapt_FVS_FIAW_NEW.log 2>&1; echo "adapt rc=$?"
ls $ADIR; cp -p $PROC/cbm_defaults_${ST}.db $OUT/; md5sum $ADIR/* $OUT/cbm_defaults_${ST}.db > $OUT/MANIFEST.md5; echo "SIT_DONE $ST"
