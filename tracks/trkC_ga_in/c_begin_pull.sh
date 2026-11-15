#!/bin/bash
# trkC_gain section 3 (c): TREE_GRM_BEGIN of the GA CEM donor states (not in the multistate donor pull), public DataMart CSV, into
# this track's restricted folder, then GRM condition tables (build_grm_cond.R unchanged) from a folder of symlinks + BEGIN.
set -uo pipefail; umask 077; source ~/perseus_me/run/env.sh >/dev/null 2>&1
B=https://apps.fs.usda.gov/fia/datamart/CSV; T=~/perseus_me/tracks/trkC_gain; R=$T/restricted/cem
for ST in FL SC NC TN AL; do D=$R/donor_tabs/$ST; mkdir -p $D; chmod 700 $R/donor_tabs $D
  for f in $(ls ~/restricted/ncasi-modeleval/fia_data_rc/donors/$ST/${ST}_*.csv); do [ -e $D/$(basename $f) ] || ln -s $f $D/$(basename $f); done
  f=${ST}_TREE_GRM_BEGIN.csv
  if [ ! -s $D/$f ]; then n=0; until nice -n 19 wget -q -S -T 120 -O $D/$f.part $B/$f 2> $D/$f.hdr && mv $D/$f.part $D/$f; do n=$((n+1)); [ $n -ge 3 ] && { echo "FAIL $f"; break; }; sleep 60; done; fi
  [ -s $D/$f ] && { echo "$(date +%F_%T) OK $f $(md5sum < $D/$f | cut -c1-32) $(grep -i Last-Modified $D/$f.hdr | tail -1)"; }
  [ -f $R/done_grm_$ST ] || { nice -n 10 Rscript ~/perseus_me/code/tools/ncasi_rc/eng/build_grm_cond.R $ST $D $R/grm > $T/logs/cem_grm2_$ST.log 2>&1 && touch $R/done_grm_$ST; echo "grm $ST rc=$?"; }
done; md5sum $R/grm/*.rds > $R/grm/MD5SUMS.txt; echo BEGIN_PULL_DONE
