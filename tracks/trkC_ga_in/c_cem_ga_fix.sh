#!/bin/bash
# 2026-11-15 fix: GA CEM extraction failed (40_convert_cem_v3.R state-code table lacks GA). Track-local converter copy with GA = 13
# (FIA STATECD), one line changed and asserted; shared code untouched. Then GA extract, slot link, whole-plot GA rescoring (primary and
# shuffle) into the same dirs after _PREV_20261115 backups, and release regeneration (TABLES.md backed up first). nice 10, one worker.
set -uo pipefail; umask 077; source ~/perseus_me/run/env.sh >/dev/null 2>&1
T=~/perseus_me/tracks/trkC_gain; R=$T/restricted/cem; G=$R/grm; S=$T/restricted/slots/cem; L=$T/logs; st_(){ echo "- $(date -Iseconds) cemfix: $*" >> $T/release/STATUS.md; }
SRC=~/perseus_me/code/engines/cem/ncasi_rc_patches/40_convert_cem_v3.R; DST=$T/tools/40_convert_cem_v3_GA.R
sed 's/IL = 17, OH = 39, IN = 18)/IL = 17, OH = 39, IN = 18, GA = 13)/' $SRC > $DST
[ "$(diff $SRC $DST | grep -c '^>')" = 1 ] || { st_ "ABORT converter patch did not apply exactly"; exit 2; }
O=$R/GA/treated
( cd $O && R_LIBS=$HOME/R/ncasi_s31_lib nice -n 10 Rscript ~/perseus_me/code/tools/v23_E/cem_pair_extract_v3.R --repo ~/perseus_me/code --converter ../tracks/trkC_gain/tools/40_convert_cem_v3_GA.R --runs . --pairs lists/pairs_in.csv --state GA --fia_dir ~/restricted/ncasi-modeleval/fia_data_rc/GA --work $G --donors GA,FL,SC,NC,TN,AL --glob 'GA_*_pair_%d' --arm cem_shared --out $O/extract 2>&1 | grep -v -E 'PLT_CN|LAT|LON' ) > $L/cem_extract_treated_GA_fix.log
[ -s $O/extract/pair_rows.csv ] || { st_ "GA extract FAILED (logs/cem_extract_treated_GA_fix.log)"; exit 3; }
touch $R/done_extract_treated_GA; st_ "GA extract done ($(tail -n +2 $O/extract/pair_rows.csv | wc -l) pair rows)"
mkdir -p $S/GA/restricted; [ -e $S/GA/restricted/treated ] || ln -s $O $S/GA/restricted/treated
SS=$T/restricted/slots; RS=$T/restricted/score
for d in score_wp_GA shufflewp_GA; do [ -d $RS/$d ] && mv $RS/$d $RS/${d}_PREV_20261115; done
nice -n 10 env HOME=$T/restricted/home python $T/tools/c_score_x2.py GA $RS/score_wp_GA $SS > $L/wp_score_GA_cemfix.log 2>&1; a=$?
nice -n 10 env HOME=$T/restricted/home python $T/tools/c_score_x2.py GA $RS/shufflewp_GA $SS --shuffle 20261114 > $L/wp_shuffle_GA_cemfix.log 2>&1; b=$?
st_ "GA rescoring with CEM rc=$a shuffle rc=$b (XX1 expected to flag only the landis row, see work/xx1_wp_GA_adjudication.json)"
cp $T/release/TABLES.md $T/release/TABLES_PREV_20261115.md
VPY=~/perseus_me/venv/bin/python; cd $T && $VPY tools/c_release.py > $L/cemfix_release.log 2>&1; $VPY tools/c_release_md.py >> $L/cemfix_release.log 2>&1; nice -n 10 $VPY tools/c_claims.py >> $L/cemfix_release.log 2>&1
st_ "release tables regenerated with GA CEM (TABLES_PREV_20261115.md kept)"; echo CEMFIX_DONE $(date -Iseconds)
