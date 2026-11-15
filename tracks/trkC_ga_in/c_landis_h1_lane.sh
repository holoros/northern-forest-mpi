#!/bin/bash
# trkC_gain PREREG section 4: LANDIS-II GA hindcast chains (P2) on firebreather (compute rule of 14 Nov: no Cardinal). Prep = Track R
# h1 frame + registered GA run block with Duration 30 (c_landis_h1_prep.sh); run = ~/opt/landis_v8/landis_run.sh (the Track R
# firebreather runner, registered image and patched DLLs); extraction = engines/landis/sites/extract_sites.R; post = Track R
# landis_h1_state.py post (multistate GA frames). Then a scorer slot: symlinks to every lane LANDIS GA file + this h1 folder.
set -uo pipefail; umask 077; source ~/perseus_me/run/env.sh >/dev/null 2>&1
export R_LD_LIBRARY_PATH=$HOME/opt/geos/lib:/usr/lib/R/lib:/usr/lib/x86_64-linux-gnu
T=~/perseus_me/tracks/trkC_gain; O=$T/restricted/landis/GA/h1; L=$T/logs
step() { n=$1; shift; [ -f $T/work/done_$n ] && { echo "skip $n"; return 0; }; echo "start $n $(date -Iseconds)"; nice -n 10 "$@" > $L/landis_$n.log 2>&1; rc=$?; echo "end $n rc=$rc $(date -Iseconds)"; [ $rc = 0 ] && touch $T/work/done_$n; return $rc; }
step prep_GA bash $T/tools/c_landis_h1_prep.sh GA $O || exit 1
step run_GA bash -c "cd $O/run && bash $HOME/opt/landis_v8/landis_run.sh $O/run > landis_stdout.txt 2>&1; rc=\$?; echo rc=\$rc >> RUN_SPEC.txt; exit \$rc" || exit 1
step extract_GA Rscript ~/perseus_me/code/engines/landis/sites/extract_sites.R $O/run $O/build $O/values.csv || exit 1
step post_GA python $T/tools/c_landis_h1.py post GA $O/values.csv $O/sites/keys.csv $O/h1_values || exit 1
S=$T/restricted/slots/landis/GA/restricted; mkdir -p $S; for f in ~/perseus_me/multistate/runs/landis/GA/restricted/*; do [ -e $S/$(basename $f) ] || ln -s $f $S/$(basename $f); done
[ -e $S/h1 ] || ln -s $O/h1_values $S/h1; ls -la $S | head; echo LANDIS_H1_DONE $(date -Iseconds)
