#!/bin/bash
# trkC_gain: overlay HOME so scorer work dirs and added run slots stay inside this track. Symlinks only; nothing in the live tree changes.
set -euo pipefail
T=/home/aaron/perseus_me/tracks/trkC_gain; O=$T/restricted/home; R=/home/aaron
mkdir -p $O; chmod 700 $O
lnall() { src=$1; dst=$2; shift 2; mkdir -p $dst; for e in $(ls -A $src); do skip=0; for x in "$@"; do [ "$e" = "$x" ] && skip=1; done; [ $skip = 1 ] && continue; [ -e $dst/$e ] || [ -L $dst/$e ] || ln -s $src/$e $dst/$e; done; }
lnall $R $O perseus_me
lnall $R/perseus_me $O/perseus_me multistate
lnall $R/perseus_me/multistate $O/perseus_me/multistate restricted runs
lnall $R/perseus_me/multistate/restricted $O/perseus_me/multistate/restricted score_work
mkdir -p $O/perseus_me/multistate/restricted/score_work; chmod 700 $O/perseus_me/multistate/restricted/score_work
lnall $R/perseus_me/multistate/runs $O/perseus_me/multistate/runs
echo overlay ok; ls -la $O/perseus_me/multistate/runs
