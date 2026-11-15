#!/bin/bash
# trkC_gain PREREG section 3: CEM GA/IN C5 port, one process at a time, nice 10, CEM_CORES 1. Steps: region treated pools (fixed copy
# of ms_prep_cem_treated.py), GRM condition tables for the GA donor states, treated runs (Track R cem_pairs_state.sh with the
# multistate input and rds paths and CEM_H_DIR = this track's pool; sed substitutions asserted), extraction (cem_pair_extract_v3.R).
set -uo pipefail; umask 077; source ~/perseus_me/run/env.sh >/dev/null 2>&1
T=~/perseus_me/tracks/trkC_gain; R=$T/restricted/cem; WT=~/perseus_me/tracks/multistate; mkdir -p $R; chmod 700 $R
export NCASI_RC_STATES=$WT/config/multistate/states_ms_inputs.yaml OMP_NUM_THREADS=1 CEM_CORES=1 CEM_LOADMAX=40 LC_ALL=C.UTF-8 CEM_EXCLUDE_SUBJECT_LOCATIONS=1
step() { n=$1; shift; [ -f $R/done_$n ] && { echo "skip $n"; return 0; }; echo "start $n $(date -Iseconds)"; nice -n 10 "$@" > $T/logs/cem_$n.log 2>&1; rc=$?; echo "end $n rc=$rc $(date -Iseconds)"; [ $rc = 0 ] && touch $R/done_$n; return $rc; }
step pool_SN python $T/tools/c_prep_cem_treated.py SN $R/pool/SN
step pool_CS python $T/tools/c_prep_cem_treated.py CS $R/pool/CS
G=$R/grm; mkdir -p $G
for s in GA FL SC NC TN AL; do D=~/restricted/ncasi-modeleval/fia_data_rc/$s; [ -d $D ] || D=~/restricted/ncasi-modeleval/fia_data_rc/donors/$s
  step grm_$s Rscript ~/perseus_me/code/tools/ncasi_rc/eng/build_grm_cond.R $s $D $G; done
[ -e $G/grm_cond_IN.rds ] || ln -s ~/perseus_me/ncasi_rc/cem_work/grm/grm_cond_IN.rds $G/grm_cond_IN.rds
P=$T/tools/c_cem_pairs_state.sh
sed -e "s#^IN=\$HOME/perseus_me/ncasi_rc/inputs/\$ST/restricted;#IN=\$HOME/perseus_me/multistate/inputs/\$ST/restricted;#" \
    -e "s#RDSDIR=\$HOME/restricted/ncasi-modeleval/fia_data_rc/cem_rds#RDSDIR=\$HOME/perseus_me/tracks/trkC_gain/restricted/cem/rds#" \
    -e "s#CEM_H_DIR=\$HOME/perseus_me/ncasi_rc/cem_work/treated_pool/\$REG#CEM_H_DIR=\$HOME/perseus_me/tracks/trkC_gain/restricted/cem/pool/\$REG#" \
    -e "s#REG=NE; case \$ST in MI|WI|MN) REG=LS;; esac#REG=NE; case \$ST in MI|WI|MN) REG=LS;; GA) REG=SN;; IN) REG=CS;; esac#" \
    $WT/tools/ncasi_rc/eng/cem_pairs_state.sh > $P
[ "$(diff $WT/tools/ncasi_rc/eng/cem_pairs_state.sh $P | grep -c '^>')" = 4 ] || { echo "ABORT path patch did not apply exactly"; exit 2; }
declare -A DON=( [GA]=GA,FL,SC,NC,TN,AL [IN]=IN )
for st in GA IN; do
  mkdir -p $R/rds; [ -e $R/rds/fia_db_GA.rds ] || ln -s ~/perseus_me/multistate/restricted/cem_rds/fia_db_GA.rds $R/rds/fia_db_GA.rds
  [ -s $R/rds/fia_db_$st.rds ] || step db_$st Rscript $T/tools/c_cem_db.R $st
  O=$R/$st/treated; mkdir -p $O
  step treated_$st bash $P $st treated $O
  step extract_treated_$st bash -c "cd $O && R_LIBS=\$HOME/R/ncasi_s31_lib Rscript ~/perseus_me/code/tools/v23_E/cem_pair_extract_v3.R --repo ~/perseus_me/code --converter engines/cem/ncasi_rc_patches/40_convert_cem_v3.R --runs . --pairs lists/pairs_in.csv --state $st --fia_dir ~/restricted/ncasi-modeleval/fia_data_rc/$st --work $G --donors ${DON[$st]} --glob '${st}_*_pair_%d' --arm cem_shared --out $O/extract 2>&1 | grep -v -E 'PLT_CN|LAT|LON'"
done
echo CEM_LANE_DONE $(date -Iseconds)
