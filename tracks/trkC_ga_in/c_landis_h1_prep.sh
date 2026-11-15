#!/bin/bash
# MS round 3 port (PREREG Amendment 4 A4.6: input path to multistate/inputs, GA case added; nothing else changed). Track R original:
# NCASI report card (task #5 setup): LANDIS-II site runs on FIA pair tree lists for any state. No scoring.
# Generalized from tools/v24_H/landis_h.sh (PERSEUS-ME v2.4 track H). Usage:
#   landis_pairs_state.sh ST SET OUT [SPLIT] [LIMIT]
#     SET   pairs_hold | pairs_all | treated   (inputs ~/perseus_me/ncasi_rc/inputs/ST/restricted)
#     SPLIT d25 (registered builder, default) | obsbio (cohort biomass split by observed tree biomass; switch)
#     LIMIT first N pairs by sorted pair_id (smoke), 0 = all
# Treated set: sites are built from the residual t1 tree list (residual-direct). The registered Track H A6 path (full t1
# community, oldest cohorts removed by species) needs the t1 cut-flag frame of prep_h.py and is not this runner.
# Registered no-climate configuration: promoted table of PARAMS_ACTIVE (MN: registry-row run copy, ProbEstablish x0.1),
# Monthly_AverageAllYears climate, 5 yr, seed 1234, no harvest. Output orientation: year 0 identity (extract_sites.R).
# Writes OUT/values.csv (site_id, year, B_gm2, ANPP_gm2yr), OUT/keys.csv, OUT/pair_rows_landis.csv (model components
# per pair, g m-2 yr-1: gross = ANPP of the step ending year 5, net = (B5 - B0)/5, mortality = gross - net) and run_meta.json.
# MS round 4 copy: PREP ONLY. Identical to ms_landis_pairs.sh up to the LANDIS call; LANDIS runs on Cardinal (registered sif md5 7c255eed,
# patched DLLs md5 equal) and ms_landis_post_r4.sh does the extraction. No other line changed.
set -uo pipefail; umask 077
ST=$1; SET=h1; OUT=$2; SPLIT=d25; LIMIT=0
source $HOME/perseus_me/run/env.sh
export R_LD_LIBRARY_PATH=$HOME/opt/geos/lib:/usr/lib/R/lib:/usr/lib/x86_64-linux-gnu
CODE=$HOME/perseus_me/code; IN=$HOME/perseus_me/multistate/inputs/$ST/restricted; PUL=$HOME/perseus_me/ncasi_rc/engines/pulled/landis
LI=$HOME/landis2_fb/states/$ST/inputs; ECOL3=$HOME/perseus_me/me_state/v24_H/restricted/inputs/ecol3.gpkg
FIA=$(PME_STATE=$ST NCASI_RC_STATES=$CODE/config/ncasi_rc/states_rc_inputs.yaml python -c "from tools.ncasi_rc import stcfg; print(stcfg.get('fia_root'))" 2>/dev/null || true)
[ -n "$FIA" ] && [ "$FIA" != None ] || FIA=$HOME/restricted/ncasi-modeleval/fia_data_rc/$ST
case $ST in
  ME) T=$PUL/ME_SppEcoregionData_fia_20261020.csv; TM=b73361188430fe67dce53bbdf9957df9; C2=$HOME/perseus_me/me_state/landis_full/work/c2/run
      SPF=$C2/species.txt; SDF=$C2/SpeciesData.csv; EF=$C2/ecoregions.txt; CLIM=$C2/climate.csv; LOOK=$HOME/perseus_me/me_state/v24_H/restricted/inputs/species_lookup_ME_s43.R ;;
  NH|VT|MI|WI) T=$PUL/${ST}_SppEcoregionData_fia_20261015.csv
      TM=$(awk -F, -v s=$ST '$1==s && $4=="promoted"{print $3}' $PUL/PARAMS_ACTIVE.csv)
      SPF=$LI/species.txt; SDF=$LI/SpeciesData.csv; EF=$LI/ecoregions_20260929.txt; CLIM=$LI/Daymet_${ST}_l3.csv; LOOK=$LI/species_lookup_${ST}.R ;;
  MN) T=$PUL/MN_SppEcoregionData_fia_pe0.1.csv; TM=fb295178181728b699893d875d468966
      SPF=$LI/species.txt; SDF=$LI/SpeciesData.csv; EF=$LI/ecoregions.txt; CLIM=$LI/PRISM_MN_l3.csv; LOOK=$LI/species_lookup.R ;;
  GA) LI=$HOME/perseus_me/multistate/restricted/pulled_cardinal/landis/GA; T=$LI/SppEcoregionData_fia.csv; TM=31490d0b5360c55dcd78aa0bf8b241ba
      SPF=$LI/species.txt; SDF=$LI/SpeciesData.csv; EF=$LI/ecoregions.txt; CLIM=$LI/PRISM_GA_l3.csv; LOOK=$LI/species_lookup.R ;;  # MS A4.6
  *) echo "unknown state $ST"; exit 2 ;;
esac
[ "$(md5sum < $T | cut -c1-32)" = "$TM" ] || { echo "registered table md5 mismatch for $ST"; exit 3; }
mkdir -p $OUT; chmod 700 $OUT; t0=$(date +%s)
python $HOME/perseus_me/tracks/trkC_gain/tools/c_landis_h1.py frame $ST $OUT || exit 1   # trkC_gain: sites = hindcast chains (Track R h1 frame)
FIT=$(awk -F, 'NR>1{print $2}' $T | sort -un | paste -sd,)
Rscript $CODE/tools/ncasi_rc/eng/build_sites_state.R --sites $OUT/sites/sites.csv --trees $OUT/sites/trees.csv --lookup $LOOK \
  --ecol3 $ECOL3 --plot_csv $FIA/${ST}_PLOT.csv --spacing 4 --cell 270 --fit $FIT --state $ST --split $SPLIT --out $OUT/build || exit 1
B=$OUT/build; RUN=$OUT/run; rm -rf $RUN; mkdir -p $RUN/output/biomass
( cd $RUN
  cp $SPF species.txt; cp $SDF SpeciesData.csv
  cp $B/ecoregions.tif $B/initial-communities.tif $B/initial-communities.csv $B/build_summary.csv .
  ECOS=$(awk -F, 'NR>1{print $1}' "$B/active_cells_by_ecoregion.csv" | tr '\n' ' ')
  awk -F, -v E="$ECOS" 'BEGIN{OFS=","; n=split(E,x," "); for(i=1;i<=n;i++) k[x[i]]=1} NR==1{print; next} ($2 in k){print}' "$T" > SppEcoregionData_fia.csv
  awk -v E="$ECOS" 'BEGIN{n=split(E,x," "); for(i=1;i<=n;i++) k[x[i]]=1} ($1!="yes") || ($2 in k)' "$EF" > ecoregions.txt
  awk -F, -v E="$ECOS" 'BEGIN{n=split(E,x," "); for(i=1;i<=n;i++) k[x[i]]=1} {sub(/\r$/, ""); if (NR==1) for(i=1;i<=NF;i++) keep[i]=(i<=3 || ($i in k))
    o=""; for(i=1;i<=NF;i++) if (keep[i]) o=o (o=="" ? "" : ",") $i; print o}' "$CLIM" > climate.csv
  printf 'LandisData "Climate Config"\nClimateTimeSeries          Monthly_AverageAllYears\nClimateFile                climate.csv\nSpinUpClimateTimeSeries    Monthly_AverageAllYears\nSpinUpClimateFile          climate.csv\nGenerateClimateOutputFiles no\nUsingFireClimate           no\n' > climate-generator.txt
  { echo 'LandisData "Biomass Succession"'; echo 'Timestep 5'; echo 'SeedingAlgorithm WardSeedDispersal'
    echo 'InitialCommunities       initial-communities.csv'; echo 'InitialCommunitiesMap    initial-communities.tif'
    echo 'ClimateConfigFile        climate-generator.txt'; echo 'MinRelativeBiomass'
    printf '        '; for e in $ECOS; do printf '%s ' "$e"; done; echo
    for pct in 20 40 60 80 100; do sh=$((pct/20)); printf '    %d  ' "$sh"; for e in $ECOS; do printf '%s%% ' "$pct"; done; echo; done
    printf 'SufficientLight\n    1   1.0  0.5  0.0  0.0  0.0  0.0\n    2   1.0  1.0  0.5  0.0  0.0  0.0\n    3   1.0  1.0  1.0  0.5  0.0  0.0\n    4   1.0  1.0  1.0  1.0  0.5  0.0\n    5   1.0  1.0  1.0  1.0  1.0  0.5\n'
    echo 'SpeciesDataFile        SpeciesData.csv'; echo 'EcoregionParameters'; echo '>>  AET (mm)'
    for e in $ECOS; do echo "$e   600"; done
    echo 'SpeciesEcoregionDataFile  SppEcoregionData_fia.csv'
    printf 'FireReductionParameters\n    1        0.0        0.5\n    2        0.0        0.75\n    3        0.5        1.0\nHarvestReductionParameters\n    PerseusBA50 0.5 0.15 0.5 0.0\n'
  } > biomass_succession.txt
  SPP=$(awk 'NF>=9 && $1 !~ /^>>/ && $1!="LandisData"{print $1}' species.txt)
  { echo 'LandisData  "Output Biomass"'; echo 'Timestep   5'; echo 'MakeTable  yes'
    first=1; for s in $SPP; do if [ $first = 1 ]; then echo "Species  $s"; first=0; else echo "         $s"; fi; done
    echo 'MapNames   output/biomass/biomass-{species}-{timestep}.tif'; echo 'DeadPools  both'
    echo 'MapNames   output/biomass/biomass-{pool}-{timestep}.tif'; } > output_biomass.txt
  CELL=$(awk -F, 'NR==2{printf "%d", $2+0.5}' build_summary.csv)
  printf "LandisData  Scenario\nDuration  30\nSpecies   ./species.txt\nEcoregions      ./ecoregions.txt\nEcoregionsMap   ./ecoregions.tif\nCellLength  %s\n\"Biomass Succession\"  biomass_succession.txt\n\"Output Biomass\"      output_biomass.txt\nRandomNumberSeed   1234\n" "$CELL" > scenario.txt
  echo "state=$ST set=$SET split=$SPLIT table=$T md5=$TM ecos=$ECOS climate=$CLIM duration=30 seed=1234 harvest=none" > RUN_SPEC.txt
  echo "prep only (trkC_gain h1): LANDIS run on Cardinal" >> RUN_SPEC.txt )
exit 0
