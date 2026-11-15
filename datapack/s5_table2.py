"""Table 2 models. Text cells transcribed from repo docs (source column names the file for each row); states from fig1_scorecard.csv
(rows with a finite MPI on the multistate basis). Blank = not stated in the docs read."""
import pandas as pd
from common import *
P = f'{HOME}/perseus_me'
S_ES = f'{P}/ncasi_rc/engines/ENGINE_SETUP.md'; S_PR = f'{P}/ncasi_rc/docs_PREREG_R1.md'; S_DS = f'{P}/tracks/multistate/docs/DESIGN.md'
S_MS = f'{P}/tracks/multistate/tracks/multistate/PREREG.md'; S_RL = f'{P}/tracks/multistate/tracks/multistate/RELEASE_MS.md'
S_CORE = f'{P}/tracks/multistate/tools/ncasi_rc/score/rc_core.py'; S_RUN = f'{P}/tracks/multistate/tools/ncasi_rc/score/rc_run.py'
S_OSM = f'{P}/tracks/multistate/docs/PERSEUS_ME_ENGINES_ADDENDUM.md'
f1 = pd.read_csv(f'{DATA}/fig1_scorecard.csv')
ran = f1[f1.MPI.notna()].groupby('model').state.apply(lambda s: ';'.join([x for x in STATES if x in set(s)]))
ens_mem = f1[f1.is_ensemble & f1.MPI.notna()].groupby('model').apply(lambda g: '; '.join(f'{r.state}: {str(r.notes).replace("members ","")}' for r in g.itertuples()), include_groups=False)
T = [
 dict(model='fvs', full_name='Forest Vegetation Simulator, default state variant (NE in ME, NH, VT; LS in MI, WI, MN)', type='individual-tree empirical (tree-list engine)',
      spatial_unit='FIA plot tree list', growth_mortality='variant individual-tree equations; pair runner reports survivor growth and mortality per pair (LS large-tree DG: Deo and Froese 2013)',
      calibration='default: variant equations, no calibration keywords; fitting data not recorded', key_reference='advanced-forestry-systems/fvs-modern (firebreather builds FVSne ba9d0fa4, FVSls 14d60c05)',
      headline=True, source=f'{S_ES}; {S_PR} item 7; {S_MS} section 1'),
 dict(model='fvs_recal', full_name='FVS default with cross-fitted recalibration', type='individual-tree empirical, post-processed',
      spatial_unit='FIA plot tree list', growth_mortality='as fvs; one gross growth scalar per fold applied in post-processing, net recomputed as gross minus mortality, no feedback into the engine',
      calibration='cross-fitted (two folds by plot, scalar from the other fold)', key_reference='', headline=False, source=f'{S_PR} item 5 and v2 items 1 and 2'),
 dict(model='fvs_ncasi_recal', full_name='FVS with the NCASI recalibration arm (NE: SDIMAX x1.3, FIXMORT 0.0100; LS: SDIMAX x1.1, FIXMORT 0.0100)', type='individual-tree empirical',
      spatial_unit='FIA plot tree list', growth_mortality='as fvs with SDIMAX and FIXMORT keywords', calibration='grid fitted on the fit half of the NCASI s41 pairs; in sample (scored on C2 and C3 only)',
      key_reference='', headline=False, source=f'{S_ES}; {P}/tracks/multistate/tools/ncasi_states/ns_score.py score_view (ncasi_arm keeps C2, C3)'),
 dict(model='osm', full_name='OSM 2.26.1, Acadian variant', type='individual-tree empirical (tree-list engine)', spatial_unit='FIA plot tree list',
      growth_mortality='Acadian variant equations, 5 yr cycles', calibration='as shipped, no local calibration; Acadian calibration base ME, NB, NS, PEI, QC', key_reference='',
      headline=True, source=f'{S_ES}; {S_OSM}; {S_PR} item 7'),
 dict(model='osm_bgi3500', full_name='OSM 2.26.1 Acadian with BGI 3500 for every plot (NH, VT sensitivity)', type='individual-tree empirical', spatial_unit='FIA plot tree list',
      growth_mortality='as osm', calibration='as shipped; outside the site input domain; labelled sensitivity, in no ensemble', key_reference='', headline=False, source=f'{S_PR} item 7 and v2 item 11'),
 dict(model='acd141', full_name='AcadianGY 14.1', type='individual-tree empirical (tree-list engine)', spatial_unit='FIA plot tree list',
      growth_mortality='', calibration='region ME coded, csi 14 m; increments used FIA plots among other data (overlap with pairs likely, unverified)',
      key_reference='Kuehne et al. 2022 (fitting data); holoros/fvs-acd engine files', headline=True, source=f'{S_ES}; {S_PR} item 7'),
 dict(model='cem', full_name='CEM, shared patch build (coarsened exact matching plot imputation)', type='plot imputation (donor matching)', spatial_unit='FIA condition/plot matched to donor FIA remeasured pairs',
      growth_mortality='change imputed from matched donor FIA remeasured pairs (strict donors, same-location pairs excluded); treated runs use untreated donors (hook 32 not inserted)',
      calibration='not recalibrated; donors are FIA remeasured pairs by construction', key_reference='holoros/fia-cem-maine', headline=True, source=f'{S_DS}; {S_ES}; {S_RL} (Maine section item 3); {S_PR} item 5'),
 dict(model='landis', full_name='LANDIS-II v8, Biomass Succession with registered no-climate tables', type='landscape process (cohort biomass succession)',
      spatial_unit='LANDIS site per pair (pair runs); raster cells for the statewide anchored run', growth_mortality='gross = ANPP to year 5, net = (B5 - B0)/5, mortality = gross - net',
      calibration='SppEcoregionData fitted from FIA plots and to FIA net growth per state (in sample for statewide totals)', key_reference='holoros/landis2HPC',
      headline=True, source=f'{S_ES}; {S_DS}; {S_PR} v4 item 3; multistate basis: C5 not read and no P5 reader ({S_RL})'),
 dict(model='cbm_gated', full_name='CBM gated (libcbm 2.10.2, FVS_FIAW_NEW yield curves)', type='carbon budget model driven by yield curves', spatial_unit='area by stratum (statewide); curve per pair',
      growth_mortality='libcbm with growth from FVS-derived curve tables; C5 by curve rule (headline ageq on the multistate basis, B1)', calibration='curve tables from the state SIT and AIDB; FVS curves run on FIA plots (state 2x2 build)',
      key_reference='', headline=True, source=f'{S_ES}; {S_DS}; {S_PR} item 3; {P}/tracks/multistate/tools/ncasi_states/ns_score.py B1'),
 dict(model='cbm_af', full_name='CBM American Forests method', type='carbon budget model driven by yield curves', spatial_unit='area by stratum (statewide); curve per pair',
      growth_mortality='Gompertz merchantable volume on age by forest type group and productivity class, Pilli et al. (2013) modification after harvest, Boudewyn conversion to biomass',
      calibration='Gompertz fits to FIA conditions; state presets (ME, MI, WI, MN only)', key_reference='Papa et al. 2023, Front. For. Glob. Change 6:1259010',
      headline=True, source=f'{S_DS}; {S_ES}'),
 dict(model='fia_curves', full_name='FIA hybrid biomass yield curves (ycx)', type='empirical yield curve', spatial_unit='area by stratum (statewide); curve per pair',
      growth_mortality='net growth from curve increments only (C2 net only; C5 net only)', calibration='fitted to FIA plots (MIN_PLOTS 30); in sample: C1 excluded, no P5',
      key_reference='', headline=True, source=f'{S_ES}; {S_CORE} (C2_QUANT, ELIG P5); {S_RL}'),
 dict(model='E_state', full_name='Equal-weight ensemble of the state members (headline)', type='ensemble', spatial_unit='per scoring unit (pair, chain, state-year)',
      growth_mortality='arithmetic mean of member predictions per unit; on this basis (B1, ens_c5_curves=ageq) curve members enter its C5 under the age-reset rule (R1 v2 item 8 had excluded them)', calibration='inherits members', key_reference='', headline=True, source=f'{S_PR} item 7 and v2 items 7 to 9; {P}/tracks/multistate/tools/ncasi_states/ns_score.py B1'),
 dict(model='E_recal', full_name='Ensemble with recalibrated members', type='ensemble', spatial_unit='per scoring unit', growth_mortality='mean of members using cross-fitted recalibrated rows',
      calibration='cross-fitted members', key_reference='', headline=False, source=f'{S_RUN} docstring; {P}/tracks/multistate/tools/ncasi_states/ns_score.py score_view'),
 dict(model='E_plot', full_name='Ensemble of plot-based engines', type='ensemble', spatial_unit='per scoring unit', growth_mortality='mean of plot-based members (fvs, osm, acd141, cem where present)',
      calibration='inherits members', key_reference='', headline=False, source=f'{S_CORE} PLOT_BASED; {S_PR} v2 item 7'),
 dict(model='E_state_median', full_name='Median ensemble of the state members', type='ensemble', spatial_unit='per scoring unit', growth_mortality='median of member predictions per unit',
      calibration='inherits members', key_reference='', headline=False, source=f'{S_RUN} docstring'),
]
d = pd.DataFrame(T)
d['states_scored_multistate'] = d.model.map(ran).fillna('')
d['ensemble_members_by_state'] = d.model.map(ens_mem).fillna('')
d = d.rename(columns={'headline': 'in_headline_scorecard', 'calibration': 'calibration_status_used_here', 'growth_mortality': 'growth_and_mortality_representation'})
d.to_csv(f'{DATA}/table2_models.csv', index=False)
record('table2_models.csv', [S_ES, S_PR, S_DS, S_MS, S_RL, S_CORE, S_RUN, S_OSM, f'{DATA}/fig1_scorecard.csv'], __file__,
       'Text transcribed from the listed docs (per-row source column); states_scored_multistate and ensemble members derived from fig1_scorecard.csv. Blank cells were not stated in the docs read.')
print(d[['model', 'states_scored_multistate']].to_string(index=False))
status(f'step5 table2_models.csv done ({len(d)} rows)')
