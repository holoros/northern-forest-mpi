import pandas as pd, json, glob
from common import *
parts=[]; metas=[]; inputs=[]
for st in STATES:
    f=f'{DATA}/_work/comp_{st}.csv'; m=f'{DATA}/_work/comp_{st}_meta.json'
    parts.append(pd.read_csv(f)); metas.append(json.load(open(m))); inputs += [f,m]
d=pd.concat(parts)
f1=pd.read_csv(f'{DATA}/fig1_scorecard.csv')[['state','model','in_headline_scorecard','is_ensemble','MPI','P1','P4_C5net']]
d=d.merge(f1,on=['state','model'],how='left',validate='many_to_one')
d=d.rename(columns={'ratio':'ratio_model_over_obs','lo95':'lo95_scorer_draws','hi95':'hi95_scorer_draws','ratio_pairboot':'ratio_model_over_obs_pairboot',
                    'lo95_pairboot':'lo95_pairboot','hi95_pairboot':'hi95_pairboot'})
d['used_in_multistate_MPI']=((d.component=='C2_uncut')&d.P1.notna())|((d.component=='C5_cut')&d.P4_C5net.notna()&(d.quantity=='net_growth'))
d.loc[d.path!='scorer_draws','used_in_multistate_MPI']=False
d['basis']=d.path.map({'scorer_draws':'multistate (rule set a, common basis)','landis_b25_post_scoring':'NCASI B25 LANDIS C5 slot written after the r4 multistate scores; not in the multistate MPI'})
d['units_note']='ratio of means of model and observed pair rates (aboveground live biomass, per-area annual rates); foliage factor (1 - phi) applied to CBM rows as in the scorer'
cols=['state','model','is_ensemble','in_headline_scorecard','basis','used_in_multistate_MPI','component','quantity','ratio_model_over_obs','lo95_scorer_draws','hi95_scorer_draws','n_pairs',
      'n_draws_finite','ratio_model_over_obs_pairboot','lo95_pairboot','hi95_pairboot','n_pairs_pairboot','flag','units_note']
d=d[cols]
d.to_csv(f'{DATA}/fig2_components.csv',index=False)
pd.DataFrame(metas).drop(columns=['missing'],errors='ignore').to_csv(f'{DATA}/fig2_components_checks.csv',index=False)
note=('Primary ratio and 95% interval: the frozen scorer code path (ns_score/rc_adapt.load_state bundle, RR.comps_single and RR.ensemble_comps, '
      'rule set a, common basis = statewide, B1 options), statewide e = ln(mean model / mean obs) on units where both exist; interval = percentile 2.5/97.5 of exp(e) over the '
      "scorer's 1,000 shared location-multiplicity draws (seed 20261103). The draws are not stored by the scorer; they were regenerated deterministically by its code. "
      'C2 = uncut common (gate) pairs; used_in_multistate_MPI: C2 rows where the row has P1 (C2 enters P1 only with C1), C5 net rows where the row has P4 (the c5_net view scores net only). C5 = cut (treated) pairs, all treated pairs (treated_subset all), three quantities where the engine reports them (c5_three quantities), fia_curves net only. '
      'LANDIS C5 masked where the r4 landis row has P4 not evaluated (B25 slot written later), reproducing the multistate basis; r4 pillars P1, P2, P4, P5 reproduced for every row (fig2_components_checks.csv). Second path for single engines: explicit pair_id join (validate one_to_one) of observed and model arrays, ratio of means, pair bootstrap 1,000 draws seed 20261114 (columns *_pairboot); '
      'point estimates asserted equal to the scorer path (<=1e-9). Shuffle: every unit array permuted consistently (run_track_r.shuffle_bundle, seed 20261114); scorer-path ratios and intervals identical (max diff in fig2_components_checks.csv).')
S2=[f'{ROOT}/tools/s2_components_v1_VT_NH_ME.py',f'{ROOT}/tools/s2_components.py']
note+=' Per-state rows written by tools/s2_components.py (VT, NH, ME with the v1 snapshot s2_components_v1_VT_NH_ME.py; MI, WI, MN with FAST_P5=1, which skips only the P5 trend-series targets and leaves P5 out of the r4 check).'
record('fig2_components.csv',inputs+S2+[f'{DATA}/fig1_scorecard.csv'],__file__,note)
record('fig2_components_checks.csv',inputs,__file__,'Per-state check values written by s2_components.py (shuffle max diffs, path agreement, unit counts, wall time).')
print(len(d)); print(d.groupby(['state','component']).size().unstack())
status(f'step2 fig2_components.csv combined ({len(d)} rows)')
