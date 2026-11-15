import pandas as pd
from common import *
src=f'{HOME}/jobs/ncasi_fig4direct_20261027/2026-10-27_init-direct-by-state_DATA.csv'
d=pd.read_csv(src)
bad=[c for c in d.columns if any(k in c.lower() for k in ('lat','lon','plt_cn','cn_','coord'))]
assert not bad, bad
d.to_csv(f'{DATA}/fig3_treemap_direct_all_variants.csv',index=False)
h=d[(d.fia_eval%10000==2401)&(d.variant=='pixel_all')].copy()
assert len(h)==6 and set(h.state)==set(STATES)
h['state']=pd.Categorical(h.state,STATES); h=h.sort_values('state')
keep=['state','fia_eval','variant','fia_tg_c','fia_se','fia_lo95','fia_hi95','fia_area_ha','fia_area_se_ha','fia_mgc_ha','fia_tg_agb','fia_n_plots',
      'treemap_tg_c','treemap_lo','treemap_hi','treemap_area_ha','treemap_mgc_ha','treemap_tg_agb','n_imputed_plots',
      'pct_diff_total','pct_diff_area','pct_diff_density','share_area_logdecomp','share_density_logdecomp','diff_tg','diff_se_indep',
      'diff_outside_95_indep','tm_outside_fia95','fia_outside_tm95','tm_invyr_med_pxw','tm_invyr_min','tm_invyr_max','tm_px_share_invyr_ge2020']
h=h[keep].rename(columns={'fia_tg_c':'fia_agc_TgC','fia_se':'fia_se_TgC','fia_lo95':'fia_lo95_TgC','fia_hi95':'fia_hi95_TgC',
   'fia_mgc_ha':'fia_density_MgC_ha','fia_tg_agb':'fia_agb_Tg','treemap_tg_c':'treemap_agc_TgC','treemap_lo':'treemap_lo95_TgC',
   'treemap_hi':'treemap_hi95_TgC','treemap_mgc_ha':'treemap_density_MgC_ha','treemap_tg_agb':'treemap_agb_Tg','diff_tg':'diff_treemap_minus_fia_TgC',
   'diff_se_indep':'diff_se_indep_TgC','pct_diff_total':'pct_diff_total_pct','pct_diff_area':'pct_diff_area_pct','pct_diff_density':'pct_diff_density_pct'})
h.to_csv(f'{DATA}/fig3_treemap_direct.csv',index=False)
note=('Headline rows: 2024 evaluation (EVALID ST2401), variant pixel_all; these reproduce the Oct 26 draft values for all six states. '
      'Source script jobs/ncasi_fig4direct_20261027/fig4_direct.R lists its input as Cardinal s49/init_direct/init_direct_by_state.csv md5 d974947ba1a8f0c98102ca8c0603c04b (not re-fetched; Cardinal not used). '
      'fia_* are FIA design-based estimates with SE; treemap_lo/hi are TreeMap bootstrap 95% CI as stored. Units: Tg C, Mg C/ha, ha. Area/density split columns as stored.')
record('fig3_treemap_direct.csv',[src],__file__,note)
record('fig3_treemap_direct_all_variants.csv',[src],__file__,'Verbatim copy of all 168 rows (6 states x 4 FIA evaluations x 7 TreeMap variants).')
print(h[['state','fia_agc_TgC','treemap_agc_TgC','treemap_lo95_TgC','treemap_hi95_TgC']].round(1).to_string(index=False))
status('step3b fig3_treemap_direct.csv done (6 rows; all variants 168 rows)')
