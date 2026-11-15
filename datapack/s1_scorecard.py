import pandas as pd, numpy as np, sys
from common import *
ms=f'{HOME}/perseus_me/tracks/multistate/tracks/multistate/multistate_scorecard.csv'
sc=pd.read_csv(ms)
out=[];inputs=[ms]
for s in STATES:
    p=f'{HOME}/perseus_me/multistate/work/r4_scores/score_{s}/ms_scores.csv'; inputs.append(p)
    r=pd.read_csv(p)
    r=r[(r.rule_set=='a')&(r.basis=='common basis')&(r.view=='c5_net')].copy()
    assert (r.state==s).all() and r.row.is_unique
    k=sc[sc.state==s].set_index('model')
    r['in_headline_scorecard']=r.row.isin(k.index)
    r['scorecard_status']=r.row.map(k.status)
    r['variant']=k.variant.iloc[0]
    # headline consistency
    for _,x in r[r.in_headline_scorecard].iterrows():
        y=k.loc[x.row]
        if y.status=='scored':
            assert abs(x.MPI-y.MPI)<1e-6 and abs(x.MPI_lo95-y.MPI_lo95)<1e-6 and x.grade==y.grade, (s,x.row)
    out.append(r)
d=pd.concat(out)
d['is_ensemble']=d.row.str.startswith('E_')
cols=['state','variant','row','role','is_ensemble','in_headline_scorecard','scorecard_status','MPI','MPI_lo95','MPI_hi95','grade','tier','coverage','pillars_evaluated']
for p in ['P1','P2','P3','P4','P5']: cols+= [p,f'{p}_lo95',f'{p}_hi95']
cols+=['rule_set','basis','view','notes']
d=d[cols].rename(columns={'row':'model','MPI_lo95':'lo95','MPI_hi95':'hi95','P4':'P4_C5net','P4_lo95':'P4_C5net_lo95','P4_hi95':'P4_C5net_hi95'})
d.to_csv(f'{DATA}/fig1_scorecard.csv',index=False)
record('fig1_scorecard.csv',inputs,__file__,'Rows: every row of r4_scores ms_scores.csv with rule_set a, basis common basis, view c5_net. in_headline_scorecard marks rows in multistate_scorecard.csv; MPI, CI and grade asserted equal to the scorecard for scored rows. Pillar points 0-100; P4 = C5 cut-pair net view.')
print(len(d)); print(d.groupby('state').size())
status(f'step1 fig1_scorecard.csv done ({len(d)} rows)')
