# Fig 1 scorecard v2: LANDIS-II cut-pair (C5, B25) results folded in for NH, VT, MI, WI, MN from the released
# NCASI states scores (same scorer, rule set a, common basis, c5_net). Maine stays on r4 (no B25 Maine run).
# Aggregates only. Writes data/fig1_scorecard_v2.csv and data/fig1_scorecard_v2_diff.csv. Does not touch v1.
import pandas as pd, numpy as np, os, hashlib, json
H=os.path.expanduser('~/perseus_me'); D=f'{H}/tracks/ncasi_ms/data'
v1=pd.read_csv(f'{D}/fig1_scorecard.csv')
f=lambda d: d[(d.rule_set=='a')&(d.basis=='common basis')&(d.view=='c5_net')].copy()
md5=lambda p: hashlib.md5(open(p,'rb').read()).hexdigest()
out=[v1[v1.state=='ME'].assign(source='multistate r4',P_rank1=np.nan)]; diff=[]; src={}
num=['MPI','lo95','hi95','coverage','P1','P2','P3','P4_C5net','P5']
for s in ['NH','VT','MI','WI','MN']:
    p=f'{H}/me_state/ncasi_{s.lower()}/release/ns_{s.lower()}_scores.csv'; w=f'{H}/ncasi_states/work/score_{s}/ns_scores.csv'
    src[s]={'release':md5(p),'work':md5(w)}
    n=f(pd.read_csv(p)); n['state']=s; assert n.row.is_unique
    n=n.rename(columns={'row':'model','MPI_lo95':'lo95','MPI_hi95':'hi95','P4':'P4_C5net','P4_lo95':'P4_C5net_lo95','P4_hi95':'P4_C5net_hi95'})
    a=v1[v1.state==s]
    m=a.merge(n,on='model',how='outer',suffixes=('_v1','_ns'),indicator=True)
    assert (m._merge=='both').all(), (s, m[m._merge!='both'].model.tolist())
    for _,r in m.iterrows():
        d={c: (r[c+'_ns']-r[c+'_v1']) if pd.notna(r[c+'_ns']) and pd.notna(r[c+'_v1']) else (np.nan if pd.isna(r[c+'_ns']) and pd.isna(r[c+'_v1']) else 'changed_na') for c in num}
        diff.append(dict(state=s,model=r.model,in_headline=r.in_headline_scorecard,**{('d_'+c):d[c] for c in num},tier_v1=r.tier_v1,tier_ns=r.tier_ns,grade_v1=r.grade_v1,grade_ns=r.grade_ns))
    k=a[['model','variant','role','is_ensemble','in_headline_scorecard','scorecard_status']].merge(n,on='model',suffixes=('','_ns'))
    k['source']='ncasi states release'
    out.append(k[[c for c in v1.columns if c in k.columns]+['P_rank1','source']])
v2=pd.concat(out,ignore_index=True)[list(v1.columns)+['P_rank1','source']]
v2.to_csv(f'{D}/fig1_scorecard_v2.csv',index=False)
dd=pd.DataFrame(diff); dd.to_csv(f'{D}/fig1_scorecard_v2_diff.csv',index=False)
json.dump(src,open(f'{D}/fig1_scorecard_v2_sources.json','w'),indent=1)
print(src)
ch=dd[dd.in_headline & dd[['d_'+c for c in num]].apply(lambda r: any((isinstance(x,str)) or (pd.notna(x) and abs(x)>1e-6) for x in r),axis=1)]
print('headline rows changed:',len(ch),'of',int(dd.in_headline.sum())); print(ch[['state','model','d_MPI','d_coverage','d_P4_C5net','tier_v1','tier_ns','grade_v1','grade_ns']].to_string(index=False))
h=v2[v2.in_headline_scorecard]
print(h[h.model.isin(['landis','E_state','fvs','cbm_gated','cem','cbm_af'])][['state','model','MPI','lo95','hi95','grade','tier','coverage','P4_C5net','P_rank1']].round(1).to_string(index=False))
