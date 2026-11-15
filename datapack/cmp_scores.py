import pandas as pd, os
os.chdir(os.path.expanduser('~/perseus_me'))
pd.set_option('display.width',250)
f=lambda d: d[(d.rule_set=='a')&(d.basis=='common basis')&(d.view=='c5_net')]
sc=pd.read_csv('tracks/multistate/tracks/multistate/multistate_scorecard.csv')
for s in ['ME','NH','VT','MI','WI','MN']:
    r=f(pd.read_csv(f'multistate/work/r4_scores/score_{s}/ms_scores.csv'))
    print(s, len(r), r.row.tolist())
    k=sc[(sc.state==s)&(sc.status=='scored')].merge(r,left_on='model',right_on='row',how='left',suffixes=('_sc','_r4'))
    k['d']=(k.MPI_sc-k.MPI_r4).round(4)
    print(k[['model','MPI_sc','MPI_r4','d']].to_string(index=False))
    if s!='ME':
        n=f(pd.read_csv(f'me_state/ncasi_{s.lower()}/release/ns_{s.lower()}_scores.csv'))
        m=r.merge(n,on='row',suffixes=('_r4','_ns'))
        m['d']=(m.MPI_r4-m.MPI_ns).round(3)
        print(m[['row','MPI_r4','MPI_ns','d']].to_string(index=False))
