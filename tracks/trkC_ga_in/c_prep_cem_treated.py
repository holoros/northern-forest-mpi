#!/usr/bin/env python3
"""MS round 3 port (PREREG A4.7; region pools SN = GA, CS = IN; inputs from multistate/inputs). Original: NCASI report card track T: CEM hook 32 inputs (treated donors for C5) per region pool.
Generalized from tools/v24_H/prep_h.py cem(): treated_donors.csv (T1,T2,res_share,rsc_class,loc), cut_t1_trees.csv,
subject_rsc.csv. Pool = the six-state treated pairs of the region (NE: ME,NH,VT; LS: MI,WI,MN); tercile cutpoints of
res_share fixed once from the pool (Maine rule A9 applied to the region pool). Cut t1 trees = t1 live trees (STATUSCD 1,
DIA>=1, DRYBIO_AG and TPA_UNADJ present) not in the residual list; checked against n_cut per pair (hard exit on mismatch).
Usage: prep_cem_treated.py REGION OUTDIR"""
import sys, os, json
from pathlib import Path
import numpy as np, pandas as pd
REG = {"SN": ["GA"], "CS": ["IN"]}[sys.argv[1]]   # MS A4.7: study-state treated pairs of the variant region inside the engine donor list
O = Path(sys.argv[2]); O.mkdir(parents=True, exist_ok=True); os.chmod(O, 0o700)   # trkC_gain fix: statement was inside the comment
RC = Path.home() / "perseus_me/multistate/inputs"; FIA = Path.home() / "restricted/ncasi-modeleval/fia_data_rc"
cs = lambda s: s.astype(str).str.strip()
Ps, cuts, res = [], [], {}
for st in REG:
    P = pd.read_csv(RC / st / "restricted/treated/treated_pairs.csv", dtype={"pair_id": str, "t1_id": str}); P["t1"] = cs(P.t1_id); P["t2"] = cs(P.pair_id)
    T = pd.read_csv(RC / st / "restricted/treated/treated_trees.csv", dtype={"pair_id": str, "tree_cn": str}); keep = set(cs(T.tree_cn))
    t1ids = set(P.t1); out = []
    for ch in pd.read_csv(FIA / st / f"{st}_TREE.csv", usecols=["CN", "PLT_CN", "STATUSCD", "DIA", "TPA_UNADJ", "DRYBIO_AG"], dtype={"CN": str, "PLT_CN": str}, chunksize=2_000_000):
        ch["PLT_CN"] = cs(ch.PLT_CN); out.append(ch[ch.PLT_CN.isin(t1ids)])
    tr = pd.concat(out); tr["CN"] = cs(tr.CN)
    t1 = tr[(tr.STATUSCD == 1) & (tr.DIA >= 1) & tr.DRYBIO_AG.notna() & tr.TPA_UNADJ.notna()].copy(); t1["cut"] = ~t1.CN.isin(keep)
    nc = t1.groupby("PLT_CN").cut.sum(); P["n_cut_rebuilt"] = P.t1.map(nc).fillna(0).astype(int)
    bad = int((P.n_cut_rebuilt != P.n_cut.astype(int)).sum())
    res[st] = dict(pairs=len(P), cut_trees=int(t1.cut.sum()), n_cut_mismatch_pairs=bad)
    if bad: sys.exit(f"HARD EXIT {st}: cut tree count mismatch on {bad} pairs")
    pl = pd.read_csv(FIA / st / f"{st}_PLOT.csv", usecols=["CN", "STATECD", "COUNTYCD", "PLOT"], dtype={"CN": str}); pl["CN"] = cs(pl.CN)
    loc = dict(zip(pl.CN, pl.STATECD.astype(str) + "_" + pl.COUNTYCD.astype(str) + "_" + pl.PLOT.astype(str)))
    P["loc"] = P.t1.map(loc); Ps.append(P); cuts.append(t1.loc[t1.cut, ["CN"]])
P = pd.concat(Ps, ignore_index=True); q = np.quantile(P.res_share.values, [1 / 3, 2 / 3])
P["rsc_class"] = np.where(P.res_share < q[0], 1, np.where(P.res_share < q[1], 2, 3))
P[["t1", "t2", "res_share", "rsc_class", "loc"]].rename(columns={"t1": "T1_PLT_CN", "t2": "T2_PLT_CN"}).to_csv(O / "treated_donors.csv", index=False)
pd.concat(cuts).to_csv(O / "cut_t1_trees.csv", index=False)
P[["t1", "rsc_class", "res_share"]].rename(columns={"t1": "PLT_CN"}).to_csv(O / "subject_rsc.csv", index=False)
for f in O.iterdir(): os.chmod(f, 0o600)
res.update(region=sys.argv[1], tercile_cutpoints=[float(x) for x in q], pool_pairs=len(P), class_counts=P.rsc_class.value_counts().sort_index().to_dict())
(O / "prep_summary.json").write_text(json.dumps(res, indent=1, default=int)); print(json.dumps(res, default=int))
