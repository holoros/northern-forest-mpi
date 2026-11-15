#!/usr/bin/env python3
"""trkA_hcm (PREREG md5 35055a79af9083edec3f5e3223552dd7) section 4: SC-L refitted on the county-blocked folds with the track SC code
(sc_logit.Design, choose, fit_logit, fit_lin via run_split, unchanged); only the fold column is replaced by the county fold map.
usage: trkA_scl_cty.py [--shuffle N] [--plot-folds]   writes restricted/layers/SCL_cty[_shufN].parquet (CN, p, s, fold, y)"""
import argparse, json, os, sys
from pathlib import Path
import numpy as np, pandas as pd
sys.path.insert(0, "/home/aaron/perseus_me/tracks/v26_SC/tools/v26_SC")
import sc_logit as S  # noqa: E402
T = Path("/home/aaron/perseus_me/tracks/trkA_hcm")
CFOLD = {3: 0, 21: 1, 31: 1, 27: 1, 23: 1, 25: 2, 9: 2, 15: 2, 19: 3, 7: 3, 11: 3, 13: 3, 29: 4, 17: 4, 5: 4, 1: 4}
ap = argparse.ArgumentParser(); ap.add_argument("--shuffle", type=int, default=None); ap.add_argument("--plot-folds", action="store_true"); a = ap.parse_args()
meta = json.load(open(S.O / "work/feature_list.json")); X = S.load(a.shuffle)
if not a.plot_folds:
    X["fold"] = X.COUNTYCD.astype(int).map(CFOLD).values; assert X.fold.notna().all()
miss = [c for c in meta["features"] if c not in meta["categorical"] and X[c].isna().any()]
X["yfit"] = X.cut.values.astype(float); folds = sorted(X.fold.unique()); out = []
for k in folds:
    tr = X.fold.values != k; te = ~tr
    d = S.run_split(X, tr, te, meta, miss, [f for f in folds if f != k], "cf"); d["fold"] = k; out.append(d); print(k, d.lam_p.iloc[0], d.lam_s.iloc[0], flush=True)
D = pd.concat(out, ignore_index=True).merge(X[["CN", "cut"]], on="CN"); D["y"] = D.cut.astype(int)
tag = ("SCL_plot" if a.plot_folds else "SCL_cty") + ("" if a.shuffle is None else f"_shuf{a.shuffle}")
o = T / "restricted/layers"; o.mkdir(exist_ok=True); D[["CN", "p", "s", "fold", "y", "lam_p", "lam_s"]].to_parquet(o / f"{tag}.parquet", index=False); os.chmod(o / f"{tag}.parquet", 0o600)
print(tag, len(D), float(D.p.mean()))
