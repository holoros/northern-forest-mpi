"""trkA second code path for the layer headlines: rank-sum AUC (pandas average ranks) on all pairs, explicit pairwise AUC on a seeded
3,000 pair subsample, Brier, E/O and paired dAUC point; independent of h_score.auc_w. usage: trkA_second_layer.py"""
import json
from pathlib import Path
import numpy as np, pandas as pd
T = Path("/home/aaron/perseus_me/tracks/trkA_hcm")
M = pd.read_csv(T / "work/eval/trkA_layer_metrics.csv").set_index("layer"); out = {}
def rank_auc(p, y):
    r = pd.Series(p).rank(method="average").values; n1 = y.sum(); n0 = len(y) - n1
    return float((r[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))
for k in ("A0", "A2", "SCL_cty", "PERM"):
    d = pd.read_parquet(T / f"restricted/layers/{k}.parquet"); p = d.p.values; y = d.y.values.astype(int)
    ra = rank_auc(p, y); i = np.random.default_rng(3000).choice(len(p), 3000, replace=False); pp, yy = p[i], y[i]
    pos, neg = pp[yy == 1], pp[yy == 0]; pw = float(((pos[:, None] > neg[None, :]).sum() + 0.5 * (pos[:, None] == neg[None, :]).sum()) / (len(pos) * len(neg)))
    br = float(np.mean((p - y) ** 2))
    out[k] = dict(rank_auc=ra, main_auc=float(M.loc[k, "AUC"]), diff=float(max(abs(ra - M.loc[k, "AUC"]), abs(br - M.loc[k, "brier"]))), pairwise_auc_sub3000=pw, brier=br, EO=float(p.sum() / y.sum()))
a = pd.read_parquet(T / "restricted/layers/A0.parquet").set_index("CN"); b = pd.read_parquet(T / "restricted/layers/SCL_cty.parquet").set_index("CN").reindex(a.index)
assert (a.y.values == b.y.values).all()
P = pd.read_csv(T / "work/eval/trkA_layer_paired.csv"); r = P[(P.layer == "A0") & (P.reference == "SCL_cty")].iloc[0]
d = rank_auc(a.p.values, a.y.values.astype(int)) - rank_auc(b.p.values, a.y.values.astype(int)); out["dAUC_A0_SCLcty"] = dict(second=d, main=float(r.dAUC), diff=abs(d - float(r.dAUC)))
mx = max(v["diff"] for v in out.values()); out["max_abs_diff"] = mx; out["tolerance"] = "5e-7: main values are read from CSVs written at 6 significant digits"; out["pass"] = bool(mx <= 5e-7)
json.dump(out, open(T / "work/g3_second_layer.json", "w"), indent=1); print(json.dumps(out))
