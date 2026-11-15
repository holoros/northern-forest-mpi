#!/usr/bin/env python3
"""trkA_hcm (PREREG md5 35055a79af9083edec3f5e3223552dd7) section 6: layer discrimination and calibration with the track H
plot-location cluster bootstrap stream (as tracks/v26_SC sc_eval.py), paired AUC differences, county-cluster sensitivity, promotion rule.
usage: trkA_eval.py OUTDIR [--shuffle N]   (aggregates only in OUTDIR)"""
import argparse, json, sys
from pathlib import Path
import numpy as np, pandas as pd
sys.path.insert(0, "/home/aaron/perseus_me/tracks/v26_H/tools/v26_H")
import h_score as HS  # noqa: E402
R = Path("/home/aaron/perseus_me/me_state"); T = Path("/home/aaron/perseus_me/tracks/trkA_hcm"); LY = T / "restricted/layers"
ap = argparse.ArgumentParser(); ap.add_argument("out"); ap.add_argument("--shuffle", type=int, default=None); a = ap.parse_args()
out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
P = pd.read_parquet(R / "v26_H/restricted/frame/pairs_h.parquet")
if a.shuffle is not None:
    P = P.sample(frac=1.0, random_state=a.shuffle)
P = P.sort_values("CN").reset_index(drop=True); y = P.cut.values.astype(float); n = len(P)
rng = np.random.default_rng(HS.SEED); locs = np.sort(P["loc"].unique()); li = np.searchsorted(locs, P["loc"].values); L = len(locs)
Wl = np.zeros((HS.NB, L), np.float32)
for b in range(HS.NB):
    Wl[b] = np.bincount(rng.integers(0, L, L), minlength=L)
Wp = Wl[:, li]
cty = np.sort(P.COUNTYCD.unique()); ci = np.searchsorted(cty, P.COUNTYCD.values); rc = np.random.default_rng(20261103)
Wc = np.stack([np.bincount(rc.integers(0, len(cty), len(cty)), minlength=len(cty)) for _ in range(1000)])[:, ci].astype(np.float32)


def calib(p, yy, w):
    x = np.log(np.clip(p, 1e-9, 1 - 1e-9) / (1 - np.clip(p, 1e-9, 1 - 1e-9))); X = np.column_stack([np.ones_like(x), x]); b = np.zeros(2)
    for _ in range(100):
        q = 1 / (1 + np.exp(-(X @ b))); W = w * q * (1 - q) + 1e-12
        st = np.linalg.solve(X.T @ (X * W[:, None]), X.T @ (w * (yy - q))); b += st
        if np.abs(st).max() < 1e-10:
            break
    return b


def ci(x):
    x = np.asarray(x, float); x = x[np.isfinite(x)]; return (float(np.percentile(x, 2.5)), float(np.percentile(x, 97.5)))


lay = {}
for f in sorted(LY.glob("*.parquet")):
    if "_shuf" in f.stem:
        continue
    d = pd.read_parquet(f).set_index("CN").reindex(P.CN); lay[f.stem] = (d.p.values.astype(float), d.y.values.astype(float) if f.stem == "PERM" else y)
sl = pd.read_parquet(R / "v26_SC/restricted/preds_sc_logit.parquet"); sl = sl[sl.split == "cf"].set_index("CN").reindex(P.CN)
sf = pd.read_csv(R / "v26_SC/restricted/preds_sc_rf.csv", dtype={"CN": str}); sf = sf[sf.split == "cf"].set_index("CN").reindex(P.CN)
hc = pd.read_parquet(R / "v26_HCM/restricted/hcm_preds.parquet"); hc = hc[hc.driver == "hcm_refit"].set_index("CN").reindex(P.CN)
lay["stored_SCL_plot"] = (sl.p.values.astype(float), y); lay["stored_SCF_plot"] = (sf.p.values.astype(float), y); lay["stored_HCMrefit_plot"] = (hc.p.values.astype(float), y)
rows, AB = [], {}
for k, (p, yy) in lay.items():
    m = np.isfinite(p); pm, ym = p[m], yy[m]
    auc = float(HS.auc_w(pm, ym, np.ones(m.sum()))[0]); ab = HS.auc_w(pm, ym, Wp[:, m]); ac = HS.auc_w(pm, ym, Wc[:, m]); AB[k] = (m, ab, ac)
    b0 = calib(pm, ym, np.ones(m.sum())); cb = np.array([calib(pm, ym, Wp[i, m].astype(float)) for i in range(HS.NB)])
    br = float(np.mean((pm - ym) ** 2)); brb = (Wp[:, m] @ (pm - ym) ** 2) / Wp[:, m].sum(1)
    rows.append(dict(layer=k, n=int(m.sum()), n_cut=int(ym.sum()), AUC=auc, AUC_lo=ci(ab)[0], AUC_hi=ci(ab)[1], AUC_cty_lo=ci(ac)[0], AUC_cty_hi=ci(ac)[1],
                     brier=br, brier_lo=ci(brb)[0], brier_hi=ci(brb)[1], cal_slope=float(b0[1]), cal_slope_lo=ci(cb[:, 1])[0], cal_slope_hi=ci(cb[:, 1])[1],
                     cal_int=float(b0[0]), E_over_O=float(pm.sum() / ym.sum()), mean_p=float(pm.mean())))
M = pd.DataFrame(rows)
pairs = [("A0", "SCL_cty"), ("A1", "SCL_cty"), ("A3", "SCL_cty"), ("A4", "SCL_cty"), ("A0", "A4"), ("A3", "A0"), ("A1", "A0"), ("A2", "A0"),
         ("SCL_cty", "stored_SCL_plot"), ("A5", "stored_SCL_plot"), ("A6", "stored_SCL_plot"), ("A5", "stored_HCMrefit_plot"), ("A0", "stored_SCL_plot")]
dif = []
for a1, a2 in pairs:
    if a1 not in lay or a2 not in lay:
        continue
    p1, p2 = lay[a1][0], lay[a2][0]; m = np.isfinite(p1) & np.isfinite(p2)
    d0 = float(HS.auc_w(p1[m], y[m], np.ones(m.sum()))[0] - HS.auc_w(p2[m], y[m], np.ones(m.sum()))[0])
    db = HS.auc_w(p1[m], y[m], Wp[:, m]) - HS.auc_w(p2[m], y[m], Wp[:, m]); dc = HS.auc_w(p1[m], y[m], Wc[:, m]) - HS.auc_w(p2[m], y[m], Wc[:, m])
    dbr = (Wp[:, m] @ ((p1[m] - y[m]) ** 2 - (p2[m] - y[m]) ** 2)) / Wp[:, m].sum(1)
    dif.append(dict(layer=a1, reference=a2, n=int(m.sum()), dAUC=d0, dAUC_lo=ci(db)[0], dAUC_hi=ci(db)[1], P_better=float(np.mean(db > 0)),
                    dAUC_cty_lo=ci(dc)[0], dAUC_cty_hi=ci(dc)[1], dBrier=float(np.mean((p1[m] - y[m]) ** 2 - (p2[m] - y[m]) ** 2)), dBrier_lo=ci(dbr)[0], dBrier_hi=ci(dbr)[1]))
Dd = pd.DataFrame(dif)
rel = []
for k in ("A0", "SCL_cty"):
    if k in lay:
        p = lay[k][0]; q = pd.qcut(pd.Series(p).rank(method="first"), 10, labels=False)
        for b_, g in pd.DataFrame(dict(q=q, p=p, y=y)).groupby("q"):
            rel.append(dict(layer=k, bin=int(b_) + 1, n=len(g), mean_p=g.p.mean(), observed=g.y.mean()))
r0 = Dd[(Dd.layer == "A0") & (Dd.reference == "SCL_cty")].iloc[0]; s0 = M.set_index("layer").loc["A0"]
promo = dict(rule="paired dAUC (A0 - SCL_cty) 95% interval above 0 AND A0 calibration slope in [0.8, 1.2] (PREREG section 6)",
             dAUC=float(r0.dAUC), dAUC_lo=float(r0.dAUC_lo), dAUC_hi=float(r0.dAUC_hi), cal_slope=float(s0.cal_slope),
             auc_gain_excludes_0=bool(r0.dAUC_lo > 0), slope_ok=bool(0.8 <= s0.cal_slope <= 1.2))
promo["promoted"] = promo["auc_gain_excludes_0"] and promo["slope_ok"]; promo["winner"] = "HCM2 (A0)" if promo["promoted"] else "SC-L county"
pm = M.set_index("layer").loc["PERM"]; promo["perm_check"] = dict(AUC=float(pm.AUC), pass_=bool(abs(pm.AUC - 0.5) <= 0.03))
M.to_csv(out / "trkA_layer_metrics.csv", index=False, float_format="%.6g"); Dd.to_csv(out / "trkA_layer_paired.csv", index=False, float_format="%.6g")
pd.DataFrame(rel).to_csv(out / "trkA_layer_reliability.csv", index=False, float_format="%.6g"); json.dump(promo, open(out / "trkA_promotion.json", "w"), indent=1)
pd.set_option("display.width", 250); print(M.round(4).to_string()); print(Dd.round(4).to_string()); print(json.dumps(promo, indent=1))
