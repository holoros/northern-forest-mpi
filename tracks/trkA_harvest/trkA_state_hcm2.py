#!/usr/bin/env python3
"""trkA_hcm (PREREG md5 35055a79af9083edec3f5e3223552dd7) section 8: HCM2 (A0 specification) refit within a state on its own county-blocked
folds (greedy rule of section 2), statewide public stumpage; OOF p and s; AUC / Brier / calibration with plot-location cluster bootstrap.
usage: trkA_state_hcm2.py ST [--shuffle N]"""
import argparse, json, os, sys
from pathlib import Path
import numpy as np, pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parent))
import trkA_hcm2 as H  # noqa: E402  (mnl_fit, mnl_prob, group_means, age3, onehot: the Maine HCM2 code)
sys.path.insert(0, "/home/aaron/perseus_me/tracks/v26_H/tools/v26_H")
import h_score as HS  # noqa: E402
T = Path("/home/aaron/perseus_me/tracks/trkA_hcm"); THR = 0.70
ap = argparse.ArgumentParser(); ap.add_argument("st"); ap.add_argument("--shuffle", type=int, default=None); a = ap.parse_args(); ST = a.st
X = pd.read_parquet(T / f"restricted/states/{ST}/frame.parquet")
if a.shuffle is not None:
    X = X.sample(frac=1.0, random_state=a.shuffle)
X = X.sort_values("CN").reset_index(drop=True); n = len(X)
c = X.COUNTYCD.astype(int).value_counts().sort_index().sort_values(ascending=False, kind="stable"); load = [0] * 5; cf = {}
for kk, v in c.items():
    f = min(range(5), key=lambda i: (load[i], i)); cf[int(kk)] = f; load[f] += int(v)
folds = X.COUNTYCD.astype(int).map(cf).values
P = pd.read_csv(H.STUMP, low_memory=False); P = P[(P.state_postal == ST) & (P.estimate == "Average") & (P.value > 0)]   # trkA-D4: value <= 0 is a missing-value sentinel (F7)
saw = P[(P.product_norm == "sawlog") & P.unit.astype(str).str.startswith("$/MBF")].groupby("nominal_year_num").value.median()
pu = P[(P.product_norm == "pulp") & P.unit.isin(["$/cord", "$/ton"])].copy(); pu["v"] = np.where(pu.unit == "$/ton", pu.value * H.TONS_PER_CORD, pu.value)
pul = pu.groupby("nominal_year_num").v.median()
yr1 = X.INVYR1.values.astype(float); t1y = set(range(int(yr1.min()), int(yr1.max()) + 1))
use_price = len(set(saw.index.astype(int)) & set(pul.index.astype(int)) & t1y) >= 5


def at(sr, y):
    ys = sr.index.values.astype(float); i = int(np.argmin(np.abs(ys - y))); return float(sr.iloc[i]), float(ys[i] != y)


if use_price:
    ps = np.array([at(saw, y)[0] for y in yr1]); pp = np.array([at(pul, y)[0] for y in yr1])
    pfb = np.array([max(at(saw, y)[1], at(pul, y)[1]) for y in yr1])
    val = X.volbf.fillna(0).clip(lower=0).values * ps / 1000.0 + X.volcf.fillna(0).clip(lower=0).values / H.FT3_PER_CORD * pp
cut = X.cut.values.astype(bool); rem = X.REMPER.values.astype(float); share_ba = X.share_ba.values.astype(float); h2 = X.h2_eligible.values.astype(bool)
ncut = X.n_cut_t1.values.astype(float); sh_saw = np.nan_to_num(X.share_saw.values.astype(float)); sh_pld = np.nan_to_num(X.share_pld.values.astype(float))
own = X.OWNCD.astype("Int64").astype(str).values; ogr = X.OWNGRPCD.astype("Int64").astype(str).values
cnt = pd.Series(own).value_counts(); own_p = np.where(pd.Series(own).map(cnt).values >= 30, "o" + own, "g" + ogr); lev_own = sorted(set(own_p))[1:]
lev_res = sorted(set(X.RESERVCD.astype(str)))[1:]
ftg = X.ftgrp1.fillna(-1).astype(int).astype(str).values; key = (pd.Series(ftg) + "|" + pd.Series(H.age3(X.stdage1.values.astype(float)))).values


def ychoice(sh, present):
    y = np.zeros(n, int); cc = cut & present; y[cc & (sh >= THR)] = 2; y[cc & (sh > 0) & (sh < THR)] = 1; y[cc & (ncut == 0)] = 1; return y


ys = ychoice(sh_saw, X.saw_t1.values > 0); yp = ychoice(sh_pld, X.pld_t1.values > 0)
per = np.column_stack([((yr1 >= 2007) & (yr1 <= 2011)), yr1 >= 2012]).astype(float)
ysc = X.yrs_since_cut.values.astype(float); ppc = X.prev_pair_cut.values.astype(float)
elev = X.ELEV.values.astype(float); elev = np.where(np.isfinite(elev), elev, np.nanmedian(elev)); rdd = X.RDDISTCD.values.astype(float); rdd = np.where(np.isfinite(rdd), rdd, np.nanmedian(rdd))
p = np.full(n, np.nan); s = np.full(n, np.nan)
for k in range(5):
    tr = folds != k; te = ~tr
    Gs = H.group_means(((X.saw_t2 - X.saw_t1) / rem).values, tr & ~cut, key, ftg); Gp = H.group_means(((X.pld_t2 - X.pld_t1) / rem).values, tr & ~cut, key, ftg)
    e = tr & h2; s_part = np.nanmean(share_ba[e][share_ba[e] < THR]); s_full = np.nanmean(share_ba[e][share_ba[e] >= THR]); med = np.nanmedian(ysc[tr]) if np.isfinite(ysc[tr]).any() else 0.0

    def design(stock, G):
        core = ([np.log(ps), np.log(pp), pfb] if use_price else []) + [X.lagbio.values, stock, G, G ** 2, elev, rdd, per, np.log(rem)]
        add = ([np.log1p(val)] if use_price else []) + [X.prev_trt10.values.astype(float), np.where(np.isfinite(ysc), ysc, med), (~np.isfinite(ysc)).astype(float),
                                                        np.nan_to_num(ppc), (~np.isfinite(ppc)).astype(float), H.onehot(own_p, lev_own), H.onehot(X.RESERVCD.astype(str).values, lev_res)]
        return np.column_stack(core + add)

    def fitpred(D, y):
        mu = D[tr].mean(0); sd = D[tr].std(0); kp = sd > 0; Z = np.column_stack([np.ones(n), (D[:, kp] - mu[kp]) / sd[kp]])
        B, it = H.mnl_fit(Z[tr], y[tr]); Vv = Z[te] @ B.T; return np.column_stack(H.mnl_prob(Vv[:, 0], Vv[:, 1]))

    def ce(Pm):
        den = Pm[:, 1] + Pm[:, 2]; return np.where(den > 0, (Pm[:, 1] * s_part + Pm[:, 2] * s_full) / np.where(den > 0, den, 1), (s_part + s_full) / 2)
    Ps = fitpred(design(X.saw_t1.values, Gs), ys); Pp = fitpred(design(X.pld_t1.values, Gp), yp)
    p[te] = np.maximum(1 - Ps[:, 0], 1 - Pp[:, 0]); s[te] = (ce(Ps) + ce(Pp)) / 2
assert np.isfinite(p).all()
tag = "" if a.shuffle is None else f"_shuf{a.shuffle}"
D = pd.DataFrame(dict(CN=X.CN.values, p=p, s=s, fold=folds, y=cut.astype(int))); D.to_parquet(T / f"restricted/states/{ST}/hcm2{tag}.parquet", index=False)
os.chmod(T / f"restricted/states/{ST}/hcm2{tag}.parquet", 0o600)
y = cut.astype(float); rng = np.random.default_rng(20261103); locs = np.sort(X["loc"].unique()); li = np.searchsorted(locs, X["loc"].values); L = len(locs)
W = np.stack([np.bincount(rng.integers(0, L, L), minlength=L) for _ in range(1000)])[:, li].astype(float)
auc = float(HS.auc_w(p, y, np.ones(n))[0]); ab = HS.auc_w(p, y, W); br = float(np.mean((p - y) ** 2)); brb = (W @ (p - y) ** 2) / W.sum(1)
xl = np.log(np.clip(p, 1e-9, 1 - 1e-9) / (1 - np.clip(p, 1e-9, 1 - 1e-9)))


def cal(w):
    Xc = np.column_stack([np.ones(n), xl]); b = np.zeros(2)
    for _ in range(100):
        q = 1 / (1 + np.exp(-(Xc @ b))); Wt = w * q * (1 - q) + 1e-12; st = np.linalg.solve(Xc.T @ (Xc * Wt[:, None]), Xc.T @ (w * (y - q))); b += st
        if np.abs(st).max() < 1e-10:
            break
    return b


cs = np.array([cal(W[i])[1] for i in range(1000)]); b0 = cal(np.ones(n))
pc = lambda v: [float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))]  # noqa: E731
res = dict(state=ST, shuffle=a.shuffle, pairs=n, cut=int(cut.sum()), counties=len(cf), fold_pairs=load, price_used=bool(use_price),
           price_years_saw=[int(saw.index.min()), int(saw.index.max())] if len(saw) else None, price_years_pulp=[int(pul.index.min()), int(pul.index.max())] if len(pul) else None,
           price_fallback_share=float(pfb.mean()) if use_price else None, AUC=auc, AUC_ci=pc(ab), brier=br, brier_ci=pc(brb), cal_slope=float(b0[1]), cal_slope_ci=pc(cs),
           E_over_O=float(p.sum() / y.sum()), mean_s=float(np.nanmean(s)))
json.dump(res, open(T / f"work/state_{ST}_hcm2{tag}.json", "w"), indent=1); print(json.dumps(res))
