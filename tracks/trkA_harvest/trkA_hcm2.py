#!/usr/bin/env python3
"""trkA_hcm (PREREG md5 35055a79af9083edec3f5e3223552dd7) sections 2, 3, 5: start-of-cycle Harvest Choice Model (HCM2) refit on the
track H Maine pair frame, cross-fitted on county-blocked folds. Writes restricted/layers/<variant>[_shufN].parquet (CN, p, s, fold).
usage: trkA_hcm2.py VARIANT [--shuffle N]   VARIANT in A0 A1 A2 A3 A4 A5 A6 PERM
No coordinate column is read (PLOT and fia_plots_with_owner.csv read with usecols that exclude LAT/LON)."""
import argparse, json, os
from pathlib import Path
import numpy as np, pandas as pd

R = Path("/home/aaron/perseus_me/me_state"); T = Path("/home/aaron/perseus_me/tracks/trkA_hcm"); RS = T / "restricted"
STUMP = Path("/home/aaron/jobs/hcs_stumpage_refit_20261001/stumpage_panel_full.csv")
OWN = Path("/home/aaron/repos-mirror/fia-cem-maine/config/fia_plots_with_owner.csv")
TONS_PER_CORD = 2.5; FT3_PER_CORD = 79.0; THR = 0.70; RIDGE = 1e-4
CFOLD = {3: 0, 21: 1, 31: 1, 27: 1, 23: 1, 25: 2, 9: 2, 15: 2, 19: 3, 7: 3, 11: 3, 13: 3, 29: 4, 17: 4, 5: 4, 1: 4}   # PREREG section 2
CNTY = {1: "ANDROSCOGGIN", 3: "AROOSTOOK", 5: "CUMBERLAND", 7: "FRANKLIN", 9: "HANCOCK", 11: "KENNEBEC", 13: "KNOX", 15: "LINCOLN", 17: "OXFORD",
        19: "PENOBSCOT", 21: "PISCATAQUIS", 23: "SAGADAHOC", 25: "SOMERSET", 27: "WALDO", 29: "WASHINGTON", 31: "YORK"}
GRP = {"KENNEBEC, KNOX, LINCOLN, AND WALDO COUNTIES": [11, 13, 15, 27], "ANDROSCOGGIN, CUMBERLAND, SAGADAHOC, AND YORK COUNTIES": [1, 5, 23, 31]}


# ---- copied unchanged from tracks/imf3/tools/imf3_harvest.py (price_tables, price_at) ----
def price_tables():
    P = pd.read_csv(STUMP, low_memory=False); P = P[(P.state_postal == "ME") & (P.estimate == "Average")].copy()
    P["yr"] = P.nominal_year_num.astype(int); P["nm"] = P.sub_region_name.str.upper().str.strip()
    saw = P[(P.product_norm == "sawlog") & (P.species == "Spruce-Fir") & (P.unit == "$/MBF")].copy()
    pul = P[(P.product_norm == "pulp") & (P.species == "Hardwood") & P.unit.isin(["$/cord", "$/ton"])].copy()
    pul["value"] = np.where(pul.unit == "$/ton", pul.value * TONS_PER_CORD, pul.value)
    out = {}
    for k, d in (("saw", saw), ("pulp", pul)):
        rows = []
        for _, r in d.iterrows():
            if r.nm in GRP:
                rows += [(c, r.yr, r.value, 1) for c in GRP[r.nm]]
            else:
                cc = [c for c, n in CNTY.items() if n == r.nm]
                if cc:
                    rows.append((cc[0], r.yr, r.value, 0))
        X = pd.DataFrame(rows, columns=["cty", "yr", "v", "grp"])
        X = X.sort_values(["cty", "yr", "grp"]).groupby(["cty", "yr"]).v.agg(lambda s: float(np.mean(s.values))).rename("v").reset_index()
        cy = X.set_index(["cty", "yr"]).v
        st = X.groupby("yr").v.median()
        out[k] = (cy, st)
    return out


def price_at(tab, cty, yr):
    cy, st = tab; y = int(min(max(yr, 2000), 2023)); fb = 0
    v = cy.get((int(cty), y), np.nan) if np.isfinite(cty) else np.nan
    if not np.isfinite(v):
        ys = st.index.values; v = float(st.iloc[int(np.argmin(np.abs(ys - y)))]); fb = 1
    return v, fb
# ---- end copy ----


def hcb_class(cns):
    """IMF3-D7 join by plot identity (UNITCD, COUNTYCD, PLOT) of t2 CN; no LAT/LON read"""
    O = pd.read_csv(OWN, usecols=["UNITCD", "COUNTYCD", "PLOT", "INVYR", "hcb_class"], dtype={"hcb_class": str})
    O = O.sort_values(["UNITCD", "COUNTYCD", "PLOT", "INVYR"]).drop_duplicates(["UNITCD", "COUNTYCD", "PLOT"], keep="last")
    om = {(int(a), int(b), int(c)): v for a, b, c, v in zip(O.UNITCD, O.COUNTYCD, O.PLOT, O.hcb_class)}
    PL = pd.read_csv("/home/aaron/fia_data/ME_PLOT.csv", usecols=["CN", "UNITCD", "COUNTYCD", "PLOT"], dtype={"CN": str}, low_memory=False)
    PL["CN"] = PL.CN.str.strip(); PL = PL.set_index("CN"); k = PL.reindex(pd.Series(cns).astype(str).str.strip().values)
    return np.array([om.get((int(a), int(b), int(c)), "missing") if np.isfinite(a) and np.isfinite(b) and np.isfinite(c) else "missing"
                     for a, b, c in zip(k.UNITCD.values.astype(float), k.COUNTYCD.values.astype(float), k.PLOT.values.astype(float))])


def load(shuffle=None):
    P = pd.read_parquet(R / "v26_H/restricted/frame/pairs_h.parquet")
    H = pd.read_parquet(R / "v26_HCM/restricted/hcm_features.parquet")
    S = pd.read_parquet(R / "v26_SC/restricted/features.parquet")[["CN", "ELEV", "RDDISTCD", "OWNCD", "OWNGRPCD", "RESERVCD", "volbf", "volcf",
                                                                    "prev_trt10", "yrs_since_cut", "prev_pair_cut"]]
    X = P.merge(H, on="CN", how="left", validate="1:1").merge(S, on="CN", how="left", validate="1:1")
    assert len(X) == len(P) == 13750
    if shuffle is not None:
        X = X.sample(frac=1.0, random_state=shuffle)
    X = X.sort_values("CN").reset_index(drop=True)      # order-safe: every fit sees rows in CN order
    X["cfold"] = X.COUNTYCD.astype(int).map(CFOLD); assert X.cfold.notna().all()
    return X


def age3(a):
    return np.where(np.isnan(a), "amiss", np.where(a < 40, "a0", np.where(a < 80, "a40", "a80")))


def group_means(val, tr, keys, fallback, n_min=20):   # track HCM hcm_predict.py rule
    d = pd.DataFrame(dict(v=val, k=keys, f=fallback)); t = d[tr]
    g = t.groupby("k").v.agg(["mean", "size"]); gf = t.groupby("f").v.agg(["mean", "size"]); glob = t.v.mean()
    m1 = d.k.map(g["mean"].where(g["size"] >= n_min)); m2 = d.f.map(gf["mean"].where(gf["size"] >= n_min))
    return m1.fillna(m2).fillna(glob).values


def mnl_fit(X, y, ridge=RIDGE, iters=100):   # track HCM hcm_predict.py mnl_fit (3 alternatives)
    n, p = X.shape; B = np.zeros((2, p)); pen = np.full(p, ridge); pen[0] = 0.0
    Y = np.stack([(y == 1), (y == 2)], 1).astype(float)
    for it in range(iters):
        V = X @ B.T; m = np.maximum(0, V.max(1)); E = np.exp(V - m[:, None]); den = np.exp(-m) + E.sum(1); Pr = E / den[:, None]
        g = np.concatenate([X.T @ (Y[:, j] - Pr[:, j]) - pen * B[j] for j in range(2)])
        H = np.zeros((2 * p, 2 * p))
        for j in range(2):
            for k in range(2):
                w = Pr[:, j] * ((j == k) - Pr[:, k])
                H[j * p:(j + 1) * p, k * p:(k + 1) * p] = -(X.T @ (X * w[:, None]))
            H[j * p:(j + 1) * p, j * p:(j + 1) * p] -= np.diag(pen)
        step = np.linalg.solve(H - 1e-10 * np.eye(2 * p), -g); B = B + step.reshape(2, p)
        if np.max(np.abs(step)) < 1e-9:
            break
    return B, it + 1


def mnl_prob(V1, V2):
    m = np.maximum(0, np.maximum(V1, V2)); e0, e1, e2 = np.exp(-m), np.exp(V1 - m), np.exp(V2 - m); s = e0 + e1 + e2
    return e0 / s, e1 / s, e2 / s


def onehot(v, levels):
    return np.column_stack([(v == l).astype(float) for l in levels]) if len(levels) else np.zeros((len(v), 0))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("variant"); ap.add_argument("--shuffle", type=int, default=None); a = ap.parse_args()
    V = a.variant; X = load(a.shuffle); n = len(X)
    folds = (X.fold.values if V in ("A5", "A6") else X.cfold.values).astype(int)
    cut = X.cut.values.astype(bool); rem = X.REMPER.values.astype(float); yr1 = X.INVYR1.values.astype(float)
    share_ba = X.share_ba.values.astype(float); h2 = X.h2_eligible.values.astype(bool); ncut = X.n_cut_t1.values.astype(float)
    sh_saw = np.nan_to_num(X.share_saw.values.astype(float)); sh_pld = np.nan_to_num(X.share_pld.values.astype(float))
    if V == "PERM":   # label permutation (seed 7): outcome block moved jointly across pairs
        pi = np.random.default_rng(7).permutation(n)
        cut, share_ba, h2, ncut, sh_saw, sh_pld = cut[pi], share_ba[pi], h2[pi], ncut[pi], sh_saw[pi], sh_pld[pi]
    tabs = price_tables(); ps, fs, pp, fp = [], [], [], []
    for c, y in zip(X.COUNTYCD.astype(float).values, yr1):
        v, f = price_at(tabs["saw"], c, y); ps.append(v); fs.append(f); v, f = price_at(tabs["pulp"], c, y); pp.append(v); fp.append(f)
    ps, pp = np.array(ps), np.array(pp); pfb = np.maximum(fs, fp).astype(float)
    val = X.volbf.fillna(0).clip(lower=0).values * ps / 1000.0 + X.volcf.fillna(0).clip(lower=0).values / FT3_PER_CORD * pp
    own = X.OWNCD.astype("Int64").astype(str).values; ogr = X.OWNGRPCD.astype("Int64").astype(str).values
    cnt = pd.Series(own).value_counts(); own_p = np.where(pd.Series(own).map(cnt).values >= 30, "o" + own, "g" + ogr)
    lev_own = sorted(set(own_p))[1:]                     # first level is the reference
    lev_res = sorted(set(X.RESERVCD.astype(str)))[1:]
    a3 = age3(X.stdage1.values.astype(float)); ftg = X.ftgrp1.astype(int).astype(str).values; key = (pd.Series(ftg) + "|" + pd.Series(a3)).values
    present_s = X.saw_t1.values > 0; present_p = X.pld_t1.values > 0

    def ychoice(sh, present):
        y = np.zeros(n, int); c = cut & present
        y[c & (sh >= THR)] = 2; y[c & (sh > 0) & (sh < THR)] = 1; y[c & (ncut == 0)] = 1
        return y
    ys = ychoice(sh_saw, present_s); yp = ychoice(sh_pld, present_p)
    yt = np.zeros(n, int); yt[cut] = 1; yt[cut & (np.nan_to_num(share_ba) >= THR)] = 2      # A1 total BA class
    per = np.column_stack([((yr1 >= 2007) & (yr1 <= 2011)), yr1 >= 2012]).astype(float)
    pub = np.isin(X.OWNGRPCD.values, [10, 20, 30]).astype(float)
    ysc = X.yrs_since_cut.values.astype(float); ppc = X.prev_pair_cut.values.astype(float)
    cdum = onehot(X.COUNTYCD.astype(int).values, sorted(CFOLD)[1:])
    hcb = hcb_class(X.CN.values) if V == "A3" else None
    p = np.full(n, np.nan); s = np.full(n, np.nan); info = {}
    for k in range(5):
        tr = folds != k; te = ~tr
        Gs = group_means(((X.saw_t2 - X.saw_t1) / rem).values, tr & ~cut, key, ftg)
        Gp = group_means(((X.pld_t2 - X.pld_t1) / rem).values, tr & ~cut, key, ftg)
        e = tr & h2; s_part = np.nanmean(share_ba[e][share_ba[e] < THR]); s_full = np.nanmean(share_ba[e][share_ba[e] >= THR])
        med_ysc = np.nanmedian(ysc[tr]) if np.isfinite(ysc[tr]).any() else 0.0

        def design(stock, G):
            core = [np.log(ps), np.log(pp), pfb, X.lagbio.values, stock, G, G ** 2, X.ELEV.values.astype(float), X.RDDISTCD.values.astype(float),
                    per, np.log(rem)]
            if V == "A4":
                return np.column_stack(core + [pub])
            add = [np.log1p(val), X.prev_trt10.values.astype(float), np.where(np.isfinite(ysc), ysc, med_ysc), (~np.isfinite(ysc)).astype(float),
                   np.nan_to_num(ppc), (~np.isfinite(ppc)).astype(float), onehot(own_p, lev_own), onehot(X.RESERVCD.astype(str).values, lev_res)]
            if V == "A3":
                add.append(onehot(hcb, sorted(set(hcb))[1:]))
            if V == "A6":
                add.append(cdum)
            return np.column_stack(core + add)

        def fitpred(D, y):
            mu = D[tr].mean(0); sd = D[tr].std(0); keep = sd > 0
            Z = np.column_stack([np.ones(n), (D[:, keep] - mu[keep]) / sd[keep]])
            B, it = mnl_fit(Z[tr], y[tr]); Vv = Z[te] @ B.T
            return np.column_stack(mnl_prob(Vv[:, 0], Vv[:, 1])), it

        def ce(P):
            den = P[:, 1] + P[:, 2]
            return np.where(den > 0, (P[:, 1] * s_part + P[:, 2] * s_full) / np.where(den > 0, den, 1), (s_part + s_full) / 2)
        if V == "A1":
            Pt, it = fitpred(design(X.lagbio.values, (Gs + Gp)), yt); p[te] = 1 - Pt[:, 0]; s[te] = ce(Pt); info[f"fold{k}_iters"] = it
        else:
            ss = X.saw_t1.values.copy(); sp = X.pld_t1.values.copy()
            Ds, Dp = design(ss, Gs), design(sp, Gp)
            if V == "A2":   # leak diagnostic: observed t2 (post-harvest) product biomass, as the paper's BioTot / BiopulpLD
                Ds = np.column_stack([Ds, X.saw_t2.values]); Dp = np.column_stack([Dp, X.pld_t2.values])
            Ps, i1 = fitpred(Ds, ys); Pp, i2 = fitpred(Dp, yp)
            p[te] = np.maximum(1 - Ps[:, 0], 1 - Pp[:, 0]); s[te] = (ce(Ps) + ce(Pp)) / 2; info[f"fold{k}_iters"] = [i1, i2]
        info[f"fold{k}_s_part_full"] = [float(s_part), float(s_full)]
    assert np.isfinite(p).all()
    out = RS / "layers"; out.mkdir(exist_ok=True); os.chmod(out, 0o700)
    tag = V + ("" if a.shuffle is None else f"_shuf{a.shuffle}")
    D = pd.DataFrame(dict(CN=X.CN.values, p=p, s=s, fold=folds, y=cut.astype(int)))
    D.to_parquet(out / f"{tag}.parquet", index=False); os.chmod(out / f"{tag}.parquet", 0o600)
    meta = dict(variant=V, shuffle=a.shuffle, n=n, folds="plot-location (track H)" if V in ("A5", "A6") else "county (PREREG section 2)",
                mean_p=float(p.mean()), mean_s=float(np.nanmean(s)), observed_cut=float(cut.mean()), price_fallback_share=float(pfb.mean()),
                own_levels=["ref"] + lev_own, choice_counts_saw=np.bincount(ys, minlength=3).tolist(), choice_counts_pulp=np.bincount(yp, minlength=3).tolist(),
                choice_counts_total=np.bincount(yt, minlength=3).tolist(), **info)
    json.dump(meta, open(T / "work" / f"hcm2_{tag}_meta.json", "w"), indent=1, default=float); print(json.dumps(meta, default=float)[:800])


if __name__ == "__main__":
    main()
