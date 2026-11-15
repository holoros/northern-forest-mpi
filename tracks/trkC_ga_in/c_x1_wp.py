#!/usr/bin/env python3
"""trkC_gain Amendment 1 whole-plot copy (curve root and tables _wp); PREREG 2.6: X1 capped variant of CBM gated for GA/IN (refine_curves2 PREREG section 3 form, section 5 X1 asymptote).
Annual recursion: gross = exp(c0 + c2 z) s^c3 G_pot, mort = gross u^kappa, net = gross - mort, u = B / Bmax, Bmax = p99 of current
C1 plot AGB stock within the gated forest type group (>= 30 plots, else state p99), s = B / P(a + t) (min positive training s when
undefined), G_pot = curve gross x (1 - phi); positivity-preserving step. lm(), step() and recurse() are copied from
tracks/refine_curves2/tools/rc2_fit.py (that module imports Maine paths at import time). theta by Levenberg-Marquardt on stacked
residuals [obs net - pred net, obs mortality - pred mortality] (Mg/ha/yr; a1 observed values the scorer uses) of uncut C2 pairs;
C2 rows are out-of-fold over the multistate pair folds; treated pairs and chains use the production fit. Writes Track R row files
(lb/ac over the cycle) and chain parquet into SLOTS/cbm_gated_x1/ST/restricted, plus work/x1_ST.json (aggregates).
usage: c_x1.py ST SLOTS [--permute SEED]"""
import argparse, json, os, sys, time
from pathlib import Path
import numpy as np, pandas as pd
ap = argparse.ArgumentParser(); ap.add_argument("st"); ap.add_argument("slots"); ap.add_argument("--permute", type=int, default=None); a = ap.parse_args()
st = a.st.upper(); H = Path.home(); T = H / "perseus_me/tracks/trkC_gain"; MS = H / "perseus_me/multistate"; I = MS / "inputs" / st / "restricted"
PUL = MS / "restricted/pulled_cardinal"; fia = H / f"restricted/ncasi-modeleval/fia_data_rc/{st}"
root = T / "restricted/curve_roots_wp" / st; w = root / "complete_j/work"; w.mkdir(parents=True, exist_ok=True)
if not (w / "cbm_curves_gated").exists():
    (w / "cbm_curves_gated").symlink_to(T / f"restricted/cbm_curves_wp/gated_{st}")
os.environ.update(PME_STATE=st, PME_STATE_ROOT=str(root), PME_FIA_ROOT=str(fia))
CODE = H / "perseus_me/code"
for p in ("", "tools/v24_H", "tools/complete", "tools/score", "tools/hindcast_all", "tools/hindcast", "engines/fvs"):
    sys.path.insert(0, str(CODE / p))
import prep_h as PH  # noqa: E402
import pairs as PG  # noqa: E402
import score as S  # noqa: E402
BJ = PH._bj(); BJ.FIA = fia; BJ.MEMB = PUL / f"yield_curves_conus/config/ycx_membership_{st}.csv"; BJ.FITS = PUL / f"yield_curves_conus/ycx_{st}_hybrid_fits.csv"
K = 0.00112085; MINCELL = 30; E = "cbm_gated"; t0 = time.time()
phi = float(pd.read_csv(MS / "work/common/phi_states.csv").set_index("state").phi[st]); fol = 1.0 - phi
ids = lambda s: pd.Series(s).astype(str).str.strip().str.replace(r"\.0$", "", regex=True)  # noqa: E731


def lm(resid, p0, iters=300):                       # rc2_fit.lm, unchanged
    p = np.array(p0, float); r = resid(p); sse = float(r @ r); lam = 1e-3
    for _ in range(iters):
        Jm = np.empty((len(r), len(p)))
        for j in range(len(p)):
            h = 1e-6 * max(1.0, abs(p[j])); q = p.copy(); q[j] += h; Jm[:, j] = (resid(q) - r) / h
        A = Jm.T @ Jm; g = -Jm.T @ r; improved = False
        while lam < 1e12:
            stp = np.linalg.solve(A + lam * np.diag(np.maximum(np.diag(A), 1e-12)), g); pn = p + stp; rn = resid(pn); sn = float(rn @ rn)
            if np.isfinite(sn) and sn < sse:
                conv = (sse - sn) <= 1e-12 * max(sse, 1.0); p, r, sse, lam, improved = pn, rn, sn, lam / 3, True
                break
            lam *= 4
        if not improved or conv:
            break
    return p, sse


def step(B, net):                                   # rc2_fit.step, unchanged
    with np.errstate(divide="ignore", invalid="ignore"):
        dec = B * np.exp(np.where(B > 0, net / np.where(B > 0, B, 1.0), 0.0))
    return np.where(net >= 0, B + net, dec)


def recurse(th, B0, Bmax, z, P_at, G_at, smin, nyr):   # rc2_fit.recurse with cbm=True
    c0, c2, c3, gam = th; kap = np.exp(gam); n = len(B0); T_ = int(np.max(nyr)); B = B0.astype(float).copy(); Ms = np.zeros(n)
    path = np.full((n, T_ + 1), np.nan); path[:, 0] = B; lev = np.exp(c0 + c2 * z)
    for t in range(T_):
        act = t < nyr; Pt = P_at(t)
        with np.errstate(divide="ignore", invalid="ignore"):
            s = np.where((Pt > 0) & (B > 0), B / Pt, np.nan)
        s = np.where(np.isfinite(s) & (s > 0), s, smin)
        u = np.where(Bmax > 0, np.maximum(B, 0) / Bmax, 0.0)
        gr = lev * s ** c3 * G_at(t); mo = gr * u ** kap; net = gr - mo
        Bn = step(B, net); B = np.where(act, Bn, B); Ms = Ms + np.where(act, mo, 0.0); path[:, t + 1] = B
    return B, Ms, path


def ftg(x):
    if pd.isna(x):
        return 999
    x = int(x)
    return 999 if x == 999 else (500 if 500 <= x <= 520 else (x // 10) * 10)


# ---------------------------------------------------------------- site (ESI at public coordinates) and plot groups
es = pd.read_csv(T / f"restricted/esi_plot_{st}.csv", dtype={"CN": str}); es["CN"] = ids(es.CN); es = es.drop_duplicates("CN").set_index("CN").esi
cd = pd.read_csv(fia / f"{st}_COND.csv", usecols=["PLT_CN", "CONDID", "CONDPROP_UNADJ", "COND_STATUS_CD", "FORTYPCD"], dtype={"PLT_CN": str})
cd["PLT_CN"] = ids(cd.PLT_CN); cd = cd[cd.COND_STATUS_CD == 1].sort_values(["PLT_CN", "CONDPROP_UNADJ", "CONDID"], ascending=[True, False, True]).drop_duplicates("PLT_CN")
grp_of = dict(zip(cd.PLT_CN, [ftg(x) for x in cd.FORTYPCD]))
# C1 frame: plots of the current evaluation ST2501 with a forest condition; plot AGB = sum DRYBIO_AG x TPA_UNADJ of live trees (Mg/ha)
code = int(pd.read_csv(fia / f"{st}_PLOT.csv", usecols=["STATECD"], nrows=1).STATECD.iloc[0]); ev = code * 10000 + 2501
pa = pd.read_csv(fia / f"{st}_POP_PLOT_STRATUM_ASSGN.csv", usecols=["PLT_CN", "EVALID"], dtype={"PLT_CN": str}); pa["PLT_CN"] = ids(pa.PLT_CN)
c1 = sorted(set(pa[pa.EVALID == ev].PLT_CN) & set(grp_of)); c1s = set(c1); bs = {}
for ch in pd.read_csv(fia / f"{st}_TREE.csv", usecols=["PLT_CN", "STATUSCD", "DIA", "TPA_UNADJ", "DRYBIO_AG"], dtype={"PLT_CN": str}, chunksize=2_000_000):
    ch["PLT_CN"] = ids(ch.PLT_CN); ch = ch[ch.PLT_CN.isin(c1s) & (ch.STATUSCD == 1) & (ch.DIA >= 1)]
    for k, v in (ch.TPA_UNADJ * ch.DRYBIO_AG).groupby(ch.PLT_CN).sum().items():
        bs[k] = bs.get(k, 0.0) + v
c1B = np.array([bs.get(p, 0.0) * K for p in c1]); c1g = np.array([grp_of[p] for p in c1])
q = pd.Series(c1B).groupby(pd.Series(c1g)); p99_grp = {int(g): float(v) for g, v in q.quantile(0.99).items() if q.size()[g] >= MINCELL}; p99_all = float(np.quantile(c1B, 0.99))
bmax = lambda g: np.array([p99_grp.get(int(x), p99_all) for x in g])  # noqa: E731


def curves(unit_ids, plts, offs):
    U, C = PH.units_for(BJ, pd.DataFrame(dict(unit=unit_ids, PLT_CN=plts, kind="pair")))
    V = BJ.unit_values(E, U, C, np.asarray(offs)); return V


# ---------------------------------------------------------------- units
P = pd.read_csv(I / "pairs_std/pairs_built.csv", dtype={"pair_id": str, "t1_id": str}); P["pair_id"] = ids(P.pair_id); P["t1"] = ids(P.t1_id)
P = P.sort_values("pair_id", kind="stable").reset_index(drop=True); assert P.pair_id.is_unique
a1 = pd.read_parquet(I / "a1/a1_pairs.parquet"); a1["pair_id"] = ids(a1.pair_id); assert a1.pair_id.is_unique
P = P.merge(a1[["pair_id", "obs_net_a1b", "obs_mort"]].rename(columns={"obs_mort": "obs_mort_a1"}), on="pair_id", how="left", validate="one_to_one")
fo = pd.read_csv(I / "folds/pair_folds.csv", dtype={"pair_id": str}); fo["pair_id"] = ids(fo.pair_id); fo = fo[fo.set == "gate"]; assert fo.pair_id.is_unique
P = P.merge(fo[["pair_id", "fold"]], on="pair_id", how="left", validate="one_to_one"); assert P.fold.notna().all() and len(P) == len(fo)
Tr = pd.read_csv(I / "treated/treated_pairs.csv", dtype={"pair_id": str, "t1_id": str}); Tr["pair_id"] = ids(Tr.pair_id); Tr["t1"] = ids(Tr.t1_id)
Tr = Tr.sort_values("pair_id", kind="stable").reset_index(drop=True); assert Tr.pair_id.is_unique
for D in (P, Tr):
    D["n"] = [PG.cycle_of(x) for x in D.REMPER]; D["grp"] = [grp_of.get(p, 999) for p in D.t1]; D["esi"] = es.reindex(D.t1).values
if a.permute is not None:
    P = P.iloc[np.random.default_rng(a.permute).permutation(len(P))].reset_index(drop=True)


def pair_inputs(D):
    nmax = int(D.n.max()); V = curves(D.pair_id.values, D.t1.values, np.arange(0, nmax + 1))
    P0 = np.full(len(D), np.nan); inc = np.full(len(D), np.nan); Gc = np.full(len(D), np.nan)
    for i, (u, n) in enumerate(zip(D.pair_id.values, D.n.values)):
        v, note = V[u]
        if v is None:
            continue
        _, b, m = v; P0[i] = b[0] * K * fol; inc[i] = (b[n] - b[0]) * K * fol; Gc[i] = ((b[n] - b[0]) + (m[n] - m[0])) * K * fol
    return P0, inc, Gc


P0, incP, GcP = pair_inputs(P); okP = np.isfinite(P0) & np.isfinite(GcP) & np.isfinite(P.obs_net_a1b.values) & np.isfinite(P.obs_mort_a1.values)
Q0, incQ, GcQ = pair_inputs(Tr); okQ = np.isfinite(Q0) & np.isfinite(GcQ)


class Cov:
    def __init__(self, tr):
        site = pd.Series(P.esi.values[tr]); g = pd.Series(P.grp.values[tr])
        self.med = site.groupby(g).median().dropna().to_dict(); self.all = float(site.median())
        f = self.fill(P.esi.values[tr], P.grp.values[tr]); self.mu, self.sd = float(np.mean(f)), float(np.std(f))
        with np.errstate(divide="ignore", invalid="ignore"):
            s = P.B1_Mgha.values[tr] / P0[tr]
        self.smin = float(np.nanmin(s[np.isfinite(s) & (s > 0)]))

    def fill(self, si, grp):
        return np.array([v if np.isfinite(v) else self.med.get(int(gg), self.all) for v, gg in zip(np.asarray(si, float), grp)])

    def z(self, si, grp):
        return (self.fill(si, grp) - self.mu) / self.sd


def predict(th, cv, B0, P0_, inc_, Gc_, n, esi, grp):
    z = cv.z(esi, grp); Bm = bmax(grp)
    Pa = lambda t: P0_ + inc_ * t / n; Ga = lambda t: Gc_ / n  # noqa: E731
    Be, Ms, _ = recurse(th, B0, Bm, z, Pa, Ga, cv.smin, n)
    return (Be - B0) / n, Ms / n


def fit(tr):
    cv = Cov(tr); m = tr & okP
    args = (P.B1_Mgha.values[m], P0[m], incP[m], GcP[m], P.n.values[m], P.esi.values[m], P.grp.values[m]); on, om = P.obs_net_a1b.values[m], P.obs_mort_a1.values[m]
    def res(th):
        pn, pm = predict(th, cv, *args); r = np.r_[on - pn, om - pm]; return np.where(np.isfinite(r), r, 0.0)
    th, sse = lm(res, [0.0, 0.0, 0.0, 0.0]); return th, sse, cv, int(m.sum())


out = {}; netP = np.full(len(P), np.nan); morP = np.full(len(P), np.nan); par = {}
for fd in sorted(P.fold.unique()):
    th, sse, cv, ntr = fit(P.fold.values != fd); te = (P.fold.values == fd) & np.isfinite(P0) & np.isfinite(GcP)
    pn, pm = predict(th, cv, P.B1_Mgha.values[te], P0[te], incP[te], GcP[te], P.n.values[te], P.esi.values[te], P.grp.values[te])
    netP[te], morP[te] = pn, pm; par[f"fold{int(fd)}"] = dict(theta=th.tolist(), sse=sse, n_train=ntr)
thp, ssep, cvp, ntrp = fit(np.ones(len(P), bool)); par["prod"] = dict(theta=thp.tolist(), sse=ssep, n_train=ntrp)
mQ = okQ
qn, qm = predict(thp, cvp, Tr.B1_Mgha.values[mQ], Q0[mQ], incQ[mQ], GcQ[mQ], Tr.n.values[mQ], Tr.esi.values[mQ], Tr.grp.values[mQ])
netQ = np.full(len(Tr), np.nan); morQ = np.full(len(Tr), np.nan); netQ[mQ], morQ[mQ] = qn, qm
O = Path(a.slots) / "cbm_gated_x1" / st / "restricted"; O.mkdir(parents=True, exist_ok=True); os.chmod(O, 0o700); tag = "" if a.permute is None else f"_perm{a.permute}"


def rows(D, net, mor):
    ok = np.isfinite(net) & np.isfinite(mor); n = D.n.values
    return pd.DataFrame(dict(pair_id=D.pair_id.values, arm="X1", ok=ok, note=np.where(ok, "x1", "no curve"), fvs_surv=np.where(ok, (net + mor) * n / K, np.nan),
                             fvs_mort=np.where(ok, mor * n / K, np.nan), eng_ingr=0.0, cycle=n, state=st))
if a.permute is None:
    rows(P, netP, morP).to_csv(O / "pairs_level__cbm_gated_x1__X1__pair_rows.csv", index=False, float_format="%.12g")
    rows(Tr, netQ, morQ).to_csv(O / "treated__cbm_gated_x1__X1__pair_rows.csv", index=False, float_format="%.12g")
# ---------------------------------------------------------------- chains (P2), production fit
Hs = S.h1_setup(H / f"perseus_me/tracks/multistate/config/multistate/hindcast/{st}.yaml", I / "hindcast")
h1 = list(Hs["h1"]); v1 = Hs["v1plt"].reindex(h1).astype(str).str.strip(); o1 = Hs["O"][Hs["O"].visit == 1].set_index("chain_id")
ch = pd.DataFrame(dict(chain_id=h1, PLT_CN=ids(v1.values).values)); ch["agb1"] = o1.agb_stock.reindex(h1).values; ch["agc1"] = o1.agc_stock.reindex(h1).values
ch["grp"] = [grp_of.get(p, 999) for p in ch.PLT_CN]; ch["esi"] = es.reindex(ch.PLT_CN).values
off = np.arange(0, 31); V = curves(ch.chain_id.values, ch.PLT_CN.values, off); keep = []; Pm, Gm = [], []
for i, u in enumerate(ch.chain_id.values):
    v, _ = V[u]
    if v is None or not np.isfinite(ch.agb1.values[i]):
        continue
    _, b, m = v; keep.append(i); Pm.append(b * K * fol); Gm.append(np.r_[np.diff(b) + np.diff(m), 0.0] * K * fol)
keep = np.array(keep); Pm = np.array(Pm); Gm = np.array(Gm); cc = ch.iloc[keep].reset_index(drop=True)
B0 = cc.agb1.values * K; zc = cvp.z(cc.esi.values, cc.grp.values)
_, _, path = recurse(thp, B0, bmax(cc.grp.values), zc, lambda t: Pm[:, t], lambda t: Gm[:, t], cvp.smin, np.full(len(cc), 30))
rat = np.where(cc.agb1.values > 0, cc.agc1.values / cc.agb1.values, 0.5)
X = pd.DataFrame([(c, o, path[i, o] / K, path[i, o] / K * rat[i]) for i, c in enumerate(cc.chain_id.values) for o in off], columns=["chain_id", "offset", "agb_full", "agc_full"])
X["agb_surv"] = np.nan; X["agc_surv"] = np.nan
if a.permute is None:
    (O / "h1").mkdir(exist_ok=True); X.to_parquet(O / "h1/cbm_gated_x1__X1__chain.parquet", index=False)
# ---------------------------------------------------------------- 100-yr realism gates on the C1 plots (no harvest)
V = curves(np.array(c1), np.array(c1), np.arange(0, 101)); kp, Pc, Gcx = [], [], []
for i, u in enumerate(c1):
    v, _ = V[u]
    if v is None:
        continue
    _, b, m = v; kp.append(i); Pc.append(b * K * fol); Gcx.append(np.r_[np.diff(b) + np.diff(m), 0.0] * K * fol)
kp = np.array(kp); Pc = np.array(Pc); Gcx = np.array(Gcx); B0c = c1B[kp]; gk = c1g[kp]
zk = cvp.z(es.reindex(np.array(c1)[kp]).values, gk)
_, _, pth = recurse(thp, B0c, bmax(gk), zk, lambda t: Pc[:, t], lambda t: Gcx[:, t], cvp.smin, np.full(len(kp), 100))
inc_first = np.mean(pth[:, 10] - pth[:, 0]) / 10; inc_last = np.mean(pth[:, 100] - pth[:, 90]) / 10
gates = dict(M1_no_negative=bool((pth >= 0).all()), M2_late_increment_below_early=bool(inc_last < inc_first), mean_y0=float(B0c.mean()), mean_y100=float(pth[:, 100].mean()),
             p99_y100=float(np.quantile(pth[:, 100], 0.99)), current_max=float(c1B.max()), current_p99=float(np.quantile(c1B, 0.99)), plots=int(len(kp)),
             M3_p99_y100_le_current_max=bool(np.quantile(pth[:, 100], 0.99) <= c1B.max()), M4_mean_y100_le_current_p99=bool(pth[:, 100].mean() <= np.quantile(c1B, 0.99)))
res = dict(state=st, phi=phi, folds=sorted(int(x) for x in P.fold.unique()), params=par, p99_grp=p99_grp, p99_all=p99_all, c1_plots=len(c1),
           pairs_ok=int(np.isfinite(netP).sum()), pairs=len(P), treated_ok=int(np.isfinite(netQ).sum()), treated=len(Tr), chains_valued=int(len(cc)), chains=len(ch),
           esi_missing_share_pairs=float(np.mean(~np.isfinite(P.esi.values))), gates=gates, permute=a.permute, wall_s=round(time.time() - t0, 1))
(T / "work" / f"x1wp_{st}{tag}.json").write_text(json.dumps(res, indent=1, default=float)); print(json.dumps({k: res[k] for k in ("params", "pairs_ok", "treated_ok", "chains_valued", "gates")}, default=float))
