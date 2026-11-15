"""Track trkB_engfix DDM fit (PREREG section 1): logit(p_annual) = a_k + b RD by DBH class, P(died) = 1 - (1 - p)^REMPER,
Bernoulli, weights TPA_UNADJ, training counties of each county-blocked fold of refine_engines2 (frames_cache.pkl) and all.
Inputs: restricted/A_rd.csv (pid, rd = FVS BSDI/BSDIMAX at cycle 1, Cardinal job A). Writes work/ddm_params.json (parameters)
and work/ddm_fit.json (counts, calibration ratios); tree data stay in memory."""
import json, pickle, sys
import numpy as np, pandas as pd
sys.path.insert(0, "/home/aaron/perseus_me/tracks/trkB_engfix/tools"); import b_lib as BL
T2 = BL.PME / "tracks/trkB_engfix"
fc = pickle.load(open(BL.PME / "tracks/refine_engines2/restricted/frames_cache.pkl", "rb")); C2 = fc["C2"]
U = C2[C2.uncut_fit.values][["pair_id", "fold"]].copy(); assert not U.pair_id.duplicated().any()
rd = pd.read_csv(T2 / "restricted/A_rd.csv", dtype={"pid": str}).set_index("pid").rd
PB = pd.read_csv(BL.PME / "pairs/me_state_std/pairs_built.csv", dtype={"pair_id": str}); PB["pid"] = PB.pair_id.str.strip()
rem = PB.set_index("pid").REMPER
T = pd.read_csv(BL.PME / "pairs/me_state_std/trees_built.csv", dtype={"pair_id": str, "tree_cn": str}, usecols=["pair_id", "tree_cn", "DIA", "TPA_UNADJ", "DRYBIO_AG"])
T["pid"] = T.pair_id.str.strip(); T["tree_cn"] = T.tree_cn.str.strip()
T = T[T.pid.isin(set(U.pair_id))]
F = pd.read_csv("/home/aaron/fia_data/ME_TREE.csv", usecols=["PREV_TRE_CN", "STATUSCD"], dtype={"PREV_TRE_CN": str})
F = F[F.PREV_TRE_CN.notna()]; F["PREV_TRE_CN"] = F.PREV_TRE_CN.str.strip()
assert not F.PREV_TRE_CN.duplicated().any(), "duplicate PREV_TRE_CN"
st = F.set_index("PREV_TRE_CN").STATUSCD
T["st"] = T.tree_cn.map(st)
n_all = len(T); T = T[T.st.isin([1, 2])].copy()
T["y"] = (T.st == 2).astype(float); T["rd"] = T.pid.map(rd).values; T["L"] = T.pid.map(rem).values
T["fold"] = T.pid.map(U.set_index("pair_id").fold).values
assert T[["rd", "L", "fold"]].notna().all().all(), T[["rd", "L", "fold"]].isna().sum().to_dict()
edges = [lo for lo, hi in BL.DCLS] + [1e9]
T["k"] = np.searchsorted(np.array(edges[1:-1]), T.DIA.values, side="right")


def nll(th, X):
    a, b = th[:7], th[7]
    p = expit(a[X.k.values] + b * X.rd.values); P = 1.0 - (1.0 - p) ** X.L.values
    P = np.clip(P, 1e-12, 1 - 1e-12); w = X.TPA_UNADJ.values
    return -np.sum(w * (X.y.values * np.log(P) + (1 - X.y.values) * np.log(1 - P))) / w.sum()


def expit(x):
    return 1.0 / (1.0 + np.exp(-x))


def grad(th, X):
    """analytic gradient of nll: d/d eta = -w (y/P - (1-y)/(1-P)) L p (1-p)^L / sum w"""
    a, b = th[:7], th[7]; k = X.k.values; r = X.rd.values; L = X.L.values; y = X.y.values; w = X.TPA_UNADJ.values
    p = expit(a[k] + b * r); P = np.clip(1.0 - (1.0 - p) ** L, 1e-12, 1 - 1e-12)
    ge = -w * (y / P - (1 - y) / (1 - P)) * L * p * (1.0 - p) ** L / w.sum()
    return np.r_[np.bincount(k, ge, minlength=7), np.sum(ge * r)]


class Res:
    pass


def fit(X):
    """Newton with a finite-difference Hessian of the analytic gradient and backtracking (numpy only; scipy absent in the venv)"""
    th = np.r_[np.full(7, -4.0), 0.0]; f = nll(th, X); it = 0
    for it in range(200):
        g = grad(th, X); h = 1e-6; H = np.zeros((8, 8))
        for j in range(8):
            e = np.zeros(8); e[j] = h; H[:, j] = (grad(th + e, X) - grad(th - e, X)) / (2 * h)
        H = 0.5 * (H + H.T)
        try:
            step = -np.linalg.solve(H + 1e-10 * np.eye(8), g)
        except np.linalg.LinAlgError:
            step = -g
        t = 1.0; f2 = nll(th + step, X)
        while t > 1e-8 and f2 > f:
            t *= 0.5; f2 = nll(th + t * step, X)
        th = th + t * step; done = abs(f - f2) < 1e-13 and np.max(np.abs(g)) < 1e-7; f = min(f, f2)
        if done:
            break
    r = Res(); r.x = th; r.fun = f; g = grad(th, X); r.success = bool(np.max(np.abs(g)) < 1e-6); r.message = f"newton it={it} max|grad|={np.max(np.abs(g)):.2e}"
    return r


def calib(th, X):
    a, b = th[:7], th[7]; p = expit(a[X.k.values] + b * X.rd.values); P = 1.0 - (1.0 - p) ** X.L.values
    w = X.TPA_UNADJ.values; wb = w * X.DRYBIO_AG.values
    return dict(stems_pred_over_obs=float((w * P).sum() / (w * X.y.values).sum()), biomass_pred_over_obs=float((wb * P).sum() / (wb * X.y.values).sum()))


par, diag = {}, {}
for f in ["all"] + list(range(5)):
    tr = T if f == "all" else T[T.fold != f]
    r = fit(tr); th = r.x
    par[str(f)] = dict(a=[float(x) for x in th[:7]], b=float(th[7]))
    d = dict(converged=bool(r.success), msg=str(r.message)[:60], nll=float(r.fun), n_trees=int(len(tr)), n_pairs=int(tr.pid.nunique()),
             deaths=int(tr.y.sum()), train=calib(th, tr))
    if f != "all":
        te = T[T.fold == f]; d["test"] = calib(th, te); d["n_test_trees"] = int(len(te))
    diag[str(f)] = d
diag["_counts"] = dict(t1_trees_uncut_fit=int(n_all), with_fate_1_or_2=int(len(T)), rd_range=[float(T.rd.min()), float(T.rd.max())],
                       trees_by_class=T.groupby("k").size().astype(int).to_dict(), deaths_by_class=T.groupby("k").y.sum().astype(int).to_dict())
json.dump(par, open(T2 / "work/ddm_params.json", "w"), indent=1); json.dump(diag, open(T2 / "work/ddm_fit.json", "w"), indent=1, default=str)
print(json.dumps(par["all"])); print(json.dumps({k: (v.get("train"), v.get("test")) for k, v in diag.items() if k != "_counts"}))
print(json.dumps(diag["_counts"], default=str))
