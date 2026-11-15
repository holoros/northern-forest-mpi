"""trkA_hcm (PREREG md5 35055a79af9083edec3f5e3223552dd7) section 8, derived from Track IMF2 section 6 (adds forecast cut status from the state HCM2): frozen Maine links applied out of state, observed cut status. usage: imf2_oos.py ST [--shuffle N]
Statewide S (domain area weights need the state inventory; not used, as track IMF D6), location-component cluster draws."""
import argparse, glob, json, os, sys
from pathlib import Path
import numpy as np, pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parent))
import imf2_lib as L
IC = L.IC
ap = argparse.ArgumentParser(); ap.add_argument("st"); ap.add_argument("--shuffle", type=int, default=None); a = ap.parse_args()
ST = a.st.upper(); RC = L.PME / "ncasi_rc"; NS = L.PME / "ncasi_states"
OUT = L.T / "work" / f"oos_{ST}" / ("main" if a.shuffle is None else f"shuf{a.shuffle}"); OUT.mkdir(parents=True, exist_ok=True)
ROWS = {  # engine -> (uncut file glob, cut file glob, arm or None)
    "osm": (f"{RC}/runs/osm/{ST}/restricted/pairs_all/pair_rows.csv", f"{RC}/runs/osm/{ST}/restricted/treated/pair_rows.csv", None),
    "fvsd": (f"{RC}/runs/fvs/{ST}/restricted/pairs_all_default/pair_rows.csv", f"{RC}/runs/fvs/{ST}/restricted/treated_default/pair_rows.csv", "default"),
    "fvsa": (f"{NS}/runs/fvs_acd_prod/{ST}/restricted/pairs_all_default/pair_rows.csv", f"{NS}/runs/fvs_acd_prod/{ST}/restricted/treated_default/pair_rows.csv", None),
    "cem": (f"{RC}/runs/cem/{ST}/restricted/pairs_all/extract/pair_rows.csv", f"{RC}/runs/cem/{ST}/restricted/treated/extract/pair_rows.csv", None),
    "cbmg": (f"{RC}/runs/cbm_gated/{ST}/restricted/pairs_level__*pair_rows.csv", f"{RC}/runs/cbm_gated/{ST}/restricted/treated_ageq__*pair_rows.csv", None),
    "cbma": (f"{RC}/runs/cbm_af/{ST}/restricted/pairs_level__*pair_rows.csv", f"{RC}/runs/cbm_af/{ST}/restricted/treated_ageq__*pair_rows.csv", None),
    "landis": ("", f"{RC}/runs/landis/{ST}/restricted/c5/rows_cut/pair_rows_nofol.csv", None),
}
base = RC / ST / "restricted" if (RC / ST / "restricted/pairs_std").exists() else RC / "inputs" / ST / "restricted"
PB = pd.read_csv(base / "pairs_std/pairs_built.csv", dtype={"pair_id": str, "t1_id": str}); PB["pair_id"] = L.strip(PB.pair_id); PB["t1_id"] = L.strip(PB.t1_id)
A1 = pd.read_parquet(base / "a1/a1_pairs.parquet"); A1.index = L.strip(A1.pair_id.astype(str))
nb = pd.Series(np.where(A1.has_grm, A1.grm_net_mgha, A1.obs_net), index=A1.index) if {"has_grm", "grm_net_mgha"} <= set(A1.columns) else None
TP = pd.read_csv(base / "treated/treated_pairs.csv", dtype={"pair_id": str, "t1_id": str}); TP["pair_id"] = L.strip(TP.pair_id); TP["t1_id"] = L.strip(TP.t1_id)
U = PB[~PB.pair_id.isin(set(TP.pair_id))].copy(); U["frame"] = "U"
U["obs_g"] = A1.obs_gg_a1.reindex(U.pair_id).values; U["obs_m"] = U.obs_mort
U["obs_n"] = nb.reindex(U.pair_id).values if nb is not None else U.obs_net
Cc = TP.copy(); Cc["frame"] = "C"; Cc["obs_g"] = Cc.obs_gg_a1; Cc["obs_m"] = Cc.obs_mort; Cc["obs_n"] = Cc.obs_net_a1b
cols = ["pair_id", "t1_id", "frame", "FORTYPCD", "STDAGE", "B1_Mgha", "obs_g", "obs_m", "obs_n"]
F = pd.concat([U[cols], Cc[cols + ["B1full_Mgha"]]], ignore_index=True)
if a.shuffle is not None:
    F = F.iloc[np.random.default_rng(a.shuffle).permutation(len(F))]
F = F.sort_values(["frame", "pair_id"]).reset_index(drop=True)
par = {}


def find(x):
    while par.get(x, x) != x:
        par[x] = par.get(par[x], par[x]); x = par[x]
    return x


for a_, b_ in sorted(zip(F.pair_id, F.t1_id)):
    ra, rb = find(a_), find(b_)
    if ra != rb:
        par[max(ra, rb)] = min(ra, rb)
F["loc"] = [find(x) for x in F.pair_id]
doms, aw, dmap, _, _ = L.domains(); F["dom"] = F.FORTYPCD.map(dmap)
cutset = set(TP.pair_id); pre = TP.set_index("pair_id").B1full_Mgha
F["lane_u"] = np.select([F.t1_id.isin(cutset).values & (F.B1_Mgha.values < pre.reindex(F.t1_id).values), F.STDAGE.values >= 100], ["L2b", "L3"], "L1")
c_pt = IC.foliage_factor()
found = {}
for e, (gu, gc, arm) in ROWS.items():
    for q in L.Q:
        F[f"{e}_{q}"] = np.nan
    for fr, gl in (("U", gu), ("C", gc)):
        hits = sorted(glob.glob(gl)) if gl else []
        found[f"{e}|{fr}"] = hits[0].replace(str(L.PME) + "/", "") if hits else None
        if not hits:
            continue
        X = IC.read_rows(hits[0], arm)
        if e in ("cbmg", "cbma", "landis"):
            X = X * c_pt
        msk = F.frame.eq(fr).values
        for q in L.Q:
            F.loc[msk, f"{e}_{q}"] = X[q].reindex(F.loc[msk, "pair_id"]).values
for e in L.ENGS:
    for q in L.Q:
        if f"{e}_{q}" not in F.columns:
            F[f"{e}_{q}"] = np.nan
fz = {tuple(k.split("|")): v for k, v in json.loads((L.T / "work/g1/imf2_links_frozen.json").read_text()).items()}


def route(idx, lanes, q):
    v = np.full(len(idx), np.nan); fb = np.zeros(len(idx), bool)
    for j, i in enumerate(idx):
        lane = lanes[j]; r = fz.get((lane, F.dom.values[i], q), fz[(lane, "lane", q)])
        for k, e in enumerate(r):
            if q == "all":
                x = F[f"{e}_n"].values[i]
            else:
                x = F[f"{e}_{q}"].values[i]
            if np.isfinite(x):
                v[j] = x; fb[j] = k > 0; break
    return v, fb


n = len(F); isC = F.frame.eq("C").values; lanes = np.where(isC, "L2a", F.lane_u.values); idx = np.arange(n)
g, fbg = route(idx, lanes, "g"); m, fbm = route(idx, lanes, "m")
M = {"IMF2_V1_oracle_frozen": {"g": g, "m": m, "n": g - m}}
# trkA forecast cut status: p from the state HCM2 (county-fold OOF), joined by pair id = t2 PLT_CN; unmatched pairs take the state base rate
HP = pd.read_parquet(L.T / f"restricted/states/{ST}/hcm2.parquet"); HP.index = L.strip(HP.CN.astype(str)); assert not HP.index.duplicated().any()
pf = HP.p.reindex(F.pair_id.values).values.astype(float); n_unmatched = int(np.isnan(pf).sum()); base_rate = float(HP.y.mean()); pf = np.where(np.isfinite(pf), pf, base_rate)
gU, _ = route(idx, F.lane_u.values, "g"); mU, _ = route(idx, F.lane_u.values, "m"); lc = np.array(["L2a"] * n); gC, _ = route(idx, lc, "g"); mC, _ = route(idx, lc, "m")
fg = pf * gC + (1 - pf) * gU; fm = pf * mC + (1 - pf) * mU
M["IMF2_V1_forecast_frozen"] = {"g": fg, "m": fm, "n": fg - fm}
# single-engine rule (V7): one engine per cell for all three quantities
sg = np.full(n, np.nan); sm = sg.copy(); sn = sg.copy()
for i in range(n):
    r = fz.get((lanes[i], F.dom.values[i], "all"), fz[(lanes[i], "lane", "all")])
    for e in r:
        if all(np.isfinite(F[f"{e}_{q}"].values[i]) for q in L.Q):
            sg[i], sm[i], sn[i] = (F[f"{e}_{q}"].values[i] for q in L.Q); break
M["IMF2_V7_oracle_frozen"] = {"g": sg, "m": sm, "n": sn}
for k, mem in L.MIX.items():
    have = [e for e in mem if np.isfinite(F[f"{e}_n"]).any()]
    M[f"{k}({len(have)}/{len(mem)} members)"] = {q: np.nanmean(np.stack([F[f"{e}_{q}"].values for e in have], 0), 0) for q in L.Q}
for e in ("osm", "fvsd"):
    M[f"{e} (single engine)"] = {q: F[f"{e}_{q}"].values for q in L.Q}
M = {k_: v_ for k_, v_ in M.items() if np.isfinite(v_["n"]).any()}   # trkA-D1: drop comparators with no rows in the state (MI, WI, MN have no OSM rows)
Wb = IC.boot_weights(F)
rows, wins = [], []
for fr in ("U", "C"):
    fin = F.frame.eq(fr).values & np.isfinite(F.obs_g.values) & np.isfinite(F.obs_m.values) & np.isfinite(F.obs_n.values)
    for mm in M.values():
        for q in L.Q:
            fin &= np.isfinite(mm[q])
    res = {}
    for k, mm in M.items():
        s = IC.score(mm, F, fin, doms, aw, Wb); res[k] = s
        sn_ = IC.score(mm, F, fin, doms, aw, Wb, qset=("n",))
        lo, hi = IC.ci(s["S_state_bt"]); lo2, hi2 = IC.ci(sn_["S_state_bt"])
        rows.append(dict(state=ST, frame=fr, method=k, pairs=int(fin.sum()), S_state=s["S_state"], S_state_lo=lo, S_state_hi=hi,
                         S_net_state=sn_["S_state"], S_net_lo=lo2, S_net_hi=hi2, ratio_g=s["ratio_g"], ratio_m=s["ratio_m"], ratio_n=s["ratio_n"]))
    for k in M:
        if k not in ("IMF2_V1_oracle_frozen", "IMF2_V1_forecast_frozen"):
            d = res["IMF2_V1_oracle_frozen"]["S_state_bt"] - res[k]["S_state_bt"]; lo, hi = IC.ci(d)
            wins.append(dict(state=ST, frame=fr, method="IMF2_V1_oracle_frozen", comparator=k, dS=res["IMF2_V1_oracle_frozen"]["S_state"] - res[k]["S_state"], dS_lo=lo, dS_hi=hi, significant_win=bool(hi < 0)))
        if k != "IMF2_V1_forecast_frozen":
            d = res["IMF2_V1_forecast_frozen"]["S_state_bt"] - res[k]["S_state_bt"]; lo, hi = IC.ci(d)
            wins.append(dict(state=ST, frame=fr, method="IMF2_V1_forecast_frozen", comparator=k, dS=res["IMF2_V1_forecast_frozen"]["S_state"] - res[k]["S_state"], dS_lo=lo, dS_hi=hi, significant_win=bool(hi < 0)))
pd.DataFrame(rows).to_csv(OUT / "oos_scores.csv", index=False, float_format="%.6g"); pd.DataFrame(wins).to_csv(OUT / "oos_wins.csv", index=False, float_format="%.6g")
json.dump(dict(state=ST, inputs_found=found, pairs={f"{a_}|{b_}": int(v) for (a_, b_), v in F.assign(l=lanes).groupby(["frame", "l"]).size().items()},
               fallback_share_g=float(fbg.mean()), fallback_share_m=float(fbm.mean()), status="forecast (state HCM2) and observed cut status", forecast_unmatched_pairs=n_unmatched, forecast_base_rate=base_rate, mean_p_U=float(pf[F.frame.eq("U").values].mean()), mean_p_C=float(pf[F.frame.eq("C").values].mean()),
               engines_missing=[e for e in ("osm_mc", "fvs_gc", "acd") if True]), open(OUT / "oos_meta.json", "w"), indent=1)
print(pd.DataFrame(rows).round(4).to_string()); print(pd.DataFrame(wins).round(4).to_string())
sys.stdout.flush(); os._exit(0)
