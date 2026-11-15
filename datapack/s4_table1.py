"""Table 1 data overview per state. C1 targets: the Track R / multistate design estimator (rc_design, EVALID ST2401 stock,
ST2403 GRM) with the scorer's shared location bootstrap (LocDraws, seed 20261103 = ns_score SEED, 1,000 draws); 95% CI =
percentile 2.5/97.5 of the target draws (the draws C1 compares against). Points asserted equal to targets_six_states.csv.
Pair sets: gate (uncut common) and treated (cut) sets of pair_folds.csv, the sets A.load_state builds the scorer bundle from.
Aggregates only (no plot keys written). usage: s4_table1.py"""
import sys, os, json
from pathlib import Path
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import *
sys.path.insert(0, f"{HOME}/perseus_me/tracks/multistate/tools/ncasi_rc/score")
import rc_design as D
SEED, B = 20261103, 1000
RC = Path(HOME) / "perseus_me/ncasi_rc"
tg6 = pd.read_csv(RC / "targets/targets_six_states.csv"); tA = pd.read_csv(RC / "targets/targets_pathA.csv")
ref = pd.read_csv(Path(HOME) / "perseus_me/ext/FIA/ENTIRE_REF_FOREST_TYPE.csv"); refg = pd.read_csv(Path(HOME) / "perseus_me/ext/FIA/ENTIRE_REF_FOREST_TYPE_GROUP.csv")
ft2g = dict(zip(ref.VALUE.astype(int), ref.TYPGRPCD.astype(int))); gname = dict(zip(refg.VALUE.astype(int), refg.MEANING))
inputs = [str(RC / "targets/targets_six_states.csv"), str(RC / "targets/targets_pathA.csv"),
          f"{HOME}/perseus_me/ext/FIA/ENTIRE_REF_FOREST_TYPE.csv", f"{HOME}/perseus_me/ext/FIA/ENTIRE_REF_FOREST_TYPE_GROUP.csv"]
rows = []
for st in STATES:
    I = RC / "inputs" / st / "restricted"
    pf = pd.read_csv(I / "folds/pair_folds.csv", dtype=str)
    pb = pd.read_csv(I / "pairs_std/pairs_built.csv", dtype={"pair_id": str}, usecols=["pair_id", "INVYR1", "INVYR2", "REMPER"])
    tp = pd.read_csv(I / "treated/treated_pairs.csv", dtype={"pair_id": str}, usecols=["pair_id", "INVYR1", "INVYR2", "REMPER", "t40", "obs_net_a1b"])
    a1 = pd.read_parquet(I / "a1/a1_pairs.parquet", columns=["pair_id", "obs_net_a1b"]); a1["pair_id"] = a1.pair_id.astype(str)
    for x in (pf, pb, tp): x["pair_id"] = x.pair_id.astype(str).str.strip().str.replace(r"\.0$", "", regex=True)
    a1["pair_id"] = a1.pair_id.str.strip().str.replace(r"\.0$", "", regex=True)
    g = pf[pf.set == "gate"].drop_duplicates("pair_id")[["pair_id"]]; t = pf[pf.set == "treated"].drop_duplicates("pair_id")[["pair_id"]]
    assert pb.pair_id.is_unique and tp.pair_id.is_unique and a1.pair_id.is_unique
    G = g.merge(pb, on="pair_id", how="left", validate="one_to_one").merge(a1, on="pair_id", how="left", validate="one_to_one")
    T = t.merge(tp, on="pair_id", how="left", validate="one_to_one")
    assert G.REMPER.notna().all() and T.REMPER.notna().all(), st
    r = dict(state=st, n_uncut_pairs_gate_set=len(G), n_uncut_pairs_obs_net_finite=int(np.isfinite(G.obs_net_a1b).sum()),
             n_cut_pairs_treated_set=len(T), n_cut_pairs_obs_net_finite=int(np.isfinite(T.obs_net_a1b).sum()),
             n_cut_pairs_without_TRTCD40=int((~T.t40.astype(str).str.lower().isin(["true", "1"])).sum()),
             remper_median_uncut_yr=float(G.REMPER.median()), remper_median_cut_yr=float(T.REMPER.median()),
             remper_q25_uncut_yr=float(G.REMPER.quantile(.25)), remper_q75_uncut_yr=float(G.REMPER.quantile(.75)),
             t1_invyr_range_uncut=f"{int(G.INVYR1.min())}-{int(G.INVYR1.max())}", t2_invyr_range_uncut=f"{int(G.INVYR2.min())}-{int(G.INVYR2.max())}",
             t1_invyr_range_cut=f"{int(T.INVYR1.min())}-{int(T.INVYR1.max())}", t2_invyr_range_cut=f"{int(T.INVYR2.min())}-{int(T.INVYR2.max())}")
    inputs += [str(I / "folds/pair_folds.csv"), str(I / "pairs_std/pairs_built.csv"), str(I / "treated/treated_pairs.csv"), str(I / "a1/a1_pairs.parquet")]
    # C1 targets with the scorer's draws
    dr = D.LocDraws(st, B, SEED); ev_cur = D.statecd(st) * 10000 + 2401; ev_grm = ev_cur + 2
    ys, dc = D.plot_stock(st, ev_cur); pt, bt = D.design_total(dc, ys, dr.mult(dc["loc"].values))
    T6 = tg6[tg6.state == st].set_index("quantity")
    def put(name, pt, bt, q, unit):
        pt *= D.LB2TG; bt = np.asarray(bt) * D.LB2TG
        ref_v = float(T6.loc[q, "ours"]); assert abs(pt / ref_v - 1) < 1e-6, (st, q, pt, ref_v)
        r[f"{name}_{unit}"] = pt; r[f"{name}_lo95_{unit}"] = float(np.percentile(bt, 2.5)); r[f"{name}_hi95_{unit}"] = float(np.percentile(bt, 97.5))
        r[f"{name}_se_draws_{unit}"] = float(np.std(bt, ddof=1))
    put("agb_stock", pt, bt, "stock", "Tg")
    yg, dg = D.plot_grm(st, ev_grm); Mg = dr.mult(dg["loc"].values)
    for q in ("gross_growth", "mortality", "net_growth"):
        p, b_ = D.design_total(dg, yg[q], Mg); put(q, p, b_, q, "Tg_yr")
    r["evalid_stock"] = ev_cur; r["evalid_grm"] = ev_grm; r["n_outside_universe_locs"] = dr.n_outside
    # forest area (recorded path A value) and forest type groups (same condition-area formula as rc_adapt._area_weights)
    fa = tA[(tA.state == st) & (tA.quantity == "forest_area_ha")]; assert len(fa) == 1
    r["forest_area_ha"] = float(fa.ours.iloc[0]); r["forest_area_evalid"] = int(fa.evalid.iloc[0])
    d = D.design(st, ev_cur)
    cd = D._rd(st, "COND", ("PLT_CN", "CONDID", "COND_STATUS_CD", "CONDPROP_UNADJ", "FORTYPCD", "PROP_BASIS"))
    cf = cd[cd.COND_STATUS_CD == 1].merge(d[["PLT_CN", "ADJ_FACTOR_SUBP", "ADJ_FACTOR_MACR", "c_h", "P2POINTCNT"]], on="PLT_CN")
    cf["a"] = cf.CONDPROP_UNADJ * np.where(cf.PROP_BASIS == "MACR", cf.ADJ_FACTOR_MACR, cf.ADJ_FACTOR_SUBP) * cf.c_h / cf.P2POINTCNT
    cf["grp"] = [ft2g.get(int(x), 999) if pd.notna(x) else 999 for x in cf.FORTYPCD]
    ag = cf.groupby("grp").a.sum(); tot_ha = float(ag.sum()) * 0.40468564224
    r["forest_area_recomputed_ha"] = tot_ha; r["forest_area_recomputed_rel_diff"] = tot_ha / r["forest_area_ha"] - 1
    sh = (ag / ag.sum()).sort_values(ascending=False)
    for i, (k, v) in enumerate(sh.head(3).items(), 1):
        r[f"ftg{i}_code"] = int(k); r[f"ftg{i}_name"] = gname.get(int(k), "unknown"); r[f"ftg{i}_area_share"] = float(v)
    rows.append(r); D._rd.cache_clear(); print(st, "ok", flush=True)
out = pd.DataFrame(rows)
out.to_csv(f"{DATA}/table1_data_overview.csv", index=False)
record("table1_data_overview.csv", sorted(set(inputs)), __file__,
       "C1 targets: rc_design path B estimator (aboveground live dry biomass DRYBIO_AG, forest land, live DIA>=1 in; GRM on ST2403), "
       "points equal targets_six_states.csv (asserted, rel 1e-6); 95% CI = percentile of 1,000 shared location bootstrap draws, seed 20261103 (scorer seed). "
       "forest_area_ha from targets_pathA.csv (path A, CONDPROP x adjustment, EXPVOL eval); forest type group shares recomputed on EVALID ST2401 with the "
       "rc_adapt._area_weights condition-area formula (acres x 0.40468564224 to ha). Pair sets: pair_folds gate (uncut) and treated (cut) sets. "
       "Restricted inputs read server side; only state aggregates written.")
status(f"step4 table1_data_overview.csv done ({len(out)} rows)")
