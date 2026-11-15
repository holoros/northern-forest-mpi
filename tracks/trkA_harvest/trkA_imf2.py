# trkA_hcm copy of tracks/imf2/tools (unedited except input paths: imf2 restricted inputs read only, outputs under trkA_hcm)
"""trkA_hcm IMF2 rerun with the trkA layers (PREREG md5 35055a79af9083edec3f5e3223552dd7 section 7), derived from Track IMF2 main run (PREREG.md md5 b539cc149e8f1e1606df0b1d46ccc065). usage: imf2_run.py TAG [--shuffle SEED]
Out-of-fold linked routing with a cross-fitted harvest-choice layer; scores C2 (P*, by domain) and C5 net (statewide) on the
registered draws; MPI with inherited no-harvest pillars; variants V1..V8 under forecast and oracle cut status."""
import argparse, json, os, sys
from pathlib import Path
import numpy as np, pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parent))
import imf2_lib as L

ap = argparse.ArgumentParser(); ap.add_argument("tag"); ap.add_argument("--shuffle", type=int, default=None); a = ap.parse_args()
OUT = L.T / "work" / a.tag; OUT.mkdir(parents=True, exist_ok=True); RO = L.RS / a.tag; RO.mkdir(parents=True, exist_ok=True)
doms, aw, _, _, _ = L.domains()
C2 = L.attach_meta(pd.read_parquet(L.RSI / "c2_frame.parquet"), "C2"); C5 = L.attach_meta(pd.read_parquet(L.RSI / "c5_frame.parquet"), "C5")
# registered draw matrices (Mp) are aligned to the FILE order of the frames (score.py / score_c5k7 order), not to sorted stripped ids
REF2 = C2.pair_id.values.copy(); REF5 = C5.pair_id.values.copy(); C2["ord"] = np.arange(len(C2)); C5["ord"] = np.arange(len(C5))
Mp2 = np.load(L.RSI / "c2_Mp.npy"); Mp5 = np.load(L.K7W / "c5_a/c5h_a_Mp.npy")
LP = L.layer_preds_trkA("" if a.shuffle is None else f"_shuf{a.shuffle}")
if a.shuffle is not None:
    rg = np.random.default_rng(a.shuffle)
    C2 = C2.iloc[rg.permutation(len(C2))]; C5 = C5.iloc[rg.permutation(len(C5))]; LP = LP.iloc[rg.permutation(len(LP))]
# order safety: the registered draw matrices are aligned to pair ids sorted as strings
C2 = C2.sort_values("ord").reset_index(drop=True); C5 = C5.sort_values("ord").reset_index(drop=True); LP = LP.sort_index()
assert (C2.pair_id.values == REF2).all() and (C5.pair_id.values == REF5).all(), "frame order differs from the draw-matrix order"
C2["cut_obs"] = (C2.frame_imf == "U_cut_overlap").astype(float); C5["cut_obs"] = 1.0
assert C2.fold.notna().all() and C5.fold.notna().all()
LAYERS = ["scl", "scf", "hcm"]; LAYERS4 = ["scl", "scf", "hcm", "hcm2"]; EXTRA = ["hcm2", "sclc"]; UNCUT = ["L1", "L2b", "L3"]
E = L.ENGS


def rank(X, q):
    """engines ranked by |ln(sum model / sum obs)| (q in g, m) or mean over g, m, n ('all'); >= 90 percent coverage"""
    res = []
    for e in E:
        if q == "all":
            ok = np.isfinite(X[f"{e}_g"].values) & np.isfinite(X[f"{e}_m"].values) & np.isfinite(X[f"{e}_n"].values)
        else:
            ok = np.isfinite(X[f"{e}_{q}"].values)
        if len(X) == 0 or ok.mean() < 0.9:
            continue
        qs = ("g", "m", "n") if q == "all" else (q,)
        err = []
        for qq in qs:
            o = X[f"obs_{qq}"].values[ok]; mm = X[f"{e}_{qq}"].values[ok]; fo = np.isfinite(o)
            err.append(abs(np.log(max(mm[fo].sum(), 1e-12) / o[fo].sum())) if o[fo].sum() > 0 else np.inf)
        res.append((float(np.mean(err)), e))
    return [e for _, e in sorted(res)]


def fit_links(tr2, tr5):
    """per lane, cell and component: ranked engine lists (training pairs only)"""
    links = {}
    for lane in UNCUT + ["L2a"]:
        X = tr5 if lane == "L2a" else tr2[(tr2.lane_u == lane) & (tr2.cut_obs == 0)]
        X = X[np.isfinite(X.obs_g) & np.isfinite(X.obs_m) & np.isfinite(X.obs_n)]
        for q in ("g", "m", "all"):
            links[(lane, "lane", q)] = rank(X, q)
            for d in doms:
                Xd = X[X.dom == d]
                if len(Xd) >= L.MIN_UNITS:
                    links[(lane, d, q)] = rank(Xd, q)
    return links


def pick(F, idx, ranking, q):
    v = np.full(len(idx), np.nan); who = np.array([""] * len(idx), dtype=object)
    for e in ranking:
        x = F[f"{e}_{q}"].values[idx]; f = np.isnan(v) & np.isfinite(x)
        v[f] = x[f]; who[f] = e
    return v, who


def path(F, idx, lanes, links, linked):
    g = np.full(len(idx), np.nan); m = g.copy(); n = g.copy(); wg = np.array([""] * len(idx), dtype=object); wm = wg.copy()
    for lane in np.unique(lanes):
        for d in doms:
            sel = (lanes == lane) & (F.dom.values[idx] == d)
            if not sel.any():
                continue
            ii = idx[sel]
            if linked:
                rg_ = links.get((lane, d, "g"), links[(lane, "lane", "g")]); rm_ = links.get((lane, d, "m"), links[(lane, "lane", "m")])
                gg, a1 = pick(F, ii, rg_, "g"); mm, a2 = pick(F, ii, rm_, "m")
                g[sel], m[sel], n[sel] = gg, mm, gg - mm; wg[sel], wm[sel] = a1, a2
            else:
                ra = links.get((lane, d, "all"), links[(lane, "lane", "all")])
                vv = {q: np.full(len(ii), np.nan) for q in ("g", "m", "n")}; who = np.array([""] * len(ii), dtype=object)
                for e in ra:
                    ok = np.isnan(vv["n"]) & np.isfinite(F[f"{e}_g"].values[ii]) & np.isfinite(F[f"{e}_m"].values[ii]) & np.isfinite(F[f"{e}_n"].values[ii])
                    for q in vv:
                        vv[q][ok] = F[f"{e}_{q}"].values[ii][ok]
                    who[ok] = e
                g[sel], m[sel], n[sel] = vv["g"], vv["m"], vv["n"]; wg[sel] = who; wm[sel] = who
    return g, m, n, wg, wm


N0 = 30.0


def int_class(rem):
    rem = np.asarray(rem, float)
    return np.where(rem < 0.25, "r<25", np.where(rem < 0.5, "r25-50", "r>=50"))


def kernel_fit(tr5):
    X = tr5[np.isfinite(tr5.B1.values) & np.isfinite(tr5.obs_g.values) & np.isfinite(tr5.obs_m.values)]
    tot = {q: X[f"obs_{q}"].sum() / X.B1.sum() for q in ("g", "m")}; est = {("_", "_"): tot}
    ic = int_class(1.0 - X.res_share.fillna(1.0).values)
    for d in sorted(X.dom.unique()):
        Xd = X[X.dom == d]; nd = len(Xd)
        rd = {q: (nd * (Xd[f"obs_{q}"].sum() / Xd.B1.sum()) + N0 * tot[q]) / (nd + N0) for q in ("g", "m")}; est[(d, "_")] = rd
        icd = ic[X.dom.values == d]
        for k in sorted(np.unique(icd)):
            Xk = Xd[icd == k]; nk = len(Xk)
            est[(d, k)] = {q: (nk * (Xk[f"obs_{q}"].sum() / Xk.B1.sum()) + N0 * rd[q]) / (nk + N0) for q in ("g", "m")}
    return est


def kernel_pred(F, idx, est, rem):
    ic = int_class(rem); g = np.full(len(idx), np.nan); m = g.copy()
    for j, (d, k, b) in enumerate(zip(F.dom.values[idx], ic, F.B1.values[idx])):
        r = est.get((d, k)) or est.get((d, "_")) or est[("_", "_")]
        g[j] = r["g"] * b; m[j] = r["m"] * b
    return g, m, g - m


# ----------------------------------------------------------------- out-of-fold build
VARS = {"V1": ("prob", "best", True, "engine"), "Rh": ("prob", "hcm2", True, "engine"), "Rs": ("prob", "sclc", True, "engine"),
        "R4": ("prob", "best4", True, "engine")}
P = {fr: {} for fr in ("C2", "C5")}
paths = {fr: {k: np.full((len(F), 3), np.nan) for k in ("unc_link", "unc_single", "cut_link", "cut_single")} for fr, F in (("C2", C2), ("C5", C5))}
kern = {fr: {k: np.full((len(F), 3), np.nan) for k in ("forecast", "oracle")} for fr, F in (("C2", C2), ("C5", C5))}
lay = {fr: {k: np.full(len(F), np.nan) for k in ("p_best", "s_best", "p_best4", "s_best4", "p_scl", "p_scf", "p_hcm", "s_scl", "s_scf", "s_hcm", "p_hcm2", "s_hcm2", "p_sclc", "s_sclc")} for fr, F in (("C2", C2), ("C5", C5))}
layer_choice, link_log, who_log = [], [], []
for f in range(5):
    tr2 = C2[C2.fold != f]; tr5 = C5[C5.fold != f]; links = fit_links(tr2, tr5)
    for (lane, cell, q), r in links.items():
        link_log.append(dict(fold=f, lane=lane, cell=cell, component=q, first=r[0] if r else "", second=r[1] if len(r) > 1 else ""))
    trL = LP[LP.fold != f]; dev = {k: L.deviance(trL[f"p_{k}"].values, trL.y.values) for k in LAYERS4}
    order = sorted(LAYERS, key=lambda k: (dev[k], k)); order4 = sorted(LAYERS4, key=lambda k: (dev[k], k))
    layer_choice.append(dict(fold=f, **{f"dev_{k}": dev[k] for k in LAYERS4}, chosen=order[0], chosen4=order4[0]))
    est = kernel_fit(tr5)
    for fr, F in (("C2", C2), ("C5", C5)):
        idx = np.where(F.fold.values == f)[0]
        lp = LP.reindex(F.pair_id.values[idx])
        for k in LAYERS + EXTRA:
            lay[fr][f"p_{k}"][idx] = lp[f"p_{k}"].values; lay[fr][f"s_{k}"][idx] = lp[f"s_{k}"].values
        pb = np.full(len(idx), np.nan); sb = pb.copy()
        for k in order:
            fill = np.isnan(pb) & np.isfinite(lp[f"p_{k}"].values); pb[fill] = lp[f"p_{k}"].values[fill]; sb[fill] = lp[f"s_{k}"].values[fill]
        base = float(trL.y.mean()); pb[np.isnan(pb)] = base
        lay[fr]["p_best"][idx] = pb; lay[fr]["s_best"][idx] = sb
        pb4 = np.full(len(idx), np.nan); sb4 = pb4.copy()
        for k in order4:
            fill = np.isnan(pb4) & np.isfinite(lp[f"p_{k}"].values); pb4[fill] = lp[f"p_{k}"].values[fill]; sb4[fill] = lp[f"s_{k}"].values[fill]
        pb4[np.isnan(pb4)] = base; lay[fr]["p_best4"][idx] = pb4; lay[fr]["s_best4"][idx] = sb4
        lanes_u = F.lane_u.values[idx]; lanes_c = np.array(["L2a"] * len(idx))
        for nm, lanes, linked in (("unc_link", lanes_u, True), ("unc_single", lanes_u, False), ("cut_link", lanes_c, True), ("cut_single", lanes_c, False)):
            g, m, n, wg, wm = path(F, idx, lanes, links, linked); paths[fr][nm][idx] = np.column_stack([g, m, n])
            if nm in ("unc_link", "cut_link"):
                who_log += [dict(frame=fr, path=nm, lane=l_, g_engine=x, m_engine=y) for l_, x, y in zip(lanes, wg, wm)]
        sfill = np.where(np.isfinite(sb), sb, 0.0)
        kern[fr]["forecast"][idx] = np.column_stack(kernel_pred(F, idx, est, sfill))
        kern[fr]["oracle"][idx] = np.column_stack(kernel_pred(F, idx, est, 1.0 - F.res_share.fillna(1.0).values[idx]))


def p_of(fr, F, mode, layer):
    if mode == "oracle":
        return F.cut_obs.values.astype(float)
    p = lay[fr]["p_best"] if layer == "best" else lay[fr][f"p_{layer}"]
    if layer not in ("best", "best4"):
        p = np.where(np.isfinite(p), p, lay[fr]["p_best"])
    return p


def build(fr, F, vk, status):
    rule, layer, linked, cutsrc = VARS[vk]
    p = p_of(fr, F, status, layer).copy()
    if status == "forecast" and rule == "hard05":
        p = (p > 0.5).astype(float)
    if status == "forecast" and rule == "topk":
        k = int(round(np.nansum(p))); o = np.lexsort((F.pair_id.values, -p)); q = np.zeros(len(p)); q[o[:k]] = 1.0; p = q
    unc = paths[fr]["unc_link" if linked else "unc_single"]
    cut = kern[fr]["forecast" if status == "forecast" else "oracle"] if cutsrc == "kernel" else paths[fr]["cut_link" if linked else "cut_single"]
    pr = {q: p * cut[:, i] + (1 - p) * unc[:, i] for i, q in enumerate(L.Q)}
    return pr, p


ST = L.stored()
V = pd.read_csv(L.R / "v26_W/release/mpi26w_a3_views_common_basis.csv")
P3 = L.p3_exact()
bo = ST["osm_mc"]
# ALIGNMENT TEST (F15 fix): comparator draws rebuilt from this frame order must equal the registered stored draws
for m_ in ("osm_mc", "E3"):
    p2_ = {q: C2[f"osm_mc_{q}"].values for q in L.Q} if m_ == "osm_mc" else L.mix_pred(C2, L.MIX[m_])
    _, b2_, _ = L.comp_S(p2_, C2, C2.pstar.values & np.isfinite(p2_["n"]), doms, aw, Mp2)
    assert np.max(np.abs(b2_ - ST[m_]["C2"][1])) < 1e-12, ("C2 draws misaligned", m_)
p5_ = {q: C5[f"osm_mc_{q}"].values for q in L.Q}
_, b5_, _ = L.comp_S(p5_, C5, np.isfinite(p5_["n"]), doms, aw, Mp5, qset=("n",), view="statewide")
assert np.max(np.abs(b5_ - ST["osm_mc"]["C5"][1])) < 1e-12, "C5 draws misaligned"
print("alignment test passed")
rows, draws, preds = [], {}, {}
for vk in VARS:
    for status in ("forecast", "oracle"):
        if status == "oracle" and VARS[vk][0] != "prob" and vk in ("V2", "V3"):
            continue
        if status == "oracle" and vk in ("Rh", "Rs", "R4"):
            continue
        pr2, p2 = build("C2", C2, vk, status); pr5, p5 = build("C5", C5, vk, status)
        cov2 = C2.pstar.values; assert np.isfinite(pr2["n"][cov2]).all(), (vk, status, int((~np.isfinite(pr2["n"][cov2])).sum()))
        cov5 = np.isfinite(pr5["n"])
        s2, b2, r2 = L.comp_S(pr2, C2, cov2, doms, aw, Mp2)
        s5, b5, r5 = L.comp_S(pr5, C5, cov5, doms, aw, Mp5, qset=("n",), view="statewide")
        s5t, b5t, r5t = L.comp_S(pr5, C5, cov5 & np.isfinite(pr5["g"]) & np.isfinite(pr5["m"]), doms, aw, Mp5)
        pil, mp, mb = L.mpi_from(bo["C1"], (s2, b2), bo["C3"], P3["osm_mc"], bo["C4"], (s5, b5))
        name = f"IMF2_{vk}_{status}"; draws[name] = mb
        lo, hi = np.percentile(mb, [2.5, 97.5])
        rows.append(dict(method=name, variant=vk, status=status, MPI=mp, MPI_lo=lo, MPI_hi=hi, **pil, S_C2=s2, S_C5net=s5, S_C5three_dom=s5t,
                         C2_pairs=int(cov2.sum()), C5_pairs=int(cov5.sum()), C2_ratio_g=r2["g"], C2_ratio_m=r2["m"], C2_ratio_n=r2["n"],
                         C5_ratio_g=r5t["g"], C5_ratio_m=r5t["m"], C5_ratio_n=r5t["n"], mean_p_C2=float(np.mean(p2)), mean_p_C5=float(np.mean(p5)),
                         inherited="C1, C3 (P2), P3, C4bfm (P5) from OSM mortcal"))
        preds[name] = (pr2, pr5, p2, p5)
for m in ("osm_mc", "E3", "E8f"):
    st = ST[m]; pil, mp, mb = L.mpi_from(st["C1"], st["C2"], st["C3"], P3[m], st["C4"], st["C5"]); draws[m] = mb
    lo, hi = np.percentile(mb, [2.5, 97.5])
    if m == "osm_mc":
        p2_ = {q: C2[f"osm_mc_{q}"].values for q in L.Q}; p5_ = {q: C5[f"osm_mc_{q}"].values for q in L.Q}
    else:
        p2_ = L.mix_pred(C2, L.MIX[m]); p5_ = L.mix_pred(C5, L.MIX[m])
    _, _, r2 = L.comp_S(p2_, C2, C2.pstar.values, doms, aw, Mp2)
    if m == "osm_mc":
        c5c = np.isfinite(p5_["g"]) & np.isfinite(p5_["m"]) & np.isfinite(p5_["n"])
    else:   # pairs every member covers (all three quantities), as the release text states
        c5c = np.logical_and.reduce([np.isfinite(C5[f"{e}_{q}"].values) for e in L.MIX[m] for q in L.Q])
    s5t, _, r5t = L.comp_S(p5_, C5, c5c, doms, aw, Mp5)
    rows.append(dict(method=m, variant="comparator", status="", MPI=mp, MPI_lo=lo, MPI_hi=hi, **pil, S_C2=st["C2"][0], S_C5net=st["C5"][0], S_C5three_dom=s5t,
                     C2_pairs=int(C2.pstar.sum()), C5_pairs=int(c5c.sum()), C2_ratio_g=r2["g"], C2_ratio_m=r2["m"], C2_ratio_n=r2["n"],
                     C5_ratio_g=r5t["g"], C5_ratio_m=r5t["m"], C5_ratio_n=r5t["n"], inherited="registered W a3 components; C5 three and ratios pair-mean (mix = per pair mean)"))
R_ = pd.DataFrame(rows)
# P(rank 1) and paired differences
pr_rows = []
for head in [r for r in R_.method if r.startswith("IMF2_")]:
    set_ = [head, "osm_mc", "E3", "E8f"]; X = np.stack([draws[k] for k in set_], 1)
    tie = np.abs(X - X.max(1, keepdims=True)) <= 1e-12; share = tie / tie.sum(1, keepdims=True)
    d = {f"Prank1_{k}": float(share[:, i].mean()) for i, k in enumerate(set_[1:], 1)}
    dd = {}
    for k in set_[1:]:
        diff = draws[head] - draws[k]; lo, hi = np.percentile(diff, [2.5, 97.5])
        dd.update({f"dMPI_vs_{k}": float(R_.set_index("method").MPI[head] - R_.set_index("method").MPI[k]), f"dMPI_vs_{k}_lo": lo, f"dMPI_vs_{k}_hi": hi,
                   f"win_vs_{k}": bool(lo > 0)})
    pr_rows.append(dict(method=head, Prank1_IMF2=float(share[:, 0].mean()), **d, **dd))
for r_ in pr_rows:   # trkA: paired difference against the released IMF2 V1 forecast (same registered draw index)
    if r_["method"] != "IMF2_V1_forecast" and "IMF2_V1_forecast" in draws:
        diff = draws[r_["method"]] - draws["IMF2_V1_forecast"]; lo, hi = np.percentile(diff, [2.5, 97.5])
        r_.update(dMPI_vs_IMF2V1f=float(R_.set_index("method").MPI[r_["method"]] - R_.set_index("method").MPI["IMF2_V1_forecast"]), dMPI_vs_IMF2V1f_lo=lo, dMPI_vs_IMF2V1f_hi=hi)
PR = pd.DataFrame(pr_rows)
R_.to_csv(OUT / "imf2_scores.csv", index=False, float_format="%.6g"); PR.to_csv(OUT / "imf2_prank.csv", index=False, float_format="%.6g")
pd.DataFrame(layer_choice).to_csv(OUT / "imf2_layer_choice.csv", index=False, float_format="%.6g")
LL = pd.DataFrame(link_log); LL.to_csv(OUT / "imf2_links_by_fold.csv", index=False)
W = pd.DataFrame(who_log); W.groupby(["frame", "path", "lane", "g_engine", "m_engine"]).size().rename("pairs").reset_index().to_csv(OUT / "imf2_waybill_counts.csv", index=False)
np.savez_compressed(RO / "draws.npz", **draws)
pp = {}
for k, (pr2, pr5, p2, p5) in preds.items():
    for q in L.Q:
        pp[f"{k}|C2|{q}"] = pr2[q]; pp[f"{k}|C5|{q}"] = pr5[q]
    pp[f"{k}|C2|p"] = p2; pp[f"{k}|C5|p"] = p5
np.savez_compressed(RO / "preds.npz", c2_ids=C2.pair_id.values.astype(str), c5_ids=C5.pair_id.values.astype(str), **pp)
# frozen Maine links (all folds) for transfer
fz = fit_links(C2, C5); json.dump({"|".join(k): v for k, v in fz.items()}, open(OUT / "imf2_links_frozen.json", "w"), indent=0)
# G4 hull check (V1 forecast and oracle): statewide gross and mortality ratios inside candidate engine hull per frame
hull = {}
for fr, F, cov in (("C2", C2, C2.pstar.values), ("C5", C5, np.ones(len(C5), bool))):
    for q in ("g", "m"):
        rr = []
        for e in E:
            v = F[f"{e}_{q}"].values; ok = cov & np.isfinite(v) & np.isfinite(F[f"obs_{q}"].values)
            if ok.mean() > 0.5 * cov.mean():
                rr.append(v[ok].mean() / F[f"obs_{q}"].values[ok].mean())
        hull[f"{fr}_{q}"] = (min(rr), max(rr))
g4 = {}
for st_ in [m_ for m_ in R_.method if m_.startswith("IMF2_")]:
    r = R_.set_index("method").loc[st_]
    for fr in ("C2", "C5"):
        for q in ("g", "m"):
            lo, hi = hull[f"{fr}_{q}"]; v = r[f"{fr}_ratio_{q}"]; g4[f"{st_}_{fr}_{q}"] = dict(value=float(v), hull_lo=lo, hull_hi=hi, inside=bool(lo <= v <= hi))
json.dump(dict(prereg_md5="35055a79af9083edec3f5e3223552dd7", shuffle=a.shuffle, G4=g4, G4_pass=all(x["inside"] for x in g4.values()),
               n=dict(C2=len(C2), Pstar=int(C2.pstar.sum()), C5=len(C5)), layer_choice=layer_choice), open(OUT / "imf2_meta.json", "w"), indent=1, default=float)
print(R_[["method", "MPI", "MPI_lo", "MPI_hi", "P1", "C2pts", "C5pts", "S_C2", "S_C5net", "C2_ratio_n", "C5_ratio_g", "C5_ratio_m", "C5_ratio_n"]].round(4).to_string())
print(PR.round(3).to_string()); print(pd.DataFrame(layer_choice).to_string())
sys.stdout.flush(); os._exit(0)
