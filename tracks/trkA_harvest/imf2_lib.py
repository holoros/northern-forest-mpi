# trkA_hcm copy of tracks/imf2/tools (unedited except input paths: imf2 restricted inputs read only, outputs under trkA_hcm)
"""Track IMF2 (PREREG.md md5 b539cc149e8f1e1606df0b1d46ccc065): loaders, routing and MPI component scoring.
Restricted pair values stay in tracks/imf2/restricted. Aggregates only leave."""
import json, sys
from pathlib import Path
import numpy as np
import pandas as pd

PME = Path("/home/aaron/perseus_me"); R = PME / "me_state"; T = PME / "tracks/trkA_hcm"; RS = T / "restricted"; RSI = PME / "tracks/imf2/restricted"  # trkA: imf2 restricted inputs read only
OV = R / "v26_W/rc_a3/ovl/me_state"; S1 = OV / "v25_K6/work/s1"; K7W = OV / "v25_K7/work"
sys.path.insert(0, str(PME / "tracks/v26_IMF/tools/v26_IMF"))
import imf_common as IC  # noqa: E402
import gate_pairs as GP  # noqa: E402  (path inserted by imf_common)
import score_c5 as D5  # noqa: E402

SEED, B = 20261103, 1000
FLOOR, MIN_UNITS = 0.01, 30
Q = ("g", "m", "n")
QL = {"g": "gross_growth", "m": "mortality", "n": "net_growth"}
ENGS = ["osm", "osm_mc", "fvsd", "fvs_gc", "fvsa", "cem", "acd", "cbmg", "cbma", "landis"]
C5RUN = {"osm": ("osm", "segmented"), "osm_mc": ("osm", "segmented_mortcal"), "fvsd": ("fvs", "default"),
         "fvs_gc": ("fvs", "default_growcal"), "fvsa": ("fvs", "acd_default"), "cem": ("cem", "engine_seeded_v2"),
         "acd": ("acd141", "acd_v141"), "cbmg": ("cbm_gated", "FVS_FIAW_NEW"), "cbma": ("cbm_af", "papa2023"),
         "landis": ("landis", "registered_tables")}
FOL = {"cbmg", "cbma", "landis"}
MIX = {"E3": ["osm", "fvsd", "fvsa"], "E8f": ["osm", "fvsd", "fvsa", "cem", "landis", "cbmg", "cbma"]}
MPIKEY = {"osm_mc": "osm|segmented_mortcal", "E3": "ensemble|E3_fvs_osm", "E8f": "ensemble|E8f"}


def pts(s):
    return 100.0 * np.maximum(0.0, 1.0 - (np.exp(s) - 1.0) / 0.5)


def ln(m, f):
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(f > 0, np.log(np.maximum(m, FLOOR * f) / f), np.nan)


def strip(s):
    return pd.Series(s).astype(str).str.strip().values


def score_cfg():
    return OV.parent.parent / "code/D/config/score/score_v1.yaml"


def domains():
    ft2grp, kept, doms, aw = D5.domains_and_weights(OV, score_cfg())
    dmap = lambda x: (lambda g: str(g) if g in kept else "other")(ft2grp.get(int(x), 999) if pd.notna(x) else 999)  # noqa: E731
    return doms, np.asarray(aw, float), dmap, ft2grp, kept


def foliage():
    import score_c5k7 as K7
    doms, aw, dmap, ft2grp, kept = domains()
    phi_pt, phi_bt = K7.phi_draws(OV, score_cfg(), B, SEED, doms, kept, ft2grp)
    return 1.0 - phi_pt, 1.0 - np.asarray(phi_bt)


# ------------------------------------------------------------------ C2 frame
def c2_frame():
    doms, aw, dmap, ft2grp, kept = domains()
    PB = pd.read_csv(OV.parent / "pairs/me_state_std/pairs_built.csv", dtype={"pair_id": str})
    x = pd.read_csv(OV.parent / "pairs/me_state_std/pair_rows.csv", dtype={"pair_id": str}); x = x[x.arm == "default"]
    Tt, _ = GP.pair_table(x, PB); Tt = Tt.sort_values("pair_id").reset_index(drop=True)
    pdom = PB.set_index("pair_id").FORTYPCD.map(dmap)
    base_ids = Tt.pair_id.values; dd = pd.Series(base_ids).map(pdom).values
    ss = np.random.SeedSequence(SEED).spawn(4); rngp = np.random.default_rng(ss[1])
    Mp = np.zeros((B, len(base_ids)), np.int32)
    for k in doms:
        ix = np.where(dd == k)[0]
        if len(ix):
            dr = rngp.integers(0, len(ix), (B, len(ix)))
            for b in range(B):
                Mp[b, ix] += np.bincount(dr[b], minlength=len(ix)).astype(np.int32)
    ids = strip(base_ids)
    A1 = pd.read_parquet(OV / "mpi2/plotlevel/a1/a1_pairs.parquet"); A1["pid"] = strip(A1.pair_id)
    A1["nb"] = np.where(A1.has_grm, A1.grm_net_mgha, A1.obs_net)
    a1g = dict(zip(A1.pid, A1.obs_gg_a1)); a1n = dict(zip(A1.pid, A1.nb))
    F = pd.DataFrame({"pair_id": ids, "dom": dd})
    F["obs_g"] = pd.Series(ids).map(a1g).values.astype(float)
    F["obs_m"] = Tt.obs_mort.values.astype(float)
    F["obs_n"] = pd.Series(ids).map(a1n).values.astype(float)
    M = pd.read_parquet(S1.parent.parent / "restricted/s1/a/c2_E9_all_members.parquet"); M.index = strip(M.pair_id)
    F["pstar"] = M.in_Pstar.reindex(ids).values.astype(bool)
    c_pt, _ = foliage()
    for e in ENGS:
        for q in Q:
            col = f"{e}_{QL[q]}"
            v = M[col].reindex(ids).values.astype(float) if col in M.columns else np.full(len(ids), np.nan)
            F[f"{e}_{q}"] = v * (c_pt if e in FOL else 1.0)
    return F, Mp, doms, aw


# ------------------------------------------------------------------ C5 frame
def c5_frame():
    doms, aw, dmap, ft2grp, kept = domains()
    cfg = json.loads((K7W / "cfg_k7_a.json").read_text())
    P = pd.read_csv(cfg["treated_pairs"], dtype={"pair_id": str, "t1_id": str}).sort_values("pair_id").reset_index(drop=True)
    Mp = np.load(K7W / "c5_a/c5h_a_Mp.npy")
    assert Mp.shape == (B, len(P)), Mp.shape
    F = pd.DataFrame({"pair_id": strip(P.pair_id), "dom": P.FORTYPCD.map(dmap).values})
    oc = D5.OBSV["a1b"]
    F["obs_g"], F["obs_m"], F["obs_n"] = (P[c].values.astype(float) for c in oc)
    c_pt, _ = foliage()
    runs = {(r["engine"], r["arm"]): r for r in cfg["runs"]}
    pos = pd.Series(np.arange(len(P)), index=P.pair_id.values)
    for e, ea in C5RUN.items():
        r = runs[ea]; x = pd.read_csv(r["rows"], dtype={"pair_id": str}); x = x[x.arm == r.get("rows_arm", ea[1])]
        Tt, _ = GP.pair_table(x, P[["pair_id", "B1_Mgha", "obs_gg", "obs_mort", "obs_net", "obs_surv"]])
        m = np.full((len(P), 3), np.nan); j = pos.reindex(Tt.pair_id.values).values.astype(int)
        m[j] = Tt[["fvs_gg", "fvs_mort_r", "fvs_net"]].values.astype(float)
        if r.get("foliage"):
            m = m * c_pt
        for i, q in enumerate(Q):
            F[f"{e}_{q}"] = m[:, i]
    return F, Mp, doms, aw


def attach_meta(F, which):
    """fold, loc and past-observation lane from the track IMF frame and track H pairs_h (restricted)"""
    H = pd.read_parquet(R / "v26_H/restricted/frame/pairs_h.parquet"); H.index = strip(H.CN)
    G = pd.read_parquet(RSI / "frame_imf.parquet"); G = G[G.frame.isin(["C"] if which == "C5" else ["U", "U_cut_overlap"])].set_index("pair_id")
    assert not G.index.duplicated().any()
    h = H.reindex(F.pair_id.values)
    F["fold"] = h.fold.values; F["loc"] = h["loc"].values
    g = G.reindex(F.pair_id.values)
    F["B1"] = g.B1_Mgha.values; F["res_share"] = g.res_share.values; F["frame_imf"] = g.frame.values
    rw = g.prev_cut.fillna(False).values.astype(bool) & (g.ba1_h.values < g.pre_cut_ba.values)
    old = g.STDAGE.values >= 100
    F["lane_u"] = np.select([rw, old], ["L2b", "L3"], "L1")
    return F


def layer_preds():
    """cross-fitted p_cut and s per CN for the three layers + track H y and fold"""
    H = pd.read_parquet(R / "v26_H/restricted/frame/pairs_h.parquet"); H.index = strip(H.CN)
    out = pd.DataFrame({"y": H.cut.astype(float).values, "fold": H.fold.values}, index=H.index)
    sl = pd.read_parquet(R / "v26_SC/restricted/preds_sc_logit.parquet"); sl = sl[sl.split == "cf"]; sl.index = strip(sl.CN)
    sf = pd.read_csv(R / "v26_SC/restricted/preds_sc_rf.csv", dtype={"CN": str}); sf = sf[sf.split == "cf"]; sf.index = strip(sf.CN)
    hc = pd.read_parquet(R / "v26_HCM/restricted/hcm_preds.parquet"); hc = hc[(hc.driver == "hcm_refit") & hc.covered.astype(bool)]; hc.index = strip(hc.CN)
    for k, d in (("scl", sl), ("scf", sf), ("hcm", hc)):
        out[f"p_{k}"] = d.p.reindex(out.index).values.astype(float); out[f"s_{k}"] = d.s.reindex(out.index).values.astype(float)
    return out


def deviance(p, y):
    p = np.clip(p, 1e-6, 1 - 1e-6); ok = np.isfinite(p) & np.isfinite(y)
    return float(-np.mean(y[ok] * np.log(p[ok]) + (1 - y[ok]) * np.log(1 - p[ok])))


# ------------------------------------------------------------------ scoring with registered draws
def dmeans(v, c, dom, doms, Mp):
    Dm = np.column_stack([(dom == k) for k in doms] + [np.ones(len(dom), bool)]).astype(float)
    Dc = Dm * c[:, None]; nc = Dc.sum(0); vq = np.where(c, np.nan_to_num(v), 0.0)
    with np.errstate(invalid="ignore", divide="ignore"):
        return (vq[:, None] * Dc).sum(0) / nc, (Mp @ (vq[:, None] * Dc)) / (Mp @ Dc), nc


def comp_S(pred, F, cov, doms, aw, Mp, qset=Q, view="domain"):
    """MPI component arithmetic: per q ratio of domain means (ln, floor, cells < 30 dropped); S by domain (area weighted,
    renormalized, mean over q) or statewide (pooled cell, mean over q). Returns S_pt, S_bt, statewide ratios."""
    dom = F.dom.values; K = len(doms); Ep, Eb, rat = [], [], {}
    for q in qset:
        mp, mb, nc = dmeans(pred[q], cov, dom, doms, Mp); op, ob, _ = dmeans(F[f"obs_{q}"].values, cov, dom, doms, Mp)
        ep = ln(mp, op); eb = ln(mb, ob)
        ep[:K][nc[:K] < MIN_UNITS] = np.nan; eb[:, :K][:, nc[:K] < MIN_UNITS] = np.nan
        Ep.append(np.abs(ep)); Eb.append(np.abs(eb)); rat[q] = float(np.exp(ep[K])); rat[f"{q}_ln_bt"] = eb[:, K]
    Ep = np.array(Ep); Eb = np.stack(Eb, 1)
    if view == "statewide":
        return float(Ep[:, K].mean()), Eb[:, :, K].mean(1), rat
    w = aw

    def wm(E):
        ww = np.broadcast_to(w, E.shape); v = np.isfinite(E)
        return np.nanmean(np.nansum(np.where(v, E * ww, 0), -1) / np.where(v, ww, 0).sum(-1), -1)
    return float(wm(Ep[None, :, :K])[0]), wm(Eb[:, :, :K]), rat


def mix_pred(F, mem):
    return {q: np.nanmean(np.stack([F[f"{e}_{q}"].values for e in mem], 0), 0) for q in Q}


def stored():
    z = np.load(S1 / "mpi23_a.npz"); d = np.load(S1 / "draws_a.npz"); c5 = np.load(K7W / "c5_a/c5h_a.npz")

    def g(k):
        for Z in (z, d):
            if k + "|pt" in Z.files:
                return float(Z[k + "|pt"][0]), np.asarray(Z[k + "|bt"], float)
        raise KeyError(k)
    out = {}
    for m, key in MPIKEY.items():
        out[m] = dict(C1=g(f"S|{key}|C1|primary_sw"), C2=g(f"S|{key}|C2|a1b_common"), C3=g(f"S|{key}|C3|primary"),
                      C4=g(f"S|{key}|C4bfm|primary_sw"),
                      C5=(float(c5[f"S|{key}|C5|a1b_netonly_statewide|pt"][0]), np.asarray(c5[f"S|{key}|C5|a1b_netonly_statewide|bt"], float)))
    return out


def mpi_from(C1, C2, C3, P3, C4, C5):
    """C* are (pt, bt) S pairs; P3 points (fixed). Returns pillars (pt) and MPI (pt, bt)."""
    P1 = (pts((0.35 * C1[0] + 0.25 * C2[0]) / 0.6), pts((0.35 * C1[1] + 0.25 * C2[1]) / 0.6))
    P2 = (pts(C3[0]), pts(C3[1])); P5 = (pts(C4[0]), pts(C4[1])); P4 = (pts(C5[0]), pts(C5[1]))
    mpt = 0.30 * P1[0] + 0.20 * P2[0] + 0.20 * P3 + 0.15 * P4[0] + 0.15 * P5[0]
    mbt = 0.30 * P1[1] + 0.20 * P2[1] + 0.20 * P3 + 0.15 * P4[1] + 0.15 * P5[1]
    return dict(P1=float(P1[0]), P2=float(P2[0]), P3=float(P3), C5pts=float(P4[0]), P5=float(P5[0]), C2pts=float(pts(C2[0]))), float(mpt), mbt


def p3_exact():
    """P3 points of each comparator recovered exactly from the stored W a3 MPI point (the release CSV is rounded to 4 dp)"""
    st = stored(); zd = np.load(K7W / "mpi_draws_B1_a2_a.npz"); out = {}
    for m, key in MPIKEY.items():
        s = st[m]; mpt = float(zd[f"common basis|c5_net|{key}|pt"][0])
        P1 = pts((0.35 * s["C1"][0] + 0.25 * s["C2"][0]) / 0.6); P2 = pts(s["C3"][0]); P5 = pts(s["C4"][0]); P4 = pts(s["C5"][0])
        out[m] = float((mpt - 0.30 * P1 - 0.20 * P2 - 0.15 * P4 - 0.15 * P5) / 0.20)
    return out


def layer_preds_trkA(shuffle_tag=""):
    """trkA: IMF2 layers plus the trkA county-fold cross-fitted layers hcm2 (HCM2 A0) and sclc (SC-L county), joined by CN"""
    out = layer_preds()
    for k, f in (("hcm2", "A0"), ("sclc", "SCL_cty")):
        d = pd.read_parquet(RS / "layers" / f"{f}{shuffle_tag}.parquet"); d.index = strip(d.CN)
        assert not d.index.duplicated().any()
        out[f"p_{k}"] = d.p.reindex(out.index).values.astype(float); out[f"s_{k}"] = d.s.reindex(out.index).values.astype(float)
    return out
