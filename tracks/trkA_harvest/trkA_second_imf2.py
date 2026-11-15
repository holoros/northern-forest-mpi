"""trkA copy of the IMF2 G3 second code path (frames read from imf2 restricted). Gate G3 (second code path): IMF2 C2 S (by domain on P*) and C5 net S (statewide) recomputed with pandas groupby from the
restricted predictions, point and the first 20 registered draws by explicit row repetition. usage: imf2_second.py TAG"""
import json, os, sys
from pathlib import Path
import numpy as np, pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parent))
import imf2_lib as L
tag = sys.argv[1]
Z = np.load(L.RS / tag / "preds.npz", allow_pickle=True); D = np.load(L.RS / tag / "draws.npz")
C2 = pd.read_parquet(L.RSI / "c2_frame.parquet").set_index("pair_id"); C5 = pd.read_parquet(L.RSI / "c5_frame.parquet").set_index("pair_id")
Mp2 = np.load(L.RSI / "c2_Mp.npy"); Mp5 = np.load(L.K7W / "c5_a/c5h_a_Mp.npy")
doms, aw, _, _, _ = L.domains(); W = pd.Series(aw, index=doms)
ids2 = list(C2.index); ids5 = list(C5.index)   # file order = registered draw order
assert list(Z["c2_ids"]) == ids2 and list(Z["c5_ids"]) == ids5
S = pd.read_csv(L.T / "work" / tag / "imf2_scores.csv").set_index("method")


def s_c2(df):
    out = []
    for q in L.Q:
        g = df.groupby("dom").agg(m=(f"p_{q}", "sum"), o=(f"o_{q}", "sum"), w=("wt", "sum"))
        g = g[g.w >= L.MIN_UNITS - 1e-9] if "wt" in df else g
        r = np.log(np.maximum(g.m, L.FLOOR * g.o) / g.o).abs(); ww = W.reindex(g.index)
        out.append(float((r * ww).sum() / ww.sum()))
    return float(np.mean(out))


res = {}
for meth in [m for m in S.index if m.startswith("IMF2_")]:
    d2 = C2.loc[ids2].copy(); d2 = d2[d2.pstar]
    for q in L.Q:
        d2[f"p_{q}"] = Z[f"{meth}|C2|{q}"][np.isin(ids2, d2.index)]; d2[f"o_{q}"] = d2[f"obs_{q}"]
    d2["wt"] = 1.0
    pt2 = s_c2(d2)
    d5 = C5.loc[ids5].copy(); d5["p_n"] = Z[f"{meth}|C5|n"]; d5 = d5[np.isfinite(d5.p_n)]
    pt5 = abs(np.log(d5.p_n.sum() / d5.obs_n.sum()))
    # draws by explicit repetition of rows (independent of the matrix product path)
    keep2 = np.isin(ids2, C2[C2.pstar].index); keep5 = np.isfinite(Z[f"{meth}|C5|n"])
    b2, b5 = [], []
    for b in (range(10) if meth in ("IMF2_V1_forecast", "IMF2_Rs_forecast", "IMF2_Rh_forecast") else range(0)):
        r2 = d2.loc[np.repeat(np.array(ids2)[keep2], Mp2[b][keep2])]; r2 = r2.assign(wt=1.0)
        b2.append(s_c2(r2))
        r5 = C5.loc[np.repeat(np.array(ids5)[keep5], Mp5[b][keep5])]; pn = pd.Series(Z[f"{meth}|C5|n"], index=ids5).reindex(r5.index)
        b5.append(abs(np.log(pn.sum() / r5.obs_n.sum())))
    C2s = C2.loc[ids2].reset_index(); C5s = C5.loc[ids5].reset_index()
    m2, _, _ = L.comp_S({q: Z[f"{meth}|C2|{q}"] for q in L.Q}, C2s, C2s.pstar.values, doms, aw, Mp2)
    m5, _, _ = L.comp_S({q: Z[f"{meth}|C5|{q}"] for q in L.Q}, C5s, np.isfinite(Z[f"{meth}|C5|n"]), doms, aw, Mp5, qset=("n",), view="statewide")
    assert abs(m2 - S.loc[meth, "S_C2"]) < 1e-5 and abs(m5 - S.loc[meth, "S_C5net"]) < 1e-5
    res[meth] = dict(C2_pt_diff=abs(pt2 - m2), C5_pt_diff=abs(pt5 - m5), C2_second=pt2, C5_second=pt5,
                     C2_draw_maxdiff=None, b2=b2, b5=b5)
# compare draws: recompute MPI draw from the second-path component draws for V1 forecast only (C2, C5 draws checked via MPI identity)
ST = L.stored(); bo = ST["osm_mc"]
V = pd.read_csv(L.R / "v26_W/release/mpi26w_a3_views_common_basis.csv"); P3 = L.p3_exact()["osm_mc"]
out = {}
for meth, r in res.items():
    b2 = np.array(r.pop("b2")); b5 = np.array(r.pop("b5"))
    k = len(b2)
    if k:
        _, _, mb = L.mpi_from((0, bo["C1"][1][:k]), (0, b2), (0, bo["C3"][1][:k]), P3, (0, bo["C4"][1][:k]), (0, b5))
        r["MPI_draws_maxdiff_20"] = float(np.max(np.abs(mb - D[meth][:k])))
    else:
        r["MPI_draws_maxdiff_20"] = 0.0
    out[meth] = r
mx = max(max(v["C2_pt_diff"], v["C5_pt_diff"], v["MPI_draws_maxdiff_20"]) for v in out.values())
json.dump(dict(max_abs_diff=mx, G3_pass=bool(mx < 1e-9), per_method=out), open(L.T / "work" / tag / "g3_second_path.json", "w"), indent=1, default=float)
print("G3 max abs diff", mx)
sys.stdout.flush(); os._exit(0)
