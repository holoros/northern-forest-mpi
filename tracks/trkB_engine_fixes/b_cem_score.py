"""round4 section 3: CEM default vs refined. C5 net (statewide) rerun from pair_rows by the registered comp_S arithmetic on the C5 frame and
the K7 draw-index matrix; C1, C2, C3, C4bfm and P3 inherited from W a3 cem|engine_seeded_v2 (P3 back-solved from the W a3 MPI as e_score).
Predictions joined to the frame by explicit pair id (asserted). Gate: the default pair_rows reproduce the W a3 C5 S and MPI.
usage: cem_score.py REFINED_PAIR_ROWS [--shuffle SEED]"""
import argparse, json, pickle, sys
import numpy as np, pandas as pd
sys.path.insert(0, "/home/aaron/perseus_me/tracks/refine_engines2/tools")
import e1_lib as E
R = E.R; L = R.L
ap = argparse.ArgumentParser(); ap.add_argument("refined"); ap.add_argument("--shuffle", type=int, default=None); a = ap.parse_args()
if a.shuffle is None:
    fc = pickle.load(open(E.RS2 / "frames_cache.pkl", "rb")); C5 = fc["C5"]
else:
    _, C5, _ = R.frames(a.shuffle)
doms, aw, _, _, _ = L.domains(); Mp5 = np.load(L.K7W / "c5_a/c5h_a_Mp.npy")
z = np.load(L.S1 / "mpi23_a.npz"); d = np.load(L.S1 / "draws_a.npz"); zd = np.load(L.K7W / "mpi_draws_B1_a2_a.npz"); c5s = np.load(L.K7W / "c5_a/c5h_a.npz"); pts = L.pts
def g(k):
    for Z in (z, d):
        if k + "|pt" in Z.files:
            return float(Z[k + "|pt"][0]), np.asarray(Z[k + "|bt"], float)
    raise KeyError(k)
key = "cem|engine_seeded_v2"
st = dict(C1=g(f"S|{key}|C1|primary_sw"), C3=g(f"S|{key}|C3|primary"), C4=g(f"S|{key}|C4bfm|primary_sw")); s2 = g(f"S|{key}|C2|a1b_common")
s5_reg = float(c5s[f"S|{key}|C5|a1b_netonly_statewide|pt"][0]); mw = float(zd[f"common basis|c5_net|{key}|pt"][0]); mw_bt = np.asarray(zd[f"common basis|c5_net|{key}|bt"], float)
P1s = pts((0.35 * st["C1"][0] + 0.25 * s2[0]) / 0.6)
P3 = float((mw - 0.30 * P1s - 0.20 * pts(st["C3"][0]) - 0.15 * pts(s5_reg) - 0.15 * pts(st["C4"][0])) / 0.20)
def score(path):
    sys.path.insert(0, "/home/aaron/perseus_me/code/workflow/scripts"); import gate_pairs as GP
    D = pd.read_csv(path, dtype={"pair_id": str}); D = D[D.ok.astype(str).str.lower().isin(["true", "1"])].copy(); D["pid"] = D.pair_id.str.strip()
    D["g"] = D.fvs_surv / D.cycle * GP.K; D["m"] = D.fvs_mort / D.cycle * GP.K; D["n"] = D.g - D.m   # e1_lib.rates arithmetic (no tag column)
    assert not D.pid.duplicated().any(); X = D.set_index("pid").reindex(C5.pair_id.values)
    assert (X.index.values == C5.pair_id.values).all()
    p5 = {q: X[q].values.astype(float) for q in R.Q}; cov5 = np.isfinite(p5["n"])
    s5, bt5, _ = L.comp_S(p5, C5, cov5, doms, aw, Mp5, qset=("n",), view="statewide")
    _, _, rat = L.comp_S(p5, C5, cov5 & np.isfinite(p5["g"]) & np.isfinite(p5["m"]), doms, aw, Mp5, view="statewide")
    P1 = (pts((0.35 * st["C1"][0] + 0.25 * s2[0]) / 0.6), pts((0.35 * st["C1"][1] + 0.25 * s2[1]) / 0.6))
    mp = (0.30 * P1[0] + 0.20 * pts(st["C3"][0]) + 0.20 * P3 + 0.15 * pts(s5) + 0.15 * pts(st["C4"][0]))
    mb = (0.30 * P1[1] + 0.20 * pts(st["C3"][1]) + 0.20 * P3 + 0.15 * pts(bt5) + 0.15 * pts(st["C4"][1]))
    return dict(S_C5net=s5, C5pts=float(pts(s5)), ratio_g=rat["g"], ratio_m=rat["m"], ratio_n=rat["n"], n_scored=int(cov5.sum()), n_frame=len(C5), mpi=float(mp)), mb
dflt, bd = score("/home/aaron/perseus_me/me_state/v25_C/restricted/values/cem_v2p_d3/pair_rows.csv")
gate = dict(S_C5_reproduced=dflt["S_C5net"], S_C5_registered=s5_reg, mpi_reproduced=dflt["mpi"], mpi_registered=mw, max_abs_draw_diff=float(np.max(np.abs(bd - mw_bt))))
gate["pass"] = bool(abs(dflt["S_C5net"] - s5_reg) < 1e-4 and abs(dflt["mpi"] - mw) < 1e-4)
out = dict(gate=gate, P3_inherited=P3, inherited_share=round(1 - 0.15, 4), default=dict(dflt, lo=float(np.percentile(bd, 2.5)), hi=float(np.percentile(bd, 97.5))))
if a.refined != "none":
    rf, rb = score(a.refined); dd = rb - bd
    out["refined"] = dict(rf, lo=float(np.percentile(rb, 2.5)), hi=float(np.percentile(rb, 97.5)))
    out["paired"] = dict(diff=rf["mpi"] - dflt["mpi"], lo=float(np.percentile(dd, 2.5)), hi=float(np.percentile(dd, 97.5)), p_better=float((dd > 0).mean()))
    if a.shuffle is None:
        np.savez("/home/aaron/perseus_me/tracks/trkB_engfix/restricted/cem_draws_v2s.npz", default=bd, refined=rb)
tag = "" if a.shuffle is None else f"_shuf{a.shuffle}"
json.dump(out, open(f"/home/aaron/perseus_me/tracks/trkB_engfix/work/cem_score_v2s{tag}.json", "w"), indent=1, default=float); print(json.dumps(out, indent=1, default=float))
