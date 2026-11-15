"""Track trkB_engfix scoring (PREREG sections 2, 3): FVS-NE default (W a3), default_fb (Cardinal rerun, same binary as ddm) and
ddm on the shared draws. Fresh pillars: C2 and C5 net from pair runs; P3 from work/p3.json; C3 from work/c3_<ver>.npz when
present (else W a3). Joins by explicit pair id (asserted). Writes work/<tag>/{scores,diffs,decomp}.csv, restricted/<tag>/draws.npz.
usage: b_score.py TAG [--shuffle SEED]"""
import argparse, json, os, pickle, sys
import numpy as np, pandas as pd
sys.path.insert(0, "/home/aaron/perseus_me/tracks/refine_engines2/tools")
import e1_lib as E
R = E.R; L = R.L
T = E.T2.parent / "trkB_engfix"; RS = T / "restricted"; W = T / "work"
ap = argparse.ArgumentParser(); ap.add_argument("tag"); ap.add_argument("--shuffle", type=int, default=None); a = ap.parse_args()
OUT = W / a.tag; OUT.mkdir(parents=True, exist_ok=True); RO = RS / a.tag; RO.mkdir(parents=True, exist_ok=True); os.chmod(RO, 0o700)
if a.shuffle is None:
    fc = pickle.load(open(E.RS2 / "frames_cache.pkl", "rb")); C2, C5 = fc["C2"], fc["C5"]
else:
    C2, C5, _ = R.frames(a.shuffle)
    rg = np.random.default_rng(a.shuffle + 1)   # the frames come back aligned to draw order; shuffle the prediction tables instead
doms, aw, _, _, _ = L.domains()
Mp2 = np.load(L.RS / "c2_Mp.npy"); Mp5 = np.load(L.K7W / "c5_a/c5h_a_Mp.npy")
z = np.load(L.S1 / "mpi23_a.npz"); d = np.load(L.S1 / "draws_a.npz"); zd = np.load(L.K7W / "mpi_draws_B1_a2_a.npz"); c5s = np.load(L.K7W / "c5_a/c5h_a.npz")
pts = L.pts


def g(k):
    for Z in (z, d):
        if k + "|pt" in Z.files:
            return float(Z[k + "|pt"][0]), np.asarray(Z[k + "|bt"], float)
    raise KeyError(k)


def preds(c2p, c5p):
    D2 = E.rates(RS / c2p / "pair_rows.csv"); D5 = E.rates(RS / c5p / "pair_rows.csv")
    assert not D2.pid.duplicated().any() and not D5.pid.duplicated().any()
    if a.shuffle is not None:
        D2 = D2.iloc[rg.permutation(len(D2))]; D5 = D5.iloc[rg.permutation(len(D5))]
    x2 = D2.set_index("pid").reindex(C2.pair_id.values); x5 = D5.set_index("pid").reindex(C5.pair_id.values)
    assert len(x2) == len(C2) and len(x5) == len(C5) and (x2.index.values == C2.pair_id.values).all() and (x5.index.values == C5.pair_id.values).all()
    return {q: x2[q].values.astype(float) for q in R.Q}, {q: x5[q].values.astype(float) for q in R.Q}


lab, (pre, key) = "FVS-NE", R.MODELS["FVS-NE"]
st = dict(C1=g(f"S|{key}|C1|primary_sw"), C3=g(f"S|{key}|C3|primary"), C4=g(f"S|{key}|C4bfm|primary_sw"))
s5st = float(c5s[f"S|{key}|C5|a1b_netonly_statewide|pt"][0]); s2st = g(f"S|{key}|C2|a1b_common")
mw = float(zd[f"common basis|c5_net|{key}|pt"][0])
P1s = pts((0.35 * st["C1"][0] + 0.25 * s2st[0]) / 0.6); P2s = pts(st["C3"][0]); P5s = pts(st["C4"][0]); P4s = pts(s5st)
P3wa3 = float((mw - 0.30 * P1s - 0.20 * P2s - 0.15 * P4s - 0.15 * P5s) / 0.20)
P3J = json.loads((W / "p3.json").read_text()) if (W / "p3.json").exists() else {}
vers = [("default", {q: C2[f"{pre}_{q}"].values.astype(float) for q in R.Q}, {q: C5[f"{pre}_{q}"].values.astype(float) for q in R.Q}, P3wa3, None)]
for ver, c2p, c5p in (("default_fb", "A_c2_default", "A_c5_default"), ("ddm", "B_c2_ddm", "B_c5_ddm"), ("ddm_all", "S1_c2_ddm_all", "S1_c5_ddm_all")):
    if (RS / c2p / "pair_rows.csv").exists() and (RS / c5p / "pair_rows.csv").exists():
        p2, p5 = preds(c2p, c5p)
        p3 = P3J.get({"ddm_all": "ddm"}.get(ver, ver), {}).get("P3", np.nan)
        c3f = W / f"c3_{ {'ddm_all': 'ddm'}.get(ver, ver)}.npz"; c3 = None
        if c3f.exists():
            cz = np.load(c3f); c3 = (float(cz["pt"]), np.asarray(cz["bt"], float))
        vers.append((ver, p2, p5, p3, c3))
rows, draws, pil_d, out, SEC = [], {}, {}, {}, []
cov2 = C2.pstar.values; cov5 = np.isfinite(C5[f"{pre}_n"].values)
for ver, p2, p5, p3, c3 in vers:
    s2, bt2, rat2 = L.comp_S(p2, C2, cov2, doms, aw, Mp2)
    s5, bt5, _ = L.comp_S(p5, C5, cov5, doms, aw, Mp5, qset=("n",), view="statewide")
    _, _, rat5 = L.comp_S(p5, C5, cov5 & np.isfinite(p5["g"]) & np.isfinite(p5["m"]), doms, aw, Mp5, view="statewide")
    C3v = c3 if c3 is not None else st["C3"]
    P1 = (pts((0.35 * st["C1"][0] + 0.25 * s2) / 0.6), pts((0.35 * st["C1"][1] + 0.25 * bt2) / 0.6))
    P2 = (pts(C3v[0]), pts(C3v[1])); P5 = (pts(st["C4"][0]), pts(st["C4"][1])); P4 = (pts(s5), pts(bt5))
    w3 = 0.0 if np.isnan(p3) else 0.20; p3v = 0.0 if np.isnan(p3) else p3; den = 0.80 + w3
    mp = (0.30 * P1[0] + 0.20 * P2[0] + w3 * p3v + 0.15 * P4[0] + 0.15 * P5[0]) / den
    mb = (0.30 * P1[1] + 0.20 * P2[1] + w3 * p3v + 0.15 * P4[1] + 0.15 * P5[1]) / den
    lo, hi = np.percentile(mb, [2.5, 97.5]); grade = "A" if mp >= 85 else "B" if mp >= 70 else "C" if mp >= 55 else "D" if mp >= 40 else "F"
    fresh = "none (W a3)" if ver == "default" else "C2, C5" + (", P3" if np.isfinite(p3) else "") + (", C3" if c3 is not None else "")
    rows.append(dict(model=lab, version=ver, overall=mp, lo=lo, hi=hi, grade=grade, P1=float(P1[0]), P2=float(P2[0]), P3=p3, C5pts=float(P4[0]),
                     P5=float(P5[0]), C2pts=float(pts(s2)), S_C2=s2, S_C5net=s5, S_C3=C3v[0], uncut_gross_ratio=rat2["g"], uncut_mort_ratio=rat2["m"],
                     uncut_net_ratio=rat2["n"], cut_gross_ratio=rat5["g"], cut_mort_ratio=rat5["m"], cut_net_ratio=rat5["n"],
                     n_C2_missing=int((cov2 & ~np.isfinite(p2["n"])).sum()), n_C5_missing=int((cov5 & ~np.isfinite(p5["n"])).sum()),
                     fresh_pillars=fresh, W_a3_overall=mw))
    draws[ver] = mb; out[ver] = mp
    pil_d[ver] = dict(P1=(0.30 / den, P1), P2=(0.20 / den, P2), P3=(w3 / den, (p3v, np.full_like(mb, p3v))), C5=(0.15 / den, P4), P5=(0.15 / den, P5))
    if ver in ("default_fb", "ddm"):
        for nm, F, pp, cv in (("c2", C2, p2, cov2), ("c5", C5, p5, cov5)):
            SEC.append(pd.DataFrame(dict(version=ver, set=nm, pair_id=F.pair_id.values, dom=F.dom.values.astype(str), cov=cv, obs_g=F.obs_g.values,
                                         obs_m=F.obs_m.values, obs_n=F.obs_n.values, pred_g=pp["g"], pred_m=pp["m"], pred_n=pp["n"])))
diffs, dec = [], []
for v1, v0 in (("ddm", "default"), ("ddm", "default_fb"), ("default_fb", "default"), ("ddm_all", "ddm")):
    if v1 in out and v0 in out:
        dd = draws[v1] - draws[v0]; lo, hi = np.percentile(dd, [2.5, 97.5])
        diffs.append(dict(contrast=f"{v1} - {v0}", d_overall=out[v1] - out[v0], d_lo=lo, d_hi=hi, p_better=float((dd > 0).mean()),
                          p_worse=float((dd < 0).mean()), lower_limit_above_0=bool(lo > 0)))
        for p in ("P1", "P2", "P3", "C5", "P5"):
            w1, (a1, b1) = pil_d[v1][p]; w0, (a0, b0) = pil_d[v0][p]
            cd = w1 * np.asarray(b1) - w0 * np.asarray(b0); l2, h2 = np.percentile(cd, [2.5, 97.5])
            dec.append(dict(contrast=f"{v1} - {v0}", pillar=p, points_v1=float(a1), points_v0=float(a0), contribution=float(w1 * a1 - w0 * a0), c_lo=l2, c_hi=h2))
S = pd.DataFrame(rows); D = pd.DataFrame(diffs); DC = pd.DataFrame(dec)
S.to_csv(OUT / "scores.csv", index=False, float_format="%.10g"); D.to_csv(OUT / "diffs.csv", index=False, float_format="%.10g"); DC.to_csv(OUT / "decomp.csv", index=False, float_format="%.10g")
np.savez_compressed(RO / "draws.npz", **draws)
if SEC:
    pd.concat(SEC).to_csv(RO / "second_in.csv", index=False, float_format="%.17g")
    json.dump(dict(doms=[str(x) for x in doms], aw=[float(x) for x in aw]), open(RO / "second_doms.json", "w"))
reg = float(np.max(np.abs(draws["default"] - np.asarray(zd[f"common basis|c5_net|{key}|bt"], float))))
json.dump(dict(alignment_regression_default_draws_maxabs=reg, shuffle=a.shuffle, P3_wa3=P3wa3), open(OUT / "meta.json", "w"), indent=1)
pd.set_option("display.width", 250)
print("alignment", reg); print(S[["version", "overall", "lo", "hi", "grade", "P1", "P2", "P3", "C2pts", "C5pts", "uncut_mort_ratio", "cut_mort_ratio", "cut_net_ratio"]].round(3).to_string())
print(D.round(3).to_string()); print(DC.round(3).to_string()); sys.stdout.flush(); os._exit(0)
