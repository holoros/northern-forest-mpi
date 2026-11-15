"""Track trkB_engfix C3 (P2) from H1 hindcast chain values by the score.py C3 arithmetic (functions imported unchanged: h1_setup,
h1_frame, BD, dom_matrix, read_inv, plot_group, fia_current, fia_grm), domains and area weights as score.py main, draw stream ss[2].
S = area-weighted mean over domains of |ln(model/FIA)|, then mean over horizon bins h10, h20 (score.py wmean). Gate G1-C3 on the
registered values. Writes work/c3_<label>.npz (pt, bt) and work/c3_check.json. usage: b_c3.py label=VALUES_DIR [...]"""
import json, sys
from pathlib import Path
import numpy as np, pandas as pd, yaml
ROOT = Path("/home/aaron/perseus_me/code"); sys.path.insert(0, str(ROOT / "tools/score"))
import score as SC  # noqa: E402
T = Path("/home/aaron/perseus_me/tracks/trkB_engfix")
cfg = yaml.safe_load((ROOT / "config/score/score_v1.yaml").read_text()); R = Path(cfg["run_root"]); B = int(cfg["bootstrap"]["draws"])
seed = int(cfg["bootstrap"]["seed"]); E = cfg["evalid"]; inv = R / "inventory"
ft = pd.read_csv(cfg["ref_forest_type"]["path"]); ft2grp = dict(zip(ft.VALUE.astype(int), ft.TYPGRPCD.astype(int)))
cur = SC.fia_current(inv, E["current"]); grm = SC.fia_grm(inv, E["grm"])
grp = SC.plot_group(cur["CD"], ft2grp); nrem = grp[grp.index.isin(grm["remeasured"])].value_counts()
kept = sorted(int(g) for g, n in nrem.items() if n >= cfg["domains"]["min_remeasured_plots"] and g != 999)
domains = [str(g) for g in kept] + ["other"]; dom = grp.map(lambda g: str(g) if g in kept else "other")
ss = np.random.SeedSequence(seed).spawn(4)
bd1 = SC.BD(inv, E["current"], B, np.random.default_rng(ss[0])); K = len(domains)
A_pt = bd1.point(SC.dom_matrix(bd1, cur["area"], dom, domains)); aw = A_pt[:K] / A_pt[:K].sum()
Hd = R / "hindcast/restricted"; H = SC.h1_setup(ROOT / "config/hindcast/perseus_me_hindcast.yaml", Hd)
_, CDh, _ = SC.read_inv(H["Bv"]); gh = SC.plot_group(CDh, ft2grp).map(lambda g: str(g) if g in kept else "other")


def c3(values_dir):
    bdc = SC.BD(H["Bv"], H["ev"], B, np.random.default_rng(ss[2]))
    pt, bt, rat, nf = [], [], {}, {}
    for hb in cfg["quantities"]["C3"]:
        fr, nfail = SC.h1_frame(H, values_dir, "fvs", "default", hb); nf[hb] = nfail
        fr["dom"] = fr.chain_id.map(H["v1plt"]).astype(str).map(gh).fillna("other"); ud = dict(zip(fr.unit, fr.dom))
        Yy = SC.dom_matrix(bdc, fr.groupby("unit").y.sum(), ud, domains); Yx = SC.dom_matrix(bdc, fr.groupby("unit").x.sum(), ud, domains)
        nch = np.array([(fr.dom == k).sum() for k in domains] + [len(fr)])
        r_pt = bdc.point(Yy) / bdc.point(Yx); r_bt = bdc.boot(Yy) / bdc.boot(Yx)
        e_pt = SC.ln(r_pt, np.ones_like(r_pt)) if hasattr(SC, "ln") else np.log(r_pt); e_bt = SC.ln(r_bt, np.ones_like(r_bt)) if hasattr(SC, "ln") else np.log(r_bt)
        e_pt[:K][nch[:K] < SC.MIN_UNITS] = np.nan; e_bt[:, :K][:, nch[:K] < SC.MIN_UNITS] = np.nan
        pt.append(np.abs(e_pt[:K])); bt.append(np.abs(e_bt[:, :K])); rat[hb] = float(r_pt[K])
    P = np.array(pt)[None]; Bt = np.stack(bt, 1)

    def wmean(x):
        w = np.broadcast_to(aw, x.shape); v = np.isfinite(x)
        return np.nanmean(np.nansum(np.where(v, x * w, 0), -1) / np.where(v, w, 0).sum(-1), -1)
    return float(wmean(P)[0]), wmean(Bt), rat, nf


z = np.load(R / "v26_W/rc_a3/ovl/me_state/v25_K6/work/s1/mpi23_a.npz"); d = np.load(R / "v26_W/rc_a3/ovl/me_state/v25_K6/work/s1/draws_a.npz")
k = "S|fvs|default|C3|primary"; Zr = z if k + "|pt" in z.files else d
ref_pt, ref_bt = float(Zr[k + "|pt"][0]), np.asarray(Zr[k + "|bt"], float)
chk = {}
s, b, rat, nf = c3(Hd / "values")
chk["registered_values"] = dict(S=s, S_wa3=ref_pt, abs_diff=abs(s - ref_pt), max_abs_draw_diff=float(np.max(np.abs(b - ref_bt))), statewide_ratio=rat, failed=nf)
chk["G1_C3_pass"] = bool(abs(s - ref_pt) < 1e-6 and np.max(np.abs(b - ref_bt)) < 1e-6)
for spec in sys.argv[1:]:
    lab, vd = spec.split("=", 1); s, b, rat, nf = c3(Path(vd))
    np.savez(T / f"work/c3_{lab}.npz", pt=s, bt=b); chk[lab] = dict(S=s, lo=float(np.percentile(b, 2.5)), hi=float(np.percentile(b, 97.5)), statewide_ratio=rat, failed=nf)
json.dump(chk, open(T / "work/c3_check.json", "w"), indent=1, default=float); print(json.dumps(chk, indent=1, default=float))
