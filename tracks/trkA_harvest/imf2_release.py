"""Track IMF2 release: aggregates, component intervals, gates G1 to G5, figure data. usage: imf2_release.py"""
import json, os, re, sys
from pathlib import Path
import numpy as np, pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parent))
import imf2_lib as L
W = L.T / "work"; REL = L.T / "release"; REL.mkdir(exist_ok=True)
doms, aw, _, _, _ = L.domains()
S = pd.read_csv(W / "main/imf2_scores.csv"); PR = pd.read_csv(W / "main/imf2_prank.csv"); meta = json.load(open(W / "main/imf2_meta.json"))
Z = np.load(L.RS / "main/preds.npz", allow_pickle=True); D = np.load(L.RS / "main/draws.npz")
C2 = pd.read_parquet(L.RS / "c2_frame.parquet")
C5 = pd.read_parquet(L.RS / "c5_frame.parquet")
Mp2 = np.load(L.RS / "c2_Mp.npy"); Mp5 = np.load(L.K7W / "c5_a/c5h_a_Mp.npy")
assert list(Z["c2_ids"]) == list(C2.pair_id) and list(Z["c5_ids"]) == list(C5.pair_id)
ST = L.stored(); P3 = L.p3_exact()
comp = {}
for m in S.method:
    if m.startswith("IMF2_"):
        p2 = {q: Z[f"{m}|C2|{q}"] for q in L.Q}; p5 = {q: Z[f"{m}|C5|{q}"] for q in L.Q}
        s2, b2, _ = L.comp_S(p2, C2, C2.pstar.values, doms, aw, Mp2)
        s5, b5, _ = L.comp_S(p5, C5, np.isfinite(p5["n"]), doms, aw, Mp5, qset=("n",), view="statewide")
        assert abs(s2 - S.set_index("method").S_C2[m]) < 1e-5 and abs(s5 - S.set_index("method").S_C5net[m]) < 1e-5
    else:
        (s2, b2), (s5, b5) = ST[m]["C2"], ST[m]["C5"]
    comp[m] = dict(C2=(s2, b2), C5=(s5, b5))
rows = []
for m in S.method:
    c = comp[m]; c2p, c2b = L.pts(c["C2"][0]), L.pts(c["C2"][1]); c5p, c5b = L.pts(c["C5"][0]), L.pts(c["C5"][1])
    fci, fcb = (12.5 * c2p + 15 * c5p) / 27.5, (12.5 * c2b + 15 * c5b) / 27.5
    comp[m]["FCI_bt"] = fcb
    rows.append(dict(method=m, C2pts_lo=np.percentile(c2b, 2.5), C2pts_hi=np.percentile(c2b, 97.5), C5pts_lo=np.percentile(c5b, 2.5), C5pts_hi=np.percentile(c5b, 97.5),
                     FCI=fci, FCI_lo=np.percentile(fcb, 2.5), FCI_hi=np.percentile(fcb, 97.5)))
S = S.merge(pd.DataFrame(rows), on="method")
# fresh-component index P(rank 1) (post hoc sensitivity, deviation IMF2-D2)
fr = []
for head in [m for m in S.method if m.startswith("IMF2_")]:
    set_ = [head, "osm_mc", "E3", "E8f"]; X = np.stack([comp[k]["FCI_bt"] for k in set_], 1)
    tie = np.abs(X - X.max(1, keepdims=True)) <= 1e-12; sh = (tie / tie.sum(1, keepdims=True)).mean(0)
    fr.append(dict(method=head, FCI_Prank1_IMF2=sh[0], FCI_Prank1_osm_mc=sh[1], FCI_Prank1_E3=sh[2], FCI_Prank1_E8f=sh[3]))
PR = PR.merge(pd.DataFrame(fr), on="method")
# ---------------- gates
zd = np.load(L.K7W / "mpi_draws_B1_a2_a.npz"); g1 = json.load(open(W / "g1_repro.json")); g1x = {}
for m, key in L.MPIKEY.items():
    st = ST[m]; _, mp, mb = L.mpi_from(st["C1"], st["C2"], st["C3"], P3[m], st["C4"], st["C5"])
    g1x[m] = dict(MPI_pt_absdiff=abs(mp - float(zd[f"common basis|c5_net|{key}|pt"][0])), MPI_bt_maxabsdiff=float(np.max(np.abs(mb - zd[f"common basis|c5_net|{key}|bt"]))),
                  C2_pt_absdiff=abs(g1[m]["C2_pt_ours"] - g1[m]["C2_pt_reg"]), C5_pt_absdiff=abs(g1[m]["C5_pt_ours"] - g1[m]["C5_pt_reg"]),
                  C2_bt_maxdiff_rebuilt=g1[m]["C2_bt_maxdiff"], C5_bt_maxdiff_rebuilt=g1[m]["C5_bt_maxdiff"])
G1 = all(v["MPI_pt_absdiff"] < 1e-6 and v["C2_pt_absdiff"] < 1e-9 and v["C5_pt_absdiff"] < 1e-9 for v in g1x.values())
s7 = pd.read_csv(W / "shuf7/imf2_scores.csv"); p7 = pd.read_csv(W / "shuf7/imf2_prank.csv")
m0 = pd.read_csv(W / "main/imf2_scores.csv")
num = [c for c in m0.columns if m0[c].dtype.kind == "f"]
g2 = float(max(np.nanmax(np.abs(m0[num].values - s7[num].values)), np.nanmax(np.abs(pd.read_csv(W / "main/imf2_prank.csv").select_dtypes("float").values - p7.select_dtypes("float").values))))
d7 = np.load(L.RS / "shuf7/draws.npz"); g2d = float(max(np.max(np.abs(D[k] - d7[k])) for k in D.files))
g3 = json.load(open(W / "main/g3_second_path.json"))
oos = {}
for tag in ("main", "shuf7"):
    f = W / f"oos_NH/{tag}/oos_scores.csv"
    if f.exists():
        oos[tag] = pd.read_csv(f)
g2oos = float(np.nanmax(np.abs(oos["main"].select_dtypes("float").values - oos["shuf7"].select_dtypes("float").values))) if len(oos) == 2 else None
# ---------------- release files
keep = ["method", "variant", "status", "MPI", "MPI_lo", "MPI_hi", "P1", "P2", "P3", "C5pts", "P5", "C2pts", "C2pts_lo", "C2pts_hi", "C5pts_lo", "C5pts_hi",
        "FCI", "FCI_lo", "FCI_hi", "S_C2", "S_C5net", "S_C5three_dom", "C2_pairs", "C5_pairs", "C2_ratio_g", "C2_ratio_m", "C2_ratio_n",
        "C5_ratio_g", "C5_ratio_m", "C5_ratio_n", "mean_p_C2", "mean_p_C5", "inherited"]
S[keep].to_csv(REL / "imf2_scores.csv", index=False, float_format="%.6g")
PR.to_csv(REL / "imf2_prank.csv", index=False, float_format="%.6g")
for f in ("imf2_layer_choice.csv", "imf2_waybill_counts.csv", "imf2_links_frozen.json"):
    (REL / f).write_bytes((W / "main" / f).read_bytes())
lk = pd.read_csv(W / "main/imf2_links_by_fold.csv"); lk.to_csv(REL / "imf2_links_by_fold.csv", index=False)
if "main" in oos:
    oos["main"].to_csv(REL / "oos_NH_scores.csv", index=False, float_format="%.6g")
    (REL / "oos_NH_wins.csv").write_bytes((W / "oos_NH/main/oos_wins.csv").read_bytes()); (REL / "oos_NH_meta.json").write_bytes((W / "oos_NH/main/oos_meta.json").read_bytes())
# G5 leak scan
bad = []
for f in REL.iterdir():
    if f.suffix in (".csv", ".json", ".md"):
        txt = f.read_text()
        head = txt.splitlines()[0] if f.suffix == ".csv" else ""
        if re.search(r"(^|,)(CN|PLT_CN|pair_id|t1_id|loc|LAT|LON|PLOT)(,|$)", head):
            bad.append(f"{f.name}: id column")
        if re.search(r"(?<![\d.])\d{10,}(?![\d.])", txt):
            bad.append(f"{f.name}: 10+ digit integer")
checks = dict(prereg_md5="b539cc149e8f1e1606df0b1d46ccc065", G1=dict(pass_=G1, detail=g1x, note="C2 and C5 rebuilt on the registered draws; E8f bt differs by the per-draw foliage factor in the rebuild only (stored draws used for comparators)"),
              G2=dict(pass_=bool(g2 < 1e-9 and g2d < 1e-9 and (g2oos is None or g2oos < 1e-9)), scores_maxdiff=g2, draws_maxdiff=g2d, oos_NH_maxdiff=g2oos),
              G3=dict(pass_=g3["G3_pass"], max_abs_diff=g3["max_abs_diff"]), G4=dict(pass_=meta["G4_pass"], detail=meta["G4"]), G5=dict(pass_=not bad, findings=bad))
json.dump(checks, open(REL / "imf2_checks.json", "w"), indent=1, default=float)
fig = S[S.method.isin(["IMF2_V1_forecast", "IMF2_V1_oracle", "osm_mc", "E3", "E8f"]) | S.method.str.endswith("_forecast")][["method", "MPI", "MPI_lo", "MPI_hi", "C2pts", "C2pts_lo", "C2pts_hi", "C5pts", "C5pts_lo", "C5pts_hi"]]
fig.to_csv(L.T / "work/figure_data.csv", index=False)
print(json.dumps({k: (v["pass_"] if isinstance(v, dict) and "pass_" in v else v) for k, v in checks.items()}, default=float))
print(S[["method", "MPI", "MPI_lo", "MPI_hi", "C2pts", "C2pts_lo", "C2pts_hi", "C5pts", "C5pts_lo", "C5pts_hi", "FCI", "FCI_lo", "FCI_hi"]].round(2).to_string())
print(PR[["method", "FCI_Prank1_IMF2", "FCI_Prank1_osm_mc", "FCI_Prank1_E3", "FCI_Prank1_E8f"]].round(3).to_string())
sys.stdout.flush(); os._exit(0)
