"""IMF2 gate G1: rebuild C2 and C5 net S and the overall MPI of OSM mortcal, E3_fvs_osm, E8f from frames and registered draws."""
import json, os, sys
from pathlib import Path
import numpy as np, pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parent))
import imf2_lib as L

C2, Mp2, doms, aw = L.c2_frame(); C5, Mp5, _, _ = L.c5_frame()
C2.to_parquet(L.RS / "c2_frame.parquet"); C5.to_parquet(L.RS / "c5_frame.parquet")
np.save(L.RS / "c2_Mp.npy", Mp2)
ST = L.stored()
V = pd.read_csv(L.R / "v26_W/release/mpi26w_a3_views_common_basis.csv")
P3 = {m: float(V[(V.engine == k.split("|")[0]) & (V.arm == k.split("|")[1])].P3.iloc[0]) for m, k in L.MPIKEY.items()}
MPIW = {m: float(V[(V.engine == k.split("|")[0]) & (V.arm == k.split("|")[1])].MPI_c5_net.iloc[0]) for m, k in L.MPIKEY.items()}
zd = np.load(L.K7W / "mpi_draws_B1_a2_a.npz")
rep = {"n_c2": int(len(C2)), "n_pstar": int(C2.pstar.sum()), "n_c5": int(len(C5))}
for m in ("osm_mc", "E3", "E8f"):
    if m == "osm_mc":
        p2 = {q: C2[f"osm_mc_{q}"].values for q in L.Q}; p5 = {q: C5[f"osm_mc_{q}"].values for q in L.Q}
        c5cov = np.isfinite(p5["n"])
    else:
        mem = L.MIX[m]; p2 = L.mix_pred(C2, mem); p5 = None
        c5cov = np.logical_and.reduce([np.isfinite(C5[f"{e}_n"].values) for e in mem])
    s2, b2, _ = L.comp_S(p2, C2, C2.pstar.values & np.isfinite(p2["n"]), doms, aw, Mp2)
    if p5 is None:   # ensemble rule: mean over members of their domain means on the common covered pairs
        dom = C5.dom.values; K = len(doms); mem = L.MIX[m]
        mpt = np.mean([L.dmeans(C5[f"{e}_n"].values, c5cov, dom, doms, Mp5)[0] for e in mem], 0)
        mbt = np.mean([L.dmeans(C5[f"{e}_n"].values, c5cov, dom, doms, Mp5)[1] for e in mem], 0)
        op, ob, _ = L.dmeans(C5.obs_n.values, c5cov, dom, doms, Mp5)
        s5, b5 = float(abs(L.ln(mpt, op)[K])), np.abs(L.ln(mbt, ob)[:, K])
    else:
        s5, b5, _ = L.comp_S(p5, C5, c5cov, doms, aw, Mp5, qset=("n",), view="statewide")
    st = ST[m]
    pil, mp, mb = L.mpi_from(st["C1"], (s2, b2), st["C3"], P3[m], st["C4"], (s5, b5))
    key = f"common basis|c5_net|{L.MPIKEY[m]}"
    rep[m] = dict(C2_pt_ours=s2, C2_pt_reg=st["C2"][0], C2_bt_maxdiff=float(np.nanmax(np.abs(b2 - st["C2"][1]))),
                  C5_pt_ours=s5, C5_pt_reg=st["C5"][0], C5_bt_maxdiff=float(np.nanmax(np.abs(b5 - st["C5"][1]))),
                  MPI_ours=mp, MPI_W=MPIW[m], MPI_bt_maxdiff_vs_K7draws=float(np.nanmax(np.abs(mb - zd[key + "|bt"]))) if key + "|bt" in zd.files else None,
                  MPI_ci=list(np.percentile(mb, [2.5, 97.5])), **pil)
(L.T / "work").mkdir(exist_ok=True)
json.dump(rep, open(L.T / "work/g1_repro.json", "w"), indent=1, default=float)
print(json.dumps(rep, indent=1, default=float))
sys.stdout.flush(); os._exit(0)
