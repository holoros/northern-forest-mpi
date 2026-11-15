"""POST HOC (not pre-registered): IMF2 V1 own components (C2, C5 net) on E3's own non-pair pillars (C1, C3, P3, C4bfm), paired draws."""
import json, sys, os, numpy as np, pandas as pd
sys.path.insert(0, "/home/aaron/perseus_me/tracks/imf2/tools"); import imf2_lib as L
doms, aw, *_ = L.domains(); ST = L.stored(); P3 = L.p3_exact()
Z = np.load(L.RS / "main/preds.npz", allow_pickle=True)
C2 = pd.read_parquet(L.RS / "c2_frame.parquet"); C5 = pd.read_parquet(L.RS / "c5_frame.parquet")
assert list(Z["c2_ids"]) == list(C2.pair_id) and list(Z["c5_ids"]) == list(C5.pair_id)
Mp2 = np.load(L.RS / "c2_Mp.npy"); Mp5 = np.load(L.K7W / "c5_a/c5h_a_Mp.npy"); e3 = ST["E3"]; out = {}
_, e3pt, e3bt = L.mpi_from(e3["C1"], e3["C2"], e3["C3"], P3["E3"], e3["C4"], e3["C5"])
for m in ("IMF2_V1_forecast", "IMF2_V1_oracle"):
    p2 = {q: Z[f"{m}|C2|{q}"] for q in L.Q}; p5 = {q: Z[f"{m}|C5|{q}"] for q in L.Q}
    c2 = L.comp_S(p2, C2, C2.pstar.values, doms, aw, Mp2)[:2]; c5 = L.comp_S(p5, C5, np.isfinite(p5["n"]), doms, aw, Mp5, qset=("n",), view="statewide")[:2]
    _, mp, mb = L.mpi_from(e3["C1"], c2, e3["C3"], P3["E3"], e3["C4"], c5); d = mb - e3bt
    out[m] = dict(MPI_on_E3_pillars=mp, ci=list(np.percentile(mb, [2.5, 97.5])), E3=e3pt, diff=mp - e3pt, diff_ci=list(np.percentile(d, [2.5, 97.5])))
json.dump(out, open(L.T / "release/imf2_posthoc_e3_pillars.json", "w"), indent=1, default=float); print(json.dumps(out, default=float)); sys.stdout.flush(); os._exit(0)
