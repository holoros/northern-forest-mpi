"""Track trkB_engfix task 2 diagnostic (PREREG section 4): where the CEM C5 net gap sits. Reads the C5 frame and the CEM pair rows
(default cem_v2p_d3, round4 refined); sums of observed and model gross, mortality and net (scorer units, e1_lib.rates arithmetic)
by removal class, joined by explicit pair id. Aggregates only."""
import json, pickle, sys
import numpy as np, pandas as pd
sys.path.insert(0, "/home/aaron/perseus_me/tracks/refine_engines2/tools"); import e1_lib as E
sys.path.insert(0, "/home/aaron/perseus_me/code/workflow/scripts"); import gate_pairs as GP
fc = pickle.load(open(E.RS2 / "frames_cache.pkl", "rb")); C5 = fc["C5"]
def rows(p):
    D = pd.read_csv(p, dtype={"pair_id": str}); D = D[D.ok.astype(str).str.lower().isin(["true", "1"])].copy(); D["pid"] = D.pair_id.str.strip()
    D["g"] = D.fvs_surv / D.cycle * GP.K; D["m"] = D.fvs_mort / D.cycle * GP.K; D["n"] = D.g - D.m
    assert not D.pid.duplicated().any(); X = D.set_index("pid").reindex(C5.pair_id.values); assert (X.index.values == C5.pair_id.values).all(); return X
out = {}
cls = E.R.icls(C5.rem.values)
for lab, p in (("default", sys.argv[1]), ("refined", sys.argv[2])):
    X = rows(p); ok = np.isfinite(X.n.values) & np.isfinite(C5.obs_n.values)
    res = {}
    for c in ["all", "r<25", "r25-50", "r>=50"]:
        s = ok & ((cls == c) if c != "all" else True)
        o = {q: float(C5[f"obs_{q}"].values[s].sum()) for q in "gmn"}; m = {q: float(X[q].values[s].sum()) for q in "gmn"}
        res[c] = dict(n=int(s.sum()), obs=o, model=m, model_over_obs={q: m[q] / o[q] for q in "gmn"},
                      net_gap=m["n"] - o["n"], gap_from_gross=m["g"] - o["g"], gap_from_mort=-(m["m"] - o["m"]))
    out[lab] = res
json.dump(out, open("work/cem_diag.json", "w"), indent=1); print(json.dumps(out, indent=1))
