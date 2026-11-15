"""trkA G2 row-shuffle comparisons (main vs shuffle seed 7) for every released stage present. usage: trkA_g2.py"""
import json, glob
from pathlib import Path
import numpy as np, pandas as pd
T = Path("/home/aaron/perseus_me/tracks/trkA_hcm"); out = {}
for v in ("A0", "SCL_cty"):
    a = pd.read_parquet(T / f"restricted/layers/{v}.parquet").set_index("CN"); b = pd.read_parquet(T / f"restricted/layers/{v}_shuf7.parquet").set_index("CN").reindex(a.index)
    out[f"layer_{v}"] = float(max(np.abs(a.p - b.p).max(), np.nanmax(np.abs(a.s - b.s))))
if (T / "restricted/shuf7/draws.npz").exists():
    a = np.load(T / "restricted/main/draws.npz"); b = np.load(T / "restricted/shuf7/draws.npz"); out["imf2_draws"] = max(float(np.max(np.abs(a[k] - b[k]))) for k in a.files)
    x = pd.read_csv(T / "work/main/imf2_scores.csv").set_index("method"); y = pd.read_csv(T / "work/shuf7/imf2_scores.csv").set_index("method")
    num = x.select_dtypes("number").columns; out["imf2_scores"] = float((x[num] - y.loc[x.index, num]).abs().max().max())
for f in ("trkA_layer_metrics.csv", "trkA_layer_paired.csv"):
    if (T / "work/eval_shuf7" / f).exists():
        x = pd.read_csv(T / "work/eval" / f); y = pd.read_csv(T / "work/eval_shuf7" / f); num = x.select_dtypes("number").columns
        out[f"eval_{f}"] = float((x[num] - y[num]).abs().max().max())
for st in ("NH", "VT", "MI", "WI", "MN"):
    a, b = T / f"work/state_{st}_hcm2.json", T / f"work/state_{st}_hcm2_shuf7.json"
    if a.exists() and b.exists():
        ja, jb = json.load(open(a)), json.load(open(b)); out[f"state_{st}"] = float(max(abs(ja[k] - jb[k]) for k in ("AUC", "brier", "cal_slope")))
    a, b = T / f"work/oos_{st}/main/oos_scores.csv", T / f"work/oos_{st}/shuf7/oos_scores.csv"
    if a.exists() and b.exists():
        x = pd.read_csv(a); y = pd.read_csv(b); num = x.select_dtypes("number").columns; out[f"oos_{st}"] = float((x[num] - y[num]).abs().max().max())
out["pass"] = all(v <= 1e-9 for k, v in out.items() if k != "pass")
json.dump(out, open(T / "work/g2_checks.json", "w"), indent=1); print(json.dumps(out))
