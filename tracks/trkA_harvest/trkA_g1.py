"""trkA G1: the copied IMF2 tools reproduce the released IMF2 run (draws, predictions, scores) and trkA_eval reproduces the released SC-L AUC."""
import json
from pathlib import Path
import numpy as np, pandas as pd
T = Path("/home/aaron/perseus_me/tracks/trkA_hcm"); I = Path("/home/aaron/perseus_me/tracks/imf2")
a = np.load(T / "restricted/g1/draws.npz"); b = np.load(I / "restricted/main/draws.npz")
pa = np.load(T / "restricted/g1/preds.npz"); pb = np.load(I / "restricted/main/preds.npz")
s = pd.read_csv(T / "work/g1/imf2_scores.csv").set_index("method"); r = pd.read_csv(I / "release/imf2_scores.csv").set_index("method")
m = pd.read_csv(T / "work/eval/trkA_layer_metrics.csv").set_index("layer").loc["stored_SCL_plot"]
# trkA main run (V1 rows) must equal the g1 run too
c = np.load(T / "restricted/main/draws.npz")
out = dict(draws_max_abs_diff=max(float(np.max(np.abs(a[k] - b[k]))) for k in b.files), draw_keys_equal=sorted(a.files) == sorted(b.files),
           preds_max_abs_diff=max(float(np.nanmax(np.abs(pa[k] - pb[k]))) for k in pb.files if not k.endswith("_ids")),
           scores_max_abs_diff_vs_release=float((s.loc[r.index, ["MPI", "MPI_lo", "MPI_hi", "S_C2", "S_C5net"]] - r[["MPI", "MPI_lo", "MPI_hi", "S_C2", "S_C5net"]]).abs().max().max()),
           main_run_V1_vs_g1=max(float(np.max(np.abs(c[k] - a[k]))) for k in ("IMF2_V1_forecast", "IMF2_V1_oracle", "osm_mc", "E3", "E8f")),
           IMF2_V1_forecast_MPI=float(s.loc["IMF2_V1_forecast", "MPI"]), SCL_AUC=float(m.AUC), SCL_AUC_lo=float(m.AUC_lo), SCL_AUC_hi=float(m.AUC_hi))
out["pass"] = out["draws_max_abs_diff"] < 1e-9 and out["preds_max_abs_diff"] < 1e-9 and out["scores_max_abs_diff_vs_release"] < 1e-9 and out["main_run_V1_vs_g1"] < 1e-9 \
    and abs(round(out["SCL_AUC"], 3) - 0.676) < 1e-9 and abs(round(out["SCL_AUC_lo"], 3) - 0.665) < 1e-9 and abs(round(out["SCL_AUC_hi"], 3) - 0.688) < 1e-9
json.dump(out, open(T / "work/g1_check.json", "w"), indent=1); print(json.dumps(out))
