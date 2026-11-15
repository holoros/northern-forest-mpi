#!/usr/bin/env python3
"""trkC_gain checks on written score files: XX1 (released ms_score_x rows reproduced by each extension run, %.6g both sides) and G2
(row-shuffled rescoring equals the primary for every row; key = rule_set, basis, view, row). Writes work/xx1_files.json, work/g2_ST.json."""
import json
from pathlib import Path
import numpy as np, pandas as pd
T = Path.home() / "perseus_me/tracks/trkC_gain"; MS = Path.home() / "perseus_me/multistate"; K = ["rule_set", "basis", "view", "row"]
def cmp(a, b):
    J = a.merge(b, on=K, how="inner", suffixes=("", "_x"), validate="one_to_one"); num = [c for c in a.columns if a[c].dtype.kind == "f" and not c.startswith("P_rank1") and c + "_x" in J]
    d = 0.0
    for c in num:
        x, y = J[c].values, J[c + "_x"].values; m = np.isfinite(x) & np.isfinite(y)
        d = max(d, float(np.abs(x[m] - y[m]).max()) if m.any() else 0.0, np.inf if (np.isfinite(x) != np.isfinite(y)).any() else 0.0)
    return dict(rows_a=len(a), rows_b=len(b), joined=len(J), cols=len(num), max_abs_diff=d)
xx = {}
for d in sorted((T / "restricted/score").glob("*_*")):
    f = d / "ms_scores.csv"
    if not f.exists() or d.name.startswith("shuffle"):
        continue
    st = d.name.split("_")[-1]; r = cmp(pd.read_csv(MS / "work/final" / f"score_{st}/ms_scores.csv"), pd.read_csv(f)); r["passed"] = r["max_abs_diff"] <= 1e-9 and r["joined"] == r["rows_a"]; xx[d.name] = r
(T / "work/xx1_files.json").write_text(json.dumps(xx, indent=1)); print("XX1", json.dumps(xx))
for st in ("GA", "IN"):
    a, b = T / f"restricted/score/score_{st}/ms_scores.csv", T / f"restricted/score/shuffle_{st}/ms_scores.csv"
    if a.exists() and b.exists():
        r = cmp(pd.read_csv(a), pd.read_csv(b)); r["passed"] = r["max_abs_diff"] <= 1e-9 and r["joined"] == r["rows_a"] == r["rows_b"]
        (T / f"work/g2_{st}.json").write_text(json.dumps(r, indent=1)); print("G2", st, r)
