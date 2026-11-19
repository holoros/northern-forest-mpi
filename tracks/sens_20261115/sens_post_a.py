#!/usr/bin/env python3
"""Task (a) post step from scores/ME/sens_scores_base.csv (switch-off run): reproduction, second path, before/after, optional v3. Aggregates only."""
import os, sys, json
import numpy as np, pandas as pd
pd.set_option("display.width", 240); pd.set_option("display.max_columns", 40)
H = os.path.expanduser("~/perseus_me"); S = f"{H}/tracks/ncasi_ms/sens_20261115"; D = f"{H}/tracks/ncasi_ms/data"
f = lambda d: d[(d.rule_set == "a") & (d.basis == "common basis") & (d["view"] == "c5_net")].copy()
reg = f(pd.read_csv(f"{H}/multistate/work/r4_scores/score_ME/ms_scores.csv")).set_index("row"); b = pd.read_csv(f"{S}/scores/ME/sens_scores_base.csv").set_index("row")
excl = list(b.index[b.notes.fillna("").str.contains("members .*landis") | (b.index == "landis")]); keep = [r for r in reg.index if r not in excl]
num = [c for c in reg.columns if c in b.columns and reg[c].dtype.kind == "f" and c != "P_rank1"]
a1, a2 = reg.loc[keep, num].astype(float).values, b.reindex(keep)[num].astype(float).values
d = np.where(np.isnan(a1) & np.isnan(a2), 0.0, np.abs(a1 - a2))
txt = bool((reg.loc[keep, ["tier", "grade", "pillars_evaluated"]].astype(str).values == b.reindex(keep)[["tier", "grade", "pillars_evaluated"]].astype(str).values).all())
W = dict(P1=0.30, P2=0.20, P3=0.20, P4=0.15, P5=0.15)
sp = b.apply(lambda r: sum(W[p] * r[p] for p in W if pd.notna(r[p])) / sum(W[p] for p in W if pd.notna(r[p])), axis=1)
# rows that contain LANDIS-II: every column other than those the C5 slot can move
fixed = [c for c in num if c.startswith(("P1", "P2", "P3", "P5"))]; dl = np.abs(reg.loc[excl, fixed].astype(float).values - b.loc[excl, fixed].astype(float).values); dl = np.where(np.isnan(reg.loc[excl, fixed].astype(float).values) & np.isnan(b.loc[excl, fixed].astype(float).values), 0.0, dl)
chk = dict(rows_registered=len(reg), rows_new=len(b), same_row_set=sorted(reg.index) == sorted(b.index), unchanged_rows=keep, unchanged_cols=len(num), unchanged_max_abs=float(np.nanmax(d)), unchanged_nan_mismatch=int(np.isnan(d).sum()),
           unchanged_text_equal=txt, replaced_rows=excl, replaced_rows_P1_P2_P3_P5_max_abs=float(np.nanmax(dl)), replaced_nan_mismatch=int(np.isnan(dl).sum()), second_path_max_abs=float((sp - b.MPI).abs().max()))
json.dump(chk, open(f"{S}/task_a_checks.json", "w"), indent=1); print(json.dumps(chk))
c5 = ["MPI", "MPI_lo95", "MPI_hi95", "grade", "tier", "coverage", "P1", "P2", "P3", "P4", "P5", "P_rank1"]
T = pd.concat([reg.loc[excl, c5].assign(when="before_r4"), b.loc[excl, c5].assign(when="after")]); T.to_csv(f"{S}/task_a_maine_before_after.csv", float_format="%.6g"); print(T.round(2).to_string())
print(b[["MPI", "tier", "coverage", "P_rank1"]].join(reg[["MPI", "tier", "P_rank1"]], rsuffix="_r4").round(3).to_string()); print(b.loc["landis", "notes"])
if "--v3" in sys.argv:
    v2 = pd.read_csv(f"{D}/fig1_scorecard_v2.csv"); v3 = v2.copy(); nb = b.rename(columns={"MPI_lo95": "lo95", "MPI_hi95": "hi95", "P4": "P4_C5net", "P4_lo95": "P4_C5net_lo95", "P4_hi95": "P4_C5net_hi95"})
    for r in excl:
        i = v3.index[(v3.state == "ME") & (v3.model == r)]; assert len(i) == 1
        for cc in v2.columns:
            if cc in nb.columns and cc not in ("state", "role"):
                v3.loc[i, cc] = nb.loc[r, cc]
        v3.loc[i, "source"] = "sens_20261115 (ns_score path; ME LANDIS-II B25 C5 run of 2026-11-15; post hoc; registered r4 row stays in v2)"
    out = f"{D}/fig1_scorecard_v3.csv"; assert list(v3.columns) == list(v2.columns) and len(v3) == len(v2) and not os.path.exists(out)
    v3.to_csv(out, index=False); print("wrote", out, "rows replaced", excl)
