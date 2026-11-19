#!/usr/bin/env python3
"""Post step (aggregates only): reproduction check, second path MPI, result tables, fig1_scorecard_v3.csv. usage: sens_post.py [--v3]"""
import os, sys, json
import numpy as np, pandas as pd
pd.set_option("display.width", 240); pd.set_option("display.max_columns", 40)
H = os.path.expanduser("~/perseus_me"); S = f"{H}/tracks/ncasi_ms/sens_20261115"; D = f"{H}/tracks/ncasi_ms/data"
f = lambda d: d[(d.rule_set == "a") & (d.basis == "common basis") & (d["view"] == "c5_net")].copy()
W = dict(P1=0.30, P2=0.20, P3=0.20, P4=0.15, P5=0.15); ST = ["ME", "NH", "VT", "MI", "WI", "MN"]
def second(r):
    ev = [p for p in W if pd.notna(r[p])]
    return sum(W[p] * r[p] for p in ev) / sum(W[p] for p in ev) if ev else np.nan
rep, B, Cc, chk = [], [], [], {}
for s in ST:
    p = f"{S}/scores/{s}/sens_scores.csv"
    if not os.path.exists(p):
        continue
    X = pd.read_csv(p); reg = f(pd.read_csv(f"{H}/multistate/work/r4_scores/score_ME/ms_scores.csv" if s == "ME" else f"{H}/me_state/ncasi_{s.lower()}/release/ns_{s.lower()}_scores.csv")).set_index("row")
    b = X[X.run == "base"].set_index("row"); num = [c for c in reg.columns if c in b.columns and reg[c].dtype.kind == "f"]
    has_l = b.notes.fillna("").str.contains("members .*landis") | (b.index == "landis")
    excl = list(b.index[has_l]) if s == "ME" else []
    cols = [c for c in num if not (s == "ME" and c == "P_rank1")]
    keep = [r for r in reg.index if r not in excl]
    a1, a2 = reg.loc[keep, cols].astype(float).values, b.reindex(keep)[cols].astype(float).values
    d = np.where(np.isnan(a1) & np.isnan(a2), 0.0, np.abs(a1 - a2)); nanmis = int(np.isnan(d).sum())
    txt = bool((reg.loc[keep, ["tier", "grade", "pillars_evaluated"]].astype(str).values == b.reindex(keep)[["tier", "grade", "pillars_evaluated"]].astype(str).values).all())
    X["MPI_second_path"] = X.apply(second, axis=1); sp = float((X.MPI_second_path - X.MPI).abs().max())
    k = json.load(open(f"{S}/scores/{s}/checks.json"))
    chk[s] = dict(repro_rows=len(keep), repro_cols=len(cols), repro_max_abs=float(np.nanmax(d)), repro_nan_mismatch=nanmis, repro_text_equal=txt, rows_missing=int(b.reindex(keep).MPI.isna().sum() - reg.loc[keep].MPI.isna().sum()),
                  excluded_rows=excl, second_path_max_abs=sp, shuffle=k.get("shuffle_base_and_loo"), shuffle_c1in=k.get("shuffle_c1in"),
                  base_repeat=k.get("base_repeat_after_switches"), loo_base=k.get("loo_base_rows_equal_base"), register_switch=k.get("register_switch"), wall_s=k.get("wall_s"))
    c = X[X.run == "c1in"].set_index("row")
    for r in ("landis", "E_state"):
        B.append(dict(state=s, row=r, reg_MPI=reg.loc[r, "MPI"], reg_lo95=reg.loc[r, "MPI_lo95"], reg_hi95=reg.loc[r, "MPI_hi95"], reg_tier=reg.loc[r, "tier"], reg_coverage=reg.loc[r, "coverage"],
                      base_MPI=b.loc[r, "MPI"], base_lo95=b.loc[r, "MPI_lo95"], base_hi95=b.loc[r, "MPI_hi95"], base_tier=b.loc[r, "tier"], base_coverage=b.loc[r, "coverage"], base_P1=b.loc[r, "P1"],
                      sens_MPI=c.loc[r, "MPI"], sens_lo95=c.loc[r, "MPI_lo95"], sens_hi95=c.loc[r, "MPI_hi95"], sens_tier=c.loc[r, "tier"], sens_coverage=c.loc[r, "coverage"], sens_P1=c.loc[r, "P1"],
                      sens_pillars=c.loc[r, "pillars_evaluated"], sens_P_rank1=c.loc[r, "P_rank1"], base_P_rank1=b.loc[r, "P_rank1"]))
    cm = c[(c.role == "member") & (c.tier == 1)]; B[-1]["sens_best_tier1_member"] = cm.MPI.idxmax() if len(cm) else ""; B[-1]["sens_best_tier1_member_MPI"] = cm.MPI.max() if len(cm) else np.nan
    L = X[(X.run == "loo") & X.row.str.startswith("E_state_wo_")]
    for _, r in L.iterrows():
        Cc.append(dict(state=s, left_out=r.row.replace("E_state_wo_", ""), reg_E_state_MPI=reg.loc["E_state", "MPI"], base_E_state_MPI=b.loc["E_state", "MPI"], base_E_state_lo95=b.loc["E_state", "MPI_lo95"], base_E_state_hi95=b.loc["E_state", "MPI_hi95"],
                       loo_MPI=r.MPI, loo_lo95=r.MPI_lo95, loo_hi95=r.MPI_hi95, loo_tier=r.tier, loo_coverage=r.coverage, pillars=r.pillars_evaluated,
                       best_tier1_member=r.best_tier1_member, best_tier1_member_MPI=r.best_tier1_member_MPI, loo_is_highest_tier1=bool(r.tier == 1 and r.MPI > r.best_tier1_member_MPI), P_rank1_loo_set=r.P_rank1_loo_set))
    if s == "ME":
        cols5 = ["MPI", "MPI_lo95", "MPI_hi95", "grade", "tier", "coverage", "P1", "P2", "P3", "P4", "P5", "P_rank1"]
        print("ME before (registered r4) / after (driver base):"); print(pd.concat([reg.loc[excl, cols5].assign(when="before"), b.loc[excl, cols5].assign(when="after")]).round(2).to_string())
        print(b.loc[["landis"], ["pillars_evaluated", "notes"]].to_string())
        if "--v3" in sys.argv:
            v2 = pd.read_csv(f"{D}/fig1_scorecard_v2.csv"); v3 = v2.copy(); ren = {"MPI_lo95": "lo95", "MPI_hi95": "hi95", "P4": "P4_C5net", "P4_lo95": "P4_C5net_lo95", "P4_hi95": "P4_C5net_hi95"}
            nb = b.rename(columns=ren)
            for r in excl:
                i = v3.index[(v3.state == "ME") & (v3.model == r)]; assert len(i) == 1
                for cc in v2.columns:
                    if cc in nb.columns and cc not in ("state", "role"):
                        v3.loc[i, cc] = nb.loc[r, cc]
                v3.loc[i, "source"] = "sens_20261115 (ns_score path, ME LANDIS-II B25 C5 run of 2026-11-15; post hoc, registered r4 row kept in v2)"
            assert list(v3.columns) == list(v2.columns) and len(v3) == len(v2); out = f"{D}/fig1_scorecard_v3.csv"; assert not os.path.exists(out)
            v3.to_csv(out, index=False); print("wrote", out, "rows replaced", excl)
json.dump(chk, open(f"{S}/sens_checks.json", "w"), indent=1, default=float)
print(json.dumps(chk, default=float))
if B:
    Bd = pd.DataFrame(B); Bd["repro_max_abs_diff_switch_off"] = Bd.state.map(lambda x: chk[x]["repro_max_abs"]); Bd.to_csv(f"{S}/landis_c1_sensitivity.csv", index=False, float_format="%.6g"); print(Bd.round(2).to_string(index=False))
if Cc:
    Cd = pd.DataFrame(Cc); Cd["repro_max_abs_diff_switch_off"] = Cd.state.map(lambda x: chk[x]["repro_max_abs"]); Cd.to_csv(f"{S}/ensemble_leave_one_out.csv", index=False, float_format="%.6g"); print(Cd.round(2).to_string(index=False))
