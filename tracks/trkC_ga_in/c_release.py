#!/usr/bin/env python3
"""trkC_gain release tables (aggregates only): before (released multistate r5, work/final/score_ST) and after (this track's scores)
for GA and IN, rule set a, common basis, c5_net view (headline) and c5_three view; C5 three-quantity ratios (model sum / observed sum,
treated pairs, all and registered subset) per row from the unit row files, joined by explicit pair_id. usage: c_release.py"""
import json, sys
from pathlib import Path
import numpy as np, pandas as pd
H = Path.home(); T = H / "perseus_me/tracks/trkC_gain"; MS = H / "perseus_me/multistate"; R = T / "release"; K = 0.00112085
ids = lambda s: pd.Series(s).astype(str).str.strip().str.replace(r"\.0$", "", regex=True)  # noqa: E731
COLS = ["row", "MPI", "MPI_lo95", "MPI_hi95", "grade", "tier", "coverage", "P1", "P2", "P3", "P4", "P5"]


def sel(f, view="c5_net"):
    if not Path(f).exists():
        return pd.DataFrame(columns=COLS)
    x = pd.read_csv(f); x = x[(x.rule_set == "a") & (x.basis == "common basis") & (x.view == view)]
    return x[COLS].set_index("row")


def c5_ratios(st, f, label):
    tp = pd.read_csv(MS / "inputs" / st / "restricted/treated/treated_pairs.csv", dtype={"pair_id": str}); tp["pair_id"] = ids(tp.pair_id); assert tp.pair_id.is_unique
    x = pd.read_csv(f, dtype={"pair_id": str}); x["pair_id"] = ids(x.pair_id); x = x[x.ok.astype(str).str.lower().isin(["true", "1"])]
    if "arm" in x.columns and x.arm.nunique() > 1:
        x = x[x.arm == x.arm.iloc[0]]
    assert x.pair_id.is_unique, f
    j = tp.merge(x[["pair_id", "fvs_surv", "fvs_mort", "cycle"]], on="pair_id", how="inner", validate="one_to_one")
    fol = (1 - float(pd.read_csv(MS / "work/common/phi_states.csv").set_index("state").phi[st])) if (label.startswith("cbm_gated (") or label.startswith("cbm_gated wholeplot")) else 1.0   # FOLIAGE rule
    g = j.fvs_surv / j.cycle * K * fol; m = j.fvs_mort / j.cycle * K * fol; out = {}
    for nm, mask in (("all", np.ones(len(j), bool)), ("registered", ~j.t40.astype(str).str.lower().isin(["true", "1"]).values)):
        jj = j[mask]; gg, mm = g[mask], m[mask]
        out[nm] = dict(n=int(len(jj)), gross=float(gg.sum() / jj.obs_gg_a1.sum()), mortality=float(mm.sum() / jj.obs_mort.sum()), net=float((gg - mm).sum() / jj.obs_net_a1b.sum()))
    return dict(state=st, row=label, file=str(f).replace(str(H), "~"), **{f"{k}_{q}": v for k, d in out.items() for q, v in d.items()})


rows, c5 = [], []
for st in ("GA", "IN"):
    b = sel(MS / "work/final" / f"score_{st}/ms_scores.csv"); b3 = sel(MS / "work/final" / f"score_{st}/ms_scores.csv", "c5_three")
    srcs = {"primary": T / f"restricted/score/score_{st}/ms_scores.csv", "cp100": T / f"restricted/score/score_cp100_{st}/ms_scores.csv",
            "cpw": T / f"restricted/score/score_cpw_{st}/ms_scores.csv", "all": T / f"restricted/score/score_all_{st}/ms_scores.csv", "landis_h1": T / f"restricted/score/score_landis_{st}/ms_scores.csv", "wholeplot": T / f"restricted/score/score_wp_{st}/ms_scores.csv"}
    A = {k: sel(v) for k, v in srcs.items()}; A3 = {k: sel(v, "c5_three") for k, v in srcs.items()}
    for row in ("fvs", "fvsb", "E_fvs3", "fia_curves", "cbm_gated", "cbm_gated_x1", "cem", "landis"):
        for ver in ("wholeplot", "primary", "cp100", "cpw", "all", "landis_h1"):
            if ver == "wholeplot" and row not in ("cbm_gated", "cbm_gated_x1"):
                continue
            if ver not in ("primary", "wholeplot") and not (row in ("cbm_gated",) and ver in ("cp100", "cpw")) and not (ver == "all" and row in ("cem", "landis")) and not (ver == "landis_h1" and row == "landis"):
                continue
            a = A[ver]; a3 = A3[ver]
            r = dict(state=st, row=row, version=ver)
            for p, D in (("before", b), ("before_c5three", b3), ("after", a), ("after_c5three", a3)):
                if row in D.index:
                    d = D.loc[row]; r.update({f"{p}_{c}": d[c] for c in COLS[1:] if c in d.index})
            rows.append(r)
    S = T / "restricted/slots"; F = MS / "runs/fvs_default" / st / "restricted/treated_default/pair_rows.csv"
    for lab, f in (("fvs", F), ("cbm_gated (level)", S / f"cbm_gated/{st}/restricted/treated_level__cbm_gated__level__pair_rows.csv"),
                   ("cbm_gated (ageq)", S / f"cbm_gated/{st}/restricted/treated_ageq__cbm_gated__ageq__pair_rows.csv"),
                   ("cbm_gated_x1", S / f"cbm_gated_x1/{st}/restricted/treated__cbm_gated_x1__X1__pair_rows.csv"),
                   ("cem", T / f"restricted/cem/{st}/treated/extract/pair_rows.csv"),
                   ("cbm_gated wholeplot (ageq)", T / f"restricted/slots_wp/cbm_gated/{st}/restricted/treated_ageq__cbm_gated__ageq__pair_rows.csv"),
                   ("cbm_gated_x1 wholeplot", T / f"restricted/slots_wp/cbm_gated_x1/{st}/restricted/treated__cbm_gated_x1__X1__pair_rows.csv")):
        if Path(f).exists():
            c5.append(c5_ratios(st, f, lab))
R.mkdir(exist_ok=True); X = pd.DataFrame(rows); X.to_csv(R / "trkC_before_after.csv", index=False, float_format="%.4f")
pd.DataFrame(c5).to_csv(R / "trkC_c5_three_ratios.csv", index=False, float_format="%.4f"); print(X.to_string()[:6000]); print(pd.DataFrame(c5).drop(columns=["file"]).to_string())
