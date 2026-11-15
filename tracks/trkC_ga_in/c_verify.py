#!/usr/bin/env python3
"""trkC_gain G3 second code path (no scorer import): (1) MPI recomputed from the written pillar points with the registered weights
(P1 30, P2 20, P3 20, P4 15, P5 15, renormalised over evaluated pillars); (2) C5 net point recomputed from the treated pair rows and
the a1 observed values: S = |ln(sum model net / sum observed net)| over all treated pairs (common basis) and the registered subset (t40 false); CBM gated rows x (1 - phi) (FOLIAGE rule), points = 100 max(0,
1 - (e^S - 1) / 0.5); (3) C2 statewide net ratio from pair rows. Joins by explicit pair_id with asserted uniqueness. Writes
work/verifier_ST.json."""
import json, sys
from pathlib import Path
import numpy as np, pandas as pd
T = Path.home() / "perseus_me/tracks/trkC_gain"; MS = Path.home() / "perseus_me/multistate"; K = 0.00112085; W = dict(P1=30, P2=20, P3=20, P4=15, P5=15)
ids = lambda s: pd.Series(s).astype(str).str.strip().str.replace(r"\.0$", "", regex=True)  # noqa: E731
pts = lambda S: 100.0 * max(0.0, 1.0 - (np.exp(S) - 1.0) / 0.5)  # noqa: E731
def rows(f):
    x = pd.read_csv(f, dtype={"pair_id": str}); x["pair_id"] = ids(x.pair_id); x = x[x.ok.astype(str).str.lower().isin(["true", "1"])]; assert x.pair_id.is_unique; return x
import os
TAG = os.environ.get("TRKC_TAG", ""); SC = "score_" + (TAG + "_" if TAG else ""); SLN = "slots" + ("_" + TAG if TAG else "")
for st in sys.argv[1:]:
    S = pd.read_csv(T / f"restricted/score/{SC}{st}/ms_scores.csv"); S = S[(S.rule_set == "a") & (S.basis == "common basis")]; out = {}
    for view in ("c5_net", "c5_three"):
        for _, r in S[S.view == view].iterrows():
            ev = [p for p in W if np.isfinite(r[p])]; m = sum(W[p] * r[p] for p in ev) / sum(W[p] for p in ev) if ev else np.nan
            out[f"{view}/{r.row}"] = dict(mpi_file=float(r.MPI), mpi_recomputed=float(m), abs_diff=float(abs(m - r.MPI)), P4_file=float(r.P4) if np.isfinite(r.P4) else None)
    tp = pd.read_csv(MS / "inputs" / st / "restricted/treated/treated_pairs.csv", dtype={"pair_id": str}); tp["pair_id"] = ids(tp.pair_id); assert tp.pair_id.is_unique
    reg = tp[~tp.t40.astype(str).str.lower().isin(["true", "1"])]
    sl = T / "restricted" / SLN; c5 = {}
    for lab, f in (("cbm_gated_ageq", sl / f"cbm_gated/{st}/restricted/treated_ageq__cbm_gated__ageq__pair_rows.csv"), ("cbm_gated_level", sl / f"cbm_gated/{st}/restricted/treated_level__cbm_gated__level__pair_rows.csv"),
                   ("cbm_gated_x1", sl / f"cbm_gated_x1/{st}/restricted/treated__cbm_gated_x1__X1__pair_rows.csv")):
        x = rows(f); fol = (1 - float(pd.read_csv(MS / "work/common/phi_states.csv").set_index("state").phi[st])) if lab.startswith("cbm_gated_") and lab != "cbm_gated_x1" else 1.0
        for sub, D in (("all", tp), ("registered", reg)):
            j = D.merge(x, on="pair_id", how="inner", validate="one_to_one"); net = (j.fvs_surv - j.fvs_mort) / j.cycle * K * fol
            Sv = abs(np.log(net.sum() / j.obs_net_a1b.sum())); c5[f"{lab}/{sub}"] = dict(n=len(j), ratio=float(net.sum() / j.obs_net_a1b.sum()), S=float(Sv), points=pts(Sv))
    out["c5_second_path"] = c5
    a1 = pd.read_parquet(MS / "inputs" / st / "restricted/a1/a1_pairs.parquet"); a1["pair_id"] = ids(a1.pair_id)
    for lab, f in (("cbm_gated", sl / f"cbm_gated/{st}/restricted/pairs_level__cbm_gated__level__pair_rows.csv"), ("cbm_gated_x1", sl / f"cbm_gated_x1/{st}/restricted/pairs_level__cbm_gated_x1__X1__pair_rows.csv")):
        x = rows(f); j = a1.merge(x, on="pair_id", how="inner", validate="one_to_one"); g = j.fvs_surv / j.cycle * K; mo = j.fvs_mort / j.cycle * K
        out[f"c2_{lab}"] = dict(n=len(j), gross=float(g.sum() / j.obs_gg_a1.sum()), mortality=float(mo.sum() / j.obs_mort.sum()), net=float((g - mo).sum() / j.obs_net_a1b.sum()))
    out["max_mpi_abs_diff"] = max(v["abs_diff"] for k, v in out.items() if "/" in k and np.isfinite(v["abs_diff"]))
    (T / f"work/verifier_{TAG + chr(95) if TAG else chr(32).strip()}{st}.json").write_text(json.dumps(out, indent=1, default=float)); print(st, json.dumps({k: out[k] for k in ("max_mpi_abs_diff", "c5_second_path")}, default=float), {k: v for k, v in out.items() if k.startswith("c2_")})
