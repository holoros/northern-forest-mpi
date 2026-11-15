"""trkA claim ledger: every number quoted in RELEASE.md is read here from the release files (source file, row, column); claims.csv."""
import json
from pathlib import Path
import pandas as pd
T = Path("/home/aaron/perseus_me/tracks/trkA_hcm"); R = T / "release"; rows = []
def add(cid, text, val, src, key, kind="number"):
    rows.append(dict(claim_id=cid, claim=text, value=val, source=f"release/{src}", key=key, type=kind, status="VERIFIED by tools/trkA_second_layer.py or trkA_second_imf2.py" if cid.startswith(("L1", "L2", "L3", "M1", "M2")) else "read from release file"))
M = pd.read_csv(R / "trkA_layer_metrics.csv").set_index("layer"); D = pd.read_csv(R / "trkA_layer_paired.csv")
for i, k in enumerate(["A0", "SCL_cty", "stored_SCL_plot", "A2", "A4", "A5", "stored_HCMrefit_plot", "PERM", "A1", "A3", "A6", "stored_SCF_plot"]):
    r = M.loc[k]
    for c in ("AUC", "AUC_lo", "AUC_hi", "brier", "brier_lo", "brier_hi", "cal_slope", "cal_slope_lo", "cal_slope_hi", "E_over_O", "AUC_cty_lo", "AUC_cty_hi"):
        add(f"L{i+1}.{c}", f"layer {k} {c}", round(float(r[c]), 4), "trkA_layer_metrics.csv", f"layer={k};{c}")
for j, r in D.iterrows():
    for c in ("dAUC", "dAUC_lo", "dAUC_hi", "dAUC_cty_lo", "dAUC_cty_hi", "dBrier", "dBrier_lo", "dBrier_hi"):
        add(f"P{j+1}.{c}", f"paired {r.layer} minus {r.reference} {c}", round(float(r[c]), 4), "trkA_layer_paired.csv", f"layer={r.layer};reference={r.reference};{c}")
pr = json.load(open(R / "trkA_promotion.json")); add("D1", "HCM2 promoted over SC-L", pr["promoted"], "trkA_promotion.json", "promoted", "decision")
S = pd.read_csv(R / "trkA_imf2_scores.csv").set_index("method"); P = pd.read_csv(R / "trkA_imf2_prank.csv").set_index("method")
for i, m in enumerate(["IMF2_V1_forecast", "IMF2_Rs_forecast", "IMF2_Rh_forecast", "IMF2_R4_forecast", "IMF2_V1_oracle", "osm_mc", "E3", "E8f"]):
    for c in ("MPI", "MPI_lo", "MPI_hi", "C2pts", "C5pts"):
        add(f"M{i+1}.{c}", f"{m} {c}", round(float(S.loc[m, c]), 2), "trkA_imf2_scores.csv", f"method={m};{c}")
for i, m in enumerate([x for x in P.index]):
    for c in [c for c in P.columns if c.startswith("dMPI") or c.startswith("Prank1")]:
        if pd.notna(P.loc[m, c]):
            add(f"Q{i+1}.{c}", f"{m} {c}", round(float(P.loc[m, c]), 3), "trkA_imf2_prank.csv", f"method={m};{c}")
St = pd.read_csv(R / "trkA_states_layer.csv").set_index("state")
for st in St.index:
    for c in ("pairs", "cut", "AUC", "AUC_lo", "AUC_hi", "brier", "cal_slope", "cal_slope_lo", "cal_slope_hi", "E_over_O"):
        add(f"S.{st}.{c}", f"{st} state HCM2 {c}", round(float(St.loc[st, c]), 4), "trkA_states_layer.csv", f"state={st};{c}")
O = pd.read_csv(R / "trkA_oos_scores.csv"); W = pd.read_csv(R / "trkA_oos_wins.csv")
for _, r in O[O.method.str.startswith("IMF2_V1") | O.method.str.startswith("E3") | O.method.str.startswith("E8f")].iterrows():
    for c in ("S_state", "S_state_lo", "S_state_hi", "pairs"):
        add(f"O.{r.state}.{r.frame}.{r.method}.{c}", f"{r.state} {r.frame} {r.method} {c}", round(float(r[c]), 4), "trkA_oos_scores.csv", f"state={r.state};frame={r.frame};method={r.method};{c}")
for _, r in W.iterrows():
    for c in ("dS", "dS_lo", "dS_hi"):
        add(f"W.{r.state}.{r.frame}.{r.method}.{r.comparator}.{c}", f"{r.state} {r.frame} {r.method} minus {r.comparator} {c}", round(float(r[c]), 4), "trkA_oos_wins.csv",
            f"state={r.state};frame={r.frame};method={r.method};comparator={r.comparator};{c}")
pd.DataFrame(rows).to_csv(R / "claims.csv", index=False); print(len(rows))
