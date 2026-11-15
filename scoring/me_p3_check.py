"""Track R realism, Maine check (aggregates only): (a) score_v2.R on the Maine v2.4 Bakuzis work files through rc_p3 reproduces the
released noR5 realism and coverage of each single engine (release bakuzis_scores.csv, and P3 in mpi23 common basis);
(b) rc_p3.p3s_score on the Maine v2.4 S lines with the Maine row set reproduces S_p3s_table P3s and evaluability;
(c) the state_p3 adapter path for ME (members of E_state) returns the same single-engine values."""
import json, os, sys
from pathlib import Path
import numpy as np
import pandas as pd
sys.path.insert(0, str(Path(__file__).parent))
import rc_p3 as P3

H = Path(os.path.expanduser("~")); MS = H / "perseus_me/me_state"; OUT = H / "perseus_me/ncasi_rc/scorer_tests/ME"; OUT.mkdir(parents=True, exist_ok=True)
res = {}
rel = P3.realism_from_scores(MS / "bakuzis/release/bakuzis_scores.csv")
W = OUT / "p3_me_work"
import shutil
shutil.rmtree(W, ignore_errors=True); shutil.copytree(MS / "bakuzis/work", W)
ours = P3.run_bakuzis(W, OUT / "p3_me_out")
cmp = []
for k, v in rel.items():
    if k in ours:
        cmp.append(dict(engine=k[0], arm=k[1], realism_rel=v["realism"], realism_rc=ours[k]["realism"], cov_rel=v["coverage"], cov_rc=ours[k]["coverage"]))
C = pd.DataFrame(cmp)
res["bakuzis_vs_release"] = dict(rows=len(C), max_abs_realism=float(np.nanmax(np.abs(C.realism_rel - C.realism_rc))), max_abs_cov=float(np.nanmax(np.abs(C.cov_rel - C.cov_rc))))
T = pd.read_csv(MS / "mpi23/release/mpi23_common_basis_scores.csv"); T = T[(T.view == "common basis") & (T.row_group == "original arm")]
bk = {("fvs", "default"): ("fvs", "default"), ("osm", "segmented"): ("osm", "continuous"), ("acd141", "acd_v141"): ("acd141", "acd_v141_csi14"),
      ("cem", "engine_seeded_v2"): ("cem", "no_harvest_baseline"), ("fvs", "acd_default"): ("fvs", "acd_default")}
res["P3_vs_mpi23"] = {f"{e}.{a}": dict(registered=float(T[(T.engine == e) & (T.arm == a)].P3.iloc[0]), rc=100 * ours[k]["realism"]) for (e, a), k in bk.items() if k in ours}
# (b)
man = pd.read_csv(MS / "v24_S/work/score/manifest.csv"); man = man[man.path.notna() & (man.path != "")]
lines = {f"{r.engine}.{r.arm}": (r.engine, Path(r.path)) for r in man.itertuples() if Path(r.path).exists()}
hz = {f"{r.engine}.{r.arm}": int(r.horizon_yr) for r in man.itertuples() if f"{r.engine}.{r.arm}" in lines}
S, drop = P3.p3s_score(lines, hz)
ref = pd.read_csv(MS / "v24_S/release/S_p3s_table.csv")
d = []
for r in ref.itertuples():
    k = f"{r.engine}.{r.arm}"
    if k in S:
        d.append(dict(row=k, P3s_rel=r.P3s, P3s_rc=100 * S[k]["realism"] if np.isfinite(S[k]["realism"]) else np.nan, ok_rel=bool(r.P3s_evaluable), ok_rc=S[k]["evaluable"]))
Dd = pd.DataFrame(d)
res["p3s_vs_release"] = dict(rows=len(Dd), max_abs_P3s=float(np.nanmax(np.abs(Dd.P3s_rel - Dd.P3s_rc))), evaluable_equal=bool((Dd.ok_rel == Dd.ok_rc).all()),
                             dropped=drop, missing_paths=int(len(man) - len(lines)))
# (c)
eng = {"fvs": "fvs", "osm": "osm", "acd141": "acd141", "cem": "cem", "landis": "landis", "cbm_gated": "cbm_gated", "cbm_af": "cbm_af", "fia_curves": "fia_curves"}
p3, p3s, ens, rep = P3.state_p3("ME", eng, {"E_state": (list(eng), "mean")}, OUT / "p3_state_me")
res["state_p3_ME"] = dict(P3={k: round(100 * v["realism"], 4) for k, v in p3.items()}, P3s={k: (round(100 * v["realism"], 4) if v["evaluable"] else None) for k, v in p3s.items()},
                          ensemble={f"{k[0]}|{k[1]}": round(100 * v["realism"], 4) for k, v in ens.items() if np.isfinite(v.get("realism", np.nan))},
                          note="P3s here uses the ME Track R row set (non-discriminating drop over this set), so it can differ from the v2.4 table by design")
(OUT / "me_p3_check.json").write_text(json.dumps(res, indent=1, default=float)); print(json.dumps(res, indent=1, default=float))
