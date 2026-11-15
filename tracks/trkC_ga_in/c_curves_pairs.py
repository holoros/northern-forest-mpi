#!/usr/bin/env python3
"""NCASI report card (task #5 setup): curve engines (CBM gated FVS_FIAW_NEW, CBM American Forests papa2023, FIA hybrid
curves) on FIA pairs for any state. No scoring: per pair model components only (fvs_surv, fvs_mort, eng_ingr in lb ac-1 over
the cycle, the format of the C2 and C5 scorers).

Rules (tools/v24_H/prep_h.py, Track A and A7):
  level  curve level x residual biomass share (treated) or share 1 (pairs), then the curve over the interval (default)
  ageq   age reset sensitivity: common age shift at which the plot curve equals the residual stock; needs --phi (state
         foliage share of CBM AG stock, from the state inventory bundle; 0 for fia_curves)
Inputs: CBM curve tables built by tools/complete/cbm_curves.py from the state's registered SIT and AIDB
(~/perseus_me/ncasi_rc/engines/cbm_curves/{gated,af}_ST; Maine: me_state/complete_j/work), FIA hybrid fits and membership
(Cardinal yield_curves_conus, pulled; Maine: v23_trackA/restricted), state FIA DataMart tables.
Usage: curves_pairs_state.py --state ST --set pairs_hold|pairs_all|treated --engines cbm_gated,cbm_af,fia_curves --out DIR
"""
import argparse, json, os, sys, time
from pathlib import Path
H = Path.home(); RC = H / "perseus_me/ncasi_rc"; ENGD = RC / "engines"
ap = argparse.ArgumentParser()
ap.add_argument("--state", required=True); ap.add_argument("--set", required=True, choices=["pairs_hold", "pairs_all", "treated"])
ap.add_argument("--engines", default="cbm_gated,cbm_af,fia_curves"); ap.add_argument("--rule", default="level", choices=["level", "ageq"])
ap.add_argument("--phi", type=float, default=None); ap.add_argument("--limit", type=int, default=0); ap.add_argument("--out", required=True)
a = ap.parse_args(); st = a.state.upper()
PUL = H / "perseus_me/multistate/restricted/pulled_cardinal"
if st == "ME":
    root = H / "perseus_me/me_state"; fits = root / "v23_trackA/restricted/ycx_ME_hybrid_fits.csv"; memb = root / "v23_trackA/restricted/ycx_membership_ME.csv"
    fia = H / "restricted/ncasi-modeleval/fia_data_rc/ME"
else:
    root = H / "perseus_me/tracks/trkC_gain/restricted/curve_roots" / st; w = root / "complete_j/work"; w.mkdir(parents=True, exist_ok=True); os.chmod(root.parent, 0o700)
    for e in ("gated", "af"):
        src = H / f"perseus_me/tracks/trkC_gain/restricted/cbm_curves/{e}_{st}"; dst = w / f"cbm_curves_{e}"
        if src.exists() and not dst.exists():
            dst.symlink_to(src)
    fits = PUL / f"yield_curves_conus/ycx_{st}_hybrid_fits.csv"; memb = PUL / f"yield_curves_conus/config/ycx_membership_{st}.csv"
    fia = H / f"restricted/ncasi-modeleval/fia_data_rc/{st}"
os.environ.update(PME_STATE=st, PME_STATE_ROOT=str(root), PME_FIA_ROOT=str(fia))
CODE = H / "perseus_me/code"; sys.path.insert(0, str(CODE)); sys.path.insert(0, str(CODE / "tools/v24_H")); sys.path.insert(0, str(CODE / "tools/complete"))
import numpy as np, pandas as pd  # noqa: E402
import prep_h as PH  # noqa: E402
import build_j  # noqa: E402
import pairs as PG  # noqa: E402
BJ = PH._bj(); BJ.FIA = fia; BJ.MEMB = memb; BJ.FITS = fits
IN = H / "perseus_me/multistate/inputs" / st / "restricted"
if a.set == "treated":
    P = pd.read_csv(IN / "treated/treated_pairs.csv", dtype={"pair_id": str, "t1_id": str}); share = (P.B1_Mgha / P.B1full_Mgha).values
else:
    P = pd.read_csv(IN / "pairs_std/pairs_built.csv", dtype={"pair_id": str, "t1_id": str})
    if a.set == "pairs_hold":
        P = P[[PG.half_of(p) == "hold" for p in P.pair_id]]
    share = np.ones(len(P))
o = np.argsort(P.pair_id.values, kind="stable"); P = P.iloc[o].reset_index(drop=True); share = share[o]
if a.limit:
    P = P.head(a.limit); share = share[:a.limit]
out = Path(a.out); out.mkdir(parents=True, exist_ok=True); os.chmod(out, 0o700); t0 = time.time()
P["t1"] = PH.cs(P.t1_id); P["cycle"] = [PG.cycle_of(x) for x in P.REMPER]
U, C = PH.units_for(BJ, pd.DataFrame(dict(unit=P.pair_id.values, PLT_CN=P.t1.values, kind="pair")))
pr = pd.DataFrame(dict(pair_id=P.pair_id.values, unit=P.pair_id.values, cycle=P.cycle.values))
meta = dict(state=st, set=a.set, rule=a.rule, pairs=len(P), fits=str(fits), membership=str(memb), scoring="none (setup runner)", engines={})
for e in a.engines.split(","):
    if e == "cbm_af" and not (root / "complete_j/work/cbm_curves_af/curve_table.csv").exists():
        meta["engines"][e] = "not run: no CBM American Forests curve table for this state (no preset)"; continue
    if a.rule == "level":
        X = PH.level_rows(BJ, e, U, C, pr, share)
    else:
        if a.phi is None and e != "fia_curves":
            raise SystemExit("--phi required for the age reset rule on CBM engines")
        X = PH.ageq_rows(BJ, e, U, C, pr, P.B1_Mgha.values / PH.KM / ((1 - a.phi) if e != "fia_curves" else 1.0))
    X["state"] = st; f = out / f"{e}__{a.rule}__pair_rows.csv"; X.to_csv(f, index=False, float_format="%.12g"); os.chmod(f, 0o600)
    meta["engines"][e] = dict(rows=len(X), ok=int(X.ok.astype(bool).sum()), notes=X.note.astype(str).str[:40].value_counts().head(4).to_dict())
meta["wall_s"] = round(time.time() - t0, 1); meta["s_per_1000_pairs"] = round(1000 * meta["wall_s"] / max(1, len(P)), 1)
(out / "run_meta.json").write_text(json.dumps(meta, indent=1, default=str)); print(json.dumps(meta, default=str))
