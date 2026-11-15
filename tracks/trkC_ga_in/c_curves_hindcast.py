#!/usr/bin/env python3
"""Track R (NCASI report card) P2 hindcast chains for the curve engines (CBM gated, CBM American Forests, FIA curves),
Maine track J chain rule (tools/complete/build_j.py apply): chain trajectory = observed visit 1 stock + k x (curve(age + t)
- curve(age)), k = observed / unexpanded visit 1 stock (k = 1 without visit 1 live trees), offsets 0..30 yr.
Writes OUT/<engine>__<arm>__chain.parquet (chain_id, offset, agb_full, agc_full, agb_surv, agc_surv) and the unanchored
version, plus run_meta.json. Plot level, restricted. No scoring. Usage: hindcast_curves_state.py ST OUT"""
import json, os, sys, time
from pathlib import Path
import numpy as np, pandas as pd
st = sys.argv[1].upper(); OUT = Path(sys.argv[2]); OUT.mkdir(parents=True, exist_ok=True); os.chmod(OUT, 0o700)
H = Path.home(); RC = H / "perseus_me/ncasi_rc"; ENGD = RC / "engines"; PUL = H / "perseus_me/multistate/restricted/pulled_cardinal"
if st == "ME":
    root = H / "perseus_me/me_state"; fits = root / "v23_trackA/restricted/ycx_ME_hybrid_fits.csv"; memb = root / "v23_trackA/restricted/ycx_membership_ME.csv"
else:
    root = H / "perseus_me/tracks/trkC_gain/restricted/curve_roots" / st; fits = PUL / f"yield_curves_conus/ycx_{st}_hybrid_fits.csv"; memb = PUL / f"yield_curves_conus/config/ycx_membership_{st}.csv"
fia = H / f"restricted/ncasi-modeleval/fia_data_rc/{st}"
os.environ.update(PME_STATE=st, PME_STATE_ROOT=str(root), PME_FIA_ROOT=str(fia))
CODE = H / "perseus_me/code"
for p in ("", "tools/v24_H", "tools/complete", "tools/score", "tools/hindcast_all", "tools/hindcast"):
    sys.path.insert(0, str(CODE / p))
import prep_h as PH  # noqa: E402
import score as S  # noqa: E402
BJ = PH._bj(); BJ.FIA = fia; BJ.MEMB = memb; BJ.FITS = fits; t0 = time.time()
Hs = S.h1_setup(H / f"perseus_me/tracks/multistate/config/multistate/hindcast/{st}.yaml", H / f"perseus_me/multistate/inputs/{st}/restricted/hindcast")
h1 = Hs["h1"]; v1 = Hs["v1plt"].reindex(h1).astype(str).str.strip(); O = Hs["O"]; o1 = O[O.visit == 1].set_index("chain_id")
Tx = Hs["Tx"]; t1 = Tx[(Tx.visit == 1) & Tx.live]
un_c = (t1.CARBON_AG.fillna(0) * t1.TPA_UNADJ).groupby(t1.chain_id).sum(); un_b = (t1.DRYBIO_AG * t1.TPA_UNADJ).groupby(t1.chain_id).sum()
ch = pd.DataFrame({"chain_id": h1, "PLT_CN": v1.values}); ch["agc1"] = o1.agc_stock.reindex(h1).values; ch["agb1"] = o1.agb_stock.reindex(h1).values
ch["kc"] = (o1.agc_stock / un_c).reindex(h1).values; ch["kb"] = (o1.agb_stock / un_b).reindex(h1).values
U, C = PH.units_for(BJ, pd.DataFrame(dict(unit=ch.chain_id.values, PLT_CN=ch.PLT_CN.values, kind="chain")))
off = np.arange(0, 31); arms = dict(cbm_gated="FVS_FIAW_NEW"); meta = dict(state=st, chains=len(ch), engines={}, scoring="none")
have_af = st == "ME" or (root / "complete_j/work/cbm_curves_af/curve_table.csv").exists()
for e, arm in arms.items():
    if e == "cbm_af" and not have_af:
        meta["engines"][e] = "not run: no preset"; continue
    V = BJ.unit_values(e, U, C, off); rows, rows_u, stc = [], [], {}
    for r in ch.itertuples():
        v, note = V[r.chain_id]; k = note.split(" for ")[0]; stc[k] = stc.get(k, 0) + 1
        if v is None:
            continue
        c, b, _ = v; kc = r.kc if np.isfinite(r.kc) else 1.0; kb = r.kb if np.isfinite(r.kb) else 1.0
        rows.append(pd.DataFrame(dict(chain_id=r.chain_id, offset=off, agb_full=r.agb1 + kb * (b - b[0]), agc_full=r.agc1 + kc * (c - c[0]), agb_surv=np.nan, agc_surv=np.nan)))
        rows_u.append(pd.DataFrame(dict(chain_id=r.chain_id, offset=off, agb_full=kb * b, agc_full=kc * c, agb_surv=np.nan, agc_surv=np.nan)))
    P = pd.concat(rows, ignore_index=True); P.to_parquet(OUT / f"{e}__{arm}__chain.parquet", index=False)
    pd.concat(rows_u, ignore_index=True).to_parquet(OUT / f"{e}__{arm}__unanchored__chain.parquet", index=False)
    meta["engines"][e] = dict(chains_valued=int(P.chain_id.nunique()), status=stc)
for f in OUT.iterdir():
    os.chmod(f, 0o600)
meta["wall_s"] = round(time.time() - t0, 1); (OUT / "run_meta.json").write_text(json.dumps(meta, indent=1, default=str)); print(json.dumps(meta, default=str))
