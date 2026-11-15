#!/usr/bin/env python3
"""Track R: LANDIS-II hindcast chains (P2) for any state, Maine track I rule (tools/landis_full/post.py h1, arm
registered_tables): sites from the chains' observed visit 1 live trees, 30 yr run, chain trajectory = LANDIS relative change
B_t / B_0 applied to the observed visit 1 stock (agb and agc); chains without a built site get 0 (Maine rule).
  frame ST OUT                 -> OUT/sites/{sites,keys,trees}.csv
  post ST VALUES KEYS OUTDIR   -> OUTDIR/landis__registered_tables__chain.parquet
No scoring; plot level, restricted."""
import json, os, sys
from pathlib import Path
import numpy as np, pandas as pd
H = Path.home(); RC = H / "perseus_me/ncasi_rc"; CODE = H / "perseus_me/code"
for p in ("", "tools/score", "tools/hindcast", "tools/hindcast_all"):
    sys.path.insert(0, str(CODE / p))
cmd, st = sys.argv[1], sys.argv[2].upper(); os.environ.update(PME_STATE=st, PME_FIA_ROOT=str(H / f"restricted/ncasi-modeleval/fia_data_rc/{st}"))
import score as S  # noqa: E402
Hs = S.h1_setup(H / f"perseus_me/tracks/multistate/config/multistate/hindcast/{st}.yaml", H / f"perseus_me/multistate/inputs/{st}/restricted/hindcast")   # trkC_gain: multistate GA frames
ids = list(Hs["h1"]); O = Hs["O"]; o1 = O[O.visit == 1].set_index("chain_id"); Tx = Hs["Tx"]; t1 = Tx[(Tx.visit == 1) & Tx.live]
if cmd == "frame":
    out = Path(sys.argv[3]) / "sites"; out.mkdir(parents=True, exist_ok=True); os.chmod(out.parent, 0o700); os.chmod(out, 0o700)
    v1 = Hs["v1plt"].reindex(ids).astype(str).str.strip()
    fia = H / f"restricted/ncasi-modeleval/fia_data_rc/{st}"
    cd = pd.read_csv(fia / f"{st}_COND.csv", usecols=["PLT_CN", "CONDID", "CONDPROP_UNADJ", "STDAGE", "COND_STATUS_CD"], dtype={"PLT_CN": str})
    cd["PLT_CN"] = cd.PLT_CN.str.strip(); cd = cd[cd.COND_STATUS_CD == 1].sort_values(["PLT_CN", "CONDPROP_UNADJ", "CONDID"], ascending=[True, False, True]).drop_duplicates("PLT_CN")
    age = dict(zip(cd.PLT_CN, cd.STDAGE))
    have = sorted(set(t1.chain_id)); Sx = pd.DataFrame(dict(key=have)); Sx["site_id"] = np.arange(1, len(Sx) + 1)
    Sx["coord_cn"] = Sx.key.map(v1); Sx["STDAGE"] = Sx.coord_cn.map(age); Sx["REMPER"] = 5.0
    sid = dict(zip(Sx.key, Sx.site_id)); T = t1[t1.chain_id.isin(sid)].copy(); T["site_id"] = T.chain_id.map(sid)
    Sx[["site_id", "key"]].to_csv(out / "keys.csv", index=False); Sx.drop(columns=["key"]).to_csv(out / "sites.csv", index=False)
    T[["site_id", "SPCD", "DIA", "TPA_UNADJ", "DRYBIO_AG"]].to_csv(out / "trees.csv", index=False)
    for f in out.iterdir():
        os.chmod(f, 0o600)
    print(json.dumps(dict(state=st, chains=len(ids), chains_with_visit1_trees=len(have), trees=len(T))))
else:
    v = pd.read_csv(sys.argv[3]); k = pd.read_csv(sys.argv[4]); od = Path(sys.argv[5]); od.mkdir(parents=True, exist_ok=True); os.chmod(od, 0o700)
    key = dict(zip(k.site_id, k.key)); v["chain_id"] = v.site_id.map(key).astype(type(ids[0]))
    b0 = v[v.year == 0].set_index("chain_id").B_gm2; P = v[["chain_id", "year", "B_gm2"]].rename(columns={"year": "offset"})
    f = P.B_gm2 / P.chain_id.map(b0); P["agb_full"] = f * P.chain_id.map(o1.agb_stock); P["agc_full"] = f * P.chain_id.map(o1.agc_stock)
    built = set(P.chain_id); empty = [c for c in ids if c not in built]
    Z = pd.DataFrame([(c, o) for c in empty for o in range(0, 35, 5)], columns=["chain_id", "offset"]); Z["agb_full"] = 0.0; Z["agc_full"] = 0.0
    P = pd.concat([P.drop(columns=["B_gm2"]), Z], ignore_index=True); P["agb_surv"] = np.nan; P["agc_surv"] = np.nan
    P = P.sort_values(["chain_id", "offset"]).reset_index(drop=True); P.to_parquet(od / "landis__registered_tables__chain.parquet", index=False)
    os.chmod(od / "landis__registered_tables__chain.parquet", 0o600)
    m = dict(state=st, chains=len(ids), chains_built=len(built), chains_zero=len(empty), b0_zero=int((b0 <= 0).sum()), offsets=sorted(P.offset.unique().tolist()))
    (od / "landis_h1_meta.json").write_text(json.dumps(m, indent=1)); print(json.dumps(m))
