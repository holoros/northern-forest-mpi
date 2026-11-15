"""Track trkB_engfix probe: FVS-NE pair runs that keep their work dir, to check (a) Event Monitor PARMS in FIXMORT
and (b) BSDI / BSDIMAX from COMPUTE. Restricted output only."""
import sys, os
sys.path.insert(0, "/home/aaron/perseus_me/tracks/refine_engines2/tools")
import fvs_pairs as FP
PG = FP.PG
PG.shutil.rmtree = lambda *a, **k: None
scratch = sys.argv[2]; os.makedirs(scratch, exist_ok=True)
P, T, coef = FP.load_set("c2"); P["pid"] = FP.strip(P.pair_id)
tg = dict(tuple(T.groupby("pair_id")))
S = P.iloc[:2]
blk = open(sys.argv[1]).read()
jobs = PG.make_jobs("ne", S, tg, coef, ["default"], lambda arm, iy, cyc: blk.replace("{Y}", str(iy)).replace("{L}", str(cyc)),
                    FP.BIN["ne"], FP.FVS_MODERN / "calibration/python", scratch)
for j in jobs:
    r = PG.run_pair(j); print({k: r.get(k) for k in ("ok", "fvs_mort", "tpa_start", "tpa_end", "note")}, "cycle", j["cycle"])
