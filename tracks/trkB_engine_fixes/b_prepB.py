"""Track trkB_engfix: job B inputs. Pair specs (test-fold parameters per pair; S1 'all'), Bakuzis keyword file ('all', 20 cycles),
hindcast chain files by visit 1 county fold (county not in the fold map -> 'all'). Fold map = refine_engines2 frames (county_folds)."""
import json, pickle, sys
from pathlib import Path
import pandas as pd
sys.path.insert(0, "/home/aaron/perseus_me/tracks/trkB_engfix/tools"); import b_lib as BL
sys.path.insert(0, "/home/aaron/perseus_me/tracks/refine_engines/tools"); import re_lib as R
T = BL.PME / "tracks/trkB_engfix"; O = T / "restricted/xferB"; O.mkdir(exist_ok=True)
P = json.loads((T / "work/ddm_params.json").read_text())
fc = pickle.load(open(BL.PME / "tracks/refine_engines2/restricted/frames_cache.pkl", "rb")); C2, C5 = fc["C2"], fc["C5"]
cc = pd.concat([pd.Series(C2.cty.values, index=C2.pair_id.values), pd.Series(C5.cty.values, index=C5.pair_id.values)]); cc = cc[~cc.index.duplicated()]
fm = R.county_folds(cc.value_counts().to_dict())
for F in (C2, C5):
    assert (pd.Series(F.cty.values).map(fm).values == F.fold.values).all()
for s, F in (("c2", C2), ("c5", C5)):
    jobs = [dict(tag=f"ddm_f{k}", kw=BL.ddm_block(P[str(k)]["a"], P[str(k)]["b"], 1), ids=F.pair_id[F.fold.values == k].tolist()) for k in range(5)]
    assert sum(len(j["ids"]) for j in jobs) == len(F)
    json.dump(dict(jobs=jobs), open(O / f"spec_{s}_ddm.json", "w"))
    json.dump(dict(jobs=[dict(tag="ddm_all", kw=BL.ddm_block(P["all"]["a"], P["all"]["b"], 1), ids=F.pair_id.tolist())]), open(O / f"spec_{s}_ddm_all.json", "w"))
(O / "kw_bk_ddm.txt").write_text(BL.ddm_block(P["all"]["a"], P["all"]["b"], 20) + "\n"); (O / "kw_bk_none.txt").write_text("\n")
json.dump(P, open(O / "ddm_params.json", "w"))
H = BL.PME / "me_state/hindcast/restricted"
V = pd.read_parquet(H / "chains.parquet"); pl = pd.read_parquet(H / "inventory_1999_2003/plots.parquet")
v1 = V[V.visit == 1].drop_duplicates("chain_id")
cty = v1.PLT_CN.map(dict(zip(pl.CN, pl.COUNTYCD)))
fold = cty.map(fm)
cnt = {}
for k in list(range(5)) + ["all"]:
    ids = v1.chain_id[fold == k] if k != "all" else v1.chain_id[fold.isna()]
    (O / f"chains_f{k}.txt").write_text("\n".join(str(x) for x in ids) + "\n"); cnt[str(k)] = int(len(ids))
json.dump(dict(chain_fold_counts=cnt, county_missing_plot=int(cty.isna().sum())), open(T / "work/chain_folds.json", "w"), indent=1)
print(cnt, int(cty.isna().sum()), {s: len(F) for s, F in (("c2", C2), ("c5", C5))})
