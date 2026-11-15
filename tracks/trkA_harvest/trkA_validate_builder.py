"""trkA section 8 builder validation: the state builder run on the Maine DataMart must reproduce pairs_h (pairs, cut, share_ba) and the
track HCM / track SC features used by HCM2. Writes work/builder_validation.json (aggregates)."""
import json
from pathlib import Path
import numpy as np, pandas as pd
R = Path("/home/aaron/perseus_me/me_state"); T = Path("/home/aaron/perseus_me/tracks/trkA_hcm")
B = pd.read_parquet(T / "restricted/states/ME/frame.parquet").set_index("CN")
P = pd.read_parquet(R / "v26_H/restricted/frame/pairs_h.parquet").set_index("CN")
H = pd.read_parquet(R / "v26_HCM/restricted/hcm_features.parquet").set_index("CN")
S = pd.read_parquet(R / "v26_SC/restricted/features.parquet").set_index("CN")
out = dict(pairs_builder=len(B), pairs_h=len(P), same_keys=bool(set(B.index) == set(P.index)))
B = B.reindex(P.index)
out["cut_mismatch"] = int((B.cut.astype(bool) != P.cut.astype(bool)).sum()); out["h2_mismatch"] = int((B.h2_eligible.astype(bool) != P.h2_eligible.astype(bool)).sum())
def md(a, b):
    a = np.asarray(a, float); b = np.asarray(b, float); both = np.isfinite(a) & np.isfinite(b)
    return dict(max_abs=float(np.max(np.abs(a[both] - b[both]))) if both.any() else None, nan_mismatch=int((np.isfinite(a) != np.isfinite(b)).sum()))
out["share_ba"] = md(B.share_ba, P.share_ba)
for c in ("lagbio", "saw_t1", "pld_t1", "saw_t2", "pld_t2", "share_saw", "share_pld"):
    out[c] = md(B[c], H[c].reindex(P.index))
for c in ("ELEV", "RDDISTCD", "OWNCD", "OWNGRPCD", "RESERVCD", "volbf", "volcf", "prev_trt10", "yrs_since_cut", "prev_pair_cut"):
    out[c] = md(B[c], S[c].reindex(P.index))
json.dump(out, open(T / "work/builder_validation.json", "w"), indent=1); print(json.dumps(out))
