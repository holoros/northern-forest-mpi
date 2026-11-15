"""trkA release assembly: copies aggregate tables into release/, builds trkA_checks.json, claims.csv and runs the G5 leak scan.
usage: trkA_release.py"""
import json, re, shutil
from pathlib import Path
import numpy as np, pandas as pd
T = Path("/home/aaron/perseus_me/tracks/trkA_hcm"); W = T / "work"; REL = T / "release"; REL.mkdir(exist_ok=True)
for f in ("trkA_layer_metrics.csv", "trkA_layer_paired.csv", "trkA_layer_reliability.csv", "trkA_promotion.json"):
    shutil.copy(W / "eval" / f, REL / f)
for f, g in (("imf2_scores.csv", "trkA_imf2_scores.csv"), ("imf2_prank.csv", "trkA_imf2_prank.csv"), ("imf2_layer_choice.csv", "trkA_imf2_layer_choice.csv"),
             ("imf2_waybill_counts.csv", "trkA_imf2_waybill_counts.csv")):
    shutil.copy(W / "main" / f, REL / g)
st = [json.load(open(W / f"state_{s}_hcm2.json")) for s in ("NH", "VT", "MI", "WI", "MN") if (W / f"state_{s}_hcm2.json").exists()]
S = pd.DataFrame([dict(state=d["state"], pairs=d["pairs"], cut=d["cut"], counties=d["counties"], price_used=d["price_used"], price_fallback_share=d["price_fallback_share"],
                       AUC=d["AUC"], AUC_lo=d["AUC_ci"][0], AUC_hi=d["AUC_ci"][1], brier=d["brier"], brier_lo=d["brier_ci"][0], brier_hi=d["brier_ci"][1],
                       cal_slope=d["cal_slope"], cal_slope_lo=d["cal_slope_ci"][0], cal_slope_hi=d["cal_slope_ci"][1], E_over_O=d["E_over_O"]) for d in st])
S.to_csv(REL / "trkA_states_layer.csv", index=False, float_format="%.6g")
os_, ow, om = [], [], {}
for s in ("NH", "VT", "MI", "WI", "MN"):
    d = W / f"oos_{s}/main"
    if (d / "oos_scores.csv").exists():
        os_.append(pd.read_csv(d / "oos_scores.csv")); ow.append(pd.read_csv(d / "oos_wins.csv")); om[s] = json.load(open(d / "oos_meta.json"))
pd.concat(os_).to_csv(REL / "trkA_oos_scores.csv", index=False, float_format="%.6g"); pd.concat(ow).to_csv(REL / "trkA_oos_wins.csv", index=False, float_format="%.6g")
json.dump({s: {k: v for k, v in m.items() if k != "inputs_found"} for s, m in om.items()}, open(REL / "trkA_oos_meta.json", "w"), indent=1)
meta = json.load(open(W / "main/imf2_meta.json")); g2 = json.load(open(W / "g2_checks.json")); bv = json.load(open(W / "builder_validation.json"))
g1 = json.load(open(W / "g1_check.json")); promo = json.load(open(W / "eval/trkA_promotion.json"))
chk = dict(prereg_md5=open(T / "PREREG.md5").read().split()[0], G1=g1, G2=g2, G4_pass=meta["G4_pass"], builder_validation=bv, label_permutation=promo["perm_check"])
# G5 leak scan
bad = []
for f in REL.iterdir():
    if f.suffix == ".csv":
        x = pd.read_csv(f)
        if any(re.search(r"(^CN$|PLT_CN|pair_id|^loc$|LAT|LON|^PLOT$)", c) for c in x.columns):
            bad.append(f"{f.name}: id column")
        for c in x.select_dtypes("number").columns:
            v = x[c].dropna().values
            if len(v) and np.any((np.abs(v) >= 1e9) & (np.mod(v, 1) == 0)):
                bad.append(f"{f.name}:{c} 10+ digit integer")
    if f.suffix in (".json", ".md"):
        if re.search(r"(?<![\d.])\d{10,}", f.read_text()):
            bad.append(f"{f.name}: 10+ digit integer")
chk["G5_leak_scan"] = dict(pass_=not bad, findings=bad)
json.dump(chk, open(REL / "trkA_checks.json", "w"), indent=1, default=str); print(json.dumps(chk, default=str)[:3000])
g3 = dict(layer=json.load(open(W / "g3_second_layer.json")), imf2=json.load(open(W / "main/g3_second_path.json")))
chk["G3_second_code_path"] = dict(layer_pass=g3["layer"]["pass"], layer_max_abs_diff=g3["layer"]["max_abs_diff"], imf2_pass=g3["imf2"]["G3_pass"], imf2_max_abs_diff=g3["imf2"]["max_abs_diff"],
                                  note="second code path written in this session; the brief's fresh-context verifier subagent could not be launched (no agent tool in this session)")
json.dump(chk, open(REL / "trkA_checks.json", "w"), indent=1, default=str)
