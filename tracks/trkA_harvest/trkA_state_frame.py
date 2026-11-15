#!/usr/bin/env python3
"""trkA_hcm (PREREG md5 35055a79af9083edec3f5e3223552dd7) section 8: state harvest frame (track H pair and H1 rules, h_observe.py) plus the
HCM2 start-of-cycle features (track HCM hcm_features.py and track SC sc_features.py rules) from a public FIADB DataMart folder.
usage: trkA_state_frame.py ST FIA_DIR OUT_PARQUET    No coordinate column is read."""
import json, os, sys
from pathlib import Path
import numpy as np, pandas as pd

BAF = 0.005454154; FT2AC_M2HA = 0.229568411; LBAC_MGHA = 0.45359237e-3 * 2.4710538
REF = "/home/aaron/perseus_me/ext/FIA/ENTIRE_REF_FOREST_TYPE.csv"


def s(d, cols):
    for c in cols:
        d[c] = d[c].astype(str).str.strip().replace({"nan": np.nan})
    return d


def rd(f, cols, keep=None, key=None, dtype=None):
    out = []
    for ch in pd.read_csv(f, usecols=cols, dtype=dtype, chunksize=2_000_000, low_memory=False):
        if keep is not None:
            ch[key] = ch[key].astype(str).str.strip(); ch = ch[ch[key].isin(keep)]
        out.append(ch)
    return pd.concat(out, ignore_index=True)


def main():
    ST, F, OUT = sys.argv[1], Path(sys.argv[2]), Path(sys.argv[3])
    pl = s(pd.read_csv(F / f"{ST}_PLOT.csv", usecols=["CN", "PREV_PLT_CN", "INVYR", "MEASYEAR", "REMPER", "STATECD", "UNITCD", "COUNTYCD", "PLOT",
                                                       "PLOT_STATUS_CD", "ELEV", "RDDISTCD"], dtype={"CN": str, "PREV_PLT_CN": str}, low_memory=False), ["CN", "PREV_PLT_CN"])
    pl["loc"] = pl.STATECD.astype(str) + "_" + pl.UNITCD.astype(str) + "_" + pl.COUNTYCD.astype(str) + "_" + pl.PLOT.astype(str)
    cd = s(pd.read_csv(F / f"{ST}_COND.csv", usecols=["PLT_CN", "CONDID", "COND_STATUS_CD", "CONDPROP_UNADJ", "OWNCD", "OWNGRPCD", "FORTYPCD", "STDAGE",
                                                       "TRTCD1", "TRTCD2", "TRTCD3", "TRTYR1", "TRTYR2", "TRTYR3", "RESERVCD"], dtype={"PLT_CN": str}, low_memory=False), ["PLT_CN"])
    ref = pd.read_csv(REF, usecols=["VALUE", "TYPGRPCD"]); grp = dict(zip(ref.VALUE, ref.TYPGRPCD))
    fc = cd[cd.COND_STATUS_CD == 1]; fprop = fc.groupby("PLT_CN").CONDPROP_UNADJ.sum()
    dom = fc.assign(_n=-fc.CONDPROP_UNADJ).sort_values(["PLT_CN", "_n", "CONDID"]).groupby("PLT_CN").head(1).set_index("PLT_CN")
    cns = set(pl.CN); t1i = pl.set_index("CN"); C = {}
    pr = pl[pl.PREV_PLT_CN.notna() & pl.PREV_PLT_CN.isin(cns) & (pl.REMPER > 0) & pl.INVYR.between(2004, 2025)].copy(); C["remeasured"] = len(pr)
    pr["INVYR1"] = pr.PREV_PLT_CN.map(t1i.INVYR); pr = pr[pr.INVYR1 >= 1999]
    pr["fprop1"] = pr.PREV_PLT_CN.map(fprop).fillna(0.0); pr["fprop2"] = pr.CN.map(fprop).fillna(0.0)
    pr = pr[(pr.fprop1 > 0) & pr.PLOT_STATUS_CD.isin([1, 2])].copy(); C["frame"] = len(pr)
    k = pr.PREV_PLT_CN.values
    for a_, b_ in (("owngrp1", "OWNGRPCD"), ("fortyp1", "FORTYPCD"), ("stdage1", "STDAGE"), ("reserv1", "RESERVCD"), ("OWNCD", "OWNCD"), ("OWNGRPCD", "OWNGRPCD"), ("RESERVCD", "RESERVCD")):
        pr[a_] = dom[b_].reindex(k).values
    pr["ftgrp1"] = pr.fortyp1.map(grp)
    pr["ELEV"] = t1i.ELEV.reindex(k).values; pr["RDDISTCD"] = t1i.RDDISTCD.reindex(k).values
    c2 = cd[cd.PLT_CN.isin(set(pr.CN))]
    pr["trt10"] = pr.CN.map(c2[["TRTCD1", "TRTCD2", "TRTCD3"]].eq(10).any(axis=1).groupby(c2.PLT_CN).any()).fillna(False).astype(bool)
    c1 = cd[cd.PLT_CN.isin(set(k))]; trt = c1[["TRTCD1", "TRTCD2", "TRTCD3"]]; yrs = c1[["TRTYR1", "TRTYR2", "TRTYR3"]]
    pr["prev_trt10"] = trt.eq(10).any(axis=1).groupby(c1.PLT_CN).any().reindex(k).fillna(False).astype(int).values
    ty = np.where(trt.eq(10).values & (yrs.values > 1900) & (yrs.values < 9000), yrs.values, np.nan).astype(float)
    last = pd.Series(np.nanmax(np.where(np.isnan(ty), -np.inf, ty), axis=1), index=c1.index).replace(-np.inf, np.nan).groupby(c1.PLT_CN).max()
    pr["yrs_since_cut"] = pr.INVYR1.values - last.reindex(k).values
    keep = set(pr.CN) | set(k)
    tr = s(rd(F / f"{ST}_TREE.csv", ["CN", "PLT_CN", "PREV_TRE_CN", "CONDID", "STATUSCD", "DIA", "TPA_UNADJ", "DRYBIO_AG", "DRYBIO_SAWLOG", "VOLCFNET", "VOLBFNET"],
                keep=keep, key="PLT_CN", dtype={"CN": str, "PLT_CN": str, "PREV_TRE_CN": str}), ["CN", "PREV_TRE_CN"])
    gc = s(rd(F / f"{ST}_TREE_GRM_COMPONENT.csv", ["TRE_CN", "PREV_TRE_CN", "PLT_CN", "MICR_COMPONENT_AL_FOREST"], keep=set(pr.CN), key="PLT_CN",
              dtype={"TRE_CN": str, "PREV_TRE_CN": str, "PLT_CN": str}), ["TRE_CN", "PREV_TRE_CN"])
    comp = gc.MICR_COMPONENT_AL_FOREST.fillna("").astype(str)
    pr["grm_cut"] = pr.CN.isin(set(gc.PLT_CN[comp.str.startswith("CUT")]))
    t2 = tr[tr.PLT_CN.isin(set(pr.CN))]; pr["st3"] = pr.CN.isin(set(t2.PLT_CN[t2.STATUSCD == 3]))
    cut_t2 = set(t2.CN[t2.STATUSCD == 3]) | set(gc.TRE_CN[comp.str.startswith("CUT")])
    cut_t1 = set(t2.PREV_TRE_CN[t2.CN.isin(cut_t2)].dropna()) | set(gc.PREV_TRE_CN[comp.str.startswith("CUT")].dropna())
    live = tr[(tr.STATUSCD == 1) & (tr.DIA >= 1) & tr.TPA_UNADJ.notna()].merge(fc[["PLT_CN", "CONDID"]].drop_duplicates(), on=["PLT_CN", "CONDID"], how="inner").copy()
    live["ba"] = BAF * live.DIA ** 2 * live.TPA_UNADJ; live["tot"] = live.DRYBIO_AG.fillna(0.0) * live.TPA_UNADJ; live["saw"] = live.DRYBIO_SAWLOG.fillna(0.0) * live.TPA_UNADJ
    live["vcf"] = live.VOLCFNET.fillna(0) * live.TPA_UNADJ; live["vbf"] = live.VOLBFNET.fillna(0) * live.TPA_UNADJ; live["cut"] = live.CN.isin(cut_t1)
    g = live.groupby("PLT_CN")[["ba", "tot", "saw", "vcf", "vbf"]].sum(); gcut = live[live.cut].groupby("PLT_CN")[["ba", "tot", "saw"]].sum()
    ncut = live[live.cut].groupby("PLT_CN").size()
    pr["ba1"] = g.ba.reindex(k).fillna(0).values; pr["ba_cut"] = gcut.ba.reindex(k).fillna(0).values; pr["n_cut_t1"] = ncut.reindex(k).fillna(0).values
    f1 = np.maximum(pr.fprop1.values, 0.1); f2 = np.maximum(pr.fprop2.values, 0.1)
    pr["lagbio"] = g.tot.reindex(k).fillna(0).values * LBAC_MGHA / f1; pr["saw_t1"] = g.saw.reindex(k).fillna(0).values * LBAC_MGHA / f1; pr["pld_t1"] = pr.lagbio - pr.saw_t1
    tot2 = g.tot.reindex(pr.CN.values).fillna(0).values * LBAC_MGHA / f2; pr["saw_t2"] = g.saw.reindex(pr.CN.values).fillna(0).values * LBAC_MGHA / f2; pr["pld_t2"] = tot2 - pr.saw_t2
    sc = gcut.saw.reindex(k).fillna(0).values * LBAC_MGHA / f1; pc = (gcut.tot - gcut.saw).reindex(k).fillna(0).values * LBAC_MGHA / f1
    pr["share_saw"] = np.where(pr.saw_t1 > 0, sc / np.where(pr.saw_t1 > 0, pr.saw_t1, 1), np.nan); pr["share_pld"] = np.where(pr.pld_t1 > 0, pc / np.where(pr.pld_t1 > 0, pr.pld_t1, 1), np.nan)
    pr["volcf"] = g.vcf.reindex(k).fillna(0).values; pr["volbf"] = g.vbf.reindex(k).fillna(0).values
    pr["cut"] = pr.grm_cut | pr.st3 | pr.trt10
    pr["h2_eligible"] = pr.cut & (pr.n_cut_t1 > 0) & (pr.ba1 > 0)
    pr["share_ba"] = np.where(pr.ba1 > 0, pr.ba_cut / pr.ba1.where(pr.ba1 > 0, 1), np.nan)
    prev = pr.set_index("CN").cut
    pr["prev_pair_cut"] = pd.Series(k).map(prev).map({True: 1.0, False: 0.0}).values
    cols = ["CN", "PREV_PLT_CN", "loc", "INVYR", "INVYR1", "MEASYEAR", "REMPER", "COUNTYCD", "UNITCD", "fprop1", "owngrp1", "ftgrp1", "stdage1", "reserv1",
            "OWNCD", "OWNGRPCD", "RESERVCD", "ELEV", "RDDISTCD", "prev_trt10", "yrs_since_cut", "prev_pair_cut", "trt10", "grm_cut", "st3", "cut", "h2_eligible",
            "n_cut_t1", "ba1", "share_ba", "lagbio", "saw_t1", "pld_t1", "saw_t2", "pld_t2", "share_saw", "share_pld", "volcf", "volbf"]
    P = pr[cols].sort_values("CN").reset_index(drop=True); OUT.parent.mkdir(parents=True, exist_ok=True); os.chmod(OUT.parent, 0o700)
    P.to_parquet(OUT, index=False); os.chmod(OUT, 0o600)
    C.update(state=ST, pairs=len(P), cut=int(P.cut.sum()), cut_rate=float(P.cut.mean()), counties=int(P.COUNTYCD.nunique()), h2=int(P.h2_eligible.sum()))
    print(json.dumps(C))


if __name__ == "__main__":
    main()
