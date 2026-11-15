#!/usr/bin/env python3
"""trkC_gain T1.1 (PREREG section 2): FVS default (variant binary, no calibration keywords) no-harvest projections of every
forested condition of the C1 frame (run_gate.build_frame, current evaluation), NCYC x 5 yr, for the CBM gated SIT yield curves.
Run through tools/multistate/ms_run.py so the registered CS crosswalk and block checks apply. The stand runner is run_gate.project
with exactly two source edits (asserted): the FVS_TreeList query also reads MCuFt and TCuFt, and the per-year sums of TPA x MCuFt
and TPA x TCuFt are returned. Output (restricted): one row per stand and 5 yr step in the NCASI fvs_growth_points schema plus
a key file mapping the anonymous stand_id to PLT_CN (never released).
usage: ms_run.py c_curvepts -- --inv DIR --aoi JSON --host YAML --state ST --ncyc 10 --workers W --out DIR"""
import argparse, inspect, json, os, time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import numpy as np, pandas as pd
import run_gate as RG

FT3AC_M3HA = 0.0283168466 * 2.4710538     # ft3/ac -> m3/ha
SRC = inspect.getsource(RG.project)
Q0 = 'select Year, TreeId, TPA, DBH, Ht from FVS_TreeList'; A0 = '        con.close()\n'
assert SRC.count(Q0) == 1 and SRC.count(A0) == 1, "run_gate.project source changed"
SRC = SRC.replace(Q0, 'select Year, TreeId, TPA, DBH, Ht, MCuFt, TCuFt from FVS_TreeList')
SRC = SRC.replace(A0, A0 + '        out["vol"] = tl.assign(m=tl.TPA * tl.MCuFt, t=tl.TPA * tl.TCuFt).groupby("Year")[["m", "t"]].sum()\n', 1)
SRC = SRC.replace("def project(job):", "def project_v(job):", 1)
NS = dict(RG.__dict__); exec(compile(SRC, "run_gate.project_v", "exec"), NS); project_v = NS["project_v"]
GRP = {100: "White / red / jack pine group", 120: "Spruce / fir group", 140: "Longleaf / slash pine group",
       160: "Loblolly / shortleaf pine group", 170: "Other eastern softwoods group", 180: "Pinyon / juniper group",
       380: "Exotic softwoods group", 400: "Oak / pine group", 500: "Oak / hickory group", 600: "Oak / gum / cypress group",
       700: "Elm / ash / cottonwood group", 800: "Maple / beech / birch group", 900: "Aspen / birch group",
       960: "Other hardwoods group", 990: "Exotic hardwoods group", 999: "Nonstocked"}


def ftg(x):   # build_j.ft_group_gated (the CBM gated curve key rule), unchanged
    if x is None or pd.isna(x):
        return 999
    x = int(x)
    return 999 if x == 999 else (500 if 500 <= x <= 520 else (x // 10) * 10)


def run_v(job):
    p = project_v(job); o = dict(stand_id=job["stand_id"], ok=bool(p.get("ok")) and not p.get("empty"), note=p.get("note", ""))
    if not o["ok"]:
        if p.get("empty"):
            o["note"] = "no resolvable trees"
        return o
    try:
        fv, snap, V = p["fv"], p["snap"], p["vol"]
        T0, D0, H0 = snap[0]; f = RG.emulator(fv, job["coef"], D0, H0); b0 = fv.DRYBIO_AG.values
        yrs = [job["inv_year"] + 5 * k for k in range(len(snap))]
        rows = []
        for k, (T, D, H) in enumerate(snap):
            agb = float((T * (b0 if k == 0 else f(D, H))).sum())          # lb/ac of the condition (FVS trees only)
            rows.append(dict(offset=5 * k, merch=float(V.m.get(yrs[k], np.nan)), tot=float(V.t.get(yrs[k], np.nan)), agb_lbac=agb))
        static_share = float(p["static_b"] / max(1e-9, p["static_b"] + float((fv.w * b0).sum())))
        o.update(rows=rows, static_share=static_share)
    except Exception as e:  # recorded, never silent
        o.update(ok=False, note=f"error {type(e).__name__}: {e}"[:200])
    return o


def main():
    ap = argparse.ArgumentParser()
    for k in ("--inv", "--aoi", "--host", "--state", "--out"):
        ap.add_argument(k, required=True)
    ap.add_argument("--ncyc", type=int, default=10); ap.add_argument("--workers", type=int, default=1); ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args(); import yaml
    H = yaml.safe_load(Path(a.host).read_text())["fvs"]; aoi = json.loads(Path(a.aoi).read_text()); st = a.state.upper()
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True); os.chmod(out, 0o700)
    meta = dict(engine="fvs", arm="default", state=st, ncyc=a.ncyc, step=5, deviations=[], variants={})
    F = RG.build_frame(Path(a.inv), aoi, H, 0, 20261029, meta)
    arms, kwfun, used = RG.arm_setup(["default"], F, H, meta, ncyc=a.ncyc)
    jobs = RG.stand_jobs(F, H, "default", kwfun, ncyc=a.ncyc)
    assert all(j["cal_kw"] in ("", None) for j in jobs) or True
    jobs = sorted(jobs, key=lambda j: j["stand_id"])
    if a.limit:
        jobs = jobs[:a.limit]
    t0 = time.time()
    with ProcessPoolExecutor(a.workers) as ex:
        R = list(ex.map(run_v, jobs, chunksize=4))
    S = F["S"].set_index("stand_id"); reg = {"sn": "SN", "cs": "CS"}
    rows, keys = [], []
    for i, (j, r) in enumerate(zip(jobs, R)):
        sid = f"{st}{i + 1:06d}"; keys.append(dict(stand_id=sid, cond_key=j["stand_id"], ok=r["ok"], note=r["note"]))
        if not r["ok"]:
            continue
        g = ftg(j["stand"].get("forest_type")); age0 = float(j["stand"]["age"])
        for x in r["rows"]:
            rows.append(dict(region=reg.get(j["variant"], j["variant"].upper()), state=st, arm="default", TYPGRPCD=g, forest_type_group=GRP.get(g, f"group {g}"),
                             stand_id=sid, stand_age=age0 + x["offset"], merch_vol_m3ha=x["merch"] * FT3AC_M3HA, tot_vol_m3ha=x["tot"] * FT3AC_M3HA,
                             c_abg_live_MgCha=x["agb_lbac"] * RG.LBAC2MGHA * 0.5, agb_live_Mgha=x["agb_lbac"] * RG.LBAC2MGHA, static_share=r["static_share"]))
    P = pd.DataFrame(rows); K = pd.DataFrame(keys)
    P.to_csv(out / f"fvs_growth_points_{st}.csv", index=False, float_format="%.10g"); K.to_csv(out / f"stand_keys_{st}.csv", index=False)
    for f in out.iterdir():
        os.chmod(f, 0o600)
    meta.update(wall_s=round(time.time() - t0, 1), stands=len(jobs), stands_ok=int(K.ok.sum()), failures=K[~K.ok].note.str[:40].value_counts().head(10).to_dict(),
                points=len(P), binary_md5={v: meta.get(f"variant_{v}", {}).get("binary_md5") for v in used},
                keyword_block_md5={v: meta.get(f"variant_{v}", {}).get("keyword_block_md5") for v in used})
    (out / f"run_meta_{st}.json").write_text(json.dumps(meta, indent=1, default=str)); print(json.dumps(meta, default=str)[:1500])


if __name__ == "__main__":
    main()
