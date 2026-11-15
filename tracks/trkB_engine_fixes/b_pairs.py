#!/usr/bin/env python3
"""Track trkB_engfix: one-cycle FVS-NE pair runs (pairs.py make_jobs and run_pair, imported unchanged) with per-group keyword
blocks; optional capture of the Event Monitor COMPUTE table (XBSDI, XBSMX) from the run's database before the work dir is removed.
spec JSON {"jobs": [{"tag": str, "kw": str, "ids": [pair_id, ...] or null}]}.
usage: b_pairs.py --set c2|c5 --spec spec.json --out DIR --bin FVSne [--workers N]"""
import argparse, json, os, sys, time, sqlite3
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parent))
import b_lib as BL
for p in (BL.PME / "code" / "engines" / "fvs", BL.PME / "code" / "workflow" / "scripts"):
    sys.path.insert(0, str(p))
import pairs as PG  # noqa: E402
_rm = PG.shutil.rmtree
_CAP = {}


def _rm_capture(path, *a, **k):
    db = Path(path) / "fvsout.db"
    if db.exists():
        try:
            con = sqlite3.connect(str(db)); t = [r[0] for r in con.execute("select name from sqlite_master where type='table'")]
            if "FVS_Compute" in t:
                x = pd.read_sql("select * from FVS_Compute", con)
                _CAP["c"] = {c: float(x[c].iloc[0]) for c in x.columns if c.upper().startswith("XB")}
            con.close()
        except Exception as e:  # recorded
            _CAP["c"] = {"err": str(e)[:80]}
    return _rm(path, *a, **k)


def run_x(job):
    _CAP.clear(); PG.shutil.rmtree = _rm_capture
    r = PG.run_pair(job); r.update({k.lower(): v for k, v in _CAP.get("c", {}).items()}); return r


def strip(x):
    return pd.Series(x).astype(str).str.strip().values


def load_set(which):
    PB = pd.read_csv(BL.PME / "pairs/me_state_std/pairs_built.csv", dtype={"pair_id": str, "t1_id": str})
    PBp = PG.prepare(PB); si_med, age_med = float(PBp.si_med.iloc[0]), float(PBp.age_med.iloc[0])
    if which == "c2":
        P = PBp[PBp.half == "hold"].drop(columns=["si_med", "age_med", "cls", "half"])
        T = pd.read_csv(BL.PME / "pairs/me_state_std/trees_built.csv", dtype={"pair_id": str, "tree_cn": str})
    else:
        R = BL.PME / "me_state/v23_D/restricted/build"
        P = pd.read_csv(R / "treated_pairs.csv", dtype={"pair_id": str, "t1_id": str})
        T = pd.read_csv(R / "treated_trees.csv", dtype={"pair_id": str, "tree_cn": str})
    P = P.copy(); P["si_med"] = si_med; P["age_med"] = age_med
    T = T[T.pair_id.isin(set(P.pair_id))]
    coef = {int(k): tuple(v) for k, v in json.loads((BL.PME / "pairs/me_state_std/emulator_built.json").read_text()).items()}
    return P, T, coef


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--set", required=True, choices=["c2", "c5"]); ap.add_argument("--spec", required=True)
    ap.add_argument("--out", required=True); ap.add_argument("--bin", required=True); ap.add_argument("--workers", type=int, default=2)
    a = ap.parse_args(); t0 = time.time()
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True); os.chmod(out, 0o700)
    scratch = os.environ.get("TMPDIR", str(out / "tmp")); os.makedirs(scratch, exist_ok=True)
    P, T, coef = load_set(a.set); P["pid"] = strip(P.pair_id)
    tg = dict(tuple(T.groupby("pair_id"))); spec = json.loads(Path(a.spec).read_text()); jobs = []
    for js in spec["jobs"]:
        ids = set(str(x).strip() for x in js["ids"]) if js.get("ids") is not None else set(P.pid)
        S = P[P.pid.isin(ids)]
        assert len(S) == len(ids), f"{js['tag']}: {len(ids) - len(S)} ids are not {a.set} pairs"
        kw = js["kw"]
        jj = PG.make_jobs("ne", S, tg, coef, ["default"], lambda arm, iy, cyc, kw=kw: kw, a.bin, BL.KF / "calibration/python", scratch)
        for j in jj:
            j["tag"] = js["tag"]
        jobs += jj
    with ProcessPoolExecutor(a.workers) as ex:
        R = pd.DataFrame(list(ex.map(run_x, jobs, chunksize=8)))
    R["tag"] = [j["tag"] for j in jobs]; R["cycle"] = [j["cycle"] for j in jobs]; R["arm"] = "default"
    PG.write_ids_csv(R, out / "pair_rows.csv"); os.chmod(out / "pair_rows.csv", 0o600)
    meta = dict(set=a.set, binary=a.bin, binary_md5=BL.md5(a.bin), spec_md5=BL.md5(a.spec), jobs=len(jobs),
                ok=int(R.ok.astype(bool).sum()), wall_s=round(time.time() - t0, 1),
                failures=R.loc[~R.ok.astype(bool), "note"].astype(str).str[:80].value_counts().head(5).to_dict())
    (out / "run_meta.json").write_text(json.dumps(meta, indent=1, default=str)); print(json.dumps(meta, default=str))


if __name__ == "__main__":
    main()
