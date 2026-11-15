#!/usr/bin/env python3
"""Track trkB_engfix: Bakuzis 100-year no-harvest FVS-NE runs (copy of refine_engines2 bk_fvs.py; project.py, run_gate.py and
bk_common.py imported unchanged). One keyword block for every stand (--kw-file, may be empty), so no domain map is needed.
usage: b_bk.py --bin PATH --label LABEL --kw-file FILE --out DIR --restricted DIR [--workers N] [--max-stands N]"""
import argparse, json, os, sys, time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import pandas as pd, yaml
sys.path.insert(0, str(Path(__file__).resolve().parent))
import b_lib as BL
ROOT = BL.PME / "code"; HERE = ROOT / "tools" / "bakuzis"
for p in (HERE, ROOT / "engines" / "fvs", ROOT / "engines" / "treelist", ROOT / "workflow" / "scripts"):
    sys.path.insert(0, str(p))
import bk_common as B  # noqa: E402
import run_gate as G  # noqa: E402
import project as PJ  # noqa: E402


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--bin", required=True); ap.add_argument("--label", required=True)
    ap.add_argument("--kw-file", required=True); ap.add_argument("--out", required=True); ap.add_argument("--restricted", required=True)
    ap.add_argument("--workers", type=int, default=2); ap.add_argument("--max-stands", type=int, default=0); a = ap.parse_args(); t0 = time.time()
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True); rdir = B.restricted(a.restricted)
    H = yaml.safe_load((ROOT / "profiles/cardinal/host.yaml").read_text())["fvs"]
    H["bin_dir"] = str(Path(a.bin).parent); H["bins"] = {"ne": Path(a.bin).name}
    H["keyfile_builder_dir"] = str(BL.KF / "calibration/python")
    H["calibrated_tree"] = {"default": str(BL.KF), "ne": str(BL.KF)}; H.pop("extra_bins", None)
    aoi = json.loads((BL.PME / "me_state/aoi/aoi.json").read_text()); inv = BL.PME / "me_state/inventory"
    meta = dict(deviations=[], variants={}); F = G.build_frame(inv, aoi, H, 0, 20261029, meta); S = F["S"]
    groups, rel = B.plot_groups(inv, S.PLT_CN.unique())
    samp = set(groups[groups.sampled].PLT_CN); Ss = S[S.PLT_CN.isin(samp)]
    (out / "groups_release.json").write_text(json.dumps(rel, indent=1)); B.fia_percentiles(F, groups).to_csv(out / "fia_percentiles.csv", index=False)
    extra = Path(a.kw_file).read_text().strip()
    keep, kwfun, used = G.arm_setup(["default"], F, H, meta, ncyc=PJ.NCYC)
    Fs = dict(F); Fs["states"] = {st: dict(d, cd=d["cd"][d["cd"].PLT_CN.isin(samp)]) for st, d in F["states"].items()}
    kf = (lambda arm, v, s, y: (kwfun(arm, v, s, y) + "\n" + extra).strip())
    jobs = [dict(j, ncyc=PJ.NCYC) for j in G.stand_jobs(Fs, H, "default", kf, ncyc=PJ.NCYC)]
    assert len({j["stand_id"] for j in jobs}) == len(jobs) == len(Ss), (len(jobs), len(Ss))
    if a.max_stands:
        jobs = jobs[: a.max_stands]
    with ProcessPoolExecutor(a.workers) as ex:
        R = list(ex.map(PJ.fvs_worker, jobs, chunksize=2))
    m = dict(label=a.label, binary=a.bin, binary_md5=BL.md5(a.bin), kw_md5=BL.md5(a.kw_file), stands=len(jobs),
             stands_ok=int(sum(r["ok"] for r in R)), failures=pd.Series([r["note"] for r in R if not r["ok"]]).value_counts().head(8).to_dict(),
             wall_s=round(time.time() - t0, 1))
    if not a.max_stands:
        V, nbad = B.plot_values(R, Ss, groups)
        f = rdir / f"plot_values__fvs__{a.label}.parquet"; V.to_parquet(f, index=False); os.chmod(f, 0o600)
        B.group_means(V, "fvs", a.label).to_csv(out / f"traj__fvs__{a.label}.csv", index=False); m["plots_with_values"] = int(V.PLT_CN.nunique())
    else:
        m["m_first"] = [r.get("m", [None])[-1] for r in R[:2]]
    (out / f"meta__fvs__{a.label}.json").write_text(json.dumps(m, indent=1, default=str)); print(json.dumps(m, default=str))


if __name__ == "__main__":
    main()
