"""Track R realism: registered P3 (full Bakuzis, tools/bakuzis/score_v2.R, noR5 headline, 0.5 coverage rule), ensemble P3 on
averaged member lines over common cells (v4 item 6, tools/ensemble/p3build.py arithmetic), shared-relation P3s (tools/v24_S
tests, D4 exclude rule registered by v3 item 6 / item 6), and the A6 frame dedup gate (one plot frame row per location).
Inputs and outputs are group-level lines and counts only (no plot rows). Usage from rc_adapt.attach_p3(bundle).
"""
import glob, json, math, os, shutil, subprocess, sys
from pathlib import Path
import numpy as np
import pandas as pd

H = Path(os.path.expanduser("~")); CODE = H / "perseus_me/code"; RUNS = H / "perseus_me/ncasi_rc/runs"; MS = H / "perseus_me/me_state"
sys.path.insert(0, str(CODE / "tools/ensemble")); sys.path.insert(0, str(CODE / "tools/v24_S"))
os.environ.setdefault("NCASI_P3S_D4", "exclude")
import p3build as PB      # noqa: E402
import s_score as SS      # noqa: E402

# per engine: (arm label used in the work dir, Bakuzis trajectory, P3s line); {ST} filled; Maine uses the v2.4 work files
TRAJ = {
    "fvs": ("default", ["{runs}/fvs/{ST}/**/traj__fvs__default.csv"], ["{runs}/fvs/{ST}/**/traj__fvs__default.csv"]),
    "osm": ("continuous", ["{runs}/osm/{ST}/**/traj__osm__continuous.csv"], ["{runs}/osm/{ST}/**/traj__osm__continuous.csv"]),
    "acd141": ("acd_v141_csi14", ["{runs}/acd141/{ST}/**/traj__acd141__acd_v141_csi14.csv"], ["{runs}/acd141/{ST}/**/traj__acd141__acd_v141_csi14.csv"]),
    "cem": ("no_harvest_baseline", ["{runs}/cem/{ST}/restricted/values_bk/traj__cem__no_harvest_baseline.csv"],
            ["{runs}/cem/{ST}/restricted/values_bk/traj__cem__no_harvest_baseline_bio.csv", "{runs}/cem/{ST}/restricted/values_bk/traj__cem__no_harvest_baseline.csv"]),
    "landis": ("plotframe_no_harvest", [], ["{runs}/p3s/{ST}/landis_bk/traj__landis__plotframe_no_harvest.csv"]),
    "cbm_gated": ("plot_lines", [], ["{runs}/p3s/{ST}/traj__cbm_gated__plot_lines.csv"]),
    "cbm_af": ("plot_lines", [], ["{runs}/p3s/{ST}/traj__cbm_af__plot_lines.csv"]),
    "fia_curves": ("plot_lines", [], ["{runs}/p3s/{ST}/traj__fia_curves__plot_lines.csv"]),
}
BK_ENGINES = ("fvs", "osm", "acd141", "cem")   # tree-list engines with A6 Bakuzis lines (ncasi_states adds its FVS-ACD rows)
ME_BK = MS / "bakuzis/work"; ME_LINES = MS / "v24_S/work/lines"
ME_P3S = {"fvs": "traj__fvs__default.csv", "osm": "traj__osm__continuous.csv", "acd141": "traj__acd141__acd_v141_csi14.csv",
          "cem": "traj__cem__no_harvest_baseline_bio.csv", "landis": "traj__landis__plotframe_no_harvest.csv",
          "cbm_gated": "traj__cbm_gated__plot_lines.csv", "cbm_af": "traj__cbm_af__plot_lines.csv", "fia_curves": "traj__fia_curves__plot_lines.csv"}


def _find(pats, st):
    for p in pats:
        g = sorted(glob.glob(p.format(runs=RUNS, ST=st), recursive=True))
        if g:
            return Path(g[0])
    return None


def a6_dedup(st):
    """A6 gate: the P3s/Bakuzis plot frame holds one row per plot location (latest measurement before the run start)"""
    f = _find(["{runs}/p3s/{ST}/plot_groups.parquet"], st)
    if f is None:
        return dict(ok=None, note="no plot frame found")
    g = pd.read_parquet(f)
    sys.path.insert(0, str(Path(__file__).parent)); import rc_design as D
    loc = g.PLT_CN.astype(str).map(D.loc_of(st))
    col = "sampled" if "sampled" in g.columns else None
    sub = g[g[col]] if col else g
    lsub = loc[sub.index]
    dup = int(lsub.duplicated().sum())
    return dict(ok=dup == 0, frame_rows=int(len(g)), sampled_rows=int(len(sub)), duplicate_locations_in_sample=dup,
                note="A6 satisfied" if dup == 0 else "A6 not applied: duplicate locations in the sampled frame; P3 and P3s not scored")


def fia_percentiles(st):
    """BAPH and SDI percentiles by forest type over the frame plots (bk_common.fia_percentiles arithmetic on the DataMart TREE)"""
    sys.path.insert(0, str(Path(__file__).parent)); import rc_design as D
    g = pd.read_parquet(RUNS / "p3s" / st / "plot_groups.parquet")
    tr = D._rd(st, "TREE", ("PLT_CN", "STATUSCD", "DIA", "TPA_UNADJ"))
    t = tr[(tr.STATUSCD == 1) & (tr.DIA >= 1) & tr.TPA_UNADJ.notna() & tr.PLT_CN.astype(str).isin(set(g.PLT_CN.astype(str)))].copy()
    HA, KBA = 2.4710538, math.pi / 40000.0
    t["u"] = t.TPA_UNADJ * HA; t["ba"] = t.u * KBA * (t.DIA * 2.54) ** 2
    X = t.groupby(t.PLT_CN.astype(str)).agg(tph=("u", "sum"), ba=("ba", "sum")).join(g.assign(PLT_CN=g.PLT_CN.astype(str)).set_index("PLT_CN")[["type", "forest_prop"]], how="inner")
    X["BAPH"] = X.ba / X.forest_prop; X["TPH"] = X.tph / X.forest_prop; X["QMD"] = np.sqrt(X.BAPH / X.TPH / KBA); X["SDI"] = X.TPH * (X.QMD / 25.4) ** 1.605
    out = [dict(type=k, n_plots=len(v), BAPH_p99=float(np.quantile(v.BAPH, 0.99)), BAPH_p999=float(np.quantile(v.BAPH, 0.999)),
                SDI_p99=float(np.quantile(v.SDI, 0.99)), SDI_p999=float(np.quantile(v.SDI, 0.999)), BAPH_median=float(v.BAPH.median()),
                SDI_median=float(v.SDI.median())) for k, v in X.groupby("type")]
    return pd.DataFrame(out)


def realism_from_scores(f):
    """noR5 realism and coverage per engine x arm (tools/mpi2/mpi2.bakuzis arithmetic)"""
    b = pd.read_csv(f); b["total"] = b.groups_evaluated + b.groups_not_evaluated; b["passed"] = b.groups_evaluated - b.groups_flagged
    out = {}
    for (e, a), g in b.groupby(["engine", "arm"]):
        r5c = g.relation.str.startswith("R5c"); r5 = g.relation.str.startswith("R5") & ~r5c; h = g[~r5 & ~r5c]; ev = h.groups_evaluated.sum()
        out[(e, a)] = dict(realism=h.passed.sum() / ev if ev else np.nan, coverage=ev / h.total.sum() if h.total.sum() else 0.0, evaluated=int(ev))
    return out


def ensemble_lines(trajs, name, how, dst):
    """averaged member lines over common cells (p3_22 arithmetic): 5 yr grid, a cell missing in any member stays missing"""
    G = {}
    for k, f in trajs.items():
        x = pd.read_csv(f); x["site"] = np.round(x.site).astype(int)
        vg = "VOL" in x and bool(np.isfinite(x.VOL).all()); g = PB.grid(x, vg); g["VOL"] = g.VOLg if vg else g.VOLd
        G[k] = g.set_index(PB.KEY + ["age"]).sort_index()
    if not G:
        return None
    ks = list(G); idx = G[ks[0]].index
    for k in ks[1:]:
        idx = idx.intersection(G[k].index)
    idx = idx.sort_values(); d = pd.DataFrame(index=idx)
    for v in ("HT", "TPH", "BAPH", "VOL"):
        X = np.stack([G[k].loc[idx, v].values.astype(float) for k in ks], 0)
        d[v] = X.mean(0) if how == "mean" else np.median(X, 0)
    d["n_plots"] = np.min(np.stack([G[k].loc[idx, "n_plots"].values for k in ks], 0), 0)
    d = d.reset_index(); PB.write(d, "ensemble", name, Path(dst) / f"traj__ensemble__{name}.csv", True)
    return dict(ensemble=name, members=",".join(ks), cells=int(len(idx)))


def run_bakuzis(work, out):
    r = subprocess.run(["Rscript", str(CODE / "tools/bakuzis/score_v2.R"), str(work), str(out)], capture_output=True, text=True)
    if r.returncode != 0 or not (Path(out) / "bakuzis_scores.csv").exists():
        raise RuntimeError("score_v2.R failed: " + r.stderr[-600:])
    return realism_from_scores(Path(out) / "bakuzis_scores.csv")


def p3s_score(lines, horizon_of=None):
    """shared relations R1 to R4 per row (s_score.tests), declared by-construction relations removed (DECL), relations no row
    flags dropped (non-discriminating, over this row set), D4 exclude rule. lines: dict row -> (engine, path). Returns dict row ->
    dict(realism, coverage, evaluable)"""
    T, meta = [], {}
    for row, (e, f) in lines.items():
        h = (horizon_of or {}).get(row, SS.HORIZON.get(e, 100))
        t = SS.tests(pd.read_csv(f), h); t.insert(0, "row", row); T.append(t); meta[row] = e
    if not T:
        return {}, []
    T = pd.concat(T, ignore_index=True); C = []
    for row, e in meta.items():
        for k in SS.RELS:
            t = T[(T.row == row) & (T.rel == k)]; decl = k in SS.DECL.get(e, set())
            ev = int(t.st.isin(["PASS", "FLAG"]).sum()); fl = int((t.st == "FLAG").sum())
            C.append(dict(row=row, rel=k, evaluated=0 if decl else ev, flagged=0 if decl else fl, ne_data=0 if decl else SS.UNIV[k] - ev, ne_declared=SS.UNIV[k] if decl else 0))
    C = pd.DataFrame(C); evr = C[C.evaluated > 0]
    drop = [k for k in SS.RELS if len(evr[evr.rel == k]) and (evr[evr.rel == k].flagged == 0).all()]
    out = {}
    for row in meta:
        c = C[(C.row == row) & ~C.rel.isin(drop)]; ev, fl, ned, nde = int(c.evaluated.sum()), int(c.flagged.sum()), int(c.ne_data.sum()), int(c.ne_declared.sum())
        cov = ev / (ev + ned) if (ev + ned) else np.nan; cov2 = ev / (ev + ned + nde) if (ev + ned + nde) else np.nan
        ok = ev > 0 and (cov if SS.D4 == "exclude" else cov2) >= 0.5
        out[row] = dict(realism=(ev - fl) / ev if ev else np.nan, coverage=cov if ok else 0.0, coverage_declared_counted=cov2, evaluable=ok)
    return out, drop


def state_p3(st, engines, ensembles, work_root):
    """engines: dict row name -> engine; ensembles: dict name -> (member row names, how). Returns (P3 dict, P3s dict, ens P3 dict,
    report). Maine uses the v2.4 work files (bakuzis/work, v24_S/work/lines)."""
    rep = dict(state=st)
    P3, P3S, ENS = {}, {}, {}
    A6DIR = None
    if st == "ME":
        rep["a6"] = dict(ok=True, note="Maine v2.4 frame (registered)")
    else:
        import rc_a6
        A6DIR = Path(work_root) / "a6_lines"; rep["a6"] = rc_a6.build(st, A6DIR); rep["a6"]["ok"] = True
        rep["a6"]["note"] = "A6 applied as a filter on existing run outputs (latest measurement per location; latest run row if the chosen row was not run)"
    W = Path(work_root) / "bk_work"; O = Path(work_root) / "bk_out"
    shutil.rmtree(W, ignore_errors=True); W.mkdir(parents=True, exist_ok=True)
    traj, lines = {}, {}
    for row, e in engines.items():
        if e not in TRAJ:
            continue
        arm = TRAJ[e][0]
        f = (ME_BK / f"traj__{e}__{arm}.csv") if st == "ME" else (A6DIR / "bk" / f"traj__{e}__{arm}.csv" if e in BK_ENGINES else None)
        if f is not None and f.exists():
            shutil.copy(f, W / f"traj__{e}__{arm}.csv"); traj[row] = W / f"traj__{e}__{arm}.csv"
        fl = (ME_LINES / ME_P3S[e]) if st == "ME" else A6DIR / "s" / ME_P3S[e]
        if fl is not None and Path(fl).exists():
            lines[row] = (e, Path(fl))
    pct = (ME_BK / "fia_percentiles.csv") if st == "ME" else None
    if pct is not None:
        shutil.copy(pct, W / "fia_percentiles.csv")
    else:
        try:
            fia_percentiles(st).to_csv(W / "fia_percentiles.csv", index=False)
        except Exception as ex:  # noqa: BLE001
            rep["fia_percentiles"] = repr(ex)
    rep["p3_members_with_traj"] = sorted(traj); rep["p3s_rows_with_lines"] = sorted(lines)
    ens_rep = []
    for name, (mem, how) in ensembles.items():
        mt = {m: traj[m] for m in mem if m in traj}
        if mt:
            r = ensemble_lines(mt, name, how, W)
            if r:
                ens_rep.append(r)
    rep["ensemble_lines"] = ens_rep
    if traj and (W / "fia_percentiles.csv").exists():
        R = run_bakuzis(W, O)
        for row, e in engines.items():
            v = R.get((e, TRAJ.get(e, ("",))[0]))
            if v is not None and row in traj:
                P3[row] = v
        for r in ens_rep:
            v = R.get(("ensemble", r["ensemble"]))
            if v is not None:
                ENS[(r["ensemble"], "P3")] = v
    if lines:
        hz = {row: SS.HORIZON.get(e, 100) for row, (e, _) in lines.items()}
        ens_lines = {}
        for r in ens_rep:
            ens_lines[r["ensemble"]] = ("ensemble", W / f"traj__ensemble__{r['ensemble']}.csv")
            hz[r["ensemble"]] = 50 if any(engines.get(m) == "cem" for m in r["members"].split(",")) else 100
        S, drop = p3s_score({**lines, **ens_lines}, hz)
        rep["p3s_dropped_non_discriminating"] = drop
        for row, v in S.items():
            if row in lines:
                P3S[row] = v
            else:
                ENS[(row, "P3s")] = v
    return P3, P3S, ENS, rep
