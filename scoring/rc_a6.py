"""Amendment A6 applied as a filter on existing run outputs: the Bakuzis/P3s plot frame is reduced to one row per plot location,
the latest measurement (INVYR, then CN) among the location's sampled frame rows. Per engine, a location whose chosen row was not
run takes that engine's latest run row for the location (counted). Group lines are then rebuilt with the registered arithmetic:
Bakuzis lines (P3) by tools/bakuzis/bk_common.group_means (tree lists) and the cem_bk_state.R rule (CEM); P3s lines by FIA
productivity class (SITECLCD, tools/v24_S/build_lines.site_class rule) with group_means (tree lists), the CEM VOL and BIO arms,
and vol_groups (LANDIS-II, curves). A no-dedup rebuild must reproduce the existing line files (rebuild check).
Group-level lines only leave this module; plot rows stay in memory."""
import json, os, sys
from pathlib import Path
import numpy as np
import pandas as pd

H = Path(os.path.expanduser("~")); CODE = H / "perseus_me/code"; RUNS = H / "perseus_me/ncasi_rc/runs"
sys.path.insert(0, str(CODE / "tools/bakuzis")); sys.path.insert(0, str(Path(__file__).parent))
import bk_common as B   # noqa: E402
import rc_design as D   # noqa: E402
GCOL = ["type", "site", "origin", "offset"]

PV = {  # engine -> (plot values file, arm, kind)
    "fvs": ("fvs/{ST}/restricted/bakuzis/plot_values__fvs__default.parquet", "default", "treelist"),
    "cem": ("cem/{ST}/restricted/values_bk/cem_plot_values.csv", "no_harvest_baseline", "cem"),
    "landis": ("p3s/{ST}/landis_bk/plot_values__landis__plotframe_no_harvest.parquet", "plotframe_no_harvest", "vol"),
    "cbm_gated": ("p3s/{ST}/plot_values__cbm_gated__plot_lines.parquet", "plot_lines", "vol"),
    "cbm_af": ("p3s/{ST}/plot_values__cbm_af__plot_lines.parquet", "plot_lines", "vol"),
    "fia_curves": ("p3s/{ST}/plot_values__fia_curves__plot_lines.parquet", "plot_lines", "vol"),
}
EXISTING = {  # existing line files used by the rebuild check (bk = Bakuzis line, s = P3s line)
    ("fvs", "bk"): "fvs/{ST}/bakuzis/traj__fvs__default.csv", ("cem", "bk"): "cem/{ST}/restricted/values_bk/traj__cem__no_harvest_baseline.csv",
    ("landis", "s"): "p3s/{ST}/landis_bk/traj__landis__plotframe_no_harvest.csv", ("cbm_gated", "s"): "p3s/{ST}/traj__cbm_gated__plot_lines.csv",
    ("cbm_af", "s"): "p3s/{ST}/traj__cbm_af__plot_lines.csv", ("fia_curves", "s"): "p3s/{ST}/traj__fia_curves__plot_lines.csv",
}


def site_class(st, plts):
    c = D._rd(st, "COND", ("PLT_CN", "CONDID", "COND_STATUS_CD", "SITECLCD", "CONDPROP_UNADJ")).copy(); c["PLT_CN"] = c.PLT_CN.astype(str).str.strip()
    c = c[c.PLT_CN.isin(set(plts)) & (c.COND_STATUS_CD == 1) & c.SITECLCD.between(1, 7)]
    a = c.groupby(["PLT_CN", "SITECLCD"]).agg(area=("CONDPROP_UNADJ", "sum"), cmin=("CONDID", "min")).reset_index()
    a = a.sort_values(["PLT_CN", "area", "cmin"], ascending=[True, False, True]).drop_duplicates("PLT_CN")
    return dict(zip(a.PLT_CN, np.where(a.SITECLCD <= 3, 3, np.where(a.SITECLCD <= 5, 2, 1))))


def chosen_rows(st):
    G = pd.read_parquet(RUNS / "p3s" / st / "plot_groups.parquet"); G["PLT_CN"] = G.PLT_CN.astype(str).str.strip()
    S = G[G.sampled].copy() if "sampled" in G.columns else G.copy()
    pl = D._rd(st, "PLOT", ("CN", "INVYR")); iv = dict(zip(pl.CN.astype(str).str.strip(), pl.INVYR))
    loc = D.loc_of(st); S["loc"] = S.PLT_CN.map(loc); S["INVYR"] = S.PLT_CN.map(iv)
    C = S.sort_values(["loc", "INVYR", "PLT_CN"]).groupby("loc").tail(1)
    return G, S, C, loc, iv


def _read(f):
    x = pd.read_parquet(f) if str(f).endswith(".parquet") else pd.read_csv(f, dtype={"PLT_CN": str})
    x["PLT_CN"] = x.PLT_CN.astype(str).str.strip(); return x


def dedup(V, chosen, loc, iv):
    """keep per location the chosen row if run, else the engine's latest run row; returns filtered V and counts"""
    runp = pd.DataFrame(dict(PLT_CN=V.PLT_CN.unique())); runp["loc"] = runp.PLT_CN.map(loc); runp["INVYR"] = runp.PLT_CN.map(iv)
    ch = set(chosen.PLT_CN); keep = []; sub = 0
    for l, g in runp.groupby("loc"):
        hit = g[g.PLT_CN.isin(ch)]
        if len(hit):
            keep.append(hit.PLT_CN.iloc[0])
        else:
            keep.append(g.sort_values(["INVYR", "PLT_CN"]).PLT_CN.iloc[-1]); sub += 1
    n_loc_chosen = int(chosen["loc"].nunique()); n_loc_run = int(runp["loc"].nunique())
    return V[V.PLT_CN.isin(set(keep))], dict(run_rows=int(len(runp)), run_locations=n_loc_run, kept_rows=len(keep),
                                             dropped_duplicate_rows=int(len(runp) - len(keep)), substituted_latest_run_row=sub,
                                             chosen_locations_not_run=int(n_loc_chosen - len(set(chosen["loc"]) & set(runp["loc"]))))


def bk_lines(V, e, arm, kind):
    if kind == "treelist":
        return B.group_means(V, e, arm)
    if kind == "cem":
        P = V.assign(offset=5 * V.cycle)
        g = P.assign(age=P.STDAGE + P.offset).groupby(GCOL).agg(age=("age", "mean"), TPH=("TPH", "mean"), BAPH=("BAPH", "mean"), VOL=("VOL", "mean"),
                                                             n_plots=("PLT_CN", "nunique")).reset_index()
        g["HT"] = np.nan; g.insert(0, "arm", arm); g.insert(0, "engine", e); return g
    return None


def s_lines(V, e, arm, kind, cls):
    V = V.copy(); V["site"] = V.PLT_CN.map(cls); V = V[V.site.notna()]
    if kind == "treelist":
        return {arm: B.group_means(V, e, arm).drop(columns=["n_ht", "QMD"])}
    if kind == "cem":
        P = V.assign(offset=5 * V.cycle); out = {}
        for a, col in (("no_harvest_baseline", "VOL"), ("no_harvest_baseline_bio", "BIO")):
            Q = P[P[col].notna()]
            g = Q.assign(age=Q.STDAGE + Q.offset).groupby(GCOL).agg(age=("age", "mean"), TPH=("TPH", "mean"), BAPH=("BAPH", "mean"), VOL=(col, "mean"),
                                                                 n_plots=("PLT_CN", "nunique")).reset_index()
            g["HT"] = np.nan; g.insert(0, "arm", a); g.insert(0, "engine", e); out[a] = g
        return out
    g = V.assign(age=V.STDAGE + V.offset).groupby(GCOL).agg(age=("age", "mean"), VOL=("VOL", "mean"), n_plots=("PLT_CN", "nunique")).reset_index()
    g["HT"] = np.nan; g["TPH"] = np.nan; g["BAPH"] = np.nan; g.insert(0, "arm", arm); g.insert(0, "engine", e)
    return {arm: g}


def _maxdiff(a, b):
    k = [c for c in ["type", "site", "origin", "offset"] if c in a and c in b]
    a = a.copy(); b = b.copy()
    for d in (a, b):
        d["site"] = np.round(d.site).astype(int)
    m = a.merge(b, on=k, how="outer", suffixes=("_a", "_b"), indicator=True)
    if (m._merge != "both").any():
        return float("inf")
    cols = [c for c in ("age", "HT", "TPH", "BAPH", "VOL", "n_plots") if c + "_a" in m and c + "_b" in m]
    return float(max(np.nanmax(np.abs(m[c + "_a"].astype(float) - m[c + "_b"].astype(float))) if m[c + "_a"].notna().any() else 0.0 for c in cols))


def build(st, out):
    """writes out/bk/traj__*.csv (P3) and out/s/traj__*.csv (P3s); returns report (counts, rebuild checks)"""
    out = Path(out); (out / "bk").mkdir(parents=True, exist_ok=True); (out / "s").mkdir(parents=True, exist_ok=True)
    G, S, C, loc, iv = chosen_rows(st)
    rep = dict(state=st, frame_rows=int(len(G)), sampled_rows=int(len(S)), sampled_locations=int(S["loc"].nunique()),
               duplicate_rows_removed_from_frame=int(len(S) - len(C)), engines={}, rebuild_check={})
    cls = site_class(st, set(S.PLT_CN))
    rep["site_class_counts"] = {int(k): int(v) for k, v in pd.Series(list(cls.values())).value_counts().sort_index().items()}
    for e, (pat, arm, kind) in PV.items():
        f = RUNS / pat.format(ST=st)
        if not f.exists():
            continue
        V = _read(f)
        if "sampled" in V.columns:
            V = V[V.sampled.astype(bool)]
        # rebuild check without dedup (proves the line arithmetic before filtering)
        for which, ex in (("bk", EXISTING.get((e, "bk"))), ("s", EXISTING.get((e, "s")))):
            if ex and (RUNS / ex.format(ST=st)).exists():
                ref = pd.read_csv(RUNS / ex.format(ST=st))
                mine = bk_lines(V, e, arm, kind) if which == "bk" else list(s_lines(V, e, arm, kind, cls).values())[0]
                rep["rebuild_check"][f"{e}.{which}"] = _maxdiff(mine, ref) if mine is not None else None
        Vd, cnt = dedup(V, C, loc, iv); rep["engines"][e] = cnt
        b = bk_lines(Vd, e, arm, kind)
        if b is not None:
            b.to_csv(out / "bk" / f"traj__{e}__{arm}.csv", index=False)
        for a, g in s_lines(Vd, e, arm, kind, cls).items():
            g.to_csv(out / "s" / f"traj__{e}__{a}.csv", index=False)
    (out / "a6_report.json").write_text(json.dumps(rep, indent=1, default=float))
    return rep
