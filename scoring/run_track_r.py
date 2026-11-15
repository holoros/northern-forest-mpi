"""Track R driver: score one state (or all ready states) under the frozen pre-registration (v5) and amendment A.
Writes aggregates only: ~/perseus_me/ncasi_rc/scores/<ST>/{report_card.csv, components.csv, membership.csv, scalars.csv,
stress.csv, five_seed.csv, missing.json, acceptance.json} and six-state tables in ~/perseus_me/ncasi_rc/scores/.
usage: run_track_r.py [--states ME,NH,...] [--draws 1000] [--force] [--check-only]
"""
import argparse, hashlib, json, os, re, subprocess, sys, time
from pathlib import Path
import numpy as np
import pandas as pd
sys.path.insert(0, str(Path(__file__).parent))
import rc_core as C
import rc_adapt as A
import rc_run as RR

H = Path(os.path.expanduser("~")); RC = H / "perseus_me/ncasi_rc"; SC = Path(os.environ.get("RC_SCORES_DIR", str(RC / "scores")))
PREREG_MD5 = "8d4964a38dc01c5cf73ce16d5ddb747b"; AMEND_MD5 = "8ed0bf0710d7b65eeb3f9f68aaf78b1e"
STATES = ["ME", "NH", "VT", "MI", "WI", "MN"]


def md5(p):
    return hashlib.md5(Path(p).read_bytes()).hexdigest()


def ready(st):
    """state is ready when no run process for it is alive and no matrix row of the run reports for this state is running or
    queued. Only table rows are read (the static head lines of STATEWIDE_REPORT.md are ignored); the CEM TreeMap cell (S8 only)
    does not block."""
    ps = subprocess.run(["ps", "-eo", "args"], capture_output=True, text=True).stdout.splitlines()
    alive = [l for l in ps if re.search(rf"(\b|_|/){st}(\b|_|/)", l) and re.search(r"run_|landis|cem|fvs|osm|treelist|curves|sw_", l)
             and "run_track_r" not in l and not re.search(r"values_tm|sw_lists_tm|\btm\b", l)]
    pend = []
    for f in (RC / "runs/RUNS_REPORT.md", RC / "runs/STATEWIDE_REPORT.md"):
        if not f.exists():
            continue
        hdr = None
        for l in f.read_text().splitlines():
            if not l.startswith("|"):
                hdr = None if not l.strip() else hdr; continue
            cells = [c.strip() for c in l.strip().strip("|").split("|")]
            low = [c.lower() for c in cells]
            if "state" in low and ("engine" in low or "component" in low):
                hdr = low; continue
            if hdr is None or set("".join(cells)) <= set("-: "):
                continue
            row = dict(zip(hdr, cells))
            if row.get("state", "").upper() != st:
                continue
            for k, v in row.items():
                if k in ("engine", "state") or ("treemap" in k and "cem" in row.get("engine", "").lower()):
                    continue
                if re.search(r"\b(running|queued|pending|waiting)\b", v, re.I):
                    pend.append(f"{f.name}: {row.get('engine', '')} {k}: {v[:60]}")
    return (not alive and not pend), alive[:5], pend[:8]


def shuffle_bundle(b, seed):
    """row-shuffle test (F14): permute every unit-level array consistently"""
    import copy
    c = copy.deepcopy(b); r = np.random.default_rng(seed)
    for key in ("pairs", "treated"):
        U = c[key]; p = r.permutation(len(U["ids"]))
        for k in ("ids", "labels", "fold", "locs", "subset_mask"):
            if k in U and U[k] is not None:
                U[k] = np.asarray(U[k])[p]
        U["M"] = U["M"][:, p]; U["obs"] = {q: np.asarray(v)[p] for q, v in U["obs"].items()}
        for E in c["engines"].values():
            comp = "C2" if key == "pairs" else "C5"
            v = E.get(comp)
            if isinstance(v, dict) and v and all(isinstance(x, np.ndarray) or x is None for x in v.values()) and set(v) <= set(C.QS):
                E[comp] = {q: np.asarray(x)[p] for q, x in v.items()}
            elif isinstance(v, dict) and ("level" in v or "ageq" in v):
                E[comp] = {rk: ({q: np.asarray(x)[p] for q, x in rv.items()} if rv else rv) for rk, rv in v.items()}
    return c


def rows_df(rows):
    return pd.DataFrame([{k: v for k, v in r.items() if k != "_bt"} for r in rows])


def dirichlet_medians(T, scheme, n=10000, seed=20261105):
    rng = np.random.default_rng(seed); w0 = np.array([C.W[p] for p in C.PILLARS]) / 100.0
    Wd = rng.dirichlet(np.ones(5) if scheme == "flat" else 20 * w0, n)
    out = {}
    for r in T.itertuples():
        p = np.array([getattr(r, q) for q in C.PILLARS], float); ok = np.isfinite(p)
        if not ok.any():
            out[r.row] = np.nan; continue
        m = (Wd[:, ok] @ p[ok]) / Wd[:, ok].sum(1); out[r.row] = float(np.median(m))
    return out


def refold(b, seed):
    """five-seed disclosure (v2 item 3): reassign location folds by the registered rule (stratified by domain of the location's
    first unit, keys sorted then permuted) with another seed"""
    import copy
    c = copy.deepcopy(b)
    L = pd.DataFrame(dict(loc=np.concatenate([c["pairs"]["locs"], c["treated"]["locs"]]),
                          dom=np.concatenate([c["pairs"]["labels"], c["treated"]["labels"]]))).drop_duplicates("loc")
    f = {}
    for d, g in L.groupby("dom"):
        ks = sorted(g["loc"]); p = np.random.default_rng(seed).permutation(len(ks))
        for i, j in enumerate(p):
            f[ks[j]] = 1 if i % 2 == 0 else 2
    for key in ("pairs", "treated"):
        c[key]["fold"] = np.array([f[l] for l in c[key]["locs"]])
    return c


def score_one(st, B, log):
    t0 = time.time(); out = SC / st; out.mkdir(parents=True, exist_ok=True)
    landis_extract(st, log)
    b = A.load_state(st, B=B); A.fill_c3(b)
    A.attach_p3(b, out / "p3_work", log)
    (out / "p3_report.json").write_text(json.dumps(b.get("p3_report", {}), indent=1, default=str))
    import rc_s8
    s8 = rc_s8.contrast(st, float(b["c1_target"]["stock"][0]) * A.LB2TG)
    s8.to_csv(out / "s8_treemap.csv", index=False, float_format="%.6g")
    rows, comp, memb, scal = RR.score_state(b)
    T = rows_df(rows); T.to_csv(out / "report_card.csv", index=False, float_format="%.6g")
    comp.to_csv(out / "components.csv", index=False, float_format="%.6g"); memb.to_csv(out / "membership.csv", index=False)
    hx = comp[comp.component == "C1"].copy() if len(comp) else comp
    if len(hx):
        rc = [c for c in hx.columns if c.startswith("ratio_") and c.endswith("_statewide")]
        hx = hx[["row"] + rc + (["hard_exit"] if "hard_exit" in hx.columns else [])]
        hx["hard_exit_bound"] = "0.2 to 5 x EVALIDator (item 10, v4 item 4)"
        hx["all_quantities_within_bound"] = hx.get("hard_exit", pd.Series([np.nan] * len(hx))).isna().values
        hx.to_csv(out / "hard_exit.csv", index=False, float_format="%.6g")
    scal.to_csv(out / "scalars.csv", index=False, float_format="%.6g")
    # ---- stress alternatives (item 11, v4 item 5, v3 item 6)
    alt = {"S2": dict(c5_quant="3q"), "S3": dict(treated_subset="registered"), "S4": dict(curve_rule="ageq"), "S5": dict(p3="P3s"),
           "S9": dict(min_chains=30), "S10": dict(view="by_domain"), "S14": dict(ens_c5_curves="level"), "S15": dict(ens_c5_curves="ageq"),
           "S18a": dict(zero=0.25), "S18b": dict(zero=1.0), "A1s": dict(p5_target="raw")}
    if st == "ME":
        alt["S13"] = dict(exclude_acd_c2p2=True)
    if st == "MN":
        alt["S16"] = dict(drop_members=("landis",))
    AV = {}
    for k, o in alt.items():
        try:
            AV[k] = rows_df(RR.score_state(b, o)[0]).set_index("row").MPI
        except Exception as ex:  # noqa: BLE001
            log(f"{st} stress {k}: {ex!r}"); AV[k] = None
    for cut, lab in ((0.25, "S19a"), (0.75, "S19b")):
        import copy
        bb = copy.deepcopy(b)
        if "cem" in bb["engines"] and np.isfinite(bb["engines"]["cem"].get("a8_share", np.nan)):
            sh = bb["engines"]["cem"]["a8_share"]; bb["engines"]["cem"]["C5_flag"] = "A8" if sh < cut else ""
            AV[lab] = rows_df(RR.score_state(bb)[0]).set_index("row").MPI
    for sch, lab in (("flat", "S11_flat"), ("centered", "S11_centered")):
        AV[lab] = pd.Series(dirichlet_medians(T, sch))
    Ti = T.set_index("row"); srows = []
    applies = lambda k, r: (k not in ("S6", "S7", "S14", "S15") or r.startswith("E_")) and (k not in ("S13",) or r.startswith("acd141") or r.startswith("E_"))  # noqa: E731
    leader = Ti[Ti.P_rank1.notna()].MPI.idxmax() if Ti.P_rank1.notna().any() else None
    lead_by = {}
    for k, v in AV.items():
        if v is None or k == "A1s":
            continue
        cand = v.reindex(Ti[Ti.P_rank1.notna()].index).dropna()
        lead_by[k] = cand.idxmax() if len(cand) else None
    for r in Ti.index:
        alts = {k: (float(v[r]) if v is not None and r in v.index and np.isfinite(v[r]) else None) for k, v in AV.items() if k != "A1s" and applies(k, r)}
        if r == "E_state":
            for e2, lab in (("E_plot", "S6"), ("E_state_median", "S7")):
                alts[lab] = float(Ti.loc[e2, "MPI"]) if e2 in Ti.index and np.isfinite(Ti.loc[e2, "MPI"]) else None
        sr = C.stress_range(float(Ti.loc[r, "MPI"]), alts) if np.isfinite(Ti.loc[r, "MPI"]) else dict(lo=np.nan, hi=np.nan, n_alternatives=0, at_favourable_extreme=False, applied=[])
        srows.append(dict(row=r, MPI=Ti.loc[r, "MPI"], range_lo=sr["lo"], range_hi=sr["hi"], n_alternatives=sr["n_alternatives"],
                          at_favourable_extreme=sr["at_favourable_extreme"], alternatives=";".join(sr["applied"]),
                          A1s_MPI=float(AV["A1s"][r]) if AV.get("A1s") is not None and r in AV["A1s"].index else np.nan,
                          **{f"MPI_{k}": (alts.get(k)) for k in sorted(set(AV) | {"S6", "S7"}) if k != "A1s"}))
    S = pd.DataFrame(srows); S["leader_registered"] = leader; S["leader_retention_share"] = C.leader_share(lead_by, leader)
    S.to_csv(out / "stress.csv", index=False, float_format="%.6g")
    # ---- five fold seeds for recalibrated rows (disclosure)
    fs = []
    for sd in (20261105, 20261106, 20261107, 20261108, 20261109):
        bb = b if sd == 20261105 else refold(b, sd)
        tt = rows_df(RR.score_state(bb)[0])
        for r in tt[tt.role.isin(["recal"]) | (tt.row == "E_recal")].itertuples():
            fs.append(dict(fold_seed=sd, row=r.row, MPI=r.MPI, registered=(sd == 20261105)))
    pd.DataFrame(fs).to_csv(out / "five_seed.csv", index=False, float_format="%.6g")
    # ---- acceptance: row shuffle, second path, missing
    T2 = rows_df(RR.score_state(shuffle_bundle(b, 11))[0]).set_index("row")
    num = [c for c in Ti.columns if Ti[c].dtype.kind == "f"]
    sh = float(np.nanmax(np.abs(Ti[num].values.astype(float) - T2.loc[Ti.index, num].values.astype(float)))) if len(Ti) else 0.0
    sp = subprocess.run(["Rscript", str(Path(__file__).parent / "second_path_mpi.R"), str(out / "report_card.csv"), str(out / "second_path.json")],
                        capture_output=True, text=True)
    acc = dict(state=st, prereg_md5=PREREG_MD5, amendment_md5=AMEND_MD5, draws=B, seed=20261105, row_shuffle_max_abs=sh,
               second_path=json.loads((out / "second_path.json").read_text()) if (out / "second_path.json").exists() else sp.stderr[-400:],
               insample_register=b.get("insample_note"), units=b["units"], n_outside_universe=b["n_outside_universe"],
               members=b["members"], wall_s=round(time.time() - t0, 1))
    (out / "missing.json").write_text(json.dumps(b["missing"], indent=1))
    (out / "acceptance.json").write_text(json.dumps(acc, indent=1, default=float))
    log(json.dumps(acc, default=float))
    return T


def landis_extract(st, log):
    """aggregate LANDIS-II statewide run outputs (tools/ncasi_rc/eng/landis_sw_extract.py, venv_fig) when the run is complete
    and no extract exists; aggregates only"""
    base = RC / "runs/landis" / st / "restricted"
    for arm in ("fia", "tm"):
        run = base / f"statewide_{arm}"; out = base / f"statewide_{arm}_extract.csv"
        tifs = list((run / "output/biomass").glob("biomass-TotalBiomass-*.tif")) if (run / "output/biomass").exists() else []
        if out.exists() or not tifs:
            continue
        r = subprocess.run([str(H / "perseus_me/venv_fig/bin/python"), str(H / "perseus_me/code/tools/ncasi_rc/eng/landis_sw_extract.py"),
                            st, arm, str(run), str(out)], capture_output=True, text=True)
        log(f"{st} landis extract {arm} rc={r.returncode} {r.stderr[-200:]}")


def leak_scan(root):
    bad = []
    for f in Path(root).rglob("*.csv"):
        txt = f.read_text(errors="ignore")
        if re.search(r"\b(PLT_CN|LAT|LON|pair_id|tree_cn|chain_id)\b", txt.splitlines()[0] if txt else "") or re.search(r"(?<![\d.])\d{13,}(?![\d.])", txt):
            bad.append(str(f))
    return bad


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--states", default=",".join(STATES)); ap.add_argument("--draws", type=int, default=1000)
    ap.add_argument("--force", action="store_true"); ap.add_argument("--check-only", action="store_true")
    a = ap.parse_args(); SC.mkdir(parents=True, exist_ok=True); logf = open(SC / "run_track_r.log", "a")
    log = lambda m: (print(m, flush=True), logf.write(f"[{time.strftime('%F %T')}] {m}\n"), logf.flush())  # noqa: E731
    assert md5(RC / "docs_PREREG_R1.md") == PREREG_MD5, "frozen pre-registration md5 mismatch"
    assert md5(RC / "docs_AMENDMENT_A.md") == AMEND_MD5, "amendment A md5 mismatch"
    done = []
    for st in a.states.split(","):
        ok, alive, pend = ready(st)
        log(f"{st} ready={ok} alive={len(alive)} pending={len(pend)}" + ("" if ok else f" {pend[:3]}"))
        if a.check_only or (not ok and not a.force):
            continue
        T = score_one(st, a.draws, log); done.append(st)
        A.D._rd.cache_clear()          # free the state's FIA tables before the next state
    if done:
        tabs = {n: [] for n in ("report_card", "components", "membership", "stress", "scalars", "five_seed", "s8_treemap")}
        for st in STATES:
            for n in tabs:
                f = SC / st / f"{n}.csv"
                if f.exists():
                    tabs[n].append(pd.read_csv(f).assign(state=st))
        for n, L in tabs.items():
            if L:
                pd.concat(L, ignore_index=True).to_csv(SC / f"six_state_{n}.csv", index=False, float_format="%.6g")
        bad = leak_scan(SC); (SC / "leak_scan.json").write_text(json.dumps(dict(files_flagged=bad), indent=1))
        log(f"six-state tables written; leak scan flagged {len(bad)}")


if __name__ == "__main__":
    main()
