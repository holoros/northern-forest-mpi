"""Track R adapters: read firebreather run outputs and inputs into the scorer bundle (rc_run.score_state).
Every missing input is recorded in bundle['missing'] with the path pattern tried; nothing is invented.
Restricted: the bundle holds plot- and pair-level arrays in memory only.
"""
import glob, json, os, sys
from pathlib import Path
import numpy as np
import pandas as pd
import rc_core as C
import rc_design as D

HOME = Path(os.path.expanduser("~")); RC = HOME / "perseus_me/ncasi_rc"; RUNS = RC / "runs"; INP = RC / "inputs"
CODE = HOME / "perseus_me/code"; ME_STATE = HOME / "perseus_me/me_state"
K = 0.00112085          # lb ac-1 -> Mg ha-1 (gate_pairs.K)
LB2TG = D.LB2TG
REG_MD5 = None          # pinned in engines/IN_SAMPLE_REGISTER.md5 (checked in load_insample)


def _ids(s):
    return pd.Series(s).astype(str).str.strip().str.replace(r"\.0$", "", regex=True)


def _first(patterns):
    for p in patterns:
        g = sorted(glob.glob(str(p), recursive=True))
        if g:
            return Path(g[0])
    return None


def pair_rows(path, arm=None):
    x = pd.read_csv(path, dtype={"pair_id": str})
    if arm is not None and "arm" in x.columns and (x.arm == arm).any():
        x = x[x.arm == arm]
    x["pair_id"] = _ids(x.pair_id).values
    ok = x.ok.astype(str).str.lower().isin(["true", "1"]) if x.ok.dtype != bool else x.ok
    x = x[ok].copy()
    g = x.fvs_surv / x.cycle * K; m = x.fvs_mort / x.cycle * K
    return pd.DataFrame({"gross_growth": g.values, "mortality": m.values, "net_growth": (g - m).values}, index=x.pair_id.values)


def align(rates, ids):
    r = rates[~rates.index.duplicated()]
    return {q: r[q].reindex(ids).values.astype(float) for q in C.QS}


def load_insample(st):
    f = Path(os.environ["RC_REGISTER"]) if os.environ.get("RC_REGISTER") else RC / "engines/IN_SAMPLE_REGISTER.csv"  # post hoc sensitivity override (2026-11-07)
    if not f.exists():
        return set(), "IN_SAMPLE_REGISTER.csv missing"
    import hashlib
    md5 = hashlib.md5(f.read_bytes()).hexdigest()
    pin = RC / "engines/IN_SAMPLE_REGISTER.md5"
    note = f"register md5 {md5}" + ("" if not pin.exists() else (" (pinned match)" if md5 in pin.read_text() else " (PINNED MD5 MISMATCH)"))
    R = pd.read_csv(f)
    R = R[(R.state.astype(str) == st) | (R.state.astype(str).isin(["all", "*"]))]
    R = R[R.arm.astype(str) != "ncasi_recal"]            # the NCASI arm's calls are applied by its role (v2 item 5)
    flag = R["call"].astype(str).str.strip().str.lower().eq("in")   # 'out_flag' is scored with a flag (AcadianGY, LANDIS ME)
    comp_map = {"C3": "P2", "hindcast": "P2"}
    return {(str(e), comp_map.get(str(c), str(c))) for e, c in zip(R[flag].engine, R[flag].component)}, note


def load_state(st, B=1000, seed=20261105, members=None, log=print):
    miss = []; dr = D.LocDraws(st, B, seed); ev_cur = D.statecd(st) * 10000 + 2401; ev_grm = ev_cur + 2
    loc = D.loc_of(st)
    I = INP / st / "restricted"
    folds = pd.read_csv(I / "folds/folds_locations_appended.csv", dtype=str).drop_duplicates("loc").set_index("loc").fold.astype(int)
    pf = pd.read_csv(I / "folds/pair_folds.csv", dtype=str); pf["pair_id"] = _ids(pf.pair_id).values
    # ---- gate pairs (common set) and A1/A1b observed
    g = pf[pf.set == "gate"].drop_duplicates("pair_id")
    a1 = pd.read_parquet(I / "a1/a1_pairs.parquet"); a1.index = _ids(a1.pair_id).values
    ids = g.pair_id.values; locs = g.t1_id.map(loc).values
    obs = {"gross_growth": a1.obs_gg_a1.reindex(ids).values.astype(float), "mortality": a1.obs_mort.reindex(ids).values.astype(float),
           "net_growth": a1.obs_net_a1b.reindex(ids).values.astype(float)}
    pairs = dict(ids=ids, labels=g.domain.values, fold=g.fold.astype(int).values, obs=obs, M=dr.mult(locs), locs=locs)
    # ---- treated
    t = pf[pf.set == "treated"].drop_duplicates("pair_id"); tp = pd.read_csv(I / "treated/treated_pairs.csv", dtype={"pair_id": str})
    tp.index = _ids(tp.pair_id).values; tid = t.pair_id.values; tl = t.t1_id.map(loc).values
    tobs = {"gross_growth": tp.obs_gg_a1.reindex(tid).values.astype(float), "mortality": tp.obs_mort.reindex(tid).values.astype(float),
            "net_growth": tp.obs_net_a1b.reindex(tid).values.astype(float)}
    reg = ~tp.t40.reindex(tid).astype(str).str.lower().isin(["true", "1"]).values
    treated = dict(ids=tid, labels=t.domain.values, fold=t.fold.astype(int).values, obs=tobs, M=dr.mult(tl), subset_mask=reg, locs=tl)
    # ---- domains and area weights (2024 forest area by domain of the plot's largest forest condition)
    domains = sorted(set(pairs["labels"]) - {"other"}) + (["other"] if "other" in set(pairs["labels"]) else [])
    aw = _area_weights(st, ev_cur, domains)
    # ---- foliage share phi (point; FIA ratio, fixed across draws: disclosed)
    ph = pd.read_csv(RUNS / "common/phi_states.csv").set_index("state").phi
    phi = (1 - float(ph[st]), np.full(B, 1 - float(ph[st])))
    # ---- C1 targets per draw on the GRM evaluation, and the model design (2024 current evaluation)
    yg, dg = D.plot_grm(st, ev_grm); Mg = dr.mult(dg["loc"].values)
    tgt = {q: D.design_total(dg, yg[q], Mg) for q in C.QS}
    ys, dc = D.plot_stock(st, ev_cur); Mc = dr.mult(dc["loc"].values)
    tgt["stock"] = D.design_total(dc, ys, Mc)
    c1_fold = dc["loc"].map(fold_map(st)).fillna(0).astype(int).values      # D-F3: statewide-extended fold file
    if (c1_fold == 0).any():
        miss.append(f"C1 fold: {int((c1_fold == 0).sum())} plots without a fold (scalar s1 applied)")
    # ---- P5 targets: FIA stock per series evaluation (+ removals, A1)
    ser = pd.read_csv(INP / "trend/trend_stock_series.csv"); ser = ser[ser.state == st].sort_values("end_invyr")
    years = [int(y) for y in ser.end_invyr]; stock_s, rem_s = {}, {}
    for ev, y in zip(ser.evalid.astype(int), years):
        v, d = D.plot_stock(st, ev); stock_s[y] = D.design_total(d, v, dr.mult(d["loc"].values))
    for ev3, y in D.evalids(st, "EXPGROW"):
        if years and years[0] < y <= years[-1]:
            try:
                v, d = D.plot_grm(st, ev3); rem_s[y] = D.design_total(d, v["removals"], dr.mult(d["loc"].values))
            except Exception as ex:  # noqa: BLE001
                miss.append(f"removals {ev3}: {ex!r}")
    import rc_comp as P
    p5_target = {"raw": stock_s, "A1": P.removals_adjusted(stock_s, rem_s, years) if rem_s else stock_s}
    if not rem_s:
        miss.append("P5 removals series empty: A1 target equals raw series")
    # ---- engines
    engines = _engines(st, pairs, treated, dc, Mc, ys, miss)
    if "fvs" in engines:                                   # D-F2: P5 recal decomposition (FVS stores gross per stand)
        comp, note = fvs_p5_comp(st, years)
        engines["fvs"]["P5_comp"] = comp
        miss.append(note)
    ins, ins_note = load_insample(st)
    bundle = dict(state=st, B=B, domains=domains, aw=aw, phi=phi, pairs=pairs, treated=treated, chains=_chains(st, dr, miss),
                  c1_target=tgt, c1_design=(dc, Mc), c1_fold=c1_fold, p5_years=years, p5_target=p5_target,
                  engines=engines, insample=ins, insample_note=ins_note, members=[m for m in (members or C.MEMBERS[st]) if m in engines],
                  missing=miss, n_outside_universe=dr.n_outside, units=dict(pairs=len(ids), treated=len(tid)))
    for m in (members or C.MEMBERS[st]):
        if m not in engines:
            miss.append(f"member {m}: no outputs found")
    return bundle


def _area_weights(st, ev, domains):
    d = D.design(st, ev)
    cd = D._rd(st, "COND", ("PLT_CN", "CONDID", "COND_STATUS_CD", "CONDPROP_UNADJ", "FORTYPCD", "PROP_BASIS"))
    cf = cd[cd.COND_STATUS_CD == 1].merge(d[["PLT_CN", "ADJ_FACTOR_SUBP", "ADJ_FACTOR_MACR", "c_h", "P2POINTCNT"]], on="PLT_CN")
    cf["a"] = cf.CONDPROP_UNADJ * np.where(cf.PROP_BASIS == "MACR", cf.ADJ_FACTOR_MACR, cf.ADJ_FACTOR_SUBP) * cf.c_h / cf.P2POINTCNT
    ref = pd.read_csv(_first([HOME / "perseus_me/ext/FIA/*REF_FOREST_TYPE.csv", HOME / "perseus_me/ext/FIA/REF_FOREST_TYPE*.csv"]))
    ft2g = dict(zip(ref.VALUE.astype(int), ref.TYPGRPCD.astype(int)))
    cf["dom"] = [str(ft2g.get(int(x), 999)) if pd.notna(x) else "999" for x in cf.FORTYPCD]
    cf["dom"] = np.where(cf.dom.isin(domains), cf.dom, "other")
    a = cf.groupby("dom").a.sum().reindex(domains).fillna(0.0).values
    return a / a.sum()


def _chains(st, dr, miss):
    """hindcast chain frames per horizon (score.h1_setup / h1_frame, Maine arithmetic), observed only; model filled per engine"""
    try:
        sys.path.insert(0, str(CODE / "tools/score")); sys.path.insert(0, str(CODE / "workflow/scripts")); sys.path.insert(0, str(CODE / "engines/fvs"))
        sys.path.insert(0, str(CODE / "tools/hindcast"))
        import score as S
        H = S.h1_setup(str(CODE / f"config/ncasi_rc/hindcast/{st}.yaml"), str(INP / st / "restricted/hindcast"))
    except Exception as ex:  # noqa: BLE001
        miss.append(f"hindcast setup: {ex!r}"); return {}
    return dict(_H=H, _S=S, _dr=dr)


def chain_frames(bundle, eng, arm, values_dir):
    """per horizon model (y) and observed (x) chain stock; fills bundle['chains'][h] observed on first call"""
    ch = bundle["chains"]
    if not ch:
        return None
    H, S, dr = ch["_H"], ch["_S"], ch["_dr"]; out = {}
    loc = D.loc_of(bundle["state"])
    for h in ("h10", "h20"):
        try:
            fr, nf = S.h1_frame(H, values_dir, eng, arm, h)
        except Exception as ex:  # noqa: BLE001
            bundle["missing"].append(f"h1_frame {eng} {arm} {h}: {ex!r}"); continue
        # D-G2 (2026-11-06): a chain can have two visits in one horizon bin (short REMPER); keep its first visit, as _ensure_obs does
        if "visit" not in fr.columns:
            fr = fr.assign(_o=range(len(fr)))
        fr = fr.sort_values(["chain_id"] + (["visit"] if "visit" in fr.columns else ["_o"])).drop_duplicates("chain_id"); _ensure_obs(bundle)
        y = pd.Series(fr.y.values, index=fr.chain_id.values).reindex(ch[h]["ids"]).values.astype(float)
        out[h] = y
    return out


def _ensure_obs(bundle):
    """observed chain frames per horizon, independent of any engine: primary chain set, offsets binned, nvp keep rule"""
    ch = bundle["chains"]; H, dr = ch["_H"], ch["_dr"]; O, OFF, h1, nvp = H["O"], H["OFF"], H["h1"], H["nvp"]; h1s = set(h1)
    loc = D.loc_of(bundle["state"])
    M = O[O.chain_id.isin(h1s)].copy(); M["offset"] = [OFF.get((c, v), np.nan) for c, v in zip(M.chain_id, M.visit)]
    M["bin"] = [H["AH"].bin_of(o, H["bins"]) for o in M.offset]
    for h in ("h10", "h20"):
        if h in ch:
            continue
        sub = M[M.bin == h]
        keep = [c for c, v in zip(sub.chain_id, sub.visit) if c in set(nvp[nvp >= v].index)]
        sub = sub[sub.chain_id.isin(set(keep))].sort_values(["chain_id", "visit"]).drop_duplicates("chain_id")  # D-G2: first visit in bin, order-stable
        fr = sub[["chain_id"]].copy(); locs = fr.chain_id.map(H["v1plt"]).astype(str).map(loc).values
        ch[h] = dict(ids=fr.chain_id.values, obs=sub.agc_stock.values.astype(float), M=dr.mult(locs), labels=_chain_domains(bundle, fr))


def chain_frames_level(bundle, path):
    """chain-level hindcast values (chain_id, offset, agb_full, agc_full, ...) -> per horizon y/x, the Maine ext23 frame_all
    arithmetic (aggregate_all.chain_level interpolation, horizon cap 30, primary chain set, nvp keep rule)"""
    ch = bundle["chains"]
    if not ch:
        return None
    H, S, dr = ch["_H"], ch["_S"], ch["_dr"]
    sys.path.insert(0, str(CODE / "tools/hindcast_all")); import aggregate_all as AA
    O, OFF, h1, nvp = H["O"], H["OFF"], H["h1"], H["nvp"]; h1s = set(h1)
    Pq = pd.read_parquet(path); out = {}; loc = D.loc_of(bundle["state"])
    _ensure_obs(bundle)
    Pr = AA.chain_level(Pq, {k: v for k, v in OFF.items() if k[0] in set(Pq.chain_id)}, 30)
    failed = h1s - set(Pq.chain_id)
    M = O.merge(Pr, on=["chain_id", "visit"], how="left", suffixes=("_o", "_p")); M = M[~M.chain_id.isin(failed)].copy()
    M["offset"] = [OFF.get((c, v), np.nan) for c, v in zip(M.chain_id, M.visit)]
    M["bin"] = [H["AH"].bin_of(o, H["bins"]) for o in M.offset]
    for h in ("h10", "h20"):
        sub = M[M.bin == h]
        keep = [c for c, v in zip(sub.chain_id, sub.visit) if c in set(nvp[nvp >= v].index) & h1s]
        sub = sub[sub.chain_id.isin(set(keep)) & sub.agc_stock_p.notna()].sort_values(["chain_id", "visit"]).drop_duplicates("chain_id")  # D-G2
        out[h] = pd.Series(sub.agc_stock_p.values, index=sub.chain_id.values).reindex(ch[h]["ids"]).values.astype(float)
    return out


def _chain_domains(bundle, fr):
    """P2 by-domain view (S10): forest-type group of the chain's visit 1 plot (largest forest condition), mapped to the state's
    domains as the pairs are (other groups -> 'other')"""
    st = bundle["state"]; H = bundle["chains"]["_H"]
    cd = D._rd(st, "COND", ("PLT_CN", "CONDID", "COND_STATUS_CD", "CONDPROP_UNADJ", "FORTYPCD"))
    c = cd[cd.COND_STATUS_CD == 1].sort_values(["PLT_CN", "CONDPROP_UNADJ", "CONDID"], ascending=[True, False, True]).drop_duplicates("PLT_CN")
    ref = pd.read_csv(_first([HOME / "perseus_me/ext/FIA/*REF_FOREST_TYPE.csv", HOME / "perseus_me/ext/FIA/REF_FOREST_TYPE*.csv"]))
    ft2g = dict(zip(ref.VALUE.astype(int), ref.TYPGRPCD.astype(int)))
    grp = {str(p): str(ft2g.get(int(f), 999)) if pd.notna(f) else "999" for p, f in zip(c.PLT_CN.astype(str), c.FORTYPCD)}
    plt = fr.chain_id.map(H["v1plt"]).astype(str)
    return np.array([g if g in bundle["domains"] else "other" for g in plt.map(grp).fillna("999")])


def _engines(st, pairs, treated, dc, Mc, ys, miss):
    E = {}
    ids, tid = pairs["ids"], treated["ids"]
    def rows(pats, arm=None):
        f = _first(pats)
        return (None, None) if f is None else (pair_rows(f, arm), f)
    sp = lambda eng: RUNS / eng / st / "restricted"  # noqa: E731
    # tree-list engines
    tl = {"fvs": ("pairs_all_default", "treated_default", "default"), "osm": ("pairs_all", "treated", None), "acd141": ("pairs_all", "treated", None)}
    for eng, (pdir, tdir, arm) in tl.items():
        if eng != "fvs" and st != "ME":
            continue
        r2, f2 = rows([sp(eng) / pdir / "pair_rows.csv"], arm); r5, f5 = rows([sp(eng) / tdir / "pair_rows.csv"], arm)
        if r2 is None and r5 is None:
            miss.append(f"{eng}: no pair rows under {sp(eng)}"); continue
        e = dict(engine=eng, arm=arm or "registered", role="member", C2=align(r2, ids) if r2 is not None else None,
                 C5=align(r5, tid) if r5 is not None else None, src=dict(C2=str(f2), C5=str(f5)))
        e["C1"] = _treelist_c1(st, eng, dc, Mc, ys, miss)
        e["P5"] = _p5_total(st, eng, miss)
        e["C3_dir"] = _first([sp(eng) / "h1"]); e["C3_arm"] = arm or "registered"
        E[eng] = e
    # OSM sensitivity rows in NH and VT (BGI 3500), not members
    if st in ("NH", "VT"):
        r2, _ = rows([sp("osm") / "pairs_all/pair_rows.csv"]); r5, _ = rows([sp("osm") / "treated/pair_rows.csv"])
        if r2 is not None:
            E["osm_bgi3500"] = dict(engine="osm", arm="bgi3500", role="sensitivity", C2=align(r2, ids), C5=align(r5, tid) if r5 is not None else None,
                                    note="sensitivity: variant default BGI 3500, outside the site input domain")
    # NCASI recalibrated FVS arm: holdout pairs from plots with no fit-half pair (v2 item 5)
    rh, fh = rows([sp("fvs") / "pairs_hold_study_recal/pair_rows.csv"])
    if rh is not None:
        hold = set(rh.index); fit_locs = {l for i, l in zip(ids, pairs["locs"]) if i not in hold}
        keep = np.array([(i in hold) and (l not in fit_locs) for i, l in zip(ids, pairs["locs"])])
        c2 = align(rh, ids); c2 = {q: np.where(keep, v, np.nan) for q, v in c2.items()}
        E["fvs_ncasi_recal"] = dict(engine="fvs", arm="ncasi_recal", role="ncasi_arm", C2=c2, C5=None,
                                    note=f"NCASI registered recalibration; holdout pairs from plots with no fit-half pair ({int(keep.sum())})")
    # CEM: pair rows written by tools/v23_E/cem_pair_extract_v3.R (run_T2_worker.sh), arm cem_shared; A8 from cem_a8.py
    r2, f2 = rows([sp("cem") / "pairs_all/extract/pair_rows.csv"], "cem_shared")
    r5, f5 = rows([sp("cem") / "treated/extract/pair_rows.csv"], "cem_shared")
    if r2 is not None or r5 is not None:
        flag, E_a8 = "", np.nan
        a8 = RUNS / "cem/cem_a8_shares.csv"
        if a8.exists():
            x = pd.read_csv(a8).set_index("state")
            if st in x.index and "treated_match_share" in x.columns and np.isfinite(x.loc[st, "treated_match_share"]):
                E_a8 = float(x.loc[st, "treated_match_share"])
                if E_a8 < 0.5:
                    flag = f"mostly untreated donors (A8, treated match share {E_a8:.2f})"
            else:
                miss.append(f"cem: no A8 row for {st} in cem_a8_shares.csv")
        else:
            miss.append("cem: runs/cem/cem_a8_shares.csv not found (run tools/ncasi_rc/eng/cem_a8.py after the treated runs)")
        E["cem"] = dict(engine="cem", arm="cem_shared", role="member", C2=align(r2, ids) if r2 is not None else None,
                        C5=align(r5, tid) if r5 is not None else None, C5_flag=flag, a8_share=E_a8, src=dict(C2=str(f2), C5=str(f5)))
        f = sorted(glob.glob(str(sp("cem") / "values_c1/plot_gate_f*.csv.gz")))
        if f:
            y = pd.concat([pd.read_csv(x, dtype={"PLT_CN": str}) for x in f]).drop_duplicates("PLT_CN").set_index("PLT_CN")
            k = dc.PLT_CN.astype(str)
            E["cem"]["C1"] = dict(kind="plot", y={"gross_growth": k.map(y.gross).fillna(0).values, "mortality": k.map(y.mort).fillna(0).values,
                                                  "net_growth": k.map(y.net).fillna(0).values}, coverage=float(k.isin(y.index).mean()))
        else:
            miss.append("cem: values_c1/plot_gate_f*.csv.gz not found")
        E["cem"]["P5"] = _p5_total(st, "cem", miss)
    else:
        miss.append("cem: no pair rows (runs/cem/<ST>/restricted/{pairs_all,treated}/extract/pair_rows.csv)")
    # LANDIS-II (pair rows already x (1 - phi), registered_tables arm)
    r2, f2 = rows([sp("landis") / "pairs_sample/rows*/pair_rows.csv", sp("landis") / "pairs*/rows*/pair_rows.csv"])
    r5, f5 = rows([sp("landis") / "c5/rows_cut/pair_rows.csv"])
    if r2 is not None or r5 is not None:
        a6 = _first([sp("landis") / "c5/a6b.json"]); a6p = json.loads(a6.read_text()).get("a6b_pass") if a6 else None
        e = dict(engine="landis", arm="registered_tables", role="member", C2=align(r2, ids) if r2 is not None else None,
                 C5=align(r5, tid) if r5 is not None else None, A6b_pass=a6p, foliage_done=True)
        c1 = _first([sp("landis") / "statewide_fia_extract.c1.json", sp("landis") / "statewide_fia*.c1.json"])
        if c1 is not None:
            j = json.loads(c1.read_text())
            e["C1"] = dict(kind="total", tot={"gross_growth": j["gross"] / LB2TG, "mortality": j["mort"] / LB2TG, "net_growth": j["net"] / LB2TG})
        else:
            miss.append("landis: statewide c1.json not found")
        E["landis"] = e
    else:
        miss.append("landis: no pair rows")
    # curve-driven engines
    sw = _first([RUNS / "statewide" / st / "statewide_curves.csv"])
    SW = pd.read_csv(sw) if sw is not None else None
    for eng in ("cbm_gated", "cbm_af", "fia_curves"):
        r2, _ = rows([sp(eng) / f"pairs_level__{eng}__level__pair_rows.csv"])
        r5l, _ = rows([sp(eng) / f"treated_level__{eng}__level__pair_rows.csv"]); r5a, _ = rows([sp(eng) / f"treated_ageq__{eng}__ageq__pair_rows.csv"])
        if r2 is None and r5l is None:
            continue
        e = dict(engine=eng, arm={"cbm_gated": "FVS_FIAW_NEW", "cbm_af": "papa2023", "fia_curves": "hybrid_fits"}[eng], role="member",
                 C2=align(r2, ids) if r2 is not None else None,
                 C5={"level": align(r5l, tid) if r5l is not None else None, "ageq": align(r5a, tid) if r5a is not None else None})
        if SW is not None:
            c = SW[(SW.engine == eng) & (SW.component == "C1")].set_index("quantity").value
            qmap = {"gross": "gross_growth", "mortality": "mortality", "net": "net_growth"}
            tot = {qmap[k]: float(v) / LB2TG for k, v in c.items() if k in qmap and np.isfinite(v)}
            if eng == "fia_curves":
                tot = {k: v for k, v in tot.items() if k == "net_growth"}          # curves carry no mortality: net only
            st_ = c.get("ag_live_stock_2024")
            if st_ is not None and np.isfinite(st_):
                tot["stock"] = float(st_) / LB2TG                                   # curve stock is model derived: exposed in C1
            e["C1"] = dict(kind="total", tot=tot)
            s5 = SW[(SW.engine == eng) & (SW.component == "P5_series") & (SW.quantity == "stock_no_harvest")]
            if len(s5) and eng != "fia_curves":
                e["P5"] = dict(kind="total", series={int(y): float(v) / LB2TG for y, v in zip(s5.year, s5.value)})
        E[eng] = e
    return E


def _treelist_c1(st, eng, dc, Mc, ys, miss):
    """FVS: run_gate stand values (score.fvs_c1 arithmetic); OSM/AcadianGY ME: Maine ci_t2 stand results (score.treelist_c1)"""
    sys.path.insert(0, str(CODE / "tools/score")); sys.path.insert(0, str(CODE / "workflow/scripts")); sys.path.insert(0, str(CODE / "engines/fvs"))
    sys.path.insert(0, str(CODE / "tools/hindcast"))
    import score as S
    k = dc.PLT_CN.astype(str)
    if eng == "fvs":
        f = _first([RUNS / "state_roots" / st / "restricted/fvs_c1/fvs_gate_stand_values.parquet"] +
                   ([ME_STATE / "score/restricted/fvs_gate/fvs_gate_stand_values.parquet"] if st == "ME" else []))
        if f is None:
            miss.append(f"fvs C1: stand values not found"); return None
        s = S.fvs_c1(f)["default"]
        return dict(kind="plot", y={q: k.map(s[q]).fillna(0).values for q in C.QS}, src=str(f))
    if st != "ME":
        return None
    f = _first([ME_STATE / f"ci_t2_engines/stand_results/stand_results__{eng}__me_state_t2_{eng}.parquet"])
    if f is None:
        miss.append(f"{eng} C1: Maine stand results not found"); return None
    R = pd.read_parquet(f); plt = R.stand.str.rsplit("_", n=1).str[0]
    stock0 = R.stock0.groupby(plt).sum()
    sc = D.design_total(dc, ys)[0] / D.design_total(dc, k.map(stock0).fillna(0).values)[0]
    g = (R.gate5_gross + R.gate10_gross).groupby(plt).sum() / 10.0 * sc
    m = (R.gate5_mort_dying_tree + R.gate10_mort_dying_tree).groupby(plt).sum() / 10.0 * sc
    return dict(kind="plot", y={"gross_growth": k.map(g).fillna(0).values, "mortality": k.map(m).fillna(0).values,
                                "net_growth": k.map(g - m).fillna(0).values}, src=str(f))


def _p5_total(st, eng, miss):
    f = _first([RUNS / eng / st / "statewide_p5/p5_series_model.csv", RUNS / eng / st / "**/p5_series_model.csv"])
    if f is None:
        miss.append(f"{eng} P5: p5_series_model.csv not found"); return None
    x = pd.read_csv(f)
    if "arm" in x.columns and (x.arm == "default").any():
        x = x[x.arm == "default"]
    col = "stock_Tg" if "stock_Tg" in x.columns else [c for c in x.columns if c.startswith("stock")][0]
    return dict(kind="total", series={int(y): float(v) / LB2TG for y, v in zip(x.year, x[col])}, src=str(f))


def fill_c3(bundle):
    """P2: per engine chain stock per horizon from the hindcast values (Maine h1_frame arithmetic); every engine/arm looked up as
    <eng>__<arm>.parquet in its h1 folder or under runs/hindcast/<ST>/"""
    st = bundle["state"]
    for n, E in bundle["engines"].items():
        if E.get("role") == "ncasi_arm":
            arm = "study_recal"; d = RUNS / "fvs" / st / "restricted/h1"
        else:
            arm = E.get("C3_arm") or E.get("arm"); d = E.get("C3_dir")
        lvl = sorted(glob.glob(str(RUNS / "hindcast" / st / "**" / f"{E['engine']}__{arm}__chain.parquet"), recursive=True)) + \
              sorted(glob.glob(str(RUNS / E["engine"] / st / "restricted/**" / f"{E['engine']}__{arm}__chain.parquet"), recursive=True))
        if lvl:                                   # chain-level values (curve, CBM, LANDIS hindcasts); anchored file, not 'unanchored'
            v = chain_frames_level(bundle, lvl[0])
            if v:
                E["C3"] = v
            continue
        cands = [d] if d else []
        cands += [Path(p).parent for p in glob.glob(str(RUNS / "hindcast" / st / "**" / f"{E['engine']}__{arm}.parquet"), recursive=True)]
        cands += [Path(p).parent for p in glob.glob(str(RUNS / E["engine"] / st / "restricted/**" / f"{E['engine']}__{arm}.parquet"), recursive=True)]
        for c in cands:
            if c is not None and (Path(c) / f"{E['engine']}__{arm}.parquet").exists():
                v = chain_frames(bundle, E["engine"], arm, str(c))
                if v:
                    E["C3"] = v
                if E["engine"] == "fvs" and E.get("role") == "member":       # D-F2: P2 recal decomposition
                    E["C3_comp"] = fvs_c3_comp(bundle, str(c), arm)
                    if E["C3_comp"] is None:
                        bundle["missing"].append("fvs P2 recal: decomposition identity failed or folds missing")
                break
        else:
            bundle["missing"].append(f"{n}: no hindcast values {E['engine']}__{arm}.parquet")
    if bundle["chains"]:
        bundle["chains"] = {h: v for h, v in bundle["chains"].items() if not h.startswith("_")}
    return bundle


def attach_p3(bundle, work_root, log=print):
    """registered P3 and P3s for every engine row and ensemble (rc_p3), A6 gate; writes nothing but group-level work files"""
    import rc_p3 as P3
    eng = {n: E["engine"] for n, E in bundle["engines"].items() if E.get("role") in ("member", "sensitivity")}
    mem = bundle["members"]; plot = [m for m in mem if bundle["engines"][m]["engine"] in C.PLOT_BASED]
    ens = {"E_state": (mem, "mean"), "E_plot": (plot, "mean"), "E_state_median": (mem, "median")}
    p3, p3s, ensp, rep = P3.state_p3(bundle["state"], eng, ens, work_root)
    for n, v in p3.items():
        bundle["engines"][n]["P3"] = v
    for n, v in p3s.items():
        bundle["engines"][n]["P3s"] = v
    ensp.update({("E_recal", k): v for (e, k), v in ensp.items() if e == "E_state"})   # recalibration does not act on realism (v3 item 7)
    bundle["ens_p3"] = ensp; bundle["p3_report"] = rep
    if rep.get("a6", {}).get("ok") is False:
        bundle["missing"].append(rep["a6"]["note"])
    return bundle


# ---------------------------------------------------------------------------------------------------- D-F2 / D-F3 helpers
def fold_map(st):
    """location -> fold: the statewide-extended fold file (folds_locations_sw.csv, v4 item 7d) where it exists, else the
    appended file. D-F3 (2026-11-06): C1 and P5 locations absent from the appended file previously got fold 0."""
    I = INP / st / "restricted/folds"
    f = I / "folds_locations_sw.csv"
    if not f.exists():
        f = I / "folds_locations_appended.csv"
    return pd.read_csv(f, dtype=str).drop_duplicates("loc").set_index("loc").fold.astype(int)


def fvs_p5_comp(st, years, check_tol=1e-6):
    """D-F2: P5 recal decomposition for FVS from the restricted stand values (tools/ncasi_rc/sw/fvs_p5.py arithmetic):
    year -> (S_y, G_f1, G_f2) in lb (the scorer's series units), G = cumulative gross since the stand's start, per fold.
    Gate: S_y rebuilt here equals p5_series_model.csv within check_tol relative, else None (reported missing)."""
    sys.path.insert(0, str(CODE / "workflow/scripts")); from fiaest import Design  # noqa: E402
    rd = RUNS / "fvs" / st / "restricted/statewide_p5/fvs_p5_stand_values.parquet"
    sm = RUNS / "fvs" / st / "statewide_p5/p5_series_model.csv"
    root = RUNS / "sw_roots" / f"{st}_p5"
    if not (rd.exists() and sm.exists() and root.exists()):
        return None, "fvs P5 recal: inputs missing"
    aoi = json.loads((root / "aoi/aoi.json").read_text())
    des = Design(root / "restricted/inventory", aoi["evalids"][st]["current"])
    P = pd.read_parquet(rd); P = P[P.arm == "default"].sort_values(["stand", "offset"]).copy()
    P["cg"] = P.groupby("stand").gross.transform(lambda x: x.fillna(0).cumsum()); P["cal"] = P.INVYR + P.offset
    M = pd.read_csv(sm); M = M[M.arm == "default"].set_index("year")
    loc = D.loc_of(st); fm = fold_map(st)
    pl_fold = P.drop_duplicates("stand").set_index("stand").PLT_CN.astype(str).map(loc).map(fm)
    if pl_fold.isna().any():
        return None, f"fvs P5 recal: {int(pl_fold.isna().sum())} stands without a fold"
    LB2TG_ = 0.45359237e-9; out = {}; worst = 0.0
    g = {k: v for k, v in P.groupby("stand")}
    for yr in years:
        if yr not in M.index:
            continue
        sc = float(M.loc[yr, "stock_scale"])
        st_v = {s: np.interp(yr, v.cal.values, v.stock.values) for s, v in g.items()}
        cg_v = {s: np.interp(yr, v.cal.values, v.cg.values) for s, v in g.items()}
        idx = list(g); plt = P.drop_duplicates("stand").set_index("stand").loc[idx].PLT_CN.values
        sv = pd.Series([st_v[s] for s in idx], index=plt).groupby(level=0).sum()
        f = pl_fold.loc[idx].values
        c1 = pd.Series([cg_v[s] if ff == 1 else 0.0 for s, ff in zip(idx, f)], index=plt).groupby(level=0).sum()
        c2 = pd.Series([cg_v[s] if ff == 2 else 0.0 for s, ff in zip(idx, f)], index=plt).groupby(level=0).sum()
        Sy = sc * des.total(sv)[0]
        rel = abs(Sy * LB2TG_ - float(M.loc[yr, "stock_Tg"])) / float(M.loc[yr, "stock_Tg"]); worst = max(worst, rel)
        out[int(yr)] = (Sy, sc * des.total(c1)[0], sc * des.total(c2)[0])
    if worst > check_tol:
        return None, f"fvs P5 recal: rebuilt series differs from p5_series_model.csv by {worst:.2e} relative"
    return out, f"fvs P5 recal decomposition rebuilt, max rel diff {worst:.1e}"


def fvs_c3_comp(bundle, values_dir, arm="default"):
    """D-F2: P2 recal decomposition for FVS from the tree-level hindcast values (aggregate_h1.projected arithmetic, carbon
    with cf1): per horizon (stock1, G, Mo, fold) aligned to the chain ids, where y = stock1 + G - Mo reproduces the registered
    projected agc_stock exactly (checked; returns None on mismatch)."""
    ch = bundle["chains"]
    if not ch or "_H" not in ch:
        return None
    H = ch["_H"]; AH = H["AH"]
    Pv = pd.read_parquet(Path(values_dir) / f"fvs__{arm}.parquet")
    W = H["W"]
    P = Pv.merge(W[["w1", "cf1"]].rename_axis("tree_cn").reset_index(), on="tree_cn", how="left")
    tpa1 = P[P.visit == 1].set_index("tree_cn").tpa
    P["k"] = P.w1 / P.tree_cn.map(tpa1); P = P.sort_values(["chain_id", "tree_cn", "visit"])
    gg = P.groupby(["chain_id", "tree_cn"], sort=False)
    P["tpa_prev"] = gg.tpa.shift(1); P["b_prev"] = gg.b.shift(1)
    P["cs"] = P.tpa * P.b * P.k * P.cf1
    P["cg"] = (P.tpa * (P.b - P.b_prev) * P.k * P.cf1).fillna(0.0)
    P["cm"] = ((P.tpa_prev - P.tpa) * P.b_prev * P.k * P.cf1).fillna(0.0)
    Q = P.groupby(["chain_id", "visit"]).agg(cs=("cs", "sum"), cg=("cg", "sum"), cm=("cm", "sum")).reset_index().sort_values(["chain_id", "visit"])
    Q["G"] = Q.groupby("chain_id").cg.cumsum(); Q["Mo"] = Q.groupby("chain_id").cm.cumsum()
    S = ch["_S"]; loc = D.loc_of(bundle["state"]); fm = fold_map(bundle["state"]); out = {}
    for h in ("h10", "h20"):
        if h not in ch:
            continue
        fr = _h1_frame_v(H, Pv, h)
        fr = fr.sort_values(["chain_id", "visit"]).drop_duplicates("chain_id")
        q = fr.merge(Q, on=["chain_id", "visit"], how="left")
        cs1 = Q[Q.visit == 1].set_index("chain_id").cs
        stock1 = q.chain_id.map(cs1)                 # independent of y: the identity check below is not an echo (F15)
        ids = ch[h]["ids"]
        if len(ids) == 0:                       # empty horizon (e.g. VT h20): nothing to decompose
            out[h] = tuple(np.array([], float) for _ in range(3)) + (np.array([], int),); continue
        al = lambda v: pd.Series(np.asarray(v, float), index=q.chain_id.values).reindex(ids).values  # noqa: E731
        y_chk = al(q.y)
        locs = pd.Series(ids).map(H["v1plt"]).astype(str).map(loc)
        fold = locs.map(fm).fillna(0).astype(int).values
        if (fold == 0).any():
            return None
        out[h] = (al(stock1), al(q.G), al(q.Mo), fold)
        # identity: stock1 + G - Mo == y
        d = np.nanmax(np.abs(out[h][0] + out[h][1] - out[h][2] - y_chk)) if np.isfinite(y_chk).any() else 0.0
        if not d < 1e-6 * max(1.0, np.nanmax(np.abs(y_chk))):
            return None
    return out


def _h1_frame_v(H, Pv, hb):
    """score.h1_frame arithmetic (projected agc_stock per chain in horizon bin hb) keeping the visit column (D-F2)"""
    AH, V, O, Tx, h1 = H["AH"], H["V"], H["O"], H["Tx"], H["h1"]
    ok = set(Pv.chain_id); Pr = AH.projected(Pv, H["W"])
    M = O.merge(Pr, on=["chain_id", "visit"], how="left", suffixes=("_o", "_p"))
    v1 = V[V.visit == 1]
    failed = set(h1) - ok - set(v1.chain_id[~v1.chain_id.isin(set(Tx[(Tx.visit == 1) & Tx.live].chain_id))])
    M = M[~M.chain_id.isin(failed)]
    M["offset"] = [H["OFF"].get((c, v), np.nan) for c, v in zip(M.chain_id, M.visit)]
    M["bin"] = [AH.bin_of(o, H["bins"]) for o in M.offset]
    sub = M[M.bin == hb]
    nvp = H["nvp"]; keep = [c for c, v in zip(sub.chain_id, sub.visit) if c in set(nvp[nvp >= v].index) & set(h1)]
    sub = sub[sub.chain_id.isin(set(keep))]
    return pd.DataFrame(dict(chain_id=sub.chain_id.values, visit=sub.visit.values, y=sub.agc_stock_p.values))
