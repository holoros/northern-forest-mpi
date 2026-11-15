"""Track R per-state scoring from a loaded bundle (rc_adapt.load_state builds it from firebreather outputs).
Rows: single engines (members and sensitivity rows), cross-fit recalibrated rows, the NCASI recalibrated FVS arm,
E_state, E_recal, E_plot and the median ensemble. Returns (rows, components table, membership table).
"""
import numpy as np
import pandas as pd
import rc_core as C
import rc_comp as P
import rc_state as R


def _fol(bundle, eng):
    return bundle["phi"] if eng in C.FOLIAGE else None


def _pairs_res(bundle, key, model, eng, opt, quants=None):
    U = bundle[key]
    quants = quants or C2Q(eng, key, opt)
    sel = U.get("subset_mask") if (key == "treated" and opt["treated_subset"] == "registered") else None
    M, lab = U["M"], U["labels"]
    obs = U["obs"]; mod = model
    if sel is not None:
        M = M[:, sel]; lab = np.asarray(lab)[sel]; obs = {q: np.asarray(v)[sel] for q, v in obs.items()}
        mod = {q: (v.subset(sel) if C.is_xf(v) else np.asarray(v)[sel]) for q, v in model.items()}
    return P.pair_component(mod, obs, M, lab, bundle["domains"], bundle["aw"], quants, foliage=_fol(bundle, eng))


def C2Q(eng, key, opt):
    if key == "pairs":
        return C.C2_QUANT.get(eng, C.QS)
    # treated: net only headline (item 3), three quantities in S2 where the engine reports them
    if opt["c5_quant"] == "net" or eng in C.CURVE_DRIVEN and eng == "fia_curves":
        return ("net_growth",)
    return C.QS


def _c1_res(bundle, E, eng):
    c = E.get("C1")
    if c is None:
        return None
    T = bundle["c1_target"]; qs = [q for q in C.QS + ("stock",) if q in (c.get("y") or c.get("tot") or {})]
    if c["kind"] == "plot":
        d, M = bundle["c1_design"]
        mt = {q: P_total(d, c["y"][q], M) for q in qs}
    else:
        mt = {q: (float(c["tot"][q]), np.full(bundle["B"], float(c["tot"][q]))) for q in qs}
    if eng in C.FOLIAGE:
        f = bundle["phi"]; mt = {q: (v[0] * f[0], v[1] * np.asarray(f[1])) for q, v in mt.items()}
    return P.c1_component(mt, {q: T[q] for q in qs}, qs)


def P_total(d, y, M):
    import rc_design as D
    return D.design_total(d, y, M)


def _c3_res(bundle, E, opt):
    c = E.get("C3")
    if not c:
        return None
    ch = bundle["chains"]; hs = [h for h in c if h in ch]
    return P.c3_component({h: c[h] for h in hs}, {h: ch[h]["obs"] for h in hs}, {h: ch[h]["M"] for h in hs},
                          {h: ch[h]["labels"] for h in hs}, bundle["domains"], bundle["aw"], min_chains=opt["min_chains"])


def _p5_res(bundle, E, opt, eng):
    c = E.get("P5")
    if c is None or eng == "fia_curves":
        return None
    yrs = bundle["p5_years"]; tgt = bundle["p5_target"]["A1" if opt["p5_target"] == "A1" else "raw"]
    if c["kind"] == "plot":
        d, M = bundle["p5_design"]
        ser = {y: P_total(d, c["series"][y], M) for y in c["series"]}
    elif c["kind"] == "total_draws":            # D-F2: recalibrated series, (pt, bt) per year with refit scalars
        ser = dict(c["series"])
    else:
        ser = {y: (float(v), np.full(bundle["B"], float(v))) for y, v in c["series"].items()}
    return P.p5_component(ser, tgt, yrs)


def comps_single(bundle, name, opt, override=None):
    E = dict(bundle["engines"][name]); E.update(override or {}); eng = E["engine"]; out = {}
    if E.get("C2") is not None and C.eligible(eng, "C2"):
        out["C2"] = _pairs_res(bundle, "pairs", E["C2"], eng, opt)
    c5 = E.get("C5")
    if c5 is not None and C.eligible(eng, "C5"):
        mod = c5.get(opt["curve_rule"], c5.get("level")) if eng in C.CURVE_DRIVEN else c5
        flag = E.get("C5_flag")
        if eng == "landis" and E.get("A6b_pass") is False:
            out["C5"] = dict(S_sw_pt=np.nan, S_dom_pt=np.nan, flag="not evaluated: A6b removal check failed")
        elif mod is not None:
            out["C5"] = _pairs_res(bundle, "treated", mod, eng, opt); out["C5"]["flag"] = flag or ""
    if E.get("C3") and C.eligible(eng, "P2"):
        out["C3"] = _c3_res(bundle, E, opt)
    if C.eligible(eng, "C1"):
        out["C1"] = _c1_res(bundle, E, eng)
    if C.eligible(eng, "P5") and not (eng == "landis"):          # A2: LANDIS-II P5 not evaluated
        out["P5"] = _p5_res(bundle, E, opt, eng)
    for k in ("P3", "P3s"):
        if E.get(k) is not None and C.eligible(eng, k):
            out[k] = E[k]
    if opt.get("start") == "treemap" and E.get("TM") is not None:   # S8: C1 stock and forward P5 from TreeMap
        out["TM"] = E["TM"]
    return out


# ---------------------------------------------------------------------------------------------------- cross-fit recal
def recal_override(bundle, name):
    """v2 items 1 to 5: scalars refit inside every draw on the C2 common pairs of each fold; every unit gets the other fold's
    scalar; post-processing, no feedback; realism unchanged (v3 item 7). D-F2 (2026-11-06): the per-draw scalars now act in
    every component (XF values), and P2 and P5 are recalibrated where the engine stores a gross / mortality split
    (C3_comp per chain, P5_comp per fold). Returns (override dict, scalar report)."""
    E = bundle["engines"][name]; eng = E["engine"]; tq = C.RECAL_TARGET[eng]; U = bundle["pairs"]
    o = np.asarray(U["obs"][tq], float); m = np.asarray(E["C2"][tq], float)
    s1p, s2p = C.crossfit_scalars(o, m, U["fold"]); s1b, s2b = C.crossfit_scalars(o, m, U["fold"], U["M"])
    rep = dict(row=name + "_recal", target=tq, scalar_fold1=s1p, scalar_fold2=s2p)
    ov = {}
    S = (s1p, s1b, s2p, s2b)

    def unitwise(mod, fold):
        return C.apply_crossfit_xf(mod["gross_growth"], mod["mortality"], fold, *S, tq)

    ov["C2"] = unitwise(E["C2"], U["fold"])
    if E.get("C5") is not None:
        ov["C5"] = unitwise(E["C5"], bundle["treated"]["fold"])
    if E.get("C1") is not None and E["C1"]["kind"] == "plot":
        y = E["C1"]["y"]; f = bundle["c1_fold"]
        ov["C1"] = dict(kind="plot", y=C.apply_crossfit_xf(y["gross_growth"], y["mortality"], f, *S, tq))
    if E.get("C3_comp"):
        c3 = {}
        for h, (stock1, g, mo, fold) in E["C3_comp"].items():
            c3[h] = C.XF(stock1 - mo, [(g, fold, *S)]) if tq == "gross_growth" else C.XF(stock1 + g, [(-mo, fold, *S)])
        ov["C3"] = c3
    else:
        ov["C3"] = None
    if E.get("P5_comp"):
        # P5_comp: dict year -> (S_y, X_f1, X_f2), X = cumulative gross (gross target) or mortality (mortality target) totals
        sgn = 1.0 if tq == "gross_growth" else -1.0; ser = {}
        for y, (Sy, X1, X2) in E["P5_comp"].items():
            pt = Sy + sgn * ((s2p - 1) * X1 + (s1p - 1) * X2)
            bt = Sy + sgn * ((np.asarray(s2b) - 1) * X1 + (np.asarray(s1b) - 1) * X2)
            ser[y] = (float(pt), bt)
        ov["P5"] = dict(kind="total_draws", series=ser)
    else:
        ov["P5"] = None
    rep.update(draw_sd_fold1=float(np.nanstd(s1b)), draw_sd_fold2=float(np.nanstd(s2b)),
               p2_recal=bool(E.get("C3_comp")), p5_recal=bool(E.get("P5_comp")))
    return ov, rep


# ---------------------------------------------------------------------------------------------------- ensembles
def ensemble_comps(bundle, members, opt, overrides, how):
    """per-component ensemble on common units over contributing eligible members, half rule (v2 item 9, v3 items 1 and 8)"""
    out, memb = {}, []
    ins = bundle["insample"]
    def contrib(comp, key):
        el = [m for m in members if C.eligible(bundle["engines"][m]["engine"], comp)]
        have = []
        for m in el:
            eng = bundle["engines"][m]["engine"]
            if (eng, comp) in ins:
                continue
            E = dict(bundle["engines"][m]); E.update(overrides.get(m, {}))
            v = E.get(key)
            if comp == "C5":
                if eng in C.CURVE_DRIVEN and opt["ens_c5_curves"] == "exclude":
                    continue
                if eng == "cem" and E.get("C5_flag"):
                    continue
                if eng == "landis" and E.get("A6b_pass") is False:
                    continue
                if v is not None and eng in C.CURVE_DRIVEN:
                    v = v.get("ageq" if opt["ens_c5_curves"] == "ageq" else "level")
            if v is not None:
                have.append((m, eng, v))
        return el, have
    # C2 and C5: per unit means
    for comp, key, U in (("C2", "C2", "pairs"), ("C5", "C5", "treated")):
        el, have = contrib(comp, key)
        if not have or not C.half_rule(have, el):
            memb.append(dict(component=comp, members=",".join(m for m, _, _ in have), eligible=len(el), evaluated=False)); continue
        qs = C.QS if comp == "C2" else (("net_growth",) if opt["c5_quant"] == "net" else C.QS)
        pred = {}
        for q in qs:
            vals = []
            for m, eng, v in have:
                a = (v[q] if C.is_xf(v[q]) else np.asarray(v[q], float)) if q in v else np.full(len(bundle[U]["labels"]), np.nan)
                if eng in C.FOLIAGE:
                    a = a.scale(bundle["phi"][0]) if C.is_xf(a) else a * bundle["phi"][0]
                vals.append(a)
            vq = [x for x in vals if np.isfinite(np.asarray(x, float)).any()]
            pred[q], _ = C.ens_mean_xf(vq, how)
        out[comp] = _pairs_res(bundle, U, pred, "ensemble", opt, qs)
        memb.append(dict(component=comp, members=",".join(m for m, _, _ in have), eligible=len(el), evaluated=True, units=out[comp]["n"]))
    # C1 and P5: mean of member totals per draw (Maine c1_ens)
    for comp in ("C1", "P5"):
        el = [m for m in members if C.eligible(bundle["engines"][m]["engine"], comp) and not (comp == "P5" and bundle["engines"][m]["engine"] == "landis")]
        res = {m: comps_single(bundle, m, opt, overrides.get(m)).get(comp) for m in el}
        have = [m for m, r in res.items() if r is not None and np.isfinite(r.get("S_sw_pt", np.nan)) and (bundle["engines"][m]["engine"], comp) not in ins]
        if not have or not C.half_rule(have, el):
            memb.append(dict(component=comp, members=",".join(have), eligible=len(el), evaluated=False)); continue
        out[comp] = _ens_totals(bundle, have, comp, opt, overrides, how)
        memb.append(dict(component=comp, members=",".join(have), eligible=len(el), evaluated=True))
    # C3 per chain means
    el = [m for m in members if C.eligible(bundle["engines"][m]["engine"], "P2")]
    have = [m for m in el if (dict(bundle["engines"][m], **overrides.get(m, {}))).get("C3") and (bundle["engines"][m]["engine"], "P2") not in ins]
    if have and C.half_rule(have, el):
        mb = {}
        for h in bundle["chains"]:
            vs = [dict(bundle["engines"][m], **overrides.get(m, {}))["C3"].get(h, np.full(len(bundle["chains"][h]["labels"]), np.nan)) for m in have]
            vs = [v if C.is_xf(v) else np.asarray(v, float) for v in vs]
            mb[h], _ = C.ens_mean_xf(vs, how)
        out["C3"] = _c3_res(bundle, {"C3": mb}, opt)
    memb.append(dict(component="C3", members=",".join(have), eligible=len(el), evaluated=bool(have and C.half_rule(have, el))))
    # P3 / P3s: averaged member lines computed upstream (v4 item 6), passed in bundle['ens_p3'][(ensemble, how)]
    return out, memb


def _ens_totals(bundle, have, comp, opt, overrides, how):
    f = np.mean if how == "mean" else np.median
    if comp == "C1":
        T = bundle["c1_target"]; qs = list(C.QS); tots = {}   # stock is not averaged into ensemble C1 (by construction for most members)
        for q in qs:
            pts_, bts = [], []
            for m in have:
                E = dict(bundle["engines"][m]); E.update(overrides.get(m, {})); c = E["C1"]; eng = E["engine"]
                if c["kind"] == "plot":
                    d, M = bundle["c1_design"]; p, b = P_total(d, c["y"][q], M)
                else:
                    if q not in c["tot"]:
                        continue
                    p = float(c["tot"][q]); b = np.full(bundle["B"], p)
                if eng in C.FOLIAGE:
                    p, b = p * bundle["phi"][0], b * np.asarray(bundle["phi"][1])
                pts_.append(p); bts.append(b)
            if pts_:
                tots[q] = (float(f(pts_)), f(np.stack(bts, 0), 0))
        return P.c1_component(tots, {q: T[q] for q in tots}, list(tots))
    yrs = bundle["p5_years"]; ser = {}
    for y in yrs:
        pts_, bts = [], []
        for m in have:
            E = dict(bundle["engines"][m]); E.update(overrides.get(m, {})); c = E["P5"]
            if c["kind"] == "plot":
                d, M = bundle["p5_design"]
                if y not in c["series"]:
                    continue
                p, b = P_total(d, c["series"][y], M)
            elif c["kind"] == "total_draws":
                if y not in c["series"]:
                    continue
                p, b = c["series"][y]
            else:
                if y not in c["series"]:
                    continue
                p = float(c["series"][y]); b = np.full(bundle["B"], p)
            pts_.append(p); bts.append(b)
        if len(pts_) == len(have):
            ser[y] = (float(f(pts_)), f(np.stack(bts, 0), 0))
    tgt = bundle["p5_target"]["A1" if opt["p5_target"] == "A1" else "raw"]
    return P.p5_component(ser, tgt, yrs)


# ---------------------------------------------------------------------------------------------------- state
def score_state(bundle, opt=None):
    o = dict(R.DEFAULT_OPT); o.update(opt or {}); B = bundle["B"]; st = bundle["state"]
    rows, comp_rows, memb_rows, scal = [], [], [], []
    names = list(bundle["engines"])
    ins_by = lambda eng: {c for (e, c) in bundle["insample"] if e == eng}  # noqa: E731
    for n in names:
        E = bundle["engines"][n]
        cs = comps_single(bundle, n, o)
        if o["exclude_acd_c2p2"] and E["engine"] == "acd141":
            cs.pop("C2", None); cs.pop("C3", None)
        insmp = ins_by(E["engine"]) if E.get("role") != "ncasi_arm" else {"C1", "P3", "C5", "P5"}
        if E.get("role") == "ncasi_arm":       # v2 item 5: C2 and P2 only
            cs = {k: v for k, v in cs.items() if k in ("C2", "C3")}
        rows.append(R.score_row(n, cs, insmp, o, B, E.get("note", "")) | dict(role=E.get("role", "member"), engine=E["engine"]))
        comp_rows += _comp_rows(n, cs)
    overrides, rec = {}, {}
    for n in names:
        E = bundle["engines"][n]
        if E.get("role") == "member" and E["engine"] in C.RECAL_TARGET and E.get("C2") is not None:
            ov, rep = recal_override(bundle, n); overrides[n] = ov; scal.append(rep)
            cs = comps_single(bundle, n, o, ov)
            rows.append(R.score_row(n + "_recal", cs, ins_by(E["engine"]), o, B, "cross-fitted recalibration; scalar does not act on realism")
                        | dict(role="recal", engine=E["engine"]))
            comp_rows += _comp_rows(n + "_recal", cs)
    members = [m for m in bundle["members"] if m not in o["drop_members"]]
    plot_members = [m for m in members if bundle["engines"][m]["engine"] in C.PLOT_BASED]
    ens_defs = [("E_state", members, {}, "mean"), ("E_recal", members, overrides, "mean"),
                ("E_plot", plot_members, {}, "mean"), ("E_state_median", members, {}, "median")]
    for en, mem, ovr, how in ens_defs:
        cs, mb = ensemble_comps(bundle, mem, o, ovr, how)
        for k in ("P3", "P3s"):
            v = bundle.get("ens_p3", {}).get((en, k))
            if v is not None:
                cs[k] = v
        rows.append(R.score_row(en, cs, set(), o, B, "members " + ",".join(mem)) | dict(role="ensemble", engine="ensemble"))
        comp_rows += _comp_rows(en, cs); memb_rows += [dict(ensemble=en, **x) for x in mb]
    comparison = {r["row"] for r in rows if r["role"] == "member" and r["row"] in bundle["members"]} | {"E_state"}
    rows = R.rank_rows(rows, comparison)
    for r in rows:
        r["state"] = st
    return rows, pd.DataFrame(comp_rows).assign(state=st), pd.DataFrame(memb_rows).assign(state=st), pd.DataFrame(scal).assign(state=st)


def _comp_rows(name, cs):
    out = []
    for k, r in cs.items():
        if not isinstance(r, dict) or "S_sw_pt" not in r:
            if isinstance(r, dict) and "realism" in r:
                out.append(dict(row=name, component=k, realism=r.get("realism"), coverage=r.get("coverage")))
            continue
        d = dict(row=name, component=k, S_statewide=r.get("S_sw_pt"), S_by_domain=r.get("S_dom_pt"), n=r.get("n"), flag=r.get("flag", ""))
        for q, e in (r.get("e_pt") or {}).items():
            d[f"ratio_{q}_statewide"] = float(np.exp(np.asarray(e)[-1])) if np.isfinite(np.asarray(e)[-1]) else np.nan
        if r.get("hard_exit"):
            d["hard_exit"] = ";".join(f"{h['quantity']}:{h['ratio']:.3f}" for h in r["hard_exit"])
        out.append(d)
    return out
