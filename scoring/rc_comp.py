"""Track R components with shared draws: C1 statewide totals, C2 and C5 pair means, C3 hindcast chains, P5 series index.
Each function takes unit-level arrays already aligned (model, observed, location multiplicities M (B, n), domain labels)
and returns dict(e_pt={q: (K+1,)}, e_bt={q: (B, K+1)}, S_sw_pt, S_sw_bt, S_dom_pt, S_dom_bt, n=units, n_by_domain).
The last column of every e array is statewide. Aggregates only.
"""
import numpy as np
import rc_core as C


def _wmeans(v, M, dmask):
    """weighted means of v over units in each domain column of dmask (n, K+1): point and (B, K+1).
    v may be an rc_core.XF (D-F2): draw b then uses its own refit cross-fit scalars."""
    xf = C.is_xf(v)
    vp = np.asarray(v, float); ok = np.isfinite(vp); vz = np.where(ok, vp, 0.0)
    D = dmask * ok[:, None]
    n = D.sum(0)
    pt = np.where(n > 0, (vz[:, None] * D).sum(0) / np.where(n > 0, n, 1), np.nan)
    num = v.M_dot(M, D) if xf else M.astype(float) @ (vz[:, None] * D); den = M.astype(float) @ D
    bt = np.where(den > 0, num / np.where(den > 0, den, 1), np.nan)
    return pt, bt, n


def dmask_of(domains, labels):
    lab = np.asarray(labels).astype(str)
    return np.column_stack([lab == d for d in domains] + [np.ones(len(lab), bool)]).astype(float)


def pair_component(model, obs, M, labels, domains, aw, quantities=C.QS, min_units=C.MIN_UNITS, foliage=None):
    """C2 / C5: e_q = ln(mean model / mean obs) per domain and statewide on units where both exist.
    model, obs: dict q -> (n,) arrays; foliage: None or (pt, (B,)) factor (1 - phi) applied to the model (FOL engines)."""
    Dm = dmask_of(domains, labels); K = len(domains); e_pt, e_bt, nn = {}, {}, {}
    for q in quantities:
        mx = model[q]; ov = np.asarray(obs[q], float)
        cov = np.isfinite(np.asarray(mx, float)) & np.isfinite(ov)
        mv = mx.masked(cov) if C.is_xf(mx) else np.where(cov, np.asarray(mx, float), np.nan); ov = np.where(cov, ov, np.nan)
        m_pt, m_bt, n = _wmeans(mv, M, Dm); o_pt, o_bt, _ = _wmeans(ov, M, Dm)
        if foliage is not None:
            m_pt = m_pt * foliage[0]; m_bt = m_bt * np.asarray(foliage[1])[:, None]
        ep, eb = C.lnr(m_pt, o_pt), C.lnr(m_bt, o_bt)
        small = n[:K] < min_units; ep[:K][small] = np.nan; eb[:, :K][:, small] = np.nan
        e_pt[q], e_bt[q], nn[q] = ep, eb, n
    return _assemble(e_pt, e_bt, aw, K, nn)


def _assemble(e_pt, e_bt, aw, K, nn):
    qs = list(e_pt)
    S_sw_pt = float(C.S_quant([e_pt[q][K] for q in qs])); S_sw_bt = C.S_quant([e_bt[q][:, K] for q in qs])
    S_dom_pt = float(C.S_domains([e_pt[q][:K] for q in qs], np.asarray(aw))) if K else np.nan
    S_dom_bt = C.S_domains([e_bt[q][:, :K] for q in qs], np.asarray(aw)) if K else np.full(len(S_sw_bt), np.nan)
    n0 = nn[qs[0]]
    return dict(e_pt=e_pt, e_bt=e_bt, S_sw_pt=S_sw_pt, S_sw_bt=S_sw_bt, S_dom_pt=S_dom_pt, S_dom_bt=S_dom_bt,
                n=int(n0[K]), n_by_domain=[int(x) for x in n0[:K]], quantities=qs)


def c1_component(model_tot, target_tot, quantities):
    """C1 statewide: model_tot / target_tot: dict q -> (pt float, bt (B,)). Statewide only (common basis, item 4).
    Hard exit (item 10): a quantity whose point ratio is outside [0.2, 5] is not scored and is returned in 'hard_exit'."""
    e_pt, e_bt, hx = {}, {}, []
    for q in quantities:
        mp, mb = model_tot[q]; tp, tb = target_tot[q]
        r = mp / tp if tp else np.nan
        if not C.hard_exit_ok(r):
            hx.append(dict(quantity=q, ratio=float(r))); continue
        e_pt[q] = np.array([C.lnr(mp, tp)]); e_bt[q] = C.lnr(np.asarray(mb, float), np.asarray(tb, float))[:, None]
    if not e_pt:
        return dict(S_sw_pt=np.nan, S_sw_bt=None, e_pt={}, e_bt={}, hard_exit=hx, quantities=[])
    qs = list(e_pt)
    return dict(S_sw_pt=float(C.S_quant([e_pt[q][0] for q in qs])), S_sw_bt=C.S_quant([e_bt[q][:, 0] for q in qs]),
                e_pt=e_pt, e_bt=e_bt, hard_exit=hx, quantities=qs)


def c3_component(model_by_h, obs_by_h, M_by_h, labels_by_h, domains, aw, min_chains=100, min_units=C.MIN_UNITS):
    """P2 hindcast: per horizon h, e_h = ln(mean projected stock / mean observed stock) on chains with both; a horizon with fewer
    than min_chains chains is not scored (item 9; S9 uses 30). S = mean over scored horizons of |e_h| (statewide and by domain)."""
    out = {}
    for h in model_by_h:
        mv = model_by_h[h]; ov = np.asarray(obs_by_h[h], float)
        if (np.isfinite(np.asarray(mv, float)) & np.isfinite(ov)).sum() < min_chains:
            continue
        out[h] = pair_component({"stock": mv}, {"stock": ov}, M_by_h[h], labels_by_h[h], domains, aw, ("stock",), min_units)
    if not out:
        return None
    hs = sorted(out)
    return dict(horizons=hs, S_sw_pt=float(np.mean([out[h]["S_sw_pt"] for h in hs])),
                S_sw_bt=np.mean(np.stack([out[h]["S_sw_bt"] for h in hs], 0), 0),
                S_dom_pt=float(np.nanmean([out[h]["S_dom_pt"] for h in hs])),
                S_dom_bt=np.nanmean(np.stack([out[h]["S_dom_bt"] for h in hs], 0), 0),
                n={h: out[h]["n"] for h in hs}, per_h=out)


def p5_component(model_series, target_series, years):
    """P5 statewide: model_series, target_series: dict year -> (pt, bt (B,)); index = value / value at the first year;
    S = mean over later years of |ln(model index / target index)| (Maine C4 'S mean |ln(index/FIA trend index)|').
    target_series is the removals-adjusted FIA stock (amendment A1) or the raw series (A1s)."""
    y0 = years[0]; ys = [y for y in years[1:] if y in model_series and y in target_series]
    if not ys:
        return None
    mi = lambda y: (model_series[y][0] / model_series[y0][0], np.asarray(model_series[y][1]) / np.asarray(model_series[y0][1]))  # noqa: E731
    ti = lambda y: (target_series[y][0] / target_series[y0][0], np.asarray(target_series[y][1]) / np.asarray(target_series[y0][1]))  # noqa: E731
    E_pt = [abs(float(C.lnr(mi(y)[0], ti(y)[0]))) for y in ys]
    E_bt = np.stack([np.abs(C.lnr(mi(y)[1], ti(y)[1])) for y in ys], 0)
    return dict(S_sw_pt=float(np.mean(E_pt)), S_sw_bt=np.mean(E_bt, 0), years=ys)


def removals_adjusted(stock_series, removal_rate_series, years):
    """A1: target(t) = FIA stock(t) + cumulative statewide removals since the first year; removal_rate_series: dict year ->
    annual removals (pt, bt) from the GRM evaluation ending that year; years without an evaluation take the linear
    interpolation of their neighbours (constant at the ends). Returns dict year -> (pt, bt)."""
    ry = sorted(removal_rate_series)
    def rate(y):
        if y in removal_rate_series:
            return removal_rate_series[y]
        lo = max([r for r in ry if r < y], default=None); hi = min([r for r in ry if r > y], default=None)
        if lo is None:
            return removal_rate_series[hi]
        if hi is None:
            return removal_rate_series[lo]
        w = (y - lo) / (hi - lo)
        return ((1 - w) * removal_rate_series[lo][0] + w * removal_rate_series[hi][0],
                (1 - w) * np.asarray(removal_rate_series[lo][1]) + w * np.asarray(removal_rate_series[hi][1]))
    y0 = years[0]; out = {}; cp = 0.0; cb = 0.0
    for y in range(y0, max(years) + 1):
        if y > y0:
            r = rate(y); cp = cp + r[0]; cb = cb + np.asarray(r[1])
        if y in stock_series:
            out[y] = (stock_series[y][0] + cp, np.asarray(stock_series[y][1]) + cb)
    return out
