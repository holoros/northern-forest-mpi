"""Track R design-based targets and the shared location bootstrap (v5 item 1, v4 item 3).

Universe: every plot in the state's 2024 evaluation (EVALID ST2401, forest and nonforest), resampled with replacement within
POP_STRATUM, stratum weights fixed. One multiplicity per plot LOCATION per draw (STATECD_UNITCD_COUNTYCD_PLOT), shared by
every component, engine and ensemble. A unit at a location outside the 2024 universe gets an independent Poisson(1)
multiplicity per draw (counted and reported). Per-plot FIA values follow targets/evalidator_pathB.py (the 1% gate path).
Restricted: per-plot values never leave the caller; functions return arrays.
"""
import os
from functools import lru_cache
from pathlib import Path
import numpy as np
import pandas as pd

FIA_ROOT = Path(os.environ.get("RC_FIA_ROOT", os.path.expanduser("~/restricted/ncasi-modeleval/fia_data_rc")))
LB2TG = 0.45359237e-9


@lru_cache(maxsize=None)
def _rd(st, t, cols):
    cols = list(cols)
    return pd.read_csv(FIA_ROOT / st / f"{st}_{t}.csv", usecols=cols,
                       dtype={c: "string" for c in cols if c.endswith("CN")}, low_memory=False)


def statecd(st):
    return int(_rd(st, "PLOT", ("STATECD",)).STATECD.iloc[0])


def loc_of(st):
    pl = _rd(st, "PLOT", ("CN", "STATECD", "UNITCD", "COUNTYCD", "PLOT"))
    key = pl.STATECD.astype(str) + "_" + pl.UNITCD.astype(str) + "_" + pl.COUNTYCD.astype(str) + "_" + pl.PLOT.astype(str)
    return pd.Series(key.values, index=pl.CN.astype(str).values)


def design(st, ev):
    """plots of EVALID ev with stratum, c_h = EXPNS x P2POINTCNT, adjustment factors, REMPER, MACRO_BREAKPOINT_DIA"""
    ps = _rd(st, "POP_STRATUM", ("CN", "EVALID", "EXPNS", "P2POINTCNT", "ADJ_FACTOR_MICR", "ADJ_FACTOR_SUBP", "ADJ_FACTOR_MACR"))
    ps = ps[ps.EVALID == ev].rename(columns={"CN": "STRATUM_CN"})
    pa = _rd(st, "POP_PLOT_STRATUM_ASSGN", ("PLT_CN", "STRATUM_CN", "EVALID"))
    pa = pa[pa.EVALID == ev][["PLT_CN", "STRATUM_CN"]]
    pl = _rd(st, "PLOT", ("CN", "MACRO_BREAKPOINT_DIA", "REMPER")).rename(columns={"CN": "PLT_CN"})
    d = pa.merge(ps, on="STRATUM_CN").merge(pl, on="PLT_CN", how="left")
    d["c_h"] = d.EXPNS * d.P2POINTCNT
    d["loc"] = d.PLT_CN.astype(str).map(loc_of(st))
    return d.reset_index(drop=True)


def evalids(st, typ):
    """EVALIDs of type 'EXPCURR' (stock, ST..01) or 'EXPGROW' (GRM, ST..03), with their year"""
    t = _rd(st, "POP_EVAL_TYP", ("EVAL_TYP", "EVAL_CN"))
    e = _rd(st, "POP_EVAL", ("CN", "EVALID", "END_INVYR")).rename(columns={"CN": "EVAL_CN"})
    x = e.merge(t[t.EVAL_TYP == typ], on="EVAL_CN")
    return sorted((int(r.EVALID), int(r.END_INVYR)) for r in x.itertuples())


def plot_stock(st, ev):
    """per-plot live AG dry biomass (lb, adjusted), forest conditions, DIA >= 1 (path B)"""
    d = design(st, ev)
    cd = _rd(st, "COND", ("PLT_CN", "CONDID", "COND_STATUS_CD")); forest = cd[cd.COND_STATUS_CD == 1][["PLT_CN", "CONDID"]]
    tr = _rd(st, "TREE", ("CN", "PLT_CN", "CONDID", "STATUSCD", "DIA", "TPA_UNADJ", "DRYBIO_AG"))
    t = tr[(tr.STATUSCD == 1) & (tr.DIA >= 1) & tr.TPA_UNADJ.notna() & tr.DRYBIO_AG.notna()].merge(forest, on=["PLT_CN", "CONDID"])
    t = t.merge(d[["PLT_CN", "ADJ_FACTOR_MICR", "ADJ_FACTOR_SUBP", "ADJ_FACTOR_MACR", "MACRO_BREAKPOINT_DIA"]], on="PLT_CN")
    big = t.MACRO_BREAKPOINT_DIA.notna() & (t.DIA >= t.MACRO_BREAKPOINT_DIA)
    adj = np.where(t.DIA < 5, t.ADJ_FACTOR_MICR, np.where(big, t.ADJ_FACTOR_MACR, t.ADJ_FACTOR_SUBP))
    y = (t.TPA_UNADJ * adj * t.DRYBIO_AG).groupby(t.PLT_CN.astype(str)).sum()
    return d.PLT_CN.astype(str).map(y).fillna(0.0).values, d


def plot_grm(st, ev):
    """per-plot annual gross, mortality, net, removals (lb yr-1, adjusted) of AG biomass on forest land (path B + removals)"""
    d = design(st, ev)
    g = _rd(st, "TREE_GRM_COMPONENT", ("TRE_CN", "PLT_CN", "MICR_COMPONENT_AL_FOREST", "MICR_SUBPTYP_GRM_AL_FOREST",
                                         "MICR_TPAGROW_UNADJ_AL_FOREST", "MICR_TPAMORT_UNADJ_AL_FOREST", "MICR_TPAREMV_UNADJ_AL_FOREST")).copy()
    g = g.rename(columns={"MICR_COMPONENT_AL_FOREST": "comp", "MICR_SUBPTYP_GRM_AL_FOREST": "sp", "MICR_TPAGROW_UNADJ_AL_FOREST": "tg",
                          "MICR_TPAMORT_UNADJ_AL_FOREST": "tm", "MICR_TPAREMV_UNADJ_AL_FOREST": "tr"})   # by name: usecols keeps file order
    be = _rd(st, "TREE_GRM_BEGIN", ("TRE_CN", "DRYBIO_AG")).rename(columns={"DRYBIO_AG": "BE"})
    mi = _rd(st, "TREE_GRM_MIDPT", ("TRE_CN", "DRYBIO_AG")).rename(columns={"DRYBIO_AG": "MI"})
    en = _rd(st, "TREE", ("CN", "DRYBIO_AG")).rename(columns={"CN": "TRE_CN", "DRYBIO_AG": "EN"})
    g = g.merge(be, on="TRE_CN", how="left").merge(mi, on="TRE_CN", how="left").merge(en, on="TRE_CN", how="left")
    g = g.merge(d[["PLT_CN", "ADJ_FACTOR_MICR", "ADJ_FACTOR_SUBP", "ADJ_FACTOR_MACR", "REMPER"]], on="PLT_CN")
    g = g[g.REMPER.notna() & (g.REMPER > 0)].fillna({"tg": 0, "tm": 0, "tr": 0, "BE": 0, "MI": 0, "EN": 0})
    adj = np.select([g.sp == 1, g.sp == 2, g.sp == 3], [g.ADJ_FACTOR_SUBP, g.ADJ_FACTOR_MICR, g.ADJ_FACTOR_MACR], 0.0)
    c = g.comp.fillna("").astype(str)
    live_end = c.isin(["SURVIVOR", "INGROWTH"]) | c.str.startswith("REVERSION")
    gone = c.str.startswith(("CUT", "DIVERSION", "MORTALITY"))
    endg = np.where(live_end, g.EN, np.where(gone, g.MI, 0.0))
    begg = np.where(c.isin(["SURVIVOR", "CUT1", "DIVERSION1", "MORTALITY1"]), g.BE, 0.0)
    yg = g.tg * adj * (endg - begg) / g.REMPER
    ym = np.where(c.str.startswith("MORTALITY"), g.tm * adj * g.MI, 0.0)
    yr = np.where(c.str.startswith("CUT"), g.tr * adj * g.MI, 0.0)
    k = g.PLT_CN.astype(str).values
    agg = lambda v: d.PLT_CN.astype(str).map(pd.Series(np.asarray(v, float)).groupby(k).sum()).fillna(0.0).values  # noqa: E731
    gross, mort, rem = agg(yg), agg(ym), agg(yr)
    return dict(gross_growth=gross, mortality=mort, net_growth=gross - mort, removals=rem), d


class LocDraws:
    """shared within-stratum bootstrap of the 2024 evaluation, indexed by location"""

    def __init__(self, st, B, seed, ev=None):
        self.st = st; self.B = B; self.ev = ev or statecd(st) * 10000 + 2401
        d = design(st, self.ev); self.d = d
        rng = np.random.default_rng(np.random.SeedSequence(seed))
        M = np.zeros((B, len(d)), np.int32)
        for _, ix in d.groupby("STRATUM_CN").indices.items():
            ix = np.asarray(ix)
            dr = rng.integers(0, len(ix), (B, len(ix)))
            for b in range(B):
                M[b, ix] += np.bincount(dr[b], minlength=len(ix)).astype(np.int32)
        self.Mloc = pd.DataFrame(M, columns=d["loc"].values).T.groupby(level=0).sum().T  # (B, nloc), duplicates summed
        self._rng_out = np.random.default_rng(np.random.SeedSequence([seed, 7]))
        self._out = {}
        self.n_outside = 0

    def mult(self, locs):
        """(B, n) multiplicities for units at the given locations"""
        locs = pd.Series(locs).astype(str).values
        cols = []
        for l in locs:
            if l in self.Mloc.columns:
                cols.append(self.Mloc[l].values)
            else:
                if l not in self._out:
                    self._out[l] = self._rng_out.poisson(1.0, self.B).astype(np.int32); self.n_outside += 1
                cols.append(self._out[l])
        return np.stack(cols, 1) if cols else np.zeros((self.B, 0), np.int32)


def design_total(d, y, M=None):
    """stratified total: sum_h c_h * mean_h(y) (EVALIDator); with M (B, n) multiplicities on d's rows -> (B,)"""
    import rc_core as _C
    yx = y if _C.is_xf(y) else None                     # D-F2: cross-fit XF, draw b uses its own refit scalars
    y = np.asarray(y, float); out_pt = 0.0; out_bt = 0.0 if M is not None else None
    for _, ix in d.groupby("STRATUM_CN").indices.items():
        ix = np.asarray(ix); c = float(d.c_h.iloc[ix[0]]); n = float(d.P2POINTCNT.iloc[ix[0]])
        out_pt += c * y[ix].sum() / n
        if M is not None:
            Mi = M[:, ix].astype(float); den = Mi.sum(1)
            num = yx.subset(ix).M_dot(Mi, np.ones((len(ix), 1)))[:, 0] if yx is not None else Mi @ y[ix]
            out_bt = out_bt + c * np.where(den > 0, num / np.where(den > 0, den, 1.0), 0.0) * (len(ix) / n)
    return out_pt, out_bt
