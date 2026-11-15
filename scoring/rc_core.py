"""Track R scorer core (NCASI report card, frozen pre-registration v5 + amendment A).

Pure arithmetic, no I/O: pillar points, P1 combination, MPI with draws, coverage, tier, grade, noise-corrected typical error,
P(rank 1), half rule, eligibility, hard exit, stress range, cross-fit scalars, ensemble means on common units.
Formulas follow tools/composite/composite.py (MPI v2.2 rules) unless a frozen Track R item says otherwise (cited inline).
"""
import numpy as np

W = dict(P1=30.0, P2=20.0, P3=20.0, P4=15.0, P5=15.0)           # item 2
PILLARS = tuple(W)
TIER_COV = 0.65; P3_MIN_COV = 0.5; ERR_ZERO = 0.50; FLOOR = 0.01; MIN_UNITS = 30
HARD_LO, HARD_HI = 0.2, 5.0                                      # item 10, v4 item 4
QS = ("gross_growth", "mortality", "net_growth")

# v3 item 1 eligibility table (component -> engines); 'P3' = full Bakuzis
ELIG = {
    "C1": {"fvs", "osm", "acd141", "cem", "landis", "cbm_gated", "cbm_af", "fia_curves"},
    "C2": {"fvs", "osm", "acd141", "cem", "landis", "cbm_gated", "cbm_af", "fia_curves"},
    "P2": {"fvs", "osm", "acd141", "cem", "landis", "cbm_gated", "cbm_af", "fia_curves"},
    "P3": {"fvs", "osm", "acd141", "cem"},
    "P3s": {"fvs", "osm", "acd141", "cem", "landis", "cbm_gated", "cbm_af", "fia_curves"},
    "C5": {"fvs", "osm", "acd141", "cem", "landis", "cbm_gated", "cbm_af", "fia_curves"},
    "P5": {"fvs", "osm", "acd141", "cem", "landis", "cbm_gated", "cbm_af"},   # curves no P5 (v4 item 7b); LANDIS n.e. by A2 at run time
}
C2_QUANT = {"fia_curves": ("net_growth",)}                      # v3 item 1: curves C2 net only
CURVE_DRIVEN = {"fia_curves", "cbm_gated", "cbm_af"}            # item 3 / v3 item 8
PLOT_BASED = {"fvs", "osm", "acd141", "cem"}                    # item 7 E_plot
AGGREGATE = {"landis", "cbm_gated", "cbm_af", "fia_curves"}
FOLIAGE = {"cbm_gated", "cbm_af"}     # model AG includes foliage: x (1 - phi), Maine FOL rule (LANDIS outputs are x (1 - phi) at extraction)
RECAL_TARGET = {"fvs": "gross_growth", "osm": "mortality", "acd141": "mortality"}   # v2 item 1

# item 7 member lineup per state (headline E_state)
MEMBERS = {
    "ME": ["fvs", "osm", "acd141", "cem", "landis", "cbm_gated", "cbm_af", "fia_curves"],
    "NH": ["fvs", "cem", "landis", "cbm_gated", "fia_curves"],
    "VT": ["fvs", "cem", "landis", "cbm_gated", "fia_curves"],
    "MI": ["fvs", "cem", "landis", "cbm_gated", "cbm_af", "fia_curves"],
    "WI": ["fvs", "cem", "landis", "cbm_gated", "cbm_af", "fia_curves"],
    "MN": ["fvs", "cem", "landis", "cbm_gated", "cbm_af", "fia_curves"],
}


def eligible(engine, comp):
    return engine in ELIG[comp]


def pts(S, zero=ERR_ZERO):
    """mean |ln ratio| -> points (composite.pts; zero = 0.5 registered, S18 uses 0.25 and 1.0)"""
    return 100.0 * np.maximum(0.0, 1.0 - (np.exp(np.asarray(S, float)) - 1.0) / zero)


def grade(m):
    return "" if not np.isfinite(m) else "A" if m >= 85 else "B" if m >= 70 else "C" if m >= 55 else "D" if m >= 40 else "F"


def lnr(m, f):
    """ln(model / FIA) with the 1% floor (score.py); nan where FIA <= 0"""
    m = np.asarray(m, float); f = np.asarray(f, float)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(f > 0, np.log(np.maximum(m, FLOOR * f) / f), np.nan)


def typical_error(S):
    return np.exp(np.asarray(S, float)) - 1.0


def te_noise_corrected(S_pt, S_bt):
    """v2 item 14: sqrt(max(0, e^2 - SE^2)), e = typical error at the point, SE = bootstrap SD of e"""
    e = float(typical_error(S_pt)); bt = typical_error(np.asarray(S_bt, float)); bt = bt[np.isfinite(bt)]
    if not np.isfinite(e) or len(bt) < 2:
        return np.nan
    return float(np.sqrt(max(0.0, e * e - np.var(bt, ddof=1))))


def p1_S(sc1, sc2):
    """P1 = (0.35 S_C1 + 0.25 S_C2) / 0.60 when C2 is out of sample, else S_C1 (item 2); None-safe. Arrays allowed."""
    if sc1 is None and sc2 is None:
        return None
    if sc2 is None:
        return sc1
    if sc1 is None:
        return None        # C1 is required for P1 (registered); C2 alone does not make P1
    return (0.35 * np.asarray(sc1) + 0.25 * np.asarray(sc2)) / 0.60


def S_quant(e_list):
    """mean |e| over quantities (statewide basis)"""
    return np.nanmean(np.abs(np.stack([np.asarray(e, float) for e in e_list], 0)), 0)


def S_domains(e_by_q, aw):
    """by-domain view: area-weighted mean over scored domains of the mean |e| over quantities; e_by_q: list of (..., K)"""
    A = np.abs(np.stack(e_by_q, 0)); m = np.nanmean(A, 0)                # (..., K)
    w = np.where(np.isfinite(m), aw, 0.0); s = w.sum(-1)
    with np.errstate(invalid="ignore"):
        return np.where(s > 0, np.nansum(np.where(np.isfinite(m), m, 0.0) * w, -1) / s, np.nan)


def mpi(P, B, weights=W):
    """P: pillar -> None or (points_pt, points_bt (B,)). Returns m_pt, m_bt, coverage, evaluated list (composite.score_all)"""
    ev = [p for p in PILLARS if P.get(p) is not None and np.isfinite(P[p][0])]
    ws = sum(weights[p] for p in ev)
    if not ws:
        return np.nan, np.full(B, np.nan), 0.0, ev
    m_pt = sum(weights[p] * P[p][0] for p in ev) / ws
    m_bt = sum(weights[p] * np.asarray(P[p][1], float) for p in ev) / ws
    return float(m_pt), m_bt, ws / sum(weights.values()), ev


def tier(cov):
    return 1 if cov >= TIER_COV - 1e-12 else 2


def p_rank1(BT):
    """BT: (rows, B) MPI draws of the comparison set; share of draws each row is first (nan-safe argmax)"""
    BT = np.asarray(BT, float)
    if BT.ndim != 2 or BT.shape[0] == 0:
        return np.array([])
    X = np.where(np.isfinite(BT), BT, -np.inf); win = np.argmax(X, 0)
    return np.array([(win == i).mean() for i in range(BT.shape[0])])


def half_rule(contributing, eligible_members):
    """v2 item 9 / v3 item 1: pillar evaluated for an ensemble only if at least half of its eligible declared members contribute"""
    n_el = len(eligible_members)
    return n_el > 0 and len(contributing) >= n_el / 2.0


def hard_exit_ok(ratio):
    """item 10: a C1 cell outside [0.2, 5] x EVALIDator is not scored (diagnosis only, v4 item 4)"""
    r = float(ratio)
    return np.isfinite(r) and HARD_LO <= r <= HARD_HI


def ens_mean(member_values, how="mean"):
    """per-unit ensemble prediction on common units: member_values list of arrays (same unit order), nan = not covered.
    Returns (pred, common mask). v2 item 9: scored only on units every contributing member covers."""
    V = np.stack([np.asarray(v, float) for v in member_values], 0)
    common = np.all(np.isfinite(V), 0)
    f = np.mean if how == "mean" else np.median
    pred = np.where(common, f(np.where(np.isfinite(V), V, 0.0), 0), np.nan)
    return pred, common


def crossfit_scalars(obs_num, mod_den, fold, M=None):
    """v2 items 1 to 3: per fold scalar = sum obs / sum model on that fold's C2 common pairs (ratio of sums), refit in every
    draw. obs_num, mod_den: (n,) targeted component; fold: (n,) in {1, 2}; M: (B, n) multiplicities or None.
    Returns s1, s2 (point floats if M is None, else (B,) arrays)."""
    o = np.asarray(obs_num, float); m = np.asarray(mod_den, float); f = np.asarray(fold)
    ok = np.isfinite(o) & np.isfinite(m)
    out = []
    for k in (1, 2):
        sel = ok & (f == k)
        if M is None:
            out.append(float(o[sel].sum() / m[sel].sum()) if m[sel].sum() > 0 else np.nan)
        else:
            Ms = M[:, sel].astype(float)
            den = Ms @ m[sel]
            out.append(np.where(den > 0, (Ms @ o[sel]) / np.where(den > 0, den, 1.0), np.nan))
    return out[0], out[1]


def apply_crossfit(g, mrt, fold, s1, s2, target):
    """every unit receives the scalar fitted on the OTHER fold (item 5); returns gross, mort, net (post-processing, v2 item 2).
    Scalars may be floats (point) or (B,) arrays -> outputs (B, n)."""
    g = np.asarray(g, float); mrt = np.asarray(mrt, float); fold = np.asarray(fold)
    if np.ndim(s1) == 0:
        s = np.where(fold == 1, s2, s1)
    else:
        s = np.where(fold[None, :] == 1, np.asarray(s2)[:, None], np.asarray(s1)[:, None])
    if target == "gross_growth":
        g2 = g * s; m2 = np.broadcast_to(mrt, g2.shape)
    else:
        m2 = mrt * s; g2 = np.broadcast_to(g, m2.shape)
    return g2, m2, g2 - m2


def stress_range(registered, alts, rank_of=None, best_rank_is_high=True):
    """v3 item 6: range over the applicable alternatives (dict name -> value or None), S11 schemes counted one each.
    Flag 'favourable extreme' = registered equals the end of the range that gives the row its best within-state rank
    (higher MPI -> better rank). Returns dict."""
    vals = {k: v for k, v in alts.items() if v is not None and np.isfinite(v)}
    allv = [registered] + list(vals.values())
    lo, hi = float(np.min(allv)), float(np.max(allv))
    fav = hi if best_rank_is_high else lo
    return dict(lo=lo, hi=hi, n_alternatives=len(vals), at_favourable_extreme=bool(np.isclose(registered, fav) and len(vals) > 0),
                applied=sorted(vals))


def leader_share(leaders_by_alt, registered_leader):
    """share of applicable alternatives that keep the headline leader (v2 item 18)"""
    L = [l for l in leaders_by_alt.values() if l is not None]
    return float(np.mean([l == registered_leader for l in L])) if L else np.nan


# ---------------------------------------------------------------------------------------------------- D-F2 per-draw cross-fit
class XF:
    """Linear cross-fit model value per unit: v = base + sum_j mult_j * s_other_j, where s_other is the scalar fitted on the
    OTHER fold (item 5), refit in every bootstrap draw (v2 item 3). terms: list of (mult, fold, s1pt, s1bt, s2pt, s2bt).
    Point value via pt() / np.asarray(); draw-wise weighted sums via M_dot(). D-F2 (2026-11-06)."""
    def __init__(self, base, terms):
        self.base = np.asarray(base, float); self.terms = [(np.asarray(m, float), np.asarray(f), float(s1p), np.asarray(s1b, float),
                                                            float(s2p), np.asarray(s2b, float)) for m, f, s1p, s1b, s2p, s2b in terms]

    def pt(self):
        v = self.base.copy()
        for m, f, s1p, _, s2p, _ in self.terms:
            v = v + m * np.where(f == 1, s2p, s1p)
        return v

    def __array__(self, dtype=None, copy=None):
        a = self.pt()
        return a.astype(dtype) if dtype is not None else a

    def __len__(self):
        return len(self.base)

    def subset(self, sel):
        return XF(self.base[sel], [(m[sel], f[sel], a, b, c, d) for m, f, a, b, c, d in self.terms])

    def masked(self, keep):
        return XF(np.where(keep, self.base, np.nan), self.terms)

    def scale(self, c):
        return XF(self.base * c, [(m * c, f, a, b, cc, d) for m, f, a, b, cc, d in self.terms])

    def M_dot(self, M, W):
        """sum_i M[b, i] * v_b[i] * W[i, k] with draw b's scalars; W (n, K) already zero where the unit is not covered"""
        M = M.astype(float); ok = np.isfinite(self.pt())
        W = W * ok[:, None]
        out = M @ (np.where(ok, self.base, 0.0)[:, None] * W)
        for m, f, s1p, s1b, s2p, s2b in self.terms:
            mz = np.where(ok, m, 0.0)
            f1 = (f == 1).astype(float)
            out = out + s2b[:, None] * (M @ ((mz * f1)[:, None] * W)) + s1b[:, None] * (M @ ((mz * (1 - f1))[:, None] * W))
        return out


def is_xf(v):
    return isinstance(v, XF)


def apply_crossfit_xf(g, mrt, fold, s1p, s1b, s2p, s2b, target):
    """D-F2: per-draw version of apply_crossfit; returns dict gross, mortality, net as XF or arrays"""
    g = np.asarray(g, float); mrt = np.asarray(mrt, float); fold = np.asarray(fold)
    if target == "gross_growth":
        gx = XF(np.zeros_like(g), [(g, fold, s1p, s1b, s2p, s2b)])
        return {"gross_growth": gx, "mortality": mrt, "net_growth": XF(-mrt, [(g, fold, s1p, s1b, s2p, s2b)])}
    mx = XF(np.zeros_like(mrt), [(mrt, fold, s1p, s1b, s2p, s2b)])
    return {"gross_growth": g, "mortality": mx, "net_growth": XF(g, [(-mrt, fold, s1p, s1b, s2p, s2b)])}


def ens_mean_xf(member_values, how="mean"):
    """ens_mean that keeps XF members linear (mean only); median uses point values (disclosed in D-F2)"""
    pts = [np.asarray(v, float) for v in member_values]
    if how != "mean" or not any(is_xf(v) for v in member_values):
        return ens_mean(pts, how)
    V = np.stack(pts, 0); common = np.all(np.isfinite(V), 0); k = len(member_values)
    base = np.zeros(V.shape[1]); terms = []
    for v in member_values:
        if is_xf(v):
            base = base + np.where(common, v.base, 0.0) / k
            terms += [(np.where(common, m, 0.0) / k, f, a, b, c, d) for m, f, a, b, c, d in v.terms]
        else:
            base = base + np.where(common, np.asarray(v, float), 0.0) / k
    base = np.where(common, base, np.nan)
    return XF(base, terms), common
