"""Unit tests on synthetic data for every Track R rule (no FIA or engine data). Run: python -m pytest -q test_rc.py"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
sys.path.insert(0, str(Path(__file__).parent))
import rc_core as C
import rc_comp as P
import rc_run as RR
import rc_design as D

RNG = np.random.default_rng(1)


def test_pts_grade_match_composite():
    root = Path(__file__).resolve().parents[3]
    sys.path.insert(0, str(root / "tools/composite"))
    try:
        import composite as CO
    except Exception:
        return
    S = np.linspace(0, 1, 101)
    assert np.allclose(C.pts(S), CO.pts(S))
    for m in (np.nan, 30, 40, 54.99, 55, 70, 85, 99):
        assert C.grade(m) == CO.grade(m)


def test_p1_and_mpi():
    assert C.p1_S(0.1, None) == 0.1
    assert np.isclose(C.p1_S(0.1, 0.2), (0.035 + 0.05) / 0.6)
    assert C.p1_S(None, 0.2) is None
    B = 4
    m, mb, cov, ev = C.mpi({"P1": (80, np.full(B, 80.)), "P2": None, "P3": (60, np.full(B, 60.)), "P4": (50, np.full(B, 50.)), "P5": None}, B)
    assert np.isclose(m, (30 * 80 + 20 * 60 + 15 * 50) / 65) and np.isclose(cov, 65 / 100) and C.tier(cov) == 1
    assert C.tier(0.6499) == 2


def test_noise_corrected_te():
    S_bt = np.log(1 + RNG.normal(0.10, 0.05, 1000))
    v = C.te_noise_corrected(np.log(1.10), S_bt)
    assert 0.08 < v < 0.10
    assert C.te_noise_corrected(np.log(1.01), np.log(1 + RNG.normal(0.01, 0.05, 1000))) == 0.0


def test_crossfit_other_fold_and_no_leak():
    n = 400; fold = np.repeat([1, 2], n // 2)
    m = RNG.uniform(1, 2, n); o = m * np.where(fold == 1, 0.8, 1.2)
    s1, s2 = C.crossfit_scalars(o, m, fold)
    assert np.isclose(s1, 0.8) and np.isclose(s2, 1.2)
    g, mo, net = C.apply_crossfit(m, m * 0.3, fold, s1, s2, "gross_growth")
    assert np.allclose(g[fold == 1], m[fold == 1] * 1.2) and np.allclose(g[fold == 2], m[fold == 2] * 0.8)
    # changing fold 1 observations must not change fold 1 units' recalibrated values (no leak)
    o2 = o.copy(); o2[fold == 1] *= 3
    t1, t2 = C.crossfit_scalars(o2, m, fold)
    g2, _, _ = C.apply_crossfit(m, m * 0.3, fold, t1, t2, "gross_growth")
    assert np.allclose(g2[fold == 1], g[fold == 1])
    # refit inside draws
    M = RNG.poisson(1, (50, n))
    b1, b2 = C.crossfit_scalars(o, m, fold, M)
    assert b1.shape == (50,) and np.allclose(b1, 0.8) and np.allclose(b2, 1.2)
    gb, mb, nb = C.apply_crossfit(m, m * 0.3, fold, b1, b2, "mortality")
    assert gb.shape == (50, n) and np.allclose(mb[:, fold == 1], (m * 0.3 * 1.2)[fold == 1])


def test_eligibility_half_rule_hard_exit():
    assert not C.eligible("fia_curves", "P5") and not C.eligible("landis", "P3") and C.eligible("landis", "P3s")
    assert C.C2_QUANT["fia_curves"] == ("net_growth",)
    assert C.half_rule([1, 2], [1, 2, 3, 4]) and not C.half_rule([1], [1, 2, 3])
    assert C.hard_exit_ok(0.2) and C.hard_exit_ok(5) and not C.hard_exit_ok(0.19) and not C.hard_exit_ok(np.nan)


def test_ens_common_units_and_rank():
    a = np.array([1, 2, np.nan, 4.]); b = np.array([3, np.nan, 5, 6.])
    p, c = C.ens_mean([a, b])
    assert list(c) == [True, False, False, True] and np.allclose(p[c], [2, 5])
    pm, _ = C.ens_mean([a, b, a * 0 + 10], "median"); assert np.isclose(pm[0], 3)
    pr = C.p_rank1(np.array([[1, 2, 3], [3, 2, 1], [0, 0, np.nan]]))
    assert np.isclose(pr.sum(), 1) and pr[2] == 0


def test_stress_range_flag():
    r = C.stress_range(80, {"S2": 75, "S3": None, "S11a": 79})
    assert r["lo"] == 75 and r["hi"] == 80 and r["at_favourable_extreme"] and r["n_alternatives"] == 2
    assert C.leader_share({"a": "E_state", "b": "fvs", "c": None}, "E_state") == 0.5


def test_pair_component_point_and_min_units():
    n = 200; lab = np.array(["800"] * 150 + ["100"] * 50); M = RNG.poisson(1, (20, n))
    obs = {q: RNG.uniform(1, 2, n) for q in C.QS}; mod = {q: 1.1 * obs[q] for q in C.QS}
    r = P.pair_component(mod, obs, M, lab, ["800", "100"], [0.7, 0.3], C.QS, min_units=60)
    assert np.isclose(r["S_sw_pt"], np.log(1.1)) and np.isnan(r["e_pt"]["net_growth"][1])     # 50 < 60 units
    assert np.isclose(r["S_dom_pt"], np.log(1.1))


def test_removals_adjusted_and_p5():
    ss = {2006: (100.0, np.full(3, 100.)), 2010: (104.0, np.full(3, 104.))}
    rr = {2008: (1.0, np.ones(3)), 2010: (2.0, np.full(3, 2.))}
    out = P.removals_adjusted(ss, rr, [2006, 2010])
    # years 2007..2010: rates 1 (before first eval, constant), 1, 1.5 (interp), 2 -> 5.5
    assert np.isclose(out[2010][0], 104 + 5.5)
    p5 = P.p5_component({2006: (50, np.full(3, 50.)), 2010: (50 * 1.095, np.full(3, 50 * 1.095))}, out, [2006, 2010])
    assert np.isclose(p5["S_sw_pt"], abs(np.log(1.095 / 1.095)), atol=1e-12)


def test_design_total_identity():
    d = pd.DataFrame(dict(STRATUM_CN=["a", "a", "b"], c_h=[10.0, 10.0, 5.0], P2POINTCNT=[2, 2, 1]))
    y = np.array([1.0, 3.0, 2.0]); M = np.ones((3, 3), int)
    p, b = D.design_total(d, y, M)
    assert np.isclose(p, 10 * 4 / 2 + 5 * 2) and np.allclose(b, p)


def synth_bundle(shuffle=None, B=40, n=600, nt=200):
    rng = np.random.default_rng(7)
    idx = np.arange(n); tdx = np.arange(nt)
    if shuffle is not None:
        idx = np.random.default_rng(shuffle).permutation(n); tdx = np.random.default_rng(shuffle + 1).permutation(nt)
    lab = np.array(["800", "100"] * (n // 2))[idx]; tlab = np.array(["800", "100"] * (nt // 2))[tdx]
    fold = np.array([1, 2] * (n // 2))[np.argsort(np.argsort(idx))] if False else (np.arange(n) % 2 + 1)[idx]
    tfold = (np.arange(nt) % 2 + 1)[tdx]
    M = rng.poisson(1, (B, n))[:, idx]; Mt = rng.poisson(1, (B, nt))[:, tdx]
    obs = {q: rng.uniform(1, 2, n)[idx] for q in ("gross_growth", "mortality")}; obs["net_growth"] = obs["gross_growth"] - 0.5 * obs["mortality"]
    tobs = {q: rng.uniform(1, 2, nt)[tdx] for q in C.QS}
    d = pd.DataFrame(dict(STRATUM_CN=np.repeat(["a", "b"], 50), c_h=np.repeat([100.0, 80.0], 50), P2POINTCNT=50))
    Md = rng.poisson(1, (B, 100))
    yplot = rng.uniform(0, 3, 100)
    tgt = {q: D.design_total(d, yplot, Md) for q in C.QS}
    def eng(name, f, role="member", c5f=None, extra=None):
        e = dict(engine=name, arm="x", role=role,
                 C2={q: obs[q] * f for q in C.QS}, C5=c5f or {q: tobs[q] * f for q in C.QS},
                 C1=dict(kind="plot", y={q: yplot * f for q in C.QS}), P3=dict(realism=0.9, coverage=0.8))
        e.update(extra or {}); return e
    engines = {"fvs": eng("fvs", 1.2), "cem": eng("cem", 0.9, extra=dict(C5_flag="mostly untreated donors (A8)")),
               "landis": eng("landis", 1.05, extra=dict(A6b_pass=False, P3=None)),
               "fia_curves": eng("fia_curves", 0.8, c5f={"level": {q: tobs[q] * 0.5 for q in C.QS}, "ageq": {q: tobs[q] * 0.9 for q in C.QS}},
                                 extra=dict(P3=None)),
               "fvs_ncasi": eng("fvs", 1.1, role="ncasi_arm")}
    return dict(state="XX", B=B, domains=["800", "100"], aw=np.array([0.6, 0.4]), phi=(0.95, np.full(B, 0.95)),
                pairs=dict(labels=lab, M=M, fold=fold, obs=obs), treated=dict(labels=tlab, M=Mt, fold=tfold, obs=tobs, subset_mask=np.ones(nt, bool)),
                chains={}, c1_target=tgt, c1_design=(d, Md), c1_fold=np.arange(100) % 2 + 1, p5_years=[2006, 2010],
                p5_target={"A1": {}, "raw": {}}, engines=engines, insample=set(), members=["fvs", "cem", "landis", "fia_curves"])


def test_state_rules_and_row_shuffle():
    rows, comp, memb, scal = RR.score_state(synth_bundle())
    T = pd.DataFrame([{k: v for k, v in r.items() if k != "_bt"} for r in rows]).set_index("row")
    assert "P4 (C5) not evaluated" in T.loc["landis", "notes"]                    # A6b failed
    m = memb[(memb.ensemble == "E_state") & (memb.component == "C5")].iloc[0]
    assert m.members == "fvs" and not m.evaluated                               # curves (v3 item 8), CEM A8, LANDIS A6b excluded; half rule
    nc = T.loc["fvs_ncasi"]
    assert np.isnan(nc.P1) and np.isnan(nc.P4)                                  # NCASI arm: C2 and P2 only (v2 item 5)
    assert "fvs_recal" in T.index and np.isclose(scal.scalar_fold1.iloc[0], 1 / 1.2) and np.isclose(scal.scalar_fold2.iloc[0], 1 / 1.2)
    assert np.isnan(T.loc["fia_curves", "P5"])
    assert set(T.dropna(subset=["P_rank1"]).index) <= {"fvs", "cem", "landis", "fia_curves", "E_state"}
    rows2, _, _, _ = RR.score_state(synth_bundle(shuffle=3))
    T2 = pd.DataFrame([{k: v for k, v in r.items() if k != "_bt"} for r in rows2]).set_index("row")
    num = [c for c in T.columns if T[c].dtype.kind == "f" and c != "P_rank1"]
    assert np.allclose(T[num].values.astype(float), T2.loc[T.index, num].values.astype(float), equal_nan=True)


def test_p3s_score_rules(tmp_path):
    try:
        import rc_p3 as P3
    except Exception:
        return                                      # needs the repo (tools/v24_S, tools/ensemble); run on firebreather
    rows = []
    for t in ("pine", "spruce_fir"):
        for o in ("young", "middle", "mature"):
            for s in (1, 2, 3):
                for a in range(20, 160, 5):
                    st = 100 * s * (1 - np.exp(-a / 40.0))
                    rows.append(dict(type=t, origin=o, site=s, age=a, HT=np.nan, BAPH=np.nan, VOL=st, n_plots=10))
    f = tmp_path / "traj.csv"; pd.DataFrame(rows).to_csv(f, index=False)
    S, drop = P3.p3s_score({"a": ("landis", f), "b": ("cbm_gated", f)})
    assert set(drop) >= {"R1", "R2"}                # relations no row flags are dropped (non-discriminating)
    assert not S["a"]["evaluable"]                  # nothing left to discriminate -> not evaluable (coverage rule)
    assert set(S) == {"a", "b"}


# ---------------------------------------------------------------- D-F2 per-draw cross-fit (XF)
def test_xf_mdot_matches_explicit():
    import numpy as np, rc_core as C
    rng = np.random.default_rng(1); n, B = 40, 7
    g = rng.uniform(1, 3, n); mo = rng.uniform(0, 1, n); fold = rng.integers(1, 3, n)
    s1b = rng.uniform(0.8, 1.2, B); s2b = rng.uniform(0.8, 1.2, B)
    X = C.XF(-mo, [(g, fold, 1.0, s1b, 1.1, s2b)])
    M = rng.poisson(1, (B, n)).astype(float); W = np.ones((n, 1))
    got = X.M_dot(M, W)[:, 0]
    exp = np.array([(M[b] * (g * np.where(fold == 1, s2b[b], s1b[b]) - mo)).sum() for b in range(B)])
    assert np.allclose(got, exp)
    assert np.allclose(np.asarray(X), g * np.where(fold == 1, 1.1, 1.0) - mo)


def test_xf_ensemble_and_design_total():
    import numpy as np, pandas as pd, rc_core as C, rc_design as D
    rng = np.random.default_rng(2); n, B = 30, 5
    g = rng.uniform(1, 3, n); fold = rng.integers(1, 3, n); s1b = rng.uniform(0.9, 1.1, B); s2b = rng.uniform(0.9, 1.1, B)
    X = C.XF(np.zeros(n), [(g, fold, 1.0, s1b, 1.0, s2b)]); other = rng.uniform(1, 3, n)
    E, common = C.ens_mean_xf([X, other])
    assert common.all() and np.allclose(np.asarray(E), (g + other) / 2)
    d = pd.DataFrame(dict(STRATUM_CN=["a"] * 15 + ["b"] * 15, c_h=[10.0] * 15 + [20.0] * 15, P2POINTCNT=[15] * 15 + [15] * 15))
    M = rng.poisson(1, (B, n)).astype(float)
    pt, bt = D.design_total(d, X, M)
    for b in range(B):
        v = g * np.where(fold == 1, s2b[b], s1b[b]); _, bb = D.design_total(d, v, M[b:b + 1])
        assert np.isclose(bt[b], bb[0])
