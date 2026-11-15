"""Track R per-state scoring from a loaded bundle (see rc_adapt.load_state). Implements the frozen rules:
pillars and MPI (item 2), C5 as P4 net only (item 3), common statewide basis (item 4; by-domain view = S10),
cross-fit recalibration refit in every draw (item 5, v2 items 1 to 5), eligibility, half rule, in-sample exclusions,
ensembles E_state / E_recal / E_plot / median (item 7, v2 items 7 to 11, v3 items 1 and 8), P(rank 1) set (v2 item 15),
noise-corrected typical error (v2 item 14), stress options (item 11, v4 item 5). Aggregates only in the returned tables.
"""
import numpy as np
import pandas as pd
import rc_core as C
import rc_comp as P

DEFAULT_OPT = dict(view="statewide", c5_quant="net", treated_subset="all", curve_rule="level", p3="P3", ens_how="mean",
                   p5_target="A1", min_chains=100, zero=C.ERR_ZERO, weights=C.W, a8_cutoff=0.5, ens_c5_curves="exclude",
                   drop_members=(), start="fia", exclude_acd_c2p2=False)


def _S(res, view):
    if res is None:
        return None
    k = "S_dom" if view == "by_domain" else "S_sw"
    pt, bt = res.get(k + "_pt"), res.get(k + "_bt")
    if pt is None or not np.isfinite(pt):
        return None
    return pt, np.asarray(bt, float)


def _pp(S, zero):
    return None if S is None else (float(C.pts(S[0], zero)), C.pts(S[1], zero), S)


def score_row(name, comps, insample, opt, B, extra_note=""):
    """comps: dict component -> result dicts (C1 statewide, C2, C3, C5, P3 (realism, coverage), P5). Returns row dict + _bt."""
    view, zero = opt["view"], opt["zero"]; notes = []
    c1 = None if "C1" in insample else _S(comps.get("C1"), "statewide")   # C1 statewide in every view (Maine K-R1); in-sample C1 excluded
    if "C1" in insample and comps.get("C1") is not None:
        notes.append("P1 C1 in sample (register), excluded")
    c2 = None if "C2" in insample else _S(comps.get("C2"), view)
    if "C2" in insample and comps.get("C2") is not None:
        notes.append("P1 C2 in sample, excluded")
    P = {}
    if c1 is not None:
        s_pt = C.p1_S(c1[0], None if c2 is None else c2[0]); s_bt = C.p1_S(c1[1], None if c2 is None else c2[1])
        P["P1"] = (float(C.pts(s_pt, zero)), C.pts(s_bt, zero), (s_pt, s_bt))
    else:
        P["P1"] = None; notes.append("P1 no C1")
    c3 = None if "P2" in insample else _S(comps.get("C3"), view)
    P["P2"] = _pp(c3, zero)
    if P["P2"] is None:
        notes.append("P2 not evaluated" + (" (in sample)" if "P2" in insample else ""))
    p3 = comps.get(opt["p3"]) or (comps.get("P3") if opt["p3"] == "P3s" else None)
    if p3 is not None and p3.get("coverage", 0) >= C.P3_MIN_COV and np.isfinite(p3.get("realism", np.nan)):
        v = 100.0 * p3["realism"]; P["P3"] = (v, np.full(B, v), (p3["realism"], None))
    else:
        P["P3"] = None; notes.append("P3 not evaluated")
    c5 = _S(comps.get("C5"), view)
    P["P4"] = _pp(c5, zero)
    if P["P4"] is None:
        notes.append("P4 (C5) not evaluated" + (": " + comps["C5"]["flag"] if comps.get("C5") and comps["C5"].get("flag") else ""))
    c4 = None if "P5" in insample else _S(comps.get("P5"), "statewide")
    P["P5"] = _pp(c4, zero)
    if P["P5"] is None:
        notes.append("P5 not evaluated")
    m_pt, m_bt, cov, ev = C.mpi({k: (v[0], v[1]) if v else None for k, v in P.items()}, B, opt["weights"])
    row = dict(row=name, MPI=m_pt, coverage=cov, tier=C.tier(cov), grade=C.grade(m_pt), pillars_evaluated=",".join(ev),
               notes="; ".join(notes + ([extra_note] if extra_note else [])), _bt=m_bt)
    lo, hi = (np.nanpercentile(m_bt, [2.5, 97.5]) if np.isfinite(m_bt).any() else (np.nan, np.nan))
    row.update(MPI_lo95=lo, MPI_hi95=hi)
    for k in C.PILLARS:
        v = P.get(k)
        row[k] = v[0] if v else np.nan
        row[k + "_lo95"], row[k + "_hi95"] = (np.nanpercentile(v[1], [2.5, 97.5]) if v else (np.nan, np.nan))
        if v and k in ("P1", "P2", "P4", "P5"):
            S = v[2]
            row[k + "_typical_error"] = float(C.typical_error(S[0]))
            row[k + "_te_noise_corrected"] = C.te_noise_corrected(S[0], S[1])
    return row


def rank_rows(rows, comparison):
    """P(rank 1) over the comparison set restricted to tier 1 rows (v2 item 15, v3 item 1)"""
    ids = [i for i, r in enumerate(rows) if r["row"] in comparison and r["tier"] == 1 and np.isfinite(r["MPI"])]
    pr = C.p_rank1(np.stack([rows[i]["_bt"] for i in ids], 0)) if ids else []
    for r in rows:
        r["P_rank1"] = np.nan
    for i, p in zip(ids, pr):
        rows[i]["P_rank1"] = float(p)
    return rows


def table(rows):
    T = pd.DataFrame([{k: v for k, v in r.items() if k != "_bt"} for r in rows])
    return T
