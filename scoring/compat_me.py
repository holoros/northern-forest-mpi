"""Maine compatibility check for the Track R scorer (aggregates only). Where the Track R rules coincide with MPI v2.3/v2.4,
the scorer arithmetic must reproduce the release values; where they differ by design, the difference is listed.
(1) MPI, coverage, tier, grade from the released pillar points (mpi23 common basis, every row);
(2) pillar points from Maine's own S draws (draws_s / mpi23_s npz): P1 = C1 statewide_basis + C2 a1b_common, P2 = C3,
    P5 = C4bfm statewide_basis, point and 95% interval, for single engine rows with all keys present;
(3) C5 points from the v2.4 C5 table S (net only, statewide and by domain; three quantity);
(4) FIA design path: rc_design per-plot totals reproduce the six-state EVALIDator targets (path A/B file).
Writes ME/compat_me.json and compat_me_rows.csv under ~/perseus_me/ncasi_rc/scorer_tests."""
import json, os, sys
from pathlib import Path
import numpy as np
import pandas as pd
sys.path.insert(0, str(Path(__file__).parent))
import rc_core as C
import rc_design as D

H = Path(os.path.expanduser("~")); MS = H / "perseus_me/me_state"; OUT = H / "perseus_me/ncasi_rc/scorer_tests"; OUT.mkdir(parents=True, exist_ok=True)
res = {}
# (1)
T = pd.read_csv(MS / "mpi23/release/mpi23_common_basis_scores.csv")
rows = []
for r in T.itertuples():
    P = {p: (getattr(r, p), np.full(2, getattr(r, p))) if np.isfinite(getattr(r, p)) else None for p in C.PILLARS}
    m, _, cov, ev = C.mpi(P, 2)
    rows.append(dict(view=r.view, row=f"{r.engine}.{r.arm}", MPI_rel=r.MPI, MPI_rc=m, cov_rel=r.coverage, cov_rc=cov, tier_rel=r.tier, tier_rc=C.tier(cov),
                     grade_rel=r.grade, grade_rc=C.grade(m)))
A = pd.DataFrame(rows)
res["agg"] = dict(rows=len(A), max_abs_MPI=float(np.nanmax(np.abs(A.MPI_rel - A.MPI_rc))), max_abs_cov=float(np.nanmax(np.abs(A.cov_rel - A.cov_rc))),
                  tier_equal=bool((A.tier_rel == A.tier_rc).all()), grade_equal=bool((A.grade_rel.fillna("") == A.grade_rc).all()))
# (2)
Z = {}
for f in ("mpi23/work/draws_s.npz", "mpi23/work/mpi23_s.npz"):
    z = np.load(MS / f); Z.update({k: z[k] for k in z.files})
def S(e, a, comp, var):
    k = f"S|{e}|{a}|{comp}|{var}"
    return (float(Z[k + "|pt"][0]), Z[k + "|bt"]) if k + "|pt" in Z else None
cb = T[(T.view == "common basis") & (T.row_group == "original arm")]
RT = pd.read_csv(MS / "mpi23/work/mpi23_s_ratios.csv")
def S_ratio(e, a, comp, var):
    x = RT[(RT.engine == e) & (RT.arm == a) & (RT.component == comp) & (RT.variant == var) & (RT.domain == "statewide")]
    return float(np.mean(np.abs(np.log(x.ratio.values)))) if len(x) else None
pr = []
for r in cb.itertuples():
    c1 = S_ratio(r.engine, r.arm, "C1", "statewide_basis"); c2 = S(r.engine, r.arm, "C2", "a1b_common")
    c3 = S(r.engine, r.arm, "C3", "primary"); c4 = S_ratio(r.engine, r.arm, "C4bfm", "statewide_basis")
    d = dict(row=f"{r.engine}.{r.arm}")
    if c1 is not None:
        ins = "C2 in sample" in str(r.notes)
        d.update(P1_rel=r.P1, P1_rc=float(C.pts(C.p1_S(c1, None if (c2 is None or ins) else c2[0]))))
    if c3 is not None and np.isfinite(r.P2):
        d.update(P2_rel=r.P2, P2_rc=float(C.pts(c3[0])))
    if c4 is not None and np.isfinite(r.P5):
        d.update(P5_rel=r.P5, P5_rc=float(C.pts(c4)))
    pr.append(d)
PR = pd.DataFrame(pr)
res["pillars_from_S"] = {c: dict(n=int(PR[c + "_rel"].notna().sum()), max_abs=float(np.nanmax(np.abs(PR[c + "_rel"] - PR[c + "_rc"]))) if c + "_rel" in PR else None)
                         for c in ("P1", "P2", "P5") if c + "_rel" in PR}
res["pillars_from_S"]["note"] = "point values; ratios table printed to 6 significant figures (rounding bound about 1e-4 points)"
# (3)
c5 = pd.read_csv(MS / "mpi24/release/mpi24_c5_table.csv")
dd = {}
for v in ("a1b_netonly_statewide", "a1b_netonly", "a1b_statewide", "a1b"):
    s, p = c5[f"{v}_S"], c5[f"{v}_points"]; ok = s.notna() & p.notna()
    dd[v] = dict(n=int(ok.sum()), max_abs=float(np.max(np.abs(C.pts(s[ok]) - p[ok]))) if ok.any() else None)
res["c5_points_from_S"] = dd
# (4)
tg = pd.read_csv(H / "perseus_me/ncasi_rc/targets/targets_six_states.csv")
dz = []
for st in ("ME", "NH", "VT", "MI", "WI", "MN"):
    sc = D.statecd(st)
    y, d = D.plot_stock(st, sc * 10000 + 2401); stock = D.design_total(d, y)[0] * D.LB2TG
    g, d3 = D.plot_grm(st, sc * 10000 + 2403)
    vals = dict(stock=stock, **{q: D.design_total(d3, g[q])[0] * D.LB2TG for q in ("gross_growth", "mortality", "net_growth")})
    for q, v in vals.items():
        ref = tg[(tg.state == st) & (tg.quantity == q)]
        col = "ours" if "ours" in tg.columns else tg.columns[3]
        if len(ref):
            dz.append(dict(state=st, quantity=q, rc=v, ref=float(ref[col].iloc[0]), rel=abs(v / float(ref[col].iloc[0]) - 1)))
    D._rd.cache_clear()
DZ = pd.DataFrame(dz)
res["design_vs_targets"] = dict(cells=len(DZ), max_rel=float(DZ.rel.max()) if len(DZ) else None)
res["differs_by_design"] = ["P4: Track R uses C5 net only (item 3); Maine registered P4 is the Spearman pillar (v2.3/v2.4 registered)",
    "recalibrated rows: Track R refits cross-fit scalars in every draw on all components (item 5, v2 items 1 to 3); Maine C1 and P3 to P5 used released v2.2 multipliers",
    "ensembles: Track R E_state lineup (item 7) and common-unit C2/C5 (v2 item 9) and curve exclusion from ensemble C5 (v3 item 8) differ from Maine E5/E9",
    "C1 targets: Track R recomputes targets per draw on the 2024 evaluation universe with a location-shared draw (v5 item 1); Maine bd1/bd3 draws",
    "P5: Track R index against removals-adjusted FIA series 2004/2006 to 2024 (A1); Maine C4bfm = 2050/2025 index of bau_fia_matched vs FIA trend"]
OUT.joinpath("ME").mkdir(exist_ok=True)
(OUT / "ME/compat_me.json").write_text(json.dumps(res, indent=1, default=float))
pd.concat([A.assign(part="aggregate"), PR.assign(part="pillars_from_S"), DZ.assign(part="design")], ignore_index=True).to_csv(OUT / "ME/compat_me_rows.csv", index=False)
print(json.dumps(res, indent=1, default=float))
