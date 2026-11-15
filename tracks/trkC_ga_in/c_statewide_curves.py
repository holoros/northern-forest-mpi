#!/usr/bin/env python3
"""Track R (NCASI report card) statewide runs for the curve engines (CBM gated, CBM American Forests, FIA hybrid curves),
prereg v4 item 3: area-by-stratum runs, never summed pair predictions. Each inventory plot's forest conditions carry the
curve from their stand age (build_j.unit_values, CONDPROP weighted, lb/ac of plot area) and are expanded with the
evaluation's EXPNS x ADJ_FACTOR_SUBP. No scoring.
  C1  stock 2024 (ST2401 plots, curve level at age); gross, mortality, net Tg/yr over t1->t2 on the ST2403 plots' t1
      conditions (horizon round(REMPER), no harvest)
  P5  annual trajectories from the earliest series evaluation (2004 MI WI MN, 2006 ME NH VT) to 2024 and from 2024 to 2050,
      no harvest and FIA-matched removals (level rule: level x (1 - h)^t, h = state removals / stock, targets file)
  TM  TreeMap 2022 start (masked state pixels, pixel PLT_CN -> that plot's conditions when the plot is in this state's
      FIADB; unmatched pixel share reported): stock and 2024->2050 trajectories
Usage: statewide_curves_state.py ST OUT
trkC_gain copy: CBM gated only; evaluations ST2501/ST2503 (multistate current); multistate trend series; removal share nan (not scored); multistate pulls."""
import json, os, sys, time
from pathlib import Path
import numpy as np, pandas as pd
st = sys.argv[1].upper(); OUT = Path(sys.argv[2]); OUT.mkdir(parents=True, exist_ok=True); os.chmod(OUT, 0o700)
H = Path.home(); RC = H / "perseus_me/ncasi_rc"; ENGD = RC / "engines"; PUL = H / "perseus_me/multistate/restricted/pulled_cardinal"
if st == "ME":
    root = H / "perseus_me/me_state"; fits = root / "v23_trackA/restricted/ycx_ME_hybrid_fits.csv"; memb = root / "v23_trackA/restricted/ycx_membership_ME.csv"
else:
    root = H / "perseus_me/tracks/trkC_gain/restricted/curve_roots" / st; fits = PUL / f"yield_curves_conus/ycx_{st}_hybrid_fits.csv"; memb = PUL / f"yield_curves_conus/config/ycx_membership_{st}.csv"
fia = H / f"restricted/ncasi-modeleval/fia_data_rc/{st}"
os.environ.update(PME_STATE=st, PME_STATE_ROOT=str(root), PME_FIA_ROOT=str(fia))
CODE = H / "perseus_me/code"; sys.path.insert(0, str(CODE)); sys.path.insert(0, str(CODE / "tools/v24_H")); sys.path.insert(0, str(CODE / "tools/complete"))
import prep_h as PH  # noqa: E402
BJ = PH._bj(); BJ.FIA = fia; BJ.MEMB = memb; BJ.FITS = fits
LB2TG = 0.45359237e-9; ACRE_HA = 0.40468564224; t0 = time.time()
rd = lambda t, cols, **k: pd.read_csv(fia / f"{st}_{t}.csv", usecols=cols, **k)  # noqa: E731
cs = lambda s: s.astype(str).str.strip()  # noqa: E731
code = int(rd("PLOT", ["STATECD"], nrows=1).STATECD.iloc[0])
pl = rd("PLOT", ["CN", "PREV_PLT_CN", "REMPER", "INVYR"], dtype={"CN": str, "PREV_PLT_CN": str}); pl["CN"] = cs(pl.CN); pl["PREV_PLT_CN"] = cs(pl.PREV_PLT_CN.fillna(""))
pa = rd("POP_PLOT_STRATUM_ASSGN", ["PLT_CN", "STRATUM_CN", "EVALID"], dtype={"PLT_CN": str, "STRATUM_CN": str}); pa["PLT_CN"] = cs(pa.PLT_CN); pa["STRATUM_CN"] = cs(pa.STRATUM_CN)
ps = rd("POP_STRATUM", ["CN", "EVALID", "EXPNS", "ADJ_FACTOR_SUBP"], dtype={"CN": str}); ps["CN"] = cs(ps.CN)
def weights(ev):
    x = pa[pa.EVALID == ev].merge(ps[ps.EVALID == ev].rename(columns={"CN": "STRATUM_CN"})[["STRATUM_CN", "EXPNS", "ADJ_FACTOR_SUBP"]], on="STRATUM_CN")
    return pd.Series((x.EXPNS * x.ADJ_FACTOR_SUBP).values, index=x.PLT_CN.values)   # acres per plot (area basis)
ENGS = ["cbm_gated"]
have_af = (root / "complete_j/work/cbm_curves_af/curve_table.csv").exists() if st != "ME" else True
if st != "ME":
    w = root / "complete_j/work"; w.mkdir(parents=True, exist_ok=True)
    for e in ("gated", "af"):
        src = H / f"perseus_me/tracks/trkC_gain/restricted/cbm_curves/{e}_{st}"; dst = w / f"cbm_curves_{e}"
        if src.exists() and not dst.exists():
            dst.symlink_to(src)
h = float("nan"); meta = dict(state=st, removal_share_per_yr=h, engines={}, scoring="none")
def run(units_plt, offs, wts, tag):
    U0 = pd.DataFrame(dict(unit=units_plt, PLT_CN=units_plt, kind="plot")); U, C = PH.units_for(BJ, U0); res = {}
    for e in ENGS:
        if e == "cbm_af" and not have_af:
            continue
        pv = BJ.unit_values(e, U, C, np.asarray(offs)); B = np.zeros(len(offs)); M = np.zeros(len(offs)); ok_a = 0.0; all_a = 0.0
        for u in units_plt:
            wgt = float(wts.get(u, 0.0)); all_a += wgt; v, note = pv[u]
            if v is None:
                continue
            _, b, m = v; ok_a += wgt; B += wgt * b
            if e != "fia_curves":
                M += wgt * m
        res[e] = dict(B=B * LB2TG, M=M * LB2TG, area_ok_share=ok_a / all_a if all_a else np.nan)
    return res
rows = []
# ---- C1 stock 2024 and P5 forward
ev1 = code * 10000 + 2501; w1 = weights(ev1); P1 = sorted(set(w1.index) & set(pl.CN))
offs_f = np.arange(0, 27); R1 = run(P1, offs_f, w1, "fwd")
for e, r in R1.items():
    rows.append(dict(engine=e, component="C1", quantity="ag_live_stock_2024", value=r["B"][0], unit="Tg", plot_area_ok_share=r["area_ok_share"], note="curve level at inventory age, ST2401 plots"))
    for k, yr in enumerate(range(2024, 2051)):
        for hv, lab in ((0.0, "no_harvest"), (h, "fia_matched_removals")):
            rows.append(dict(engine=e, component="P5_forward", quantity=f"stock_{lab}", year=yr, value=r["B"][k] * (1 - hv) ** k, unit="Tg", plot_area_ok_share=r["area_ok_share"]))
# ---- C1 rates on the ST2403 plots' t1 conditions
ev3 = code * 10000 + 2503; w3 = weights(ev3); t2 = pl[pl.CN.isin(set(w3.index)) & (pl.PREV_PLT_CN != "")]
t2 = t2[t2.PREV_PLT_CN.isin(set(pl.CN))]; t2["n"] = [max(1, int(round(r))) for r in t2.REMPER.fillna(5)]
G = Gm = N = 0.0; aok = aall = 0.0; per = {e: dict(G=0.0, M=0.0, ok=0.0) for e in ENGS}
for n, g in t2.groupby("n"):
    wt = pd.Series(w3.reindex(g.CN).values, index=g.PREV_PLT_CN.values)
    U0 = pd.DataFrame(dict(unit=g.PREV_PLT_CN.values, PLT_CN=g.PREV_PLT_CN.values, kind="plot")); U, C = PH.units_for(BJ, U0)
    for e in ENGS:
        if e == "cbm_af" and not have_af:
            continue
        pv = BJ.unit_values(e, U, C, np.array([0, n]))
        for u in g.PREV_PLT_CN.values:
            v, _ = pv[u]
            if v is None:
                continue
            _, b, m = v; wgt = float(wt.get(u, 0.0)); per[e]["ok"] += wgt
            net = (b[1] - b[0]) / n; mort = 0.0 if e == "fia_curves" else (m[1] - m[0]) / n
            per[e]["G"] += wgt * (net + mort); per[e]["M"] += wgt * mort
area3 = float(w3.reindex(t2.CN).sum()); area_all3 = float(w3.sum())
for e, d in per.items():
    if e == "cbm_af" and not have_af:
        continue
    for q, v in (("gross", d["G"]), ("mortality", d["M"] if e != "fia_curves" else np.nan), ("net", d["G"] - d["M"])):
        rows.append(dict(engine=e, component="C1", quantity=q, value=v * LB2TG, unit="Tg/yr", plot_area_ok_share=d["ok"] / area3 if area3 else np.nan,
                         note=f"t1 conditions of ST2403 remeasured plots, no harvest; remeasured plot area share of eval {area3 / area_all3:.4f}" + ("; curves carry no mortality: gross = net" if e == "fia_curves" else "")))
# ---- P5 hindcast series from the earliest series evaluation
ser = pd.read_csv(H / "perseus_me/multistate/inputs/trend/trend_stock_series.csv"); ser = ser[ser.state == st].sort_values("end_invyr"); ev0 = int(ser.evalid.iloc[0]); y0 = int(ser.end_invyr.iloc[0])
w0 = weights(ev0); P0 = sorted(set(w0.index) & set(pl.CN)); offs_b = np.arange(0, 2024 - y0 + 1); R0 = run(P0, offs_b, w0, "back")
for e, r in R0.items():
    for k in range(len(offs_b)):
        for hv, lab in ((0.0, "no_harvest"), (h, "fia_matched_removals")):
            rows.append(dict(engine=e, component="P5_series", quantity=f"stock_{lab}", year=y0 + k, value=r["B"][k] * (1 - hv) ** k, unit="Tg", plot_area_ok_share=r["area_ok_share"], note=f"start evalid {ev0}"))
# ---- TreeMap start
try:
    import struct
    VAT = H / "landis2_fb/TREEMAP/TM2022/TreeMap2022_CONUS.tif.vat.dbf"; b = VAT.read_bytes(); n_, hl, rl = struct.unpack("<IHH", b[4:12]); fl = []; o = 32
    while b[o] != 0x0D:
        fl.append((b[o:o + 11].split(b"\0")[0].decode(), b[o + 16])); o += 32
    off = {}; pos = 1
    for nm, ln in fl:
        off[nm] = (pos, ln); pos += ln
    rec = np.frombuffer(b[hl:hl + n_ * rl], dtype=f"S{rl}")
    get = lambda f: [r[off[f][0]:off[f][0] + off[f][1]].decode(errors="ignore").strip() for r in rec]  # noqa: E731
    V = pd.DataFrame(dict(Value=pd.to_numeric(get("Value")), PLT_CN=[x[:-2] if x.endswith(".0") else x for x in get("PLT_CN")]))
    pxf = OUT.parent.parent / "treemap_px" / f"{st}_px.csv"     # Value, px (masked state pixels), written by treemap_px.py (venv_fig, rasterio)
    cnt = pd.read_csv(pxf).set_index("Value").px
    V = V[V.Value.isin(cnt.index)].copy(); V["px"] = V.Value.map(cnt); V["acres"] = V.px * 900.0 / 4046.8564224
    instate = set(pl.CN); V["in"] = V.PLT_CN.isin(instate); tot = V.acres.sum(); share_in = V.acres[V["in"]].sum() / tot
    wt = V[V["in"]].groupby("PLT_CN").acres.sum(); Pt = sorted(wt.index); Rt = run(Pt, offs_f, wt, "tm")
    for e, r in Rt.items():
        for k, yr in enumerate(range(2024, 2051)):
            for hv, lab in ((0.0, "no_harvest"), (h, "fia_matched_removals")):
                rows.append(dict(engine=e, component="TreeMap_start", quantity=f"stock_{lab}", year=yr, value=r["B"][k] * (1 - hv) ** k, unit="Tg",
                                 plot_area_ok_share=r["area_ok_share"], note=f"in-state imputed plot pixel area share {share_in:.4f}; not rescaled"))
    meta["treemap"] = dict(masked_acres=float(tot), in_state_plot_share=float(share_in))
except Exception as ex:  # noqa: BLE001
    meta["treemap"] = f"not run: {ex!r}"
X = pd.DataFrame(rows); X.insert(0, "state", st); X.to_csv(OUT / "statewide_curves.csv", index=False, float_format="%.10g")
meta.update(wall_s=round(time.time() - t0, 1), evalids=dict(stock=ev1, grm=ev3, series_start=ev0)); (OUT / "statewide_meta.json").write_text(json.dumps(meta, indent=1, default=str)); print(json.dumps(meta, default=str))
