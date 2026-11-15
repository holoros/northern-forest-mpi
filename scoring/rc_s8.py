"""S8, the TreeMap initialization contrast (item 8, v2 item 23, amendment A5). Per engine and state: 2024 AG live stock from the
TreeMap start against the FIA-start stock and the EVALIDator stock (C1 stock), and the forward P5 shape, the 2024 to 2050 stock
index from the TreeMap start against the same index from the FIA start (|ln| of their ratio). The mapped share of masked TreeMap
pixels that land on in-state FIA plots is reported (A5). Aggregates only. Engines: LANDIS-II, CBM gated, CBM AF, FIA curves, CEM.
"""
import glob, json, os, re
from pathlib import Path
import numpy as np
import pandas as pd

H = Path(os.path.expanduser("~")); RUNS = H / "perseus_me/ncasi_rc/runs"


def _curves(st):
    f = RUNS / "statewide" / st / "statewide_curves.csv"
    if not f.exists():
        return {}
    X = pd.read_csv(f); out = {}
    share = np.nan
    m = RUNS / "statewide" / st / "statewide_meta.json"
    if m.exists():
        t = json.loads(m.read_text()).get("treemap")
        share = t.get("in_state_plot_share", np.nan) if isinstance(t, dict) else np.nan
    for e in X.engine.unique():
        tm = X[(X.engine == e) & (X.component == "TreeMap_start") & (X.quantity == "stock_no_harvest")].set_index("year").value
        fs = X[(X.engine == e) & (X.component == "P5_forward") & (X.quantity == "stock_no_harvest")].set_index("year").value
        if len(tm) and len(fs):
            out[e] = dict(tm=tm, fia=fs, mapped_share=share)
    return out


def _landis(st):
    base = RUNS / "landis" / st / "restricted"
    def ser(tag):
        g = sorted(glob.glob(str(base / f"statewide_{tag}*.csv"))) + sorted(glob.glob(str(base / f"statewide_{tag}/*.csv")))
        for f in g:
            x = pd.read_csv(f)
            if {"year", "stock_Tg"} <= set(x.columns):
                return x.set_index("year").stock_Tg
        return None
    tm, fs = ser("tm"), ser("fia")
    if tm is None or fs is None:
        return {}
    off = lambda s: pd.Series(s.values, index=2024 + s.index.astype(int))  # noqa: E731   LANDIS years are run offsets from 2024
    return {"landis": dict(tm=off(tm), fia=off(fs), mapped_share=np.nan)}


def _cem(st):
    d = RUNS / "cem" / st / "restricted"
    def ser(sub):
        fs = sorted(glob.glob(str(d / sub / "series_f*.csv.gz")))
        if not fs:
            return None
        X = pd.concat([pd.read_csv(f) for f in fs])
        yc = [c for c in X.columns if c.lower() in ("year", "yr")]; sc = [c for c in X.columns if c.lower().startswith("stock")]
        if not yc or not sc:
            return None
        return X.groupby(yc[0])[sc[0]].sum()
    tm, fs = ser("values_tm"), ser("values_p5nh")
    if tm is None or fs is None:
        return {}
    meta = sorted(glob.glob(str(d / "sw_lists_tm/lists_meta.json")))
    share = json.loads(Path(meta[0]).read_text()).get("in_state_pixel_share", np.nan) if meta else np.nan
    return {"cem": dict(tm=tm, fia=fs, mapped_share=share)}


def contrast(st, evalidator_stock_tg, y0=2024, y1=2050):
    rows = []
    for e, v in {**_curves(st), **_landis(st), **_cem(st)}.items():
        tm, fs = v["tm"], v["fia"]
        def at(s, y):
            return float(s.loc[y]) if y in s.index else np.nan
        tm0, fs0 = at(tm, y0), at(fs, y0)
        yy = min(y1, int(max(tm.index)), int(max(fs.index)))
        itm, ifs = at(tm, yy) / tm0, at(fs, yy) / fs0
        rows.append(dict(state=st, engine=e, stock_2024_treemap_Tg=tm0, stock_2024_fia_start_Tg=fs0, stock_evalidator_Tg=evalidator_stock_tg,
                         ratio_treemap_to_evalidator=tm0 / evalidator_stock_tg, ratio_fia_start_to_evalidator=fs0 / evalidator_stock_tg,
                         starting_stock_gap=tm0 / fs0 - 1.0, index_end_year=yy, index_treemap=itm, index_fia_start=ifs,
                         p5_shape_abs_ln=abs(np.log(itm / ifs)) if itm > 0 and ifs > 0 else np.nan, mapped_pixel_share=v["mapped_share"]))
    return pd.DataFrame(rows)
