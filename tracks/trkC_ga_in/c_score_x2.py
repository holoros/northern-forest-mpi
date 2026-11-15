#!/usr/bin/env python3
"""trkC_gain scorer extension (PREREG section 1 and 2.5). Imports the frozen tools/multistate/ms_score_x.py (f775cf47, which imports
ms_score.py 25bf7b51) unchanged and, for GA and IN, scores the ms_score_x row set plus the rows of this track whose outputs exist:
cbm_gated (Track R curve layout under SLOTS/cbm_gated/ST, statewide C1/P5 under SLOTS/statewide/ST), cem (SLOTS/cem/ST) and
cbm_gated_x1 (injected like fvsc; foliage already removed in its rows, so it is not in the FOLIAGE set; C5 level = ageq = the X1
cut-pair rows; chains SLOTS/cbm_gated_x1/ST/restricted/h1). Check XX1: every value of the released ms_score_x output
(work/final/score_ST/ms_scores.csv) is reproduced within 1e-9 for the rows present in both (P(rank 1) columns excepted).
Optional --shuffle SEED permutes the order of every unit row file of this track's slots (G2: new rows only). usage: c_score_x2.py ST OUTDIR SLOTS [--shuffle SEED]"""
import argparse, hashlib, importlib.util, json, sys, time
from pathlib import Path
import numpy as np, pandas as pd
WT = Path("/home/aaron/perseus_me/tracks/multistate")
_sp = importlib.util.spec_from_file_location("ms_score_x", WT / "tools/multistate/ms_score_x.py")
MSX = importlib.util.module_from_spec(_sp); _sp.loader.exec_module(MSX)
MSS, NS, A, C = MSX.MSS, MSX.NS, MSX.A, MSX.C
FROZEN_X = "f775cf4710c06f1673489d53424a7675"; X1ARM = "X1"
md5 = lambda p: hashlib.md5(Path(p).read_bytes()).hexdigest()  # noqa: E731


def add_slots2(root, st, slots, log):
    root = Path(root); slots = Path(slots); added = []
    for eng in ("cbm_gated", "cbm_gated_x1", "statewide"):
        src = slots / eng / st
        if src.exists():
            (root / eng).mkdir(exist_ok=True); d = root / eng / st
            if not d.exists():
                d.symlink_to(src); added.append(eng)
    src = slots / "cem" / st
    if src.exists():
        d = root / "cem" / st
        if d.is_symlink() or d.exists():
            d.unlink() if d.is_symlink() else None
        (root / "cem").mkdir(exist_ok=True)
        if not d.exists():
            d.symlink_to(src); added.append("cem")
    src = slots / "landis" / st
    if src.exists():                                           # LANDIS GA: lane rows (symlinks) + this track's h1 chains
        d = root / "landis" / st
        if d.is_symlink():
            d.unlink()
        if not d.exists():
            d.symlink_to(src); added.append("landis")
    log(f"trkC slots linked: {added}")
    return added


def inject_x1(st, b, slots, log):
    d = Path(slots) / "cbm_gated_x1" / st / "restricted"
    f2, f5 = d / "pairs_level__cbm_gated_x1__X1__pair_rows.csv", d / "treated__cbm_gated_x1__X1__pair_rows.csv"
    if not f2.exists():
        return
    E = dict(engine="cbm_gated_x1", arm=X1ARM, role="member", note="refine_curves2 X1 bounded recursion on the GA/IN CBM gated curves (trkC_gain 2.6); foliage removed in the rows")
    E["C2"] = A.align(A.pair_rows(f2), b["pairs"]["ids"])
    if f5.exists():
        r5 = A.align(A.pair_rows(f5), b["treated"]["ids"]); E["C5"] = {"level": r5, "ageq": r5}
    bx = dict(b); bx["engines"] = {"cbm_gated_x1": E}; bx["chains"] = A._chains(st, A.D.LocDraws(st, MSS.B, MSS.SEED), b["missing"]); A.fill_c3(bx)
    if "C3" not in E:
        b["missing"].append("cbm_gated_x1 P2: chain file not read")
    b["engines"]["cbm_gated_x1"] = E; log(f"cbm_gated_x1 injected: {sorted(k for k in E if k[:1] in 'CP')}")


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("st"); ap.add_argument("out"); ap.add_argument("slots"); ap.add_argument("--shuffle", type=int, default=None)
    a = ap.parse_args(); st, out = a.st.upper(), Path(a.out)
    assert md5(WT / "tools/multistate/ms_score_x.py") == FROZEN_X and md5(WT / "tools/multistate/ms_score.py") == MSX.FROZEN, "frozen scorer md5 changed"
    assert st in MSS.NEW, "GA and IN only"
    out.mkdir(parents=True, exist_ok=True); t0 = time.time(); msgs = []
    log = lambda s: (msgs.append(str(s)[:300]), print(str(s)[:300]))  # noqa: E731
    if a.shuffle is not None:                                  # G2: every unit row file read in permuted order
        rng = np.random.default_rng(a.shuffle); _rc = pd.read_csv; _rp = pd.read_parquet
        def rc(*k, **kw):
            x = _rc(*k, **kw)
            return x.iloc[rng.permutation(len(x))].reset_index(drop=True) if isinstance(x, pd.DataFrame) and len(x) > 1 and "pair_id" in x.columns and "/restricted/slots/" in __import__("os").path.realpath(str(k[0] if k else kw.get("filepath_or_buffer", ""))) else x
        def rp(*k, **kw):
            x = _rp(*k, **kw)
            return x.iloc[rng.permutation(len(x))].reset_index(drop=True) if isinstance(x, pd.DataFrame) and len(x) > 1 and "chain_id" in x.columns and "/restricted/slots/" in __import__("os").path.realpath(str(k[0] if k else kw.get("path", ""))) else x
        pd.read_csv = rc; pd.read_parquet = rp
    work = MSS.MS / "restricted/score_work" / (f"{st}_x2" + ("" if a.shuffle is None else f"_s{a.shuffle}"))
    work.mkdir(parents=True, exist_ok=True)
    for k in ("C2", "P2"):
        C.ELIG[k] = set(C.ELIG[k]) | {"fvsc", "cbm_gated_x1"}
    C.ELIG["C5"] = set(C.ELIG["C5"]) | {"cbm_gated_x1"}
    C.C2_QUANT["fvsc"] = ("net_growth",)
    C.CURVE_DRIVEN = set(C.CURVE_DRIVEN) | {"cbm_gated_x1"}; C.AGGREGATE = set(C.AGGREGATE) | {"cbm_gated_x1"}
    patch, evc = MSS.new_state_patch(st)
    C.MEMBERS.setdefault(st, ["fvs"] + (["cem"] if (MSS.MSR / "cem" / st).exists() else []))
    acc = dict(state=st, scorer="c_score_x2.py", scorer_md5=md5(__file__), ms_score_x_md5=FROZEN_X, frozen_ms_score_md5=MSX.FROZEN, evalid_current=evc,
               seed=MSS.SEED, draws=MSS.B, shuffle=a.shuffle)
    with patch:
        ovb = MSX.add_slots(MSS.overlay_new(st, work / "ov_base", "fvs_default"), st); linked = add_slots2(ovb, st, a.slots, log)
        b = NS.load(st, ovb, MSS.SEED, log)
        if MSS.fvsb_present(st):
            E, m = MSS.load_fvsb(st, MSS.overlay_new(st, work / "ov_fvsb", "fvs_bayes"), log); E["role"] = "member"; b["engines"]["fvsb"] = E; b["missing"] += m
        NS.set_runs(ovb); MSX.inject_fvsc(st, b, log); inject_x1(st, b, a.slots, log)
        sets = MSS.sets_new(st, b); added = [x for x in ("fia_curves", "landis", "fvsc", "cbm_gated", "cem", "cbm_gated_x1") if x in b["engines"]]
        for S in sets.values():
            S["singles"] = list(dict.fromkeys(S["singles"] + added)); S["label"] += "; plus " + ", ".join(added) + " (A4.5, trkC_gain)" if added else ""
        acc["p3_report"] = MSS.attach_p3_ms(st, b, sets, work / "p3", log)
        rows, memb = MSS.score_all(st, b, sets); T = MSS.tab(rows)
    try:   # Aaron 14 Nov decision 4: save the MPI and pillar bootstrap draws of every row (keys rule_set|basis|view|row[|part])
        D, man = {}, []
        for r in rows:
            k0 = "|".join(str(r.get(x)) for x in ("rule_set", "basis", "view", "row")); bt = r.get("_bt")
            items = bt.items() if isinstance(bt, dict) else ([("MPI", bt)] if bt is not None else [])
            for nm, v in items:
                v = np.asarray(v, dtype=float) if v is not None else None
                if v is not None and v.size:
                    D[f"{k0}|{nm}"] = v; man.append(dict(key=f"{k0}|{nm}", n=int(v.size), md5=hashlib.md5(v.tobytes()).hexdigest()))
        np.savez_compressed(out / "draws.npz", **D)
        (out / "draws_manifest.json").write_text(json.dumps(dict(seed=MSS.SEED, draws=MSS.B, arrays=len(D), entries=man), indent=1))
    except Exception as e:  # recorded, never silent
        (out / "draws_manifest.json").write_text(json.dumps(dict(error=repr(e)[:300], seed=MSS.SEED, draws=MSS.B)))
    ref = pd.read_csv(MSS.MS / "work/final" / f"score_{st}" / "ms_scores.csv")
    key = ["rule_set", "basis", "view", "row"]; num = [c for c in ref.columns if ref[c].dtype.kind == "f" and not c.startswith("P_rank1")]
    import io   # D-C8: the released file holds %.6g values, so compare the extension table after the same CSV round trip
    Tc = pd.read_csv(io.StringIO(T.to_csv(index=False, float_format="%.6g")))
    J = ref.merge(Tc, on=key, how="inner", suffixes=("", "_x"), validate="one_to_one"); dmax = 0.0
    for c in num:
        if c + "_x" not in J:
            continue
        x, y = J[c].astype(float).values, J[c + "_x"].astype(float).values
        both = np.isfinite(x) & np.isfinite(y); bad = np.isfinite(x) != np.isfinite(y)
        dmax = max(dmax, float(np.abs(x[both] - y[both]).max()) if both.any() else 0.0, np.inf if bad.any() else 0.0)
    acc["XX1"] = dict(rows_ref=len(ref), rows_joined=len(J), cols=len(num), max_abs_diff=dmax, passed=bool(dmax <= 1e-9 and len(J) == len(ref)))
    acc["linked"] = linked; acc["added_rows"] = added; acc["missing"] = b["missing"]; acc["wall_s"] = round(time.time() - t0, 1)
    T.to_csv(out / "ms_scores.csv", index=False, float_format="%.6g"); pd.DataFrame(memb).to_csv(out / "ms_membership.csv", index=False)
    (out / "checks.json").write_text(json.dumps(acc, indent=1, default=float)); (out / "log.txt").write_text("\n".join(msgs))
    print(json.dumps({k: acc[k] for k in ("XX1", "added_rows", "wall_s")}, default=float))
    if not acc["XX1"]["passed"]:
        sys.exit("XX1 failed")


if __name__ == "__main__":
    main()
