#!/usr/bin/env python3
"""Sensitivity driver of 15 Nov 2026 (see PREREG.md in this folder). Imports the track scorer ns_score.py unchanged and calls its
own functions in the order of ns_score.main; scratch under this folder; rule set a, common basis, c5_net only. Aggregates only.
usage: sens_driver.py ST [--shuffle]"""
import copy, hashlib, importlib.util, json, os, sys, time
from pathlib import Path
import numpy as np
import pandas as pd
H = Path.home(); WT = H / "perseus_me/tracks/multistate"; S = H / "perseus_me/tracks/ncasi_ms/sens_20261115"
_sp = importlib.util.spec_from_file_location("ns_score", WT / "tools/ncasi_states/ns_score.py")
NS = importlib.util.module_from_spec(_sp); _sp.loader.exec_module(NS)
C, A, RR, R, TRK = NS.C, NS.A, NS.RR, NS.R, NS.TRK
REG = S / "IN_SAMPLE_REGISTER_landisC1in.csv"
VIEW, BASIS, K = NS.VIEWS["c5_net"], "statewide", "a"


def tab(rows, run):
    T = pd.DataFrame([{x: y for x, y in r.items() if x != "_bt"} for r in rows]); T.insert(0, "run", run); return T


def maxdiff(r1, r2):
    d = [abs(x["MPI"] - y["MPI"]) for x, y in zip(r1, r2) if np.isfinite(x["MPI"]) and np.isfinite(y["MPI"])]
    assert [x["row"] for x in r1] == [y["row"] for y in r2]
    return dict(max_abs=float(max(d)) if d else 0.0, rows_compared=len(d))


def main():
    st = sys.argv[1].upper(); shuf = "--shuffle" in sys.argv
    assert not os.environ.get("RC_REGISTER"), "RC_REGISTER must be unset at start"
    work = S / "restricted/score_work" / st; work.mkdir(parents=True, exist_ok=True); os.chmod(work, 0o700)
    out = S / "scores" / st; out.mkdir(parents=True, exist_ok=True); t0 = time.time(); tm = {}
    log = lambda s: print(str(s)[:300], flush=True)  # noqa: E731
    acc = dict(state=st, seed=NS.SEED, draws=NS.B, ns_score_md5=hashlib.md5(Path(NS.__file__).read_bytes()).hexdigest())
    ovb = NS.overlay(st, work / "ov_base", None, True); b = NS.load(st, ovb, NS.SEED, log); tm["load"] = time.time() - t0
    for row, arm in NS.ARMS.items():      # as ns_score.main
        if st in NS.NE and (NS.NSR / arm / st / "restricted/pairs_all_default/pair_rows.csv").exists():
            bx = NS.load(st, NS.overlay(st, work / f"ov_{row}", arm, True), NS.SEED, log)
            E = dict(bx["engines"]["fvs"]); E["role"] = "member"; E["arm"] = arm; b["engines"][row] = E
            b["missing"] += [f"{row}: {m}" for m in bx["missing"] if "fvs" in m.lower()]
    NS.set_runs(ovb); sets = NS.rule_sets(st, b); tm["arms"] = time.time() - t0
    mem = list(b["members"]); aug = copy.deepcopy(sets)
    for m in mem:                           # task (c): leave-one-out ensembles, scorer tuple format (name, members, recal, how, ranked)
        aug[K]["ens"].append((f"E_state_wo_{m}", [x for x in mem if x != m], False, "mean", False))
    acc["p3_report"] = NS.attach_p3_all(st, b, aug, work / "p3", log); tm["p3"] = time.time() - t0
    ins0 = set(b["insample"])
    # ---- switch off: registered path
    r0, m0 = NS.score_view(b, sets[K], K, VIEW, BASIS); tm["base"] = time.time() - t0
    T = [tab(r0, "base")]; M = [pd.DataFrame(m0).assign(run="base")]
    T[0].to_csv(out / "sens_scores_base.csv", index=False, float_format="%.6g")
    # ---- task (c)
    r1, m1 = NS.score_view(b, aug[K], K, VIEW, BASIS); tm["loo"] = time.time() - t0
    acc["loo_base_rows_equal_base"] = maxdiff(r0, r1[:len(r0)]) if [x["row"] for x in r1[:len(r0)]] == [x["row"] for x in r0] else "row order differs"
    t1 = tab(r1, "loo"); pr, top = {}, {}
    for m in mem:
        nm = f"E_state_wo_{m}"; sub = [dict(r) for r in r1]; R.rank_rows(sub, set(mem) | {nm})
        pr[nm] = next(x["P_rank1"] for x in sub if x["row"] == nm)
        t1m = [x for x in r1 if x["row"] in mem and x["tier"] == 1 and np.isfinite(x["MPI"])]
        best = max(t1m, key=lambda x: x["MPI"]) if t1m else None
        top[nm] = (best["row"], best["MPI"]) if best else ("", np.nan)
    t1["P_rank1_loo_set"] = t1.row.map(pr); t1["best_tier1_member"] = t1.row.map(lambda x: top.get(x, ("", np.nan))[0])
    t1["best_tier1_member_MPI"] = t1.row.map(lambda x: top.get(x, ("", np.nan))[1])
    T.append(t1); M.append(pd.DataFrame(m1).assign(run="loo"))
    # ---- task (b): scorer's own register override
    os.environ["RC_REGISTER"] = str(REG); ins2, note2 = A.load_insample(st); del os.environ["RC_REGISTER"]
    acc["register_switch"] = dict(added=sorted(map(list, ins2 - ins0)), removed=sorted(map(list, ins0 - ins2)), note=note2)
    b2 = dict(b); b2["insample"] = ins2
    r2, m2 = NS.score_view(b2, sets[K], K, VIEW, BASIS); tm["c1in"] = time.time() - t0
    T.append(tab(r2, "c1in")); M.append(pd.DataFrame(m2).assign(run="c1in"))
    r0b, _ = NS.score_view(b, sets[K], K, VIEW, BASIS); acc["base_repeat_after_switches"] = maxdiff(r0, r0b); tm["base2"] = time.time() - t0
    if shuf:                                # scorer hook run_track_r.shuffle_bundle
        sb = TRK.shuffle_bundle(b, 11); rs, _ = NS.score_view(sb, aug[K], K, VIEW, BASIS); acc["shuffle_base_and_loo"] = maxdiff(r1, rs)
        sb["insample"] = ins2; rs2, _ = NS.score_view(sb, sets[K], K, VIEW, BASIS); acc["shuffle_c1in"] = maxdiff(r2, rs2); tm["shuffle"] = time.time() - t0
    X = pd.concat(T, ignore_index=True); X.insert(0, "state", st)
    X.to_csv(out / "sens_scores.csv", index=False, float_format="%.6g"); pd.concat(M, ignore_index=True).to_csv(out / "sens_membership.csv", index=False)
    acc["missing"] = b["missing"]; acc["members"] = mem; acc["timing_s"] = {k: round(v, 1) for k, v in tm.items()}; acc["wall_s"] = round(time.time() - t0, 1)
    (out / "checks.json").write_text(json.dumps(acc, indent=1, default=float))
    print("SENS_DONE", st, json.dumps({k: acc.get(k) for k in ("loo_base_rows_equal_base", "register_switch", "base_repeat_after_switches", "shuffle_base_and_loo", "shuffle_c1in", "timing_s")}, default=float), flush=True)


if __name__ == "__main__":
    main()
