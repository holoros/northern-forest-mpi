#!/usr/bin/env python3
"""trkC_gain claim ledger: every headline number of RELEASE.md with its source file (aggregates only). usage: c_claims.py"""
import csv, json
from pathlib import Path
T = Path.home() / "perseus_me/tracks/trkC_gain"; R = T / "release"; rows = []
def add(cid, claim, value, src, status="VERIFIED by file read"):
    rows.append(dict(claim_id=cid, claim=claim, value=value, source=src.replace(str(Path.home()), "~"), status=status))
g1 = json.loads((T / "work/g1_IN.json").read_text()); add("G1", "IN released rows reproduced (rows; max abs diff)", f"{g1['rows_rerun']}; {g1['max_abs_diff']}", str(T / "work/g1_IN.json"))
for f in sorted((T / "work").glob("xx1_files.json")) + sorted((T / "work").glob("g2_*.json")) + sorted((T / "work").glob("xhost_sit.json")) + sorted((T / "work").glob("verifier*.json")):
    add(f.stem, f.stem, json.dumps(json.loads(f.read_text()))[:300], str(f))
for st in ("GA", "IN"):
    x = json.loads((T / f"work/x1_{st}.json").read_text()); xp = json.loads((T / f"work/x1_{st}_perm20261114.json").read_text())
    add(f"X1_{st}_theta", f"{st} X1 production theta (c0, c2, c3, gamma)", ";".join(f"{v:.4f}" for v in x["params"]["prod"]["theta"]), str(T / f"work/x1_{st}.json"))
    add(f"X1_{st}_perm", f"{st} X1 permuted refit max theta diff", f"{max(abs(a-b) for k in x['params'] for a,b in zip(x['params'][k]['theta'], xp['params'][k]['theta'])):.2e}", str(T / f"work/x1_{st}_perm20261114.json"))
    gt = x["gates"]; add(f"X1_{st}_gates", f"{st} X1 100-yr gates M1 M2 M3 M4; mean y0, y100 Mg/ha", f"{gt['M1_no_negative']} {gt['M2_late_increment_below_early']} {gt['M3_p99_y100_le_current_max']} {gt['M4_mean_y100_le_current_p99']}; {gt['mean_y0']:.1f}, {gt['mean_y100']:.1f}", str(T / f"work/x1_{st}.json"))
for r in csv.DictReader(open(R / "trkC_before_after.csv")):
    for p in ("before", "after"):
        if r.get(f"{p}_MPI"):
            add(f"{r['state']}_{r['row']}_{r['version']}_{p}", f"{r['state']} {r['row']} ({r['version']}) {p} MPI [95%] grade tier", f"{float(r[p+'_MPI']):.1f} [{float(r[p+'_MPI_lo95']):.1f}, {float(r[p+'_MPI_hi95']):.1f}] {r[p+'_grade']} {r[p+'_tier']}",
                str(R / "trkC_before_after.csv") + " (from " + ("multistate/work/final/score_ST" if p == "before" else "restricted/score") + ")")
for r in csv.DictReader(open(R / "trkC_c5_three_ratios.csv")):
    add(f"C5_{r['state']}_{r['row']}", f"{r['state']} {r['row']} C5 model/observed gross, mortality, net (all treated pairs)", f"{float(r['all_gross']):.3f}, {float(r['all_mortality']):.3f}, {float(r['all_net']):.3f} (n {r['all_n']})", str(R / "trkC_c5_three_ratios.csv"))
with open(T / "claims.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["claim_id", "claim", "value", "source", "status"]); w.writeheader(); w.writerows(rows)
print(len(rows), "claims")
