#!/usr/bin/env python3
"""trkC_gain: release/TABLES.md generated from release/trkC_before_after.csv and trkC_c5_three_ratios.csv (aggregates only)."""
from pathlib import Path
import pandas as pd, numpy as np
T = Path.home() / "perseus_me/tracks/trkC_gain"; R = T / "release"
X = pd.read_csv(R / "trkC_before_after.csv"); C = pd.read_csv(R / "trkC_c5_three_ratios.csv")
f = lambda r, p: "" if not np.isfinite(r.get(f"{p}_MPI", np.nan)) else f"{r[p+'_grade']}{'*' if r[p+'_tier'] == 2 else ''} {r[p+'_MPI']:.1f} [{r[p+'_MPI_lo95']:.1f}, {r[p+'_MPI_hi95']:.1f}]"  # noqa: E731
pp = lambda r, p: ", ".join(f"{k} {r[p+'_'+k]:.1f}" for k in ("P1", "P2", "P3", "P4", "P5") if np.isfinite(r.get(f"{p}_{k}", np.nan)))  # noqa: E731
L = ["| state | row | version | before (released r5) | after (this track) | after, c5_three view | pillars after |", "|---|---|---|---|---|---|---|"]
for _, r in X.iterrows():
    a = f(r, "after"); b = f(r, "before")
    if not a and not b:
        continue
    L.append(f"| {r.state} | {r.row} | {r.version} | {b or 'not run'} | {a or 'not run'} | {f(r, 'after_c5three')} | {pp(r, 'after')} |")
M = ["| state | row | n treated | gross | mortality | net | registered n | reg. gross | reg. mortality | reg. net |", "|---|---|---|---|---|---|---|---|---|---|"]
for _, r in C.iterrows():
    M.append(f"| {r.state} | {r.row} | {r.all_n} | {r.all_gross:.3f} | {r.all_mortality:.3f} | {r.all_net:.3f} | {r.registered_n} | {r.registered_gross:.3f} | {r.registered_mortality:.3f} | {r.registered_net:.3f} |")
(R / "TABLES.md").write_text("## Scorecard rows, rule set a, common basis, c5_net view (headline), MPI [95% bootstrap interval]; * tier 2\n\n" + "\n".join(L) +
    "\n\n## C5 three-quantity view: model / observed sums on treated pairs (CBM gated rows x (1 - phi))\n\n" + "\n".join(M) + "\n")
print((R / "TABLES.md").read_text())
