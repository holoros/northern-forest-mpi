"""Figure 2: model/observed ratios of gross growth, mortality and net growth on
uncut (C2) and cut (C5) remeasured FIA pairs, by state and model, 95% intervals
from the scorer bootstrap draws (1,000 draws, seed 20261103).
Source: data/fig2_components.csv (md5 9be2979b53bd1198dbffa43d949f6c13)."""
import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

import os
SRC = sys.argv[1] if len(sys.argv) > 1 else "data/fig2_components.csv"
OUT = (sys.argv[2] if len(sys.argv) > 2 else "figures/out") + "/"
os.makedirs(OUT, exist_ok=True)
INK = "#1F2A24"; MUTE = "#5E6862"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8, "text.color": INK,
                     "axes.labelcolor": INK, "xtick.color": INK, "ytick.color": INK,
                     "axes.edgecolor": "#8A928C"})
MM = 1 / 25.4

d = pd.read_csv(SRC)
d = d[d.in_headline_scorecard].copy()
assert not d.duplicated(["state", "model", "component", "quantity"]).any()

states = ["ME", "NH", "VT", "MI", "WI", "MN"]
scol = {"ME": "#1B4F72", "NH": "#2E86C1", "VT": "#48C9B0",
        "MI": "#B9770E", "WI": "#E67E22", "MN": "#C0392B"}
models = [("fvs", "FVS default"), ("osm", "OSM"), ("acd141", "AcadianGY 14.1"), ("cem", "CEM"),
          ("landis", "LANDIS-II"), ("cbm_gated", "CBM-CFS3, FVS curves"), ("cbm_af", "CBM-CFS3, Gompertz curves"),
          ("fia_curves", "FIA yield curves"), ("E_state", "Ensemble")]
comps = [("C2_uncut", "Uncut pairs"), ("C5_cut", "Cut pairs")]
qs = [("gross_growth", "Gross growth"), ("mortality", "Mortality"), ("net_growth", "Net growth")]
off = dict(zip(states, np.linspace(-0.33, 0.33, 6)))

fig, axes = plt.subplots(2, 3, figsize=(190 * MM, 170 * MM), dpi=300, sharey=True, sharex=True)
ny = len(models)
for r, (c, cl) in enumerate(comps):
    for k, (q, ql) in enumerate(qs):
        ax = axes[r, k]
        for j in range(ny):
            if j % 2 == 0:
                ax.axhspan(ny - 1 - j - 0.5, ny - 1 - j + 0.5, color="#F3F1EC", lw=0, zorder=0)
        ax.axvline(1, color=INK, lw=0.8, zorder=1)
        sub = d[(d.component == c) & (d.quantity == q)]
        for j, (m, _) in enumerate(models):
            for s in states:
                row = sub[(sub.model == m) & (sub.state == s)]
                if row.empty or not np.isfinite(row.ratio_model_over_obs.iloc[0]):
                    continue
                row = row.iloc[0]
                y = ny - 1 - j - off[s]
                hollow = str(row.basis).startswith("NCASI")
                ax.plot([row.lo95_scorer_draws, row.hi95_scorer_draws], [y, y], color=scol[s], lw=1.0, zorder=2)
                ax.plot(row.ratio_model_over_obs, y, "o", ms=3.6, mec=scol[s], mew=0.9,
                        mfc=scol[s], zorder=3)
        ax.set_xscale("log")
        ax.set_xlim(0.05, 4)
        ax.set_xticks([0.1, 0.25, 0.5, 1, 2, 4])
        ax.set_xticklabels(["0.1", "0.25", "0.5", "1", "2", "4"])
        ax.set_ylim(-0.6, ny - 0.4)
        ax.set_yticks(range(ny)); ax.set_yticklabels([m[1] for m in models][::-1])
        ax.tick_params(length=2)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        letter = "abcdef"[r * 3 + k]
        ax.set_title(f"({letter}) {ql}, {cl.lower()}", fontsize=7.5, loc="left", fontweight="bold")
        ax.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
h = [Line2D([], [], marker="o", ls="-", lw=1, color=scol[s], ms=4, label=s) for s in states]
fig.legend(handles=h, loc="lower center", ncol=6, frameon=False, fontsize=7.5, bbox_to_anchor=(0.5, 0.0))
fig.supxlabel("Ratio of predicted to observed rate (log scale)", fontsize=8.5, y=0.045)
fig.tight_layout(rect=(0, 0.055, 1, 1), w_pad=1.6)
fig.savefig(OUT + "Figure_2_components.png"); fig.savefig(OUT + "Figure_2_components.pdf")
print("ok", len(d))
