import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import sys
from matplotlib.colors import LinearSegmentedColormap, BoundaryNorm
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

INK = "#1F2A24"; MUTE = "#5E6862"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "text.color": INK,
                     "axes.labelcolor": INK, "xtick.color": INK, "ytick.color": INK,
                     "axes.edgecolor": "#8A928C"})
MM = 1 / 25.4
import os
OUT = (sys.argv[2] if len(sys.argv) > 2 else "figures/out") + "/"
os.makedirs(OUT, exist_ok=True)

# ---------------- Figure 1: scorecard heatmap ----------------
states = ["ME", "NH", "VT", "MI", "WI", "MN"]
models = [("fvs", "FVS\ndefault"), ("osm", "OSM"), ("acd141", "AcadianGY\n14.1"), ("cem", "CEM"),
          ("landis", "LANDIS-II"), ("cbm_gated", "CBM\ngated"), ("cbm_af", "CBM\nAF"),
          ("fia_curves", "FIA yield\ncurves"), ("E_state", "State\nensemble")]
# Scores read from the data pack (fig1_scorecard_v2.csv: LANDIS-II cut-pair results folded in for NH, VT, MI, WI, MN)
import pandas as pd, sys
DATA = sys.argv[1] if len(sys.argv) > 1 else 'data/fig1_scorecard_v2.csv'
_d = pd.read_csv(DATA); _d = _d[_d.in_headline_scorecard == True]
S = {}
for _, _r in _d.iterrows():
    S.setdefault(_r.state, {})[_r.model] = (float(_r.MPI), _r.grade, int(_r.tier))
cmap = LinearSegmentedColormap.from_list("g", ["#B23A2E", "#E08A3C", "#E9C46A", "#8DBF7A", "#2E7D4F"])
norm = BoundaryNorm([0, 40, 55, 70, 85, 100], cmap.N)
fig, ax = plt.subplots(figsize=(190 * MM, 95 * MM), dpi=300)
for i, s in enumerate(states):
    for j, (m, _) in enumerate(models):
        v = S[s].get(m); x, y = j, len(states) - 1 - i
        if v is None:
            ax.add_patch(plt.Rectangle((x - .46, y - .42), .92, .84, color="#ECEAE3", lw=0)); continue
        mpi, gr, tier = v
        ax.add_patch(plt.Rectangle((x - .46, y - .42), .92, .84, color=cmap(norm(mpi)), lw=0,
                                   alpha=1.0 if tier == 1 else 0.55))
        txt = f"{gr}{'*' if tier == 2 else ''}\n{int(mpi + 0.5)}"
        ax.text(x, y, txt, ha="center", va="center", fontsize=8.5, fontweight="bold", linespacing=1.1,
                color="white" if (tier == 1 and (mpi < 40 or mpi >= 70)) else INK)
ax.set_xticks(range(len(models))); ax.set_xticklabels([m[1] for m in models], fontsize=7.5)
ax.xaxis.tick_top()
ax.set_yticks(range(len(states))); ax.set_yticklabels(states[::-1], fontsize=9, fontweight="bold")
ax.set_xlim(-.6, len(models) - .4); ax.set_ylim(-.6, len(states) - .4)
ax.axhline(2.5, color=INK, lw=0.8)
ax.text(len(models) - .35, 4.0, "Northeast", rotation=270, va="center", fontsize=8, color=MUTE)
ax.text(len(models) - .35, 1.0, "Lake States", rotation=270, va="center", fontsize=8, color=MUTE)
for sp in ax.spines.values(): sp.set_visible(False)
ax.tick_params(length=0)
h = [Patch(color=cmap(norm(v)), label=l) for v, l in [(90, "A (85+)"), (77, "B (70 to 85)"), (62, "C (55 to 70)"),
                                                       (47, "D (40 to 55)"), (20, "F (<40)")]]
h.append(Patch(color="#ECEAE3", label="Not run"))
ax.legend(handles=h, loc="upper center", bbox_to_anchor=(0.5, -0.02), ncol=6, frameon=False, fontsize=7.5)
fig.tight_layout()
fig.savefig(OUT + "Figure_1_scorecard.png"); fig.savefig(OUT + "Figure_1_scorecard.pdf"); plt.close(fig)

# ---------------- Figure 3: initialization ----------------
D = [  # state, total, area, density, fia, tm, tm_lo, tm_hi
 ("ME", -16.42, -3.89, -13.04, 374.0, 312.6, 300.1, 326.0),
 ("NH", -17.44, -2.43, -15.39, 153.2, 126.4, 121.3, 131.8),
 ("VT", -15.26, 0.07, -15.32, 142.9, 121.1, 116.4, 126.0),
 ("MI", -15.49, 4.02, -18.75, 512.5, 433.1, 419.8, 447.1),
 ("WI", -2.46, 4.78, -6.91, 396.2, 386.5, 377.2, 396.4),
 ("MN", 1.42, 3.26, -1.78, 308.8, 313.2, 302.9, 324.4)]
P2 = [("ME", 22.8, 28.7), ("NH", -7.2, -14.5), ("VT", -2.6, 0.1), ("MI", -9.4, -21.4),
      ("WI", -3.7, -23.1), ("MN", 11.6, -1.3), ("GA", -44.7, -50.5)]
C_AREA = "#3A8FB7"; C_DENS = "#D9822B"; C_TOT = "#1F2A24"
fig, (a, b) = plt.subplots(1, 2, figsize=(190 * MM, 92 * MM), dpi=300, gridspec_kw={"width_ratios": [1, 1.1]})
y = np.arange(len(D))[::-1]
for k, (st, tot, ar, de, fia, tm, lo, hi) in enumerate(D):
    yy = y[k]
    a.barh(yy + 0.17, ar, height=0.3, color=C_AREA)
    a.barh(yy - 0.17, de, height=0.3, color=C_DENS)
    lo_p, hi_p = (lo / fia - 1) * 100, (hi / fia - 1) * 100
    a.plot([lo_p, hi_p], [yy, yy], color=C_TOT, lw=1.2, zorder=3)
    a.scatter([tot], [yy], marker="D", s=26, color=C_TOT, zorder=4)
    a.text(10.5, yy, f"{tot:+.1f}%", va="center", ha="right", fontsize=7.5, fontweight="bold")
a.axvline(0, color="#8A928C", lw=0.8)
a.set_yticks(y); a.set_yticklabels([d[0] for d in D], fontweight="bold")
a.set_xlim(-26, 12); a.set_xlabel("TreeMap minus FIA, % of FIA stock")
a.set_title("(a) Direct contrast, 2024 inventory", loc="left", fontsize=9, fontweight="bold")
a.legend(handles=[Patch(color=C_AREA, label="Area part"), Patch(color=C_DENS, label="Density part"),
                  Line2D([], [], color=C_TOT, marker="D", ms=5, label="Total, TreeMap 95% CI")],
         loc="upper center", bbox_to_anchor=(0.5, -0.17), ncol=3, frameon=False, fontsize=7)
for sp in ["top", "right"]: a.spines[sp].set_visible(False)
y2 = np.arange(len(P2))[::-1]
c25, c125 = "#9BC3DD", "#2B6A99"
for k, (st, v25, v125) in enumerate(P2):
    b.barh(y2[k] + 0.17, v25, height=0.3, color=c25)
    b.barh(y2[k] - 0.17, v125, height=0.3, color=c125)
b.axvline(0, color="#8A928C", lw=0.8)
b.axhline(0.5, color="#C9CEC9", lw=0.8, ls=":")
b.set_yticks(y2); b.set_yticklabels([p[0] for p in P2], fontweight="bold")
b.get_yticklabels()[-1].set_color(MUTE)
b.set_xlabel("TreeMap start minus FIA start, % of FIA start")
b.set_title("(b) Same yield curve model from each start", loc="left", fontsize=9, fontweight="bold")
b.legend(handles=[Patch(color=c25, label="2025"), Patch(color=c125, label="2125")], loc="upper center",
         bbox_to_anchor=(0.5, -0.17), ncol=2, frameon=False, fontsize=7)
for sp in ["top", "right"]: b.spines[sp].set_visible(False)
fig.tight_layout()
fig.savefig(OUT + "Figure_3_initialization.png"); fig.savefig(OUT + "Figure_3_initialization.pdf"); plt.close(fig)
print("ok")
