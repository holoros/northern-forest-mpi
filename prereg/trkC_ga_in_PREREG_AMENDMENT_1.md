# trkC_gain PREREG Amendment 1: whole-plot curve basis (filed 14 November 2026, before any whole-plot value exists)

Base PREREG.md md5 38940bf8ed7c463b7f0857f1eebadddc. Decision: Aaron Weiskittel, 14 November 2026 (relayed by the coordinator),
responding to release/DECISIONS.md. Known at filing: every score in release/ (default curves on the CONDPROP >= 0.5 rule, X1, the
cp100 and cpw stress variants, LANDIS GA P2). No whole-plot growth point, curve or score has been computed.

A1.1 Rule (replaces D-C2 as primary). The CBM gated growth points are built on a whole-plot basis, Aaron's standing US Forest Service
FIA convention: the subplot fixed area is the sample; no split or filter by condition. One FVS stand per forested plot of the C1 frame
(run_gate.build_frame, current evaluation): the tree list is every live tree of every forested condition of the plot, each at its
plot-level per acre expansion w = TPA_UNADJ x subplot/macroplot/microplot adjustment (run_gate's w; no division by CONDPROP_UNADJ),
so area outside forest conditions counts as unstocked. Stand attributes (site index, stand age, forest type, slope, aspect, elevation)
come from the plot's largest forested condition (CONDPROP_UNADJ, then CONDID). Forest type group of the plot = build_j.ft_group_gated of
that condition's FORTYPCD. Everything else as PREREG 2.1 (FVS default, 10 x 5 yr, no harvest, no regeneration, MCuFt response).

A1.2 Why. The D-C2 rule (CONDPROP_UNADJ >= 0.5) was chosen after seeing the first fit; the condition basis inflated stand density on
small conditions (TPA_UNADJ / CONDPROP). The whole-plot basis removes that inflation by construction without a post hoc threshold.

A1.3 Downstream, unchanged in form: registered NCASI builder (no CONDPROP filter), ECOREGION "NA" (D-C1), 2x2 FVS x FIAW x NEW SIT cell
on firebreather (D-C10 pipeline), libcbm curve tables, Track R unit rows (C2, C5 level and ageq, P2), statewide C1/P5 (evaluations
2501/2503), X1 refit (PREREG 2.6, D-C3), scoring by tools/c_score_x2.py with G2 (new rows' unit files permuted), XX1 on files, second
code path (tools/c_verify.py), MPI and pillar bootstrap draws saved per row (draws.npz + draws_manifest.json; decision 4).

A1.4 Reporting. Whole-plot rows are the primary CBM gated and X1 values for GA and IN once scored. The CONDPROP >= 0.5, = 1 and weighted
rows stay as stress variants beside them; the range over all four is reported. Default rows are the cells; X1 is a labelled
refinement only (decision 1). CEM GA/IN is scored on P4 only (tier 2) when its lane ends, with the stalled multistate C2 runs recorded
(decision 3). Compute: firebreather only, one worker, $0 at OSC.
