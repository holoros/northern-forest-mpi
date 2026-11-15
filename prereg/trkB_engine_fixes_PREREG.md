# PERSEUS-ME Track B (trkB_engfix): structural engine fixes, pre-registration

Filed 14 November 2026 on firebreather before any new engine run, fit or score of this track. Inputs read only: tracks/
refine_engines, refine_engines2, refine_curves, refine_curves2, round4_fixes, imf2, verify_20261113*, verify_20261114, the FVS
source of build 9c99ff1 (vls/morts.f90, ne/varmrt.f90) and two keyword probes (restricted/probe*, no score). Writes only under
tracks/trkB_engfix/ (restricted/ mode 700). Engine runs go to OSC Cardinal batch jobs; firebreather does file prep, fits and
scoring. No push, email or deposit. Later changes are dated entries in DEVIATIONS.md citing this file's md5. Gates batched.

## 0. Diagnosis that motivates task 1 (from source, before any run)
In the FVS-NE mortality routine (vls/morts.f90 of build 9c99ff1) MORTMULT is applied only when the stand is below the lower
self-thinning bound: when T > PMSDIL x TMD the rate is the density rate RN and the code sets X = 1 (line "IF(RIP .EQ. RN) X=1.0").
So in stands above 55 percent of maximum SDI the MORTMULT keyword has no effect, which explains why uncut mortality stayed near
0.73 of observed at the largest MORTMULT (refine_engines2). The fix must replace the mortality submodel, not scale it.

## 1. FVS-NE structural mortality (primary analysis, label "ddm")
Model. Annual tree mortality probability p: logit(p) = a_k + b RD, k = DBH class at the start of the cycle (inches) [0,3), [3,5),
[5,7), [7,9), [9,12), [12,15), [15,999); RD = BSDI / BSDIMAX, the stand SDI and maximum SDI the FVS Event Monitor reports at the
start of the cycle. Eight parameters, no hand-set constant.
Data. Training trees: live t1 trees (trees_built.csv) of C2 uncut-fit pairs (lanes L1, L3, not the cut overlap; re_lib.frames)
in the training counties of each county-blocked fold (re_lib.county_folds, K = 5, the folds of refine_engines2). Fate from FIA
ME_TREE: the t2 record with PREV_TRE_CN = t1 tree CN; STATUSCD 2 = died, 1 = survived; STATUSCD 3 (cut), 0 or no t2 record are
excluded. RD per pair = FVS BSDI/BSDIMAX at cycle 1 from a COMPUTE record in the default (no other keyword) pair run.
Likelihood: P(died) = 1 - (1 - p)^REMPER, Bernoulli, weights TPA_UNADJ. Fit by maximum likelihood (scipy BFGS) per fold and on
all uncut-fit pairs ("all").
Engine implementation. FIXMORT in replace mode (field 5 = 0), all species (field 1 = 0), one record per DBH class, rate given by
Event Monitor PARMS: 1 - (1 + EXP(a_k + b*BSDI/BSDIMAX))**(YEAR - CENDYEAR - 1), i.e. 1 - (1 - p)^cycle length, scheduled at every
cycle (cycle-number dates). Because PARMS is evaluated at the start of each cycle from the current stand, mortality responds to
the projected density in later cycles and feeds back into growth. FIXMORT runs after the native density and BAMAX logic, so the
fitted model replaces the native mortality for every tree record. No MORTMULT, SDIMAX, post-cut multiplier or cut-cell term.
Fold use. C2 and C5 pairs take the parameters of the fold that holds their county (test fold), so no pair is scored with
parameters fitted on its own county. Hindcast chains take the fold of the visit 1 county (county not in the fold map: "all", flagged).
Bakuzis runs use "all" (no observed outcome is scored in P3).

## 2. FVS-NE pillars rerun
Binary: firebreather build 9c99ff1 FVSne (md5 ba9d0fa463c3bcc8d2e7e57d84310965), copied to Cardinal with an md5 check. Each pillar
is run twice on Cardinal, default_fb (no added keyword) and ddm, so every contrast is within one binary.
- C2 (3,636 holdout pairs, P*, by domain) and C5 net (1,287 partly cut pairs, statewide): one-cycle pair runs (pairs.py run_pair
  arithmetic, unchanged), scored by the registered comp_S arithmetic on the shared draws (refine_engines2 e_score.py logic).
- P3: Bakuzis 100-year no-harvest runs of refine_engines2 bk_fvs.py (sampled plots, same frame), scored by p3_score.py (score_v2.R).
- C3 (P2): H1 hindcast chains (engines/fvs/run_hindcast.py run_chain, unchanged) with the ddm block, scored by the score.py C3
  arithmetic (h1_frame, domain design bootstrap with the registered seed stream). Condition G1-C3: my scorer applied to the
  registered values (me_state/hindcast/restricted/values/fvs__default.parquet) reproduces the W a3 C3 S point within 1e-6 and its
  draws; if it does not, C3 is reported as a within-track contrast (default_fb vs ddm) and the MPI keeps the W a3 C3 (stated).
- C1 and C4bfm (P5): rerun only if the gate stand runner and the y50 bau runner can be driven with the ddm block inside this
  session's compute; otherwise INHERITED from W a3 and the inherited MPI share is reported.
MPI = 0.30 P1 + 0.20 P2 + 0.20 P3 + 0.15 C5 + 0.15 P5, P1 = pts((0.35 S_C1 + 0.25 S_C2)/0.60), pts(S) = 100 max(0, 1 - (e^S - 1)/0.5),
grades A >= 85, B >= 70, C >= 55, D >= 40 (Track W a3 basis, cut-plot view), 1,000 shared draws seed 20261103.

## 3. Gates (hard; a failed gate is reported, nothing is tuned to pass it)
G1 default_fb pair rerun on Cardinal equals the registered default pair values (statewide sums, max abs pair difference) and
reproduces MPI 43.0; G1-P3 default_fb P3 equals refine_engines2 default_fb 94.5. G2 row shuffle (seed 7) of the C2 and C5 frames
leaves every headline unchanged. G3 a fresh-context verifier recomputes the headlines by a different code path. Joins of draws and
predictions to frames are by explicit pair id with asserted length, uniqueness and order (frames are aligned to the draw matrix order).
F16: parameters come only from the likelihood on training fates. F19 realism: ddm P3 >= default_fb P3 - 5. F5: statewide model/obs
ratios inside 0.5 to 2.0. Improvement: paired difference vs default with 95% lower limit > 0. Every paired difference is decomposed
by pillar (weight x points change).

## 4. CEM (task 2)
First verify in code whether treated donors still mix removals with mortality. If removals are already separated (GRM CUT and
DIVERSION components vs MORTALITY, TPAREMV vs TPAMORT, residual donor time 1), no CEM rerun is made, because the requested change
would be a no-op; the C5 net gap is then decomposed from the round4 pair values (gross, mortality, net by removal class) to say
where it sits. Any rerun found necessary is stated as a decision, not made after seeing scores.

## 5. OSM and OSM mortcal (task 3)
Check whether the OSM 2.26.1 binary exposes an in-engine hook for the refinement (diameter increment and mortality). If only a
diameter increment hook exists or none is documented, no engine rerun is made and the OSM rows stay labeled "post-run adjustment".

## 6. Alternatives (stress, conditional on compute)
S1 ddm with the "all" parameters on C2 and C5 (in sample) to size the fold effect. S2 DBH as ln(D) and D instead of classes (not run
unless S1 shows class edges matter). Neither replaces the primary.
