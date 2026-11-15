# PERSEUS multistate Track C (trkC_gain): CBM gated, CEM and LANDIS-II for Georgia and Indiana, pre-registration

Filed 14 November 2026 on firebreather, folder ~/perseus_me/tracks/trkC_gain/{tools, work, restricted (700), logs, launchers, release},
before any new GA or IN value of CBM gated, CEM or LANDIS-II P2 exists. Request: Aaron Weiskittel via the Track C brief (14 Nov 2026).
Basis: the multistate report card (tracks/multistate/tracks/multistate/PREREG.md, md5 at filing ac8d34bbcdc13ce2ce63cba1fc424593,
Amendments 1 to 6; RELEASE_MS.md round 5). Frozen scorers used unchanged: ms_score.py 25bf7b515c067c4ba13debe717b21b0a,
ms_score_x.py f775cf4710c06f1673489d53424a7675. Later changes are dated entries in DEVIATIONS.md beside this file.
Never touched: fvs-conus code, cells and tracks/fvsconus_audit; PID 21329 and its outputs; LANDIS dotnet and FVSne processes;
Cardinal array 16715636; formal_*, s49/cem_formal, *_restricted_quarantine*, existing release folders, multistate/runs (read only).
No GitHub push, email or Zenodo. FIA true coordinates never appear in any output; releases are aggregates only.

## 0. Known at filing (disclosure)

All released multistate GA and IN values (RELEASE_MS.md), verify_20261114/findings.md, refine_curves2 (RELEASE r4: X1 CBM gated
74.02 B on the Maine a3 basis) and round4_fixes were read. A 3-stand GA smoke of the growth-point runner (section 2.1) was run to
test plumbing (restricted/smoke_pts; no curve was fitted from it). The multistate lane is running CEM GA C2 common pairs
(ms_cem_round3.sh, not this track); its outputs are read when complete.

## 1. Gates (apply to every new headline)

G1 reproduction: before any new score, the frozen ms_score_x.py is rerun for IN under this track's overlay HOME (symlinks only;
tools/mk_overlay.sh) and every released IN value (work/final/score_IN/ms_scores.csv) must match within 1e-9 (work/g1_IN.json).
G2 row shuffle: each new scored row is rescored with its unit rows permuted (seed 20261114); MPI and pillar points must match
within 1e-9. G3 independent verifier: a fresh-context subagent recomputes each new headline (MPI and its pillars, C5 ratios) by a
different code path from the row files. Joins: every join of unit rows is by explicit pair_id or chain_id with asserted one-to-one
keys and asserted counts. C5 three-quantity view: gross, mortality and net model/observed ratios on the C5 treated pairs are
reported beside every net headline. Extension check XX1: the new extension (tools/c_score_x2.py) imports ms_score_x.py unchanged and
must reproduce every value ms_score_x.py writes for the same state within 1e-9 (P(rank 1) columns excepted, as in X1).

## 2. T1 CBM gated GA (SN) and IN (CS)

2.1 Growth points. FVS default: the variant binary (FVSsn 9deedc02, FVScs ae9cc343), no calibration keywords, no harvest, no
regeneration keywords, NOTRIPLE, 10 cycles of 5 yr (50 yr, the NCASI window), on every forested condition of the multistate C1 frame
(run_gate.build_frame on state_roots/ST inventory, current evaluation) [ASSUMPTION: the full C1 frame instead of the NCASI 300-stand
sample]. Runner tools/c_curvepts.py = run_gate.project with two asserted source edits (FVS_TreeList also reads MCuFt, TCuFt; per-year
sums of TPA x MCuFt, TPA x TCuFt returned). Response merch_vol_m3ha = sum(TPA x MCuFt) x 0.0699725 (ft3/ac to m3/ha). Stand age =
FIA STDAGE (run_gate median fill) + offset. Forest type group = build_j.ft_group_gated(FORTYPCD), the gated curve key rule. Run on
OSC Cardinal (batch) in the read-only ms_r4 mirror; launcher and md5 manifest stored in launchers/.
2.2 Yield curves. The registered NCASI builder (jobs/ncasi_yield_20260917/build_sit_yield_from_fvs_DRAFT_2026-09-17.R) copied with
only the input path, the output folder and the region labels changed (region SN = GA, CS = IN; the pooled fallback becomes the two
states pooled): Chapman-Richards on merch volume, b0 bounded in [obs max, 2 x obs max], MIN_STANDS 10, state then region then pooled
fallback, ages 0 to 200, arm label default. Output sit_yield_ST_default.csv in the GCBM schema.
2.3 SIT and curve tables. The registered 2x2 driver (cbm_states/shared/tools/run_2x2_curves_inventory_2026-09-22.sh) cell
FVS x FIAW x NEW only, run on Cardinal against a scratch copy of cbm_states (ROOT redirected; the shared states/ST folders are
not written), pixel_attributes_fia_20261113.csv (A5.4), AIDB cbm_defaults_ST.db. Then tools/complete/cbm_curves.py cbm_gated
(unchanged; second-path check.json) gives curve_table.csv.
2.4 Unit rows. Copies of the multistate curve runners (ms_curves_pairs.py, ms_curves_hindcast.py) with the curve root redirected to
this track and the hindcast arm set to cbm_gated FVS_FIAW_NEW (sed with asserted substitution counts): C2 pairs_all level, C5
treated level and ageq (phi GA 0.0304319, IN 0.0175016, multistate/work/common/phi_states.csv), P2 chains. C1 and P5: a copy of
Track R statewide_curves_state.py with evaluation ids 2501/2503 (the multistate current evaluation), the multistate trend series and
the removal-share line made optional (only the no-harvest series is scored); TreeMap block not run.
2.5 Scoring. Row cbm_gated (arm FVS_FIAW_NEW, member) through tools/c_score_x2.py on the multistate basis (rule set a, common basis,
cut-plot view headline; C5 rule as registered in ns_score B1).
2.6 X1 capped variant (refine_curves2 PREREG section 3 form and section 5 X1, as registered for Maine): annual recursion
gross = exp(c0 + c2 z) s^c3 G_pot, mort = gross u^kappa, net = gross (1 - u^kappa), u = B/Bmax, Bmax = p99 of current C1 plot AGB
stock within the forest type group (training data; groups under 30 plots take the state p99), s = B / P(a + t), G_pot = curve
gross x (1 - phi), positivity-preserving step; theta = (c0, c2, c3, gamma) by Levenberg-Marquardt on stacked net and mortality
residuals (Mg/ha/yr) of the uncut C2 pairs; five folds = the multistate pair folds (folds/pair_folds.csv); C2 rows are out-of-fold;
C5 cut pairs and P2 chains use the production fit (all uncut C2 pairs) [ASSUMPTION; Maine used folds for all units]. Site z: ESI
v6r NA 1 km at the public (fuzzed) DataMart PLOT coordinates, standardized on the training data, missing filled with the training
median of the group [ASSUMPTION: GA/IN true coordinates are not used; precedent: climate track, 1 km grids]. C1 and P5 of X1: not
evaluated (inherited would mix forms). Row cbm_gated_x1, scored as 2.5. Gates M1 (no negative stock) and M2 (year 91 to 100
increment below year 1 to 10) on a 100-yr projection of the C1 plots, reported.

## 3. T2 CEM GA and IN

Registered path on the multistate basis = Track R shared CEM build (tools/ncasi_rc/eng/cem_pairs_state.sh, strict donors, same-
location exclusion, treated set with hook 32 and a region treated pool, then tools/v23_E/cem_pair_extract_v3.R and eng/cem_a8.py),
as for the six states. The Maine round4 refinement (cem_r4_refine.R) belongs to the Maine a3 path, not the multistate basis, and is
not applied [DECISION, listed for Aaron]. Ports: (a) treated pools SN = GA, CS = IN by tools/multistate/ms_prep_cem_treated.py with
its NameError fixed (a comment swallowed the O = Path(...) statement; copy in tools/, no logic change); (b) treated runs through
cem_pairs_state.sh with the multistate input and rds paths (as ms_launch_cem_GA_IN.sh) and CEM_H_DIR pointing at the region pool;
(c) GRM condition tables (eng/build_grm_cond.R) for the donor states; (d) extractor with donors GA,FL,SC,NC,TN,AL (GA) and IN (IN),
glob ST_*_pair_%d. C2 common pairs come from the multistate lane's run (read only). C1, P2, P5 CEM items: not in this filing
(reported not evaluated). One CEM process at a time on firebreather, nice 10, CEM_CORES 1.

## 4. T3 LANDIS-II Georgia

Existing: C2 pair rows (D-MS6). P1 needs a statewide C1 extract from a registered statewide run, which does not exist for GA:
not evaluated. P2: hindcast chains by the Track R rule (eng/landis_h1_state.sh/.py: sites from visit 1 live trees, 30 yr registered
no-climate run, B_t/B_0 applied to observed visit 1 stock), GA registered table md5 31490d0b, built on firebreather and run on
Cardinal in the registered image (as A5.1). C5: not evaluated (residual-direct build, A6b check cannot run). P5: no LANDIS P5
reader in Track R: not evaluated. If the P2 runs finish, LANDIS GA is scored as a tier 2 row; otherwise reported running.

## 5. T4 LANDIS-II Indiana

Inventory only: list each missing input (species parameters, initial communities, ecoregions, climate, biomass succession
calibration) with the file evidence (path listing or its absence) on firebreather and Cardinal.

## 6. Reporting

Release tables in release/ (aggregates only), claims.csv with every number and its source file, leak scan of release/.
Grades and tiers as the multistate scorer (A >= 85, B >= 70, C >= 55, D >= 40; tier 2 below coverage 0.65).
