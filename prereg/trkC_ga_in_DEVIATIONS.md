# trkC_gain deviations

Dated entries after filing of PREREG.md (md5 in PREREG.md5).

- 2026-11-14 08:30 AST, G1 passed before any new score: ms_score_x.py IN rerun under the overlay HOME reproduced all 144 released IN
  rows, max abs diff 0 (work/g1_IN.json; wall 761 s).
- D-C1 (plumbing, 2026-11-14 08:40, before any score): the GCBM builder (build_gcbm_input_state.py, section 6) matches state-mean
  curves only when ECOREGION_L3_CODE is the string "NA"; the registered six-state inputs in ncasi_cem_yield/sit_yield_in carry "NA",
  while the firebreather copy of the builder writes an empty field (data.table fwrite). The first SIT build therefore resolved every
  GA and IN curve by the state-pool fallback (78 of 78 GA, 55 of 55 IN), one pooled curve for all groups. Fix: the empty field is
  written as "NA" (as in sit_yield_in); first build kept as sit2x2_ST_v0 (not used).
- D-C2 (data, 2026-11-14 08:45, after the first curve fit and before any score): the full-frame growth points include conditions
  with small CONDPROP_UNADJ whose stand density (TPA_UNADJ / CONDPROP, run_gate convention) is inflated; their start volumes reach
  23,221 m3/ha (GA loblolly) and pin the bounded b0 (lower bound = stratum max). Of GA stands above 1,500 m3/ha, 11 of 11 have
  CONDPROP < 0.5 (IN 1 of 1); at CONDPROP >= 0.5 the GA maximum is 1,084 and the p99 745. Rule adopted: growth points from
  conditions with CONDPROP_UNADJ >= 0.5 only. Alternatives not taken: CONDPROP = 1 only; weighting by CONDPROP. Disclosed: this rule
  was chosen after seeing the first fit (coefficients in out/sit_v0 on Cardinal).
- D-C3 (2026-11-14 09:20, before any X1 fit): the multistate pair fold file (inputs/ST/restricted/folds/pair_folds.csv, set gate) has
  two folds (1, 2), not five; X1 C2 rows are out-of-fold over these two folds (PREREG 2.6 said five).
- D-C4 (compute, 2026-11-14 ~09:15, coordinator relaying Aaron): firebreather only, $0 at OSC, at most one worker for this track; no
  further Cardinal jobs. Jobs already run on Cardinal before the rule (16747858 growth points, 16747909 first SIT, 16748064 SIT with
  D-C1/D-C2, 16748373 unit rows) were pulled to restricted/cardinal_pull by tar over ssh with a 142-file md5 manifest (all match;
  launchers/cardinal_PULL_MANIFEST.md5, launchers/cardinal_mcp_job_*.slurm). On firebreather the CEM lane is paused (SIGSTOP of its
  process group) whenever another step of this track runs, so one worker is active at a time.
- D-C5 (plumbing, 2026-11-14 09:20, before any score): the ported statewide script failed on Cardinal (job 16748373) because an
  inline comment added by the port swallowed the rest of two statements (P1 and meta); comments moved; rerun on firebreather.
- D-C6 (stress, coordinator 2026-11-14): the D-C2 rule (CONDPROP_UNADJ >= 0.5) stays primary and logged as chosen after seeing the
  first fit; two stress variants are run and scored like the primary: cp100 (CONDPROP_UNADJ = 1 only) and cpw (all conditions,
  nlsLM weights = CONDPROP_UNADJ; builder edits applied in memory with asserted counts, c_fit_sit.R weighted mode). Before the
  variants, a cross-host check runs the SIT pipeline on firebreather (venv_cbm via engines/pycbm, cbm_states files pulled from
  Cardinal, md5 manifest) on the D-C2 curves and compares the adapted SIT files with the Cardinal ones (work/xhost_sit.json).
- D-C7 (LANDIS GA P2, compute rule): PREREG section 4 said the LANDIS runs go to Cardinal; under D-C4 they run on firebreather with the
  Track R runner ~/opt/landis_v8/landis_run.sh (registered image unpacked, patched DLLs bound), queued after the CBM lane.
- D-C8 (check plumbing, 2026-11-14 09:45, after the first IN score): XX1 compared the in-memory extension table with the released
  ms_scores.csv, which ms_score_x.py writes at %.6g, so it reported max diff 5.0e-5 (rounding only). Recomputed on the written files
  (both %.6g): 144 of 144 rows joined, 0 cells differ by more than 1e-9. c_score_x2.py now compares after the same CSV round trip;
  the IN and GA primary runs and the shuffle runs that loaded the earlier code exit 1 after writing their outputs, and their XX1 is
  re-evaluated from files in work/xx1_files.json.
- D-C9 (G2 design, 2026-11-14 10:15): the first G2 rescoring permuted every CSV/parquet with a pair_id or chain_id column, including
  the frozen frame files (pairs, folds, a1, chains). Common-basis point values were identical for every row (IN), but bootstrap
  intervals and some by-domain values of rows this track did not change (fvsb, E_fvs3; IN by-domain fvsb MPI 68.38 vs 68.61) moved,
  because the bootstrap location order follows the frame file order. That run is kept as restricted/score/shufflewide_IN (finding for
  the multistate owner: interval reproducibility depends on frame file order; the registered N4 test shuffles engine rows only). G2 as
  registered (new rows' unit files permuted, frames untouched) is rerun as restricted/score/shuffle_ST by tools/c_master2.sh.
- D-C10 (2026-11-14 11:15, before any stress score was used): the first firebreather stress lane was invalid. build_gcbm_input_state.py
  failed on firebreather (GCBM template db at a Cardinal path), the pulled root still held the Cardinal D-C2 build, so the "cross-host
  check" passed vacuously and the cp100/cpw SIT files equalled the primary (their scores equalled the primary to print precision).
  Fixes: GCBM template db (md5 fb37f046) and FIA REF_FOREST_TYPE.csv (md5 d279a898) pulled from Cardinal; CBM_MAINE_ROOT and
  FIA_REF_FOREST_TYPE point at them; c_sit_2x2_fb.sh deletes the state gcbm db and bundle before each build and exits on any failed
  step. After the fix a GA rebuild from scratch reproduces the Cardinal adapted growth_and_yield.csv and inventory.csv exactly. The
  invalid run is kept as logs/stress_lane_v1_invalid.log; its outputs were deleted. Stress v2 runs after the LANDIS GA h1 lane.
  Note: the firebreather libcbm curve table from the same SIT differs in md5 from the Cardinal one (cbm_curves.py on two hosts);
  the resulting pair rows differ below print precision of the scores.
- D-C11 (check scope, 2026-11-14 13:00): with the LANDIS h1 slot the released LANDIS GA row (all values missing, NS) becomes scored, so
  XX1 reports a difference by design for that row; excluding the landis row, the run equals the released file (144 rows) and the
  primary run (192 rows), max diff 0 (work/xx1_landis_GA.json).
- D-C12 (2026-11-14 ~14:10): c_score_x2.py now saves every row's bootstrap draws (_bt of ms_score.score_all) to draws.npz with a manifest (decision 4); no scoring logic changed. Amendment 1 (whole-plot basis) filed and frozen before any whole-plot value; c_curvepts_wp.py, c_x1_wp.py and c_wp_lane.sh implement it.
