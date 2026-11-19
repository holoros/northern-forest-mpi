# Sensitivity runs of 15 November 2026 for the NCASI Northern Forest scorecard: registration

Filed before any score of these runs is computed. Requested by Aaron Weiskittel. The registered and released scores
(multistate r4 for Maine, NCASI states releases for NH, VT, MI, WI, MN, and data/fig1_scorecard_v2.csv) are not replaced,
edited or moved by anything here. Everything below is a post hoc sensitivity or an added run and is labelled so.

## Code path common to the three tasks
A driver, sens_20261115/sens_driver.py, imports the track scorer unchanged (tools/ncasi_states/ns_score.py md5 77cdcc95,
which imports tools/ncasi_rc/score rc_core, rc_adapt, rc_run, rc_state, rc_p3, rc_a6, run_track_r) and calls its own
functions in the order of ns_score.main: overlay, load (seed 20261103, 1,000 draws), FVS-ACD arm rows for NH and VT,
rule_sets, attach_p3_all, score_view. The driver differs from ns_score.main only in: (1) the scratch folder is
sens_20261115/restricted/score_work/<ST> (ns_score and ms_score delete and rebuild their own scratch folders, which this
work may not touch); (2) only rule set a, common basis, view c5_net is written (the manuscript headline filter);
(3) the N0 Track R reproduction and the ms_score fvsb row and M1 block are not rerun (they do not enter rule set a rows;
the reproduction check below takes their place); (4) the switches of tasks (b) and (c). Reason for a driver and not a
full copied scorer run: the full scorers took 0.15 to 3.8 hours per state when released.
Maine: the registered Maine rows come from ms_score.py, whose rule set a rows equal the ns_score path within 1e-9 (its M1 check).
Outputs are aggregates only (state by row tables). Plot-level intermediates stay under sens_20261115/restricted (mode 700).

## Task (a) Maine LANDIS-II cut-pair (C5, B25 weight) run and Maine rescoring
Run: cd ~/perseus_me/tracks/ncasi_states; nice -n 10 bash tools/ncasi_states/ns_landis_b25.sh ME
~/perseus_me/ncasi_states/runs/landis_b25/ME/restricted/c5 (the pipeline function lb25, same script md5 6c5bddd5, same
environment variables NS_ROOT and OMP_NUM_THREADS=1; TMPDIR is set to sens_20261115/restricted/tmp instead of the
pipeline tmp folder). The output path does not exist beforehand. No pipeline done marker is written.
Result: a6b.json aggregate fields; then the Maine rows of rule set a, common basis, c5_net from the driver with the
overlay LANDIS-II C5 slot pointing at that run (ns_score.overlay with landis_b25=True, as for the other five states).
Before and after: LANDIS-II and E_state (MPI, 95% interval, grade, tier, coverage, P1 to P5, P_rank1); every other Maine
row must be unchanged except rows that contain LANDIS-II C5 by construction (ensembles with LANDIS-II as a member).
Output file: data/fig1_scorecard_v3.csv = v2 with the Maine landis row and the Maine ensembles that include LANDIS-II
replaced; v1 and v2 untouched. Note: once the ME run folder exists, any later rerun of ns_score or ms_score for Maine will
read it (their overlay points there already); this is disclosed in the report.

## Task (c) leave-one-member-out ensembles
For each state and each member m of the registered E_state lineup (bundle members, rule set a), an added ensemble tuple
("E_state_wo_<m>", members without m, no recalibration, mean, unranked) is appended to rule set a before attach_p3_all, so
the scorer's own ensemble_comps, half rule and ensemble P3 are used. Result: MPI, 95% interval, tier, coverage per
leave-one-out ensemble; and whether it is the highest tier 1 row among the member rows of that state (point MPI), plus
P(rank 1) from the scorer's rank_rows over members and that ensemble. Maine uses the task (a) LANDIS-II C5 run if (a) succeeded.

## Task (b) LANDIS-II statewide totals (C1) treated as in sample
Mechanism: the scorer's own in-sample register, read through its RC_REGISTER override (rc_adapt.load_insample). A copy of
ncasi_rc/engines/IN_SAMPLE_REGISTER.csv is written to sens_20261115 with the call of engine landis, component C1 set to
"in" for the six states (the call fia_curves C1 already has); nothing else changes. Under the registered rule
(rc_state.score_row, rc_core.p1_S) an in-sample C1 is excluded and P1 is then not evaluated ("C1 is required for P1; C2
alone does not make P1"); this is what the fia_curves row gets. The same rule is applied to landis as coded; P1 is not
re-defined as C2 only, because that would change scoring logic. The ensemble C1 drops LANDIS-II through ensemble_comps.
Result: LANDIS-II and E_state MPI, 95% interval, tier, coverage, registered and sensitivity, by state.

## Checks
1. Reproduction: with all switches off the driver must reproduce the registered rows of that state (rule set a, common
   basis, c5_net; ns release for NH, VT, MI, WI, MN; r4 for Maine apart from rows that contain the new LANDIS-II C5) to the
   6 significant digits of the released files. Maximum absolute difference reported. Switched results are not trusted otherwise.
2. Second code path: each new MPI recomputed from its pillar points with weights 0.30, 0.20, 0.20, 0.15, 0.15 renormalised
   over evaluated pillars, by separate code on the written csv. Maximum difference reported.
3. Row shuffle: the scorer hook run_track_r.shuffle_bundle (seed 11), one state per task, maximum MPI change reported.
## Rules
One compute process at a time, nice 10. No existing file is modified. Registered scores are not replaced.
