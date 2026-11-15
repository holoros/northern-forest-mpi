# SOURCES: NCASI Northern Forest manuscript data pack

Every output file in this folder, the script that wrote it (with md5) and its input files (with md5). Paths relative to the home folder on firebreather. Restricted inputs were read server side only; every output is a state-level aggregate. No plot coordinates, plot CNs or PLT_CN are written.

## ME_basis_crosswalk.csv

- output md5: 0b202641f6fefccdd8707e7578c796e3
- script: perseus_me/tracks/ncasi_ms/tools/s3_copies.py (md5 6700ac00904ea953bebcd295a20ffd97)
- written: 2026-11-14T14:34:37
- notes: Verbatim copy of the multistate release crosswalk (multistate vs Maine a3 basis).
- inputs:
  - perseus_me/tracks/multistate/tracks/multistate/ME_basis_crosswalk.csv (md5 0b202641f6fefccdd8707e7578c796e3)

## fig1_scorecard.csv

- output md5: fdde7ce298cec5cdb31b180086d98cba
- script: perseus_me/tracks/ncasi_ms/tools/s1_scorecard.py (md5 79984cd2fc1bb623f3e931c140ab8eee)
- written: 2026-11-14T14:24:30
- notes: Rows: every row of r4_scores ms_scores.csv with rule_set a, basis common basis, view c5_net. in_headline_scorecard marks rows in multistate_scorecard.csv; MPI, CI and grade asserted equal to the scorecard for scored rows. Pillar points 0-100; P4 = C5 cut-pair net view.
- inputs:
  - perseus_me/tracks/multistate/tracks/multistate/multistate_scorecard.csv (md5 7c5855bb6f393c9f63a3941e09b78cc6)
  - perseus_me/multistate/work/r4_scores/score_ME/ms_scores.csv (md5 4b01a1073e08a0726bec01567e4700cf)
  - perseus_me/multistate/work/r4_scores/score_NH/ms_scores.csv (md5 cf22285ce54ec700e034e16c875c8ec8)
  - perseus_me/multistate/work/r4_scores/score_VT/ms_scores.csv (md5 a69ab3020ad92ddc02bc296458e6a2d0)
  - perseus_me/multistate/work/r4_scores/score_MI/ms_scores.csv (md5 7e1efdb1172c1139eb70d3f428e86f57)
  - perseus_me/multistate/work/r4_scores/score_WI/ms_scores.csv (md5 c96e07e8ab69bd9a58889d1e01d0e225)
  - perseus_me/multistate/work/r4_scores/score_MN/ms_scores.csv (md5 214568067f0381b0b34c149465595fca)

## fig2_components.csv

- output md5: 9be2979b53bd1198dbffa43d949f6c13
- script: perseus_me/tracks/ncasi_ms/tools/s6_fig2_combine.py (md5 262839fae83cd8c6f14769b414ae188e)
- written: 2026-11-15T02:27:04
- notes: Primary ratio and 95% interval: the frozen scorer code path (ns_score/rc_adapt.load_state bundle, RR.comps_single and RR.ensemble_comps, rule set a, common basis = statewide, B1 options), statewide e = ln(mean model / mean obs) on units where both exist; interval = percentile 2.5/97.5 of exp(e) over the scorer's 1,000 shared location-multiplicity draws (seed 20261103). The draws are not stored by the scorer; they were regenerated deterministically by its code. C2 = uncut common (gate) pairs; used_in_multistate_MPI: C2 rows where the row has P1 (C2 enters P1 only with C1), C5 net rows where the row has P4 (the c5_net view scores net only). C5 = cut (treated) pairs, all treated pairs (treated_subset all), three quantities where the engine reports them (c5_three quantities), fia_curves net only. LANDIS C5 masked where the r4 landis row has P4 not evaluated (B25 slot written later), reproducing the multistate basis; r4 pillars P1, P2, P4, P5 reproduced for every row (fig2_components_checks.csv). Second path for single engines: explicit pair_id join (validate one_to_one) of observed and model arrays, ratio of means, pair bootstrap 1,000 draws seed 20261114 (columns *_pairboot); point estimates asserted equal to the scorer path (<=1e-9). Shuffle: every unit array permuted consistently (run_track_r.shuffle_bundle, seed 20261114); scorer-path ratios and intervals identical (max diff in fig2_components_checks.csv). Per-state rows written by tools/s2_components.py (VT, NH, ME with the v1 snapshot s2_components_v1_VT_NH_ME.py; MI, WI, MN with FAST_P5=1, which skips only the P5 trend-series targets and leaves P5 out of the r4 check).
- inputs:
  - perseus_me/tracks/ncasi_ms/data/_work/comp_ME.csv (md5 05a68caf294900b66aa88944e70ab5e6)
  - perseus_me/tracks/ncasi_ms/data/_work/comp_ME_meta.json (md5 934ca771f5bdaf515fc092d6aa109725)
  - perseus_me/tracks/ncasi_ms/data/_work/comp_NH.csv (md5 c6c5ba7c313f8b2b14656ec8f64e5b8c)
  - perseus_me/tracks/ncasi_ms/data/_work/comp_NH_meta.json (md5 b3d0afcbd70b4c2830a6f73af94ea17e)
  - perseus_me/tracks/ncasi_ms/data/_work/comp_VT.csv (md5 60d9834eb61bd9399cb0ee98ab1e8324)
  - perseus_me/tracks/ncasi_ms/data/_work/comp_VT_meta.json (md5 c82a3f1e3bd89029d1c973470d4ca99e)
  - perseus_me/tracks/ncasi_ms/data/_work/comp_MI.csv (md5 337aabdd2c6eeddc71702190ae0ac23f)
  - perseus_me/tracks/ncasi_ms/data/_work/comp_MI_meta.json (md5 70a6b616ff66e0e1b749d833d7791231)
  - perseus_me/tracks/ncasi_ms/data/_work/comp_WI.csv (md5 32d4201e3e0eee46ff6696af419f339c)
  - perseus_me/tracks/ncasi_ms/data/_work/comp_WI_meta.json (md5 f0ce116f1aafd884038b426ec608b697)
  - perseus_me/tracks/ncasi_ms/data/_work/comp_MN.csv (md5 ef9e7be0499d0caead479b0f8e1343b0)
  - perseus_me/tracks/ncasi_ms/data/_work/comp_MN_meta.json (md5 85742f68db6f06e15d4be5dd5cad9a5c)
  - perseus_me/tracks/ncasi_ms/tools/s2_components_v1_VT_NH_ME.py (md5 cd5ae524ccb71e843c9add23196e9321)
  - perseus_me/tracks/ncasi_ms/tools/s2_components.py (md5 d2d3818276876f9cf2480a4e8bca23da)
  - perseus_me/tracks/ncasi_ms/data/fig1_scorecard.csv (md5 fdde7ce298cec5cdb31b180086d98cba)

## fig2_components_checks.csv

- output md5: a172619a2fb01c40e2f2d1f4838e8ff9
- script: perseus_me/tracks/ncasi_ms/tools/s6_fig2_combine.py (md5 262839fae83cd8c6f14769b414ae188e)
- written: 2026-11-15T02:27:04
- notes: Per-state check values written by s2_components.py (shuffle max diffs, path agreement, unit counts, wall time).
- inputs:
  - perseus_me/tracks/ncasi_ms/data/_work/comp_ME.csv (md5 05a68caf294900b66aa88944e70ab5e6)
  - perseus_me/tracks/ncasi_ms/data/_work/comp_ME_meta.json (md5 934ca771f5bdaf515fc092d6aa109725)
  - perseus_me/tracks/ncasi_ms/data/_work/comp_NH.csv (md5 c6c5ba7c313f8b2b14656ec8f64e5b8c)
  - perseus_me/tracks/ncasi_ms/data/_work/comp_NH_meta.json (md5 b3d0afcbd70b4c2830a6f73af94ea17e)
  - perseus_me/tracks/ncasi_ms/data/_work/comp_VT.csv (md5 60d9834eb61bd9399cb0ee98ab1e8324)
  - perseus_me/tracks/ncasi_ms/data/_work/comp_VT_meta.json (md5 c82a3f1e3bd89029d1c973470d4ca99e)
  - perseus_me/tracks/ncasi_ms/data/_work/comp_MI.csv (md5 337aabdd2c6eeddc71702190ae0ac23f)
  - perseus_me/tracks/ncasi_ms/data/_work/comp_MI_meta.json (md5 70a6b616ff66e0e1b749d833d7791231)
  - perseus_me/tracks/ncasi_ms/data/_work/comp_WI.csv (md5 32d4201e3e0eee46ff6696af419f339c)
  - perseus_me/tracks/ncasi_ms/data/_work/comp_WI_meta.json (md5 f0ce116f1aafd884038b426ec608b697)
  - perseus_me/tracks/ncasi_ms/data/_work/comp_MN.csv (md5 ef9e7be0499d0caead479b0f8e1343b0)
  - perseus_me/tracks/ncasi_ms/data/_work/comp_MN_meta.json (md5 85742f68db6f06e15d4be5dd5cad9a5c)

## fig3_treemap_direct.csv

- output md5: 68f310d90f653b61fae489b880b37cd1
- script: perseus_me/tracks/ncasi_ms/tools/s3b_treemap_direct.py (md5 cd50b4d951358102d3035beef8ce6394)
- written: 2026-11-14T14:41:43
- notes: Headline rows: 2024 evaluation (EVALID ST2401), variant pixel_all; these reproduce the Oct 26 draft values for all six states. Source script jobs/ncasi_fig4direct_20261027/fig4_direct.R lists its input as Cardinal s49/init_direct/init_direct_by_state.csv md5 d974947ba1a8f0c98102ca8c0603c04b (not re-fetched; Cardinal not used). fia_* are FIA design-based estimates with SE; treemap_lo/hi are TreeMap bootstrap 95% CI as stored. Units: Tg C, Mg C/ha, ha. Area/density split columns as stored.
- inputs:
  - jobs/ncasi_fig4direct_20261027/2026-10-27_init-direct-by-state_DATA.csv (md5 d974947ba1a8f0c98102ca8c0603c04b)

## fig3_treemap_direct_all_variants.csv

- output md5: f5d847fb3a7dd64e5a965e662a27dcbd
- script: perseus_me/tracks/ncasi_ms/tools/s3b_treemap_direct.py (md5 cd50b4d951358102d3035beef8ce6394)
- written: 2026-11-14T14:41:43
- notes: Verbatim copy of all 168 rows (6 states x 4 FIA evaluations x 7 TreeMap variants).
- inputs:
  - jobs/ncasi_fig4direct_20261027/2026-10-27_init-direct-by-state_DATA.csv (md5 d974947ba1a8f0c98102ca8c0603c04b)

## fig3_treemap_phase2.csv

- output md5: cfd9fd43f878788383530bd3cecf8106
- script: perseus_me/tracks/ncasi_ms/tools/s3_copies.py (md5 6700ac00904ea953bebcd295a20ffd97)
- written: 2026-11-14T14:34:37
- notes: Verbatim copy (already clean, aggregate model-pair contrast).
- inputs:
  - jobs/ncasi_p2_20260908/figs_notitle/fia-vs-treemap-divergence_DATA_2026-09-07.csv (md5 cfd9fd43f878788383530bd3cecf8106)

## table1_data_overview.csv

- output md5: 3fb5a81261c99ac826972cb6be98e2be
- script: perseus_me/tracks/ncasi_ms/tools/s4_table1.py (md5 2c00a77fdb9f772f79e396ea4b577c45)
- written: 2026-11-14T22:14:29
- notes: C1 targets: rc_design path B estimator (aboveground live dry biomass DRYBIO_AG, forest land, live DIA>=1 in; GRM on ST2403), points equal targets_six_states.csv (asserted, rel 1e-6); 95% CI = percentile of 1,000 shared location bootstrap draws, seed 20261103 (scorer seed). forest_area_ha from targets_pathA.csv (path A, CONDPROP x adjustment, EXPVOL eval); forest type group shares recomputed on EVALID ST2401 with the rc_adapt._area_weights condition-area formula (acres x 0.40468564224 to ha). Pair sets: pair_folds gate (uncut) and treated (cut) sets. Restricted inputs read server side; only state aggregates written.
- inputs:
  - perseus_me/ext/FIA/ENTIRE_REF_FOREST_TYPE.csv (md5 d279a8987440ab661dc1820f2357233a)
  - perseus_me/ext/FIA/ENTIRE_REF_FOREST_TYPE_GROUP.csv (md5 6bb024997bbe9426532b8cdc73ea4e50)
  - perseus_me/ncasi_rc/inputs/ME/restricted/a1/a1_pairs.parquet (md5 7ae5cce88e64f02413e7a22c4c871d26)
  - perseus_me/ncasi_rc/inputs/ME/restricted/folds/pair_folds.csv (md5 3a40ac4d6ec2e00263052345c2f14147)
  - perseus_me/ncasi_rc/inputs/ME/restricted/pairs_std/pairs_built.csv (md5 71a182b897279e92d3db132efc39f2ad)
  - perseus_me/ncasi_rc/inputs/ME/restricted/treated/treated_pairs.csv (md5 a2a9f077831a6b81332ce4f5c0016108)
  - perseus_me/ncasi_rc/inputs/MI/restricted/a1/a1_pairs.parquet (md5 09482f534dd6d120a8a479c797bab86e)
  - perseus_me/ncasi_rc/inputs/MI/restricted/folds/pair_folds.csv (md5 eb0b5cad891650e94fa47e0a6f9548e8)
  - perseus_me/ncasi_rc/inputs/MI/restricted/pairs_std/pairs_built.csv (md5 c9677cc9b7eae5c1e644488d4b950c15)
  - perseus_me/ncasi_rc/inputs/MI/restricted/treated/treated_pairs.csv (md5 0381440b363f8a01098c38241ff551d5)
  - perseus_me/ncasi_rc/inputs/MN/restricted/a1/a1_pairs.parquet (md5 eeef8441e6c2eaea89696d610029302a)
  - perseus_me/ncasi_rc/inputs/MN/restricted/folds/pair_folds.csv (md5 b673cdd9f322dd17a5886c73733920dc)
  - perseus_me/ncasi_rc/inputs/MN/restricted/pairs_std/pairs_built.csv (md5 f3807c82889d9363bcf89a54b03fe490)
  - perseus_me/ncasi_rc/inputs/MN/restricted/treated/treated_pairs.csv (md5 22455bbfc30ce721d8d1de2e75f8e66f)
  - perseus_me/ncasi_rc/inputs/NH/restricted/a1/a1_pairs.parquet (md5 f34e6e4950c3e256b6666b9984a037f3)
  - perseus_me/ncasi_rc/inputs/NH/restricted/folds/pair_folds.csv (md5 02865fa46df445b64e4034d83bf69966)
  - perseus_me/ncasi_rc/inputs/NH/restricted/pairs_std/pairs_built.csv (md5 5cef6abade566b5aafede9f4b9d7cbe2)
  - perseus_me/ncasi_rc/inputs/NH/restricted/treated/treated_pairs.csv (md5 5b8f8a157e8a7f40a1cd54a4256abe7c)
  - perseus_me/ncasi_rc/inputs/VT/restricted/a1/a1_pairs.parquet (md5 9c106c0f11724ddb2f48f222c5d09323)
  - perseus_me/ncasi_rc/inputs/VT/restricted/folds/pair_folds.csv (md5 fb9d541f5716bef9a04c7e811d799eb4)
  - perseus_me/ncasi_rc/inputs/VT/restricted/pairs_std/pairs_built.csv (md5 b6a54e91f7ad2ff3ea82de1173e469e9)
  - perseus_me/ncasi_rc/inputs/VT/restricted/treated/treated_pairs.csv (md5 f9de4027cb6fc2aeee06aaad1fc733fc)
  - perseus_me/ncasi_rc/inputs/WI/restricted/a1/a1_pairs.parquet (md5 3b0afa82161850d5fd41f7f2b3ec7c64)
  - perseus_me/ncasi_rc/inputs/WI/restricted/folds/pair_folds.csv (md5 a177cfdbb7cce396a0da1a2d9977293c)
  - perseus_me/ncasi_rc/inputs/WI/restricted/pairs_std/pairs_built.csv (md5 1fae6d32b46c32a2a252ae9489b0873d)
  - perseus_me/ncasi_rc/inputs/WI/restricted/treated/treated_pairs.csv (md5 5a418611ea0196e17b2c937ba394d084)
  - perseus_me/ncasi_rc/targets/targets_pathA.csv (md5 d66e67e4aad2027329dbaf8e12be80bb)
  - perseus_me/ncasi_rc/targets/targets_six_states.csv (md5 a82b75f13e409402c848078154c91ee3)

## table2_models.csv

- output md5: de5678b96e7e0d0f3e28b6050428c61d
- script: perseus_me/tracks/ncasi_ms/tools/s5_table2.py (md5 4ccecde60b92e538adc85fca3488c4a2)
- written: 2026-11-14T14:48:43
- notes: Text transcribed from the listed docs (per-row source column); states_scored_multistate and ensemble members derived from fig1_scorecard.csv. Blank cells were not stated in the docs read.
- inputs:
  - perseus_me/ncasi_rc/engines/ENGINE_SETUP.md (md5 27533a64c1c3d753e338b80798530bd3)
  - perseus_me/ncasi_rc/docs_PREREG_R1.md (md5 8d4964a38dc01c5cf73ce16d5ddb747b)
  - perseus_me/tracks/multistate/docs/DESIGN.md (md5 d84c0d41c1b96c344356eb5bde6f409d)
  - perseus_me/tracks/multistate/tracks/multistate/PREREG.md (md5 bf3cffeedbe73f358f154270b0ea7c78)
  - perseus_me/tracks/multistate/tracks/multistate/RELEASE_MS.md (md5 583d77ab580a65356790d20dda6bc345)
  - perseus_me/tracks/multistate/tools/ncasi_rc/score/rc_core.py (md5 aaa3d9048af4a2c30cbd60e0cb8ea771)
  - perseus_me/tracks/multistate/tools/ncasi_rc/score/rc_run.py (md5 3e2b5bb726f2857d4428b7c1ea1e7d74)
  - perseus_me/tracks/multistate/docs/PERSEUS_ME_ENGINES_ADDENDUM.md (md5 d0a694b1120fc7e12c8b1cf7981ab3b4)
  - perseus_me/tracks/ncasi_ms/data/fig1_scorecard.csv (md5 fdde7ce298cec5cdb31b180086d98cba)

## Helper scripts

- perseus_me/tracks/ncasi_ms/tools/common.py (md5 553bd24797dafbb7910e7b62a494a24e)
- perseus_me/tracks/ncasi_ms/tools/run_queue.sh (md5 1afea438fbf2271fcc5f9847c399f8bb)
