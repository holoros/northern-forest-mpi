# PERSEUS-ME Track A (trkA_hcm): start-of-cycle Harvest Choice Model refit, county-blocked scoring, IMF2 rerun, NCASI port (pre-registration)

Registered 14 November 2026 on firebreather before any trkA prediction or score existed. Folder ~/perseus_me/tracks/trkA_hcm
(restricted/ mode 700; release/ aggregates only). No git binary on firebreather: the freeze is the md5 in PREREG.md5. Later changes go
to DEVIATIONS.md. Never written: formal_*, s49/cem_formal, *_restricted_quarantine*, tracks/fvsconus_audit, any existing release folder,
tracks/imf2 (tools are copied, not edited). FIA coordinates (true or public) are never read; no LAT/LON column is loaded.
Compute: one worker process on firebreather at nice 10.

## 0. Inputs inspected before registration (no trkA result existed)
- IMF2 release (PREREG md5 b539cc149e8f1e1606df0b1d46ccc065): V1 forecast MPI 80.49 [77.22, 83.24], oracle 86.39; layer SC-L
  cross-fitted AUC 0.676 [0.665, 0.688] on 5 plot-location folds (track H `fold`); HCM refit 0.634. Draw/frame alignment by explicit
  pair id with assertions (IMF2-D5).
- Track HCM (PREREG_HCM md5 b52391e14eb6ef0a27a9caaa72c278dd) refit used t2 MEASYEAR price windows and an expected no-harvest BioTot.
- IMF3 K1 to K3 (MFS county stumpage, Harris-Caputo-Butler ownership, stand value added to SC-L) gave no AUC gain over SC-L
  (all paired intervals include 0). Prior expectation: an HCM-form layer is unlikely to beat SC-L; it is still refitted and scored.
- Zhao, Daigneault and Weiskittel (2020) Forest Policy and Economics 118:102251 (read 14 Nov 2026, Drive id 1ga_jSPYLbYptnW8HXK2HtHxTbAnFOtkY):
  MNL per product (sawlog, pulplog), alternatives none (0 percent removal), partial (1 to 70 percent), full (70 to 100 percent); eq. 10
  covariates PriceSaw, PricePulp (county), LagBio, BioTot / BiopulpLD, PostGrowth and its square, Mills (50 km), LandValue, County,
  HighwayDist, Conservation, Elevation, Year; periods 2002-2006, 2007-2011, 2012-2016.

## 1. Leakage rule (F17)
BioTot and BiopulpLD in the paper are current-period aboveground biomass, which includes the harvest outcome. They are NOT used. Every
predictor is measured at the pair's t1 visit or earlier (prices: years up to t1). A diagnostic variant A2 (section 5) re-adds the observed
t2 biomass to show the size of the leak; it is never eligible for promotion or for IMF2.

## 2. Frame, target, folds (Maine)
- Frame: track H pair frame me_state/v26_H/restricted/frame/pairs_h.parquet (13,750 pairs), target `cut` (H1 rule), intensity
  share_ba on h2_eligible pairs. Product choices exactly as track HCM section 4 (share of t1 product biomass cut; partial if 0 < share < 0.70
  or cut with no identified cut tree and product present; full if >= 0.70).
- County-blocked folds (fold map, fixed here; greedy assignment of counties by descending pair count to the least-loaded fold):
  fold 0: 3 (Aroostook); fold 1: 21, 31, 27, 23; fold 2: 25, 9, 15; fold 3: 19, 7, 11, 13; fold 4: 29, 17, 5, 1
  (pairs 2,994 / 2,723 / 2,681 / 2,708 / 2,644). COUNTYCD is the FIPS county code of pairs_h.
- Every trkA layer prediction is out of fold on these county folds; inner tuning, standardisation, imputation and PostGrowth / intensity
  means use training folds only.

## 3. HCM2 specification (primary)
Per product MNL (reference = none), ridge 1e-4 on standardised slopes, intercept unpenalised (track HCM convention, near MLE).
Design (identical in both product equations except the product stock term):
- ln county stumpage at t1: spruce-fir sawlog $/MBF and hardwood pulpwood $/cord, Maine Forest Service county rows of
  jobs/hcs_stumpage_refit_20261001/stumpage_panel_full.csv at INVYR1, statewide median fallback with indicator (the IMF3 price_at rule,
  copied unchanged).
- LagBio (t1 live aboveground biomass, Mg/ha, track HCM hcm_features lagbio); product stock at t1 (saw_t1 in the sawlog equation,
  pld_t1 in the pulplog equation) in place of BioTot / BiopulpLD.
- PostGrowth and PostGrowth^2: training-fold mean annual product biomass change on uncut pairs by type group x age class (track HCM rule).
- Elevation (t1 ELEV), RDDISTCD (public road distance class; substitute for HighwayDist), t1 period dummies (INVYR1 2007-2011, >= 2012;
  reference <= 2006), ln REMPER.
- Ownership: OWNCD one-hot (levels with < 30 pairs pooled into their OWNGRPCD), RESERVCD. No public conservation easement layer is on
  firebreather (searched PAD-US, conserved lands): FIA ownership only, logged. Mills, LandValue: no public layer on firebreather, dropped.
- County dummies are dropped (a held-out county has no estimable effect); county variation enters through county prices.
- Additions: ln(1 + stand value), stand value = volbf / 1000 x saw price + volcf / 79 x pulp price ($/acre, track SC volumes, 79 ft3 per
  cord as IMF3 K3); prior cut history: prev_trt10, yrs_since_cut (training median fill + missing indicator), prev_pair_cut (1, 0,
  not observable indicator).
- Collapse (track HCM section 4 rule): P(cut) = max(1 - P_none_saw, 1 - P_none_pulp); intensity s = mean over products of
  (P_partial s_part + P_full s_full) / (P_partial + P_full) with training-fold mean observed share_ba below / at or above 0.70
  (multinomial expectation; no Tobit or beta model).

## 4. Comparators (layer scoring)
- SC-L county: track SC sc_logit.py code (Design, choose, fit_logit, run_split, unchanged) with `fold` replaced by the county fold map
  (inner lambda CV over the other county folds). This is the like-for-like comparator.
- SC-L stored (plot folds, as released) and HCM refit stored: reported for reference only.

## 5. Stress test (pre-registered alternatives, all county-fold out of fold)
A1 single MNL on total basal area removal class (none, partial < 0.70, full >= 0.70), P(cut) = 1 - P_none. A2 leak diagnostic: A0 plus
observed t2 product biomass (never eligible). A3 A0 plus Harris-Caputo-Butler hcb_class (IMF3-D7 plot identity join). A4 core ZDW
start-of-cycle terms only (no stand value, no cut history, no OWNCD; OWNGRPCD public indicator as Conservation). A5 HCM2 on the track H
plot-location folds (comparison with the released SC-L 0.676 on its own folds). A6 A5 plus county dummies.

## 6. Layer metrics
On all 13,750 pairs: AUC (Mann-Whitney with ties), Brier, calibration slope and intercept (logistic of y on logit p), E/O.
Intervals: 1,000 plot-location cluster bootstrap draws, track H stream (h_score SEED, NB; as sc_eval.py), predictions fixed.
Paired AUC difference HCM2 - SC-L county on the same draws. Sensitivity: county-cluster bootstrap (16 clusters, seed 20261103).
Promotion rule (primary): HCM2 is promoted over SC-L iff the 95 percent interval of the paired AUC difference lies above 0 AND the
HCM2 calibration slope point estimate lies in [0.8, 1.2]. Checks: label permutation (seed 7, AUC within 0.5 +/- 0.03), single-feature
screen (flag univariate AUC > 0.85), feature list contains no t2 or outcome column.

## 7. IMF2 rerun (copy of tracks/imf2/tools into tools/, paths repointed, imf2 restricted frames and draws read only)
- G1: the copy reproduces IMF2 V1 forecast and oracle MPI points and all 1,000 draws of IMF2 V1, OSM mortcal, E3, E8f to 1e-9 against
  tracks/imf2/restricted/main/draws.npz.
- Runs (all linked, probabilistic routing, same registered draw index, same comparators):
  R_win: layer = the winner of section 6 under county-blocked cross-fitting (HCM2 if promoted, else SC-L county) used alone.
  R_h: HCM2 alone. R_s: SC-L county alone. R_4: IMF2 deviance rule over {SC-L, SC-F, HCM refit, HCM2}.
- Headline "after" = R_win forecast. Reported: MPI [95%], C2 and C5 net points, paired MPI differences vs IMF2 V1 forecast (released,
  SC-L plot folds), OSM mortcal, E3 (and E8f), P(rank 1). Win rule: paired interval above 0. Inherited pillars stay inherited and labelled.

## 8. NCASI states (NH, VT, MI, WI, MN)
- State harvest frame built by tools/trkA_state_frame.py from the public FIADB DataMart copy
  ~/restricted/ncasi-modeleval/fia_data_rc/<ST> with the track H frame and H1 rules (t2 2004 to 2025, t1 1999+, t1 forest), validated
  by running it on the Maine DataMart (me_state/v26_H/restricted/fiadb_raw): pair count and cut labels must equal pairs_h exactly.
- State county folds by the section 2 greedy rule within the state. HCM2 refit within each state (same specification). Prices: statewide
  annual median of sawlog $/MBF (any log rule) and pulpwood $/cord ($/ton x 2.5) from stumpage_panel_full.csv, nearest year outside the
  covered years with an indicator; a state with fewer than 5 price years overlapping its t1 years is refit without price and stand value
  (NH has no rows; expected). No outbound fetch.
- IMF2 out of state with forecast cut status: frozen Maine links (rebuilt by the copied run), state pair frames
  ncasi_rc/inputs/<ST>/restricted (NH: ncasi_rc/NH/restricted), engine rows as imf2_oos.py, p from the state HCM2 OOF predictions joined
  by pair id (unmatched pairs take the state training base rate, counted). Statewide S on U and C, forecast and oracle, against single
  engines and the E3 / E8f mixes with the members present, imf2_oos location-component draws. Descriptive; no promotion rule.

## 9. Full chain (task 5)
C1 can be routed without a stepper only if per-stand gross and mortality values exist for the linked engines on the C1 frame; if so the
frozen links route them (no harvest in the C1 gate) and C1 is scored by the registered C1 arithmetic. P2 (hindcast), P3 (Bakuzis) and
P5 (trend) require a coupled multi-engine stepper (IMF3 found none; composing single-engine chains as a proxy is not allowed). They stay
inherited from OSM mortcal unless a stepper exists; if not, a design note and launcher stub are written and the item is reported open.

## 10. Gates (hard exits; nothing promoted on failure)
G1 reproduction (section 7). G2 row shuffle (seed 7) of every input frame: every headline value identical to 1e-9. G3 fresh-context
verifier subagent recomputes each headline by a different code path. G4 hull bound as IMF2. G5 leak scan of release/ (no CN, pair id,
PLOT, loc, LAT, LON column; no 10+ digit integer). No calibration or post-run adjustment is called an engine refinement (F12); no
hand-set constant is tuned (F16).
