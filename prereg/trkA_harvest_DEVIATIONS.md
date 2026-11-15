# trkA_hcm deviations from PREREG.md (md5 35055a79af9083edec3f5e3223552dd7)

- trkA-D1 (code, 14 Nov, before any MI/WI/MN out-of-state score existed). tools/trkA_oos.py (copied from imf2_oos.py) intersected the scored
  pair set over every comparator, so a comparator with no rows in a state (OSM single engine in MI, WI, MN) would have emptied the frame.
  Comparators with no finite value in the state are now dropped before scoring. NH and VT results are unaffected (no comparator was empty;
  NH and VT ran before the patch and the shuffle reruns ran after it).
- trkA-D2 (implementation, state builder). In the state frames, missing t1 ELEV or RDDISTCD is filled with the state median and a missing
  type group is coded -1 (no such rows exist in the Maine frame, where the builder reproduces pairs_h). Not a rule change.
- trkA-D3 (scope, PREREG section 9). Per-stand C1 component values do exist on firebreather for the linked engines (mortcal stand_results for
  OSM, OSM mortcal, AcadianGY; fvs_gate_stand_values for FVS-NE default, growcal and FVS-ACD; CEM plot gate). The C1 statistic, however, is
  computed inside the T4 scorer (tools/score/score.py main: design bootstrap bd1, domain targets F_pt and F_bt) and is reachable only as a
  score.py extension hook over the full W a3 rescoring chain. That rerun was not started in this session (shared load, one worker budget), so
  C1 is still inherited from OSM mortcal. The route and the launcher plan are in FULLCHAIN_DESIGN.md. P2, P3 and P5 remain blocked as
  registered (no coupled multi-engine stepper exists).
- trkA-D4 (data, F7, found 14 Nov while the MN fit was running and before any MN score existed). stumpage_panel_full.csv codes missing
  state prices as 0: 77 percent of MN sawlog rows and 50 percent of MN pulpwood rows are 0 (statewide medians 0 from 2006 on), VT 0.3 percent
  and WI 0.1 percent of pulpwood rows; ME and MI have none. Rows with value <= 0 are now dropped before the statewide medians (the Maine
  county rule is unaffected: no ME row is <= 0). All five state layer fits and out-of-state runs were rerun with the corrected prices; the
  earlier NH, VT, MI and WI outputs are superseded (NH has no price rows and is numerically unchanged).
