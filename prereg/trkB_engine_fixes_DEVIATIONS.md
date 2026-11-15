# trkB_engfix deviations (PREREG md5 8378f9e1eaff5ff9e7d74ef45807e686)
D1 (14 Nov 2026, compute, before any score): engine runs went to OSC Cardinal as the brief then directed: job 16747881 (failed in
14 s: keyfile builder reads kf/src-converted/vls/spctrn.f90, a symlink not followed in the bundle; patched with md5 check), job
16748037 (default pair runs with COMPUTE, 18 s) and job 16748372 (ddm pairs, S1, hindcast, Bakuzis, 80 s), 48 cores each: about
1.5 core-hours in total. Inputs were a coordinate-free bundle (no LAT/LON columns; md5 manifest identical on both hosts). After
Aaron's 14 Nov rule (firebreather only, $0 at OSC) nothing more was sent to Cardinal; scoring, fits and checks ran on
firebreather at nice 19, one process. Transfers passed through the connector host's /tmp/trkB_xfer (bundle and pair tables; no
coordinates); those copies are still there.
D2 (14 Nov 2026, method as run, before any score): scipy is not in the firebreather venv, so the PREREG fit (scipy BFGS) was run as
Newton with the analytic gradient and a finite-difference Hessian (numpy). Same likelihood, weights and data; all six fits
converged (max |gradient| below 2e-16).
D3 (14 Nov 2026, scope): the C2 pair runs cover all 3,790 holdout pairs of the frame; scoring uses the registered P* coverage
(3,636), unchanged. 133 pair ids occur in both C2 (cut overlap) and C5; RD for the fit comes from the C2 run (only uncut-fit C2
pairs enter the fit).
D4 (14 Nov 2026, scope): C1 and C4bfm (P5) were not rerun (PREREG section 2 conditional): the gate stand runner and the y50 bau
runner need new scoring plumbing (design totals and the bfm index) beyond this session; they are INHERITED from W a3 (P1 keeps
W a3 C1 at weight 0.175; P5 weight 0.15). Inherited MPI share for ddm: 32.5 percent.
D5 (14 Nov 2026, task 2): PREREG section 4 check found removals already separated in the CEM treated donors (GRM MORTALITY vs
CUT/DIVERSION components since v25_C D3; round4 R-a T2 := T1 + gross - mortality; donor time 1 residual by the treated_c5 hook), so
no CEM rerun was made. The decomposition and the score of an existing registered sensitivity run (v24_H cem_v2s, density cap and
disturbance off) are reported as diagnostics; a rerun of the current configuration with cap and disturbance off is a decision
for Aaron (post hoc, not made here).
D6 (14 Nov 2026, task 3): no in-engine mortality hook was found in OSM 2.26.1 (strings of OSM.ConsoleApp show a DBH increment
"Genetic Gain" multiplier and a DBHI/HTI .osmcal calibrator, no documented mortality modifier); no OSM rerun; rows stay
"post-run adjustment".
D7 (14 Nov 2026, G3): this agent has no subagent tool, so the fresh-context verifier was not launched; independent second paths
were run instead (base R recomputation of C2 and C5 S; C3 scorer reproduces the registered W a3 C3 S and draws exactly).
