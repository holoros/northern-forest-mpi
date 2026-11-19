# northern-forest-mpi

Replication code and state-level results for the Northern Forest multi-model evaluation: eight forest carbon projection models and an equal-weight ensemble scored against US Forest Service FIA remeasurements in Maine, New Hampshire, Vermont, Michigan, Wisconsin and Minnesota with a composite Model Performance Index (MPI, 0 to 100).

The manuscript is in preparation (target: Forest Ecology and Management). The general area-of-interest workflow that runs the engines lives in a separate repository, [holoros/multimodel-carbon](https://github.com/holoros/multimodel-carbon); this repository holds only the code behind the paper's scores, figures and tables.

## Contents

| Folder | What it holds |
|:-|:-|
| `scoring/` | The MPI scorer: pillars P1 to P5 (statewide totals and uncut pairs, hindcast, Bakuzis realism, cut pairs, FIA trend), bootstrap intervals (1,000 draws, seed 20261103), the second-path R check (`second_path_mpi.R`) and unit tests (`test_rc.py`) |
| `datapack/` | Builders for the manuscript data pack: scorecard (Fig. 1), growth component ratios on uncut and cut pairs (Fig. 2), TreeMap versus FIA initialization contrasts (Fig. 3), Table 1 and Table 2; `run_queue3.sh` gives the run order |
| `tracks/trkA_harvest/` | Harvest choice layer (refit of Zhao, Daigneault and Weiskittel 2020 with county-blocked folds) and the linked IMF2 framework |
| `tracks/trkB_engine_fixes/` | FVS-NE in-engine density-dependent mortality and the CEM removals split |
| `tracks/trkC_ga_in/` | Georgia and Indiana extension: FVS SN and CS curve tables for CBM-CFS3 with FVS curves (code key `cbm_gated`, whole-plot basis), CEM cut-pair port, LANDIS-II Georgia |
| `prereg/` | Pre-registrations (frozen by md5 before scoring), the Track C amendment and the logged deviations |
| `data/` | State-level outputs used for the figures and tables (no plot-level records) |
| `figures/` | Scripts that draw Figures 1 to 3 from `data/` |

## Model names in the code and in the paper

The scripts, pre-registrations and logs use short internal keys. The paper uses descriptive names.

| Code key or term | Name in the paper | What it is |
|:-|:-|:-|
| `fvs` | FVS default | Forest Vegetation Simulator, Northeast (ME, NH, VT) or Lake States (MI, WI, MN) variant, no calibration |
| `osm` | OSM | Open Stand Model 2.26.1, Acadian variant (Maine only) |
| `acd141` | AcadianGY 14.1 | Acadian individual-tree model with refit increment equations (Kuehne et al. 2022; Maine only) |
| `cem` | CEM | Plot imputation by coarsened exact matching to remeasured FIA donor plots |
| `landis` | LANDIS-II | LANDIS-II v8 with Biomass Succession |
| `cbm_gated` ("CBM gated") | CBM-CFS3 with FVS curves | CBM-CFS3 (libcbm 2.10.2) driven by Chapman-Richards volume curves fitted to calibrated FVS projections of FIA conditions (FVS_FIAW_NEW cell of the registered gate plan; Georgia and Indiana used curves from default FVS) |
| `cbm_af` ("CBM AF", "American Forests method") | CBM-CFS3 with Gompertz curves | CBM-CFS3 driven by Gompertz volume curves fitted to FIA plots by forest type group and productivity class, following Papa et al. (2023); presets for ME, MI, WI and MN only |
| `fia_curves` | FIA yield curves | Empirical biomass yield curves fitted to the same FIA plots (in sample, net growth only) |
| `E_state` ("state ensemble") | Ensemble | Equal-weight mean of the members available in each state |
| tier 1 / tier 2 | fully scored / partially scored | Coverage of the pillar weight at least 0.65 / below 0.65 |
| mirror (AF mirror) | reconstruction | Reconstruction of the published Lake States CBM-CFS3 configurations (Papa et al. 2026a to c) |

## Reproduce the figures

```bash
pip install pandas numpy matplotlib
python figures/make_fig1_fig3.py data/fig1_scorecard_v3.csv figures/out
python figures/make_fig2.py data/fig2_components.csv figures/out
```

`data/fig1_scorecard_v3.csv` is the manuscript scorecard: `v2` plus the Maine LANDIS-II cut-plot score and the Maine ensemble rescored with it. Figure 3 values are transcribed in `make_fig1_fig3.py` from `data/fig3_treemap_direct.csv` and `data/fig3_treemap_phase2.csv`.

## Rerun the scoring

The scoring and track scripts were run on a single Linux server (Python 3.9, R 4.x with data.table) and still carry that server's absolute paths (`~/perseus_me/...`); set them to your own layout before running. Engine runs (FVS, OSM, AcadianGY, CEM, LANDIS-II, CBM-CFS3 through libcbm) come from multimodel-carbon and the engine repositories it names.

## What is not here

- FIA plot coordinates. Public FIA coordinates are perturbed; some steps (CEM matching covariates and raster extractions) used true plot locations under a US Forest Service data agreement. Those inputs and every plot-level intermediate stay on the original server. Rerunning those steps with public coordinates is possible but will not reproduce the scores exactly.
- Plot-level model outputs, engine binaries and large intermediate files.

## License and citation

Code is released under the MIT License. Please cite the manuscript once published; until then cite this repository using `CITATION.cff`.

## Funding

NCASI Model Evaluation Phase II and the USDA NIFA PERSEUS project.
