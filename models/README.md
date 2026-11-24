# Models

One folder per model, `<circuit>_<frame>` (see the table in the root `README.md`). Every folder holds
the same copy of the pipeline scripts, its rotation files, the age compilations and the outputs used in
the article. A script works on the files of the folder it is run from.

## Pipeline steps

| Script | Step | Folders |
|---|---|---|
| `02_rot_engine.py` | Rotation engine: reads a `.rot` file with covariances and composes finite rotations; imported by the other steps | all |
| `03b_hotspot_fit.py` | Fit of the hotspot positions to the dated seamounts (writes `hotspot_fit*.json` and the `_hsb` rotation files) | all |
| `04_mc_collisions.py` | Monte Carlo of the track–trench collision points | all |
| `05_vc_uncertainty.py` | Analytical uncertainty of the convergence velocity from Quiero et al. (2022) | `Q22_WK08-A` |
| `07_deposit_tracks.py` | Assignment of deposits and igneous samples to the nearest hotspot track | all |
| `08_track_ellipses.py` | Uncertainty ellipses of the present-day position of the tracks | all |
| `09_vc_uncertainty_circuit.py` | Convergence-velocity uncertainty from a Monte Carlo of the plate circuit | `*_WK08-A` |
| `10_recon_ellipses.py` | Uncertainty ellipses of the reconstructed tracks per epoch | all |
| `11_recon_model_tracks.py` | Nominal tracks of the f_PAC = 0.5 and 0.36 models per epoch | all |
| `12_deposit_track_test.py` | Test of the deposit–track association across models | `Q22_WK08-A` |
| `collision_primitives.py` | Geodesic primitives and the track–trench intersection; imported by step 04 and, through it, by the later steps | all |

`WD25_WK08-A` also holds `parse_wd25_si.py`, which reads Tables S2.2–S2.4 of Wilson and DeMets (2025)
into `wd25/` (not distributed; see `THIRD_PARTY.md`), and `build_rot_wd25.py`, which
`rotation_models/build_rotation_models.py` calls to write the WD25 plate circuit.

Steps 04, 07, 08, 10, 11 and 12 are run by `analysis/run_models.sh`, steps 05 and 09 by
`analysis/run_chain.sh`, and step 03b by `rotation_models/fit_hotspots.py`. `run_models.sh` sets for each folder the rotation file, the hotspot fit,
the prior on the hotspot motion and the frame covariance. To run one step by hand, copy that environment
from the `MODELOS` table and the `step()` function of `run_models.sh`.

## Files in each folder

| File | Content |
|---|---|
| `data/ages_*.csv` | Age compilations of the Easter–Nazca, Juan Fernández, Taltal–San Félix and Foundation chains; the source of each age is in the `#` header and in the rows |
| `data/wk08a_covariances.csv` | Finite rotations and covariances of WK08-A (Wessel and Kroenke, 2008, Table 3) |
| `data/fz_challenger_gsfml.csv` | Trace of the Challenger fracture zone (GSFML) |
| `rot/*_hsb.rot`, `rot/*_hsb_gauss_pred_age_A.rot` | Rotation files of the model with the fitted hotspots (w as confidence and as prediction term) |
| `hotspot_fit.json`, `hotspot_fit_gauss_pred_age_E.json` | Hotspot fits read by steps 04 to 11 |
| `hotspot_fit_s9_*.json` | Fits of the sensitivity tests of Supplementary Section S4 (nine main-text models) |
| `mc_collision_bands_total_hsb_cc_g13.csv` | Collision bands of step 04 (all sources) |
| `track_*`, `recon_*` | Outputs of steps 08, 10 and 11 |
| `deposit_tracks_*.csv` | Output of step 07 |
| `vc_sigma_*.csv` | Outputs of steps 05 and 09 |
| `deposit_track_test_by_model.csv` | Output of step 12 (`Q22_WK08-A`) |

The collision bands restored to the present-day margin, by variant and model, are in `results/`.
