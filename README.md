# Data and code: uncertainty budget of subducting hotspot tracks off western South America

Rotation models, geochronological compilations, deposit and igneous-sample databases, and the Python
code behind the article

> Bello-González, J.P., Lara, L.E., Contreras-Reyes, E., Quiero Silva, F. Formal uncertainty budget for
> subducting hotspot-track reconstructions and trench-collision points along the Cenozoic Andean margin
> of Chile and Peru. Submitted to Earth and Planetary Science Letters.

Every number, table and figure of the article and its Supplementary Material can be regenerated from the
files in this repository and the third-party datasets listed in [`THIRD_PARTY.md`](THIRD_PARTY.md).
`analysis/manuscript_numbers.py` recomputes the 309 numbers quoted in the text from the saved outputs and
compares each one with the value printed in the article.

## Layout

| Folder | Content |
|---|---|
| [`inputs/`](inputs/) | Starting rotation model, deposit and igneous-sample databases, restored trench, GPlates validation tracks |
| [`rotation_models/`](rotation_models/) | Construction of the 12 rotation files (plate circuit × hotspot frame) and of the hotspot fits |
| [`models/`](models/) | One folder per model: Monte Carlo pipeline (steps 02–12), age compilations, rotation files and outputs |
| [`analysis/`](analysis/) | Orchestrators, covariance and time-scale tools, sensitivity tests, summaries, `manuscript_numbers.py` |
| [`tables/`](tables/) | Restoration of the collision bands to the present-day margin |
| [`results/`](results/) | Collision bands on the present-day margin, the numbers behind the model comparison (Fig. 3) and the deposit–track distances |
| [`figures/`](figures/) | Scripts of every figure of the article and the Supplementary Material |
| [`supplementary/`](supplementary/) | Supplementary Tables S1 and S2 and the spreadsheet of supplementary data |
| [`tools/hellinger_chang/`](tools/hellinger_chang/) | Python port of the Hellinger–Chang finite-rotation toolkit, with its validation |
| `external/` | Third-party datasets, downloaded by the user (not distributed; see `THIRD_PARTY.md`) |

Each folder has its own README with the role of every file.

## Models

Nine models enter the main text (three plate circuits × three hotspot frames). The three GMHRF-A'
models are reported in the Supplementary Material; their folders are labelled `GMHRF-A`.

| Circuit \ frame | WK08-A | GMHRF-B | GWG20 | GMHRF-A' (suppl.) |
|---|---|---|---|---|
| Q22, Quiero et al. (2022) | `Q22_WK08-A` | `Q22_GMHRF-B` | `Q22_GWG20` | `Q22_GMHRF-A` |
| WD25, Wilson and DeMets (2025) | `WD25_WK08-A` | `WD25_GMHRF-B` | `WD25_GWG20` | `WD25_GMHRF-A` |
| SG12, Somoza and Ghidella (2012) | `SG12_WK08-A` | `SG12_GMHRF-B` | `SG12_GWG20` | `SG12_GMHRF-A` |

All node ages are on the Geologic Time Scale 2020.

## Requirements

Python 3 with the packages in `requirements.txt` (the versions used for the article) and `bash`
(Git Bash on Windows). The scripts call the Windows launcher `py`; `run_all.sh` maps it to `python3`
where `py` does not exist. Every script reads the repository root from the environment variable
`PAPER_ROOT`.

```bash
pip install -r requirements.txt
export PAPER_ROOT=/path/to/this/repository
```

## Quick check (minutes)

With the outputs included here, and without the third-party data:

```bash
cd analysis && py manuscript_numbers.py        # expected last line: 309/309 match
```

## Full reproduction (hours)

Place the third-party inputs as described in `THIRD_PARTY.md`, then

```bash
bash run_all.sh          # checks the inputs and lists the steps without running them
bash run_all.sh --run    # runs everything
```

The dry run also lists the intermediate outputs of the chain (Monte Carlo parts, sensitivity runs)
as MISSING, because the repository ships only the final outputs, and flags as OLD the outputs whose
file date is earlier than the hotspot fit of their folder; after a fresh download these flags come
from the file dates and do not indicate an error.

The steps, in order:

1. `models/WD25_WK08-A/parse_wd25_si.py` reads Tables S2.2–S2.4 of Wilson and DeMets (2025) from their
   Supporting Information.
2. `rotation_models/build_rotation_models.py` builds the 12 rotation files and compares them row by row
   with `models/*/rot/`.
3. `rotation_models/fit_hotspots.py` refits the hotspot positions of the 12 models and compares them with
   `models/*/hotspot_fit.json`.
4. `analysis/run_chain.sh` runs the Monte Carlo propagation (step 04, N = 1000 per model and spreading
   partition, about 25 minutes per job), the deposit–track test, the sensitivity tests, the figures and
   tables, and ends with `manuscript_numbers.py`.

Steps 2 and 3 reproduce the files used in the article exactly (76 covariances of the GWG20 models differ
in the fifth significant figure, from the rounding of the published values).

## File-name tags

| Tag | Meaning |
|---|---|
| `_hsb` | Hotspot positions from the fit of `03b_hotspot_fit.py`, with the geological dispersion w used as a confidence term |
| `_hsb_pred` | The same fit with w used as a prediction term |
| `_cc` | Circuit covariances of variant b (`analysis/circuit_covariance.py`) |
| `_g13` | Rotation files with the nodes of Granot et al. (2013) |
| `pf`, `f_PAC` | Spreading partition of the Nazca–Pacific ridge (0.36, 0.44, 0.5) |
| `total`, `circuit`, `hotspot`, `wk08` | Monte Carlo variant: all sources, plate circuit only, hotspot fit only, hotspot frame only |

## Data of this work

The deposit and igneous-sample databases update the compilation published with Bello-González et al.
(2026), doi:10.1016/j.pce.2025.104202. The rotation files combine published rotations and covariances
(sources in each row and in `rotation_models/README.md`) with covariances fitted or interpolated in this
work, as described in Supplementary Section S1.

## Licenses

Code: MIT (`LICENSE`). Data produced by the authors: CC BY 4.0 (`LICENSE-DATA.md`). Third-party datasets
keep their own licenses and are not included.
