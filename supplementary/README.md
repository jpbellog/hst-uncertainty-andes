# Supplementary tables

| File | Content |
|---|---|
| `TableS1_hotspot_fits.csv`, `tableS1.tex` | Table S1: hotspot fits of the twelve models |
| `TableS2_collision_bands.csv` | Table S2: collision bands of every model on the present-day margin |
| `TablesS_data.xlsx` | Spreadsheet of supplementary data (sheets README, S1_hotspot_fits, S2_collision_bands, S3_geochronology, S4_hotspot_seeds, S5_convergence) |

`make_supplementary_tables.py` writes the two tables and `tableS1.tex`; `make_supplementary_excel.py`
writes the spreadsheet. Both read the saved outputs of `models/` and `results/`.
