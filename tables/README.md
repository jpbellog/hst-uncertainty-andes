# Tables

| Script | Role |
|---|---|
| `make_present_bands.py` | Maps the latitude percentiles of the step 04 collision bands from the margin restored to each age onto the present-day margin |
| `margin_map.py` | Node correspondence of `inputs/margin/Margin_SAM_0.1.csv` between age t and the present; imported by `make_present_bands.py` and some figures |

`analysis/run_chain.sh` runs `make_present_bands.py` for every model with `PAPER_TABLES_DIR` set to
`results/`, where the present-day bands are kept. Without that variable the output is written here.
`PAPER_PIPE` selects the model folder and `PAPER_BANDS_SUFFIX` the suffix of the output names.
