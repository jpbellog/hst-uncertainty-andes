# Notes for Claude Code (and other coding agents)

This repository reproduces the results of an article on hotspot-track uncertainty. Read `README.md`
first; each folder has a README describing its files.

## Environment
- Set `PAPER_ROOT` to the repository root before running anything; every script reads it.
- The scripts call `py` (Windows launcher). On Linux or macOS define `py() { python3 "$@"; }` or run
  through `run_all.sh`, which does it.
- Set `PYTHONIOENCODING=utf-8` on Windows.
- Shell scripts are bash.

## Order of work
1. Quick check, no third-party data needed: `cd analysis && py manuscript_numbers.py`. The last line
   must report 309/309.
2. For a full run, place the third-party files listed in `THIRD_PARTY.md` under `external/` and run
   `bash run_all.sh` (dry run: checks inputs and prints the steps), then `bash run_all.sh --run`.
   In the dry run, MISSING marks intermediate outputs that are not shipped and OLD depends on file
   dates; neither is an error after a fresh download.
3. `rotation_models/build_rotation_models.py` and `rotation_models/fit_hotspots.py` compare what they
   build with the files used in the article and print the differences; they are the first checks after
   any change to inputs or code.

## Conventions
- Model folders are `models/<circuit>_<frame>`. The pipeline scripts are the same in every folder,
  and each folder holds the steps its model needs (step 05 and step 12 only in `Q22_WK08-A`). Each copy
  works on the files of its own folder; run a step from inside that folder with the
  environment set by `analysis/run_models.sh` (rotation file, hotspot fit, prior and frame covariance).
- Step 04 (Monte Carlo) is the long one: about 25 minutes per model, variant and spreading partition
  at N = 1000. Pass a small N as first argument (e.g. `py 04_mc_collisions.py 20 total 0.44`) and
  `PAPER_SUFFIX=_test` to test without overwriting the saved outputs.
- Outputs are compared byte by byte or number by number; do not reformat output files.
- Some variable names and some strings written to output files are in Spanish. They are identifiers
  or data read by other scripts; do not translate them.
