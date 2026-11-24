#!/bin/bash
# Runs the full chain that reproduces every number, table and figure of the article.
# See README.md. Third-party inputs must be in place first (THIRD_PARTY.md).
#
#   bash run_all.sh            dry run: checks inputs and lists the commands of the chain
#   bash run_all.sh --run      runs everything (several hours; the Monte Carlo step dominates)
set -u
export PAPER_ROOT="${PAPER_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
export PYTHONIOENCODING=utf-8
# the scripts call the Windows launcher "py"; elsewhere map it to python3
if ! command -v py > /dev/null 2>&1; then py() { python3 "$@"; }; export -f py; fi
P="$PAPER_ROOT"; A="$P/analysis"

echo "== PAPER_ROOT=$P"
missing=0
for f in "external/wilson2025/wilson2025_SI.pdf" \
         "external/doubrovine2012/Supporting_Information/jgrb17279-sup-0004-ts03.txt" \
         "external/doubrovine2012/Supporting_Information/jgrb17279-sup-0008-ds01.txt" \
         "external/gaastra2022/Supporting_Information/2021tc006772-sup-0003-table si-s02.csv" \
         "external/quiero2022/MOESM5.txt" \
         "external/muller2019/Coastlines_USGS_MedRes_AREPS_polyline.shp" \
         "external/vanhinsbergen2026/Utrecht_Paleogeography_15Dec25.rot" \
         "external/wagenaar2025/Rotation_files.zip" \
         "external/etopo1/ETOPO1_Bed_g_geotiff.tif"; do
  [ -f "$P/$f" ] || { echo "   missing third-party input: $f"; missing=1; }
done
[ $missing -eq 1 ] && echo "   (see THIRD_PARTY.md)"

if [ "${1:-}" != "--run" ]; then
  echo "== dry run"
  echo "   1. py models/WD25_WK08-A/parse_wd25_si.py      Wilson and DeMets (2025) tables from their SI"
  echo "   2. py rotation_models/build_rotation_models.py           12 rotation files, compared with models/*/rot"
  echo "   3. py rotation_models/fit_hotspots.py             12 hotspot fits, compared with models/*/hotspot_fit.json"
  echo "   4. bash $A/run_chain.sh --run         Monte Carlo, deposits, sensitivities, figures, tables, manuscript_numbers.py"
  bash "$A/run_chain.sh"
  exit 0
fi
[ $missing -eq 1 ] && { echo "== stopped: third-party inputs missing"; exit 1; }
( cd "$P/models/WD25_WK08-A" && py parse_wd25_si.py ) || exit 1
( cd "$P/rotation_models" && py build_rotation_models.py ) || exit 1
( cd "$P/rotation_models" && py fit_hotspots.py ) || exit 1
bash "$A/run_chain.sh" --run
