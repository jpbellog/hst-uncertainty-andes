#!/bin/bash
# Andean shortening with the NINE models of the text (g13): step 07 with the Utrecht positions,
# with the Utrecht trench and with our trench (own_trench).
P="${PAPER_ROOT:?define PAPER_ROOT}"; A="$P/analysis"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PAPER_CIRCUIT_COV=b
MODELOS="Q22_WK08-A|Global_Model_R_BG_2025_hsb.rot|none|wk08
Q22_GMHRF-B|Global_Model_R_BG_2025_GMHRF_gts12_hsbgmhrf.rot|gmhrf|rot
Q22_GWG20|Global_Model_R_BG_2025_GWG20_hsb.rot|none|rot
WD25_WK08-A|Global_Model_R_BG_2025_WD25_hsb.rot|none|wk08
WD25_GMHRF-B|Global_Model_R_BG_2025_WD25_GMHRF_gts12_hsbgmhrf.rot|gmhrf|rot
WD25_GWG20|Global_Model_R_BG_2025_WD25_GWG20_hsb.rot|none|rot
SG12_WK08-A|Global_Model_R_BG_2025_SG12_hsb.rot|none|wk08
SG12_GMHRF-B|Global_Model_R_BG_2025_SG12_GMHRF_gts12_hsbgmhrf.rot|gmhrf|rot
SG12_GWG20|Global_Model_R_BG_2025_SG12_GWG20_hsb.rot|none|rot"
uno() {
  local d="$1" v="$2" line rc prior cov
  line=$(echo "$MODELOS" | grep "^$d|"); IFS="|" read -r _ rc prior cov <<< "$line"
  PAPER_ROT_FILE="$P/models/$d/rot/$rc" PAPER_HS_FIT="$P/models/$d/hotspot_fit.json" PAPER_HS_PRIOR="$prior" \
    PAPER_FRAME_COV="$cov" py "$A/step07_utrecht.py" "$P/models/$d" $v \
    > "$A/outputs/andean_shortening/ds07_$d$v.log" 2>&1
  echo "$(date +%T) $d $v exit $?"
}
export -f uno; export P A MODELOS
J=""
for line in $MODELOS; do d="${line%%|*}"; J+="$d"$'\n'"$d own_trench"$'\n'; done
printf "%s" "$J" | xargs -P "${NP:-3}" -L 1 bash -c 'uno "$@"' _
echo "end $(date '+%F %T')"
