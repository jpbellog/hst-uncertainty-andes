#!/bin/bash
# Foundation and the Selkirk microplate. The microplate extends from the failed rift
# at 124 45'W to the eastern boundary at 121 40'W (O'Connor et al. 1998, EPSL 161, section 2).
# Samples used in the fit inside that segment: SO100-33DS-1 (12.7 Ma) and SO100-38DS-1 (11.8 Ma);
# SO100-41DS-1 (11.8 Ma) lies ~0.1 degrees east of the eastern boundary.
#   (a) _s9_sinSelkirk : without SO100-33DS-1 and SO100-38DS-1
#   (b) _s9_sinSelkirkB: also without SO100-41DS-1
# Base-case 03b (no seed, no forced drift, the BIC chooses). The .rot files go to rot/sens9/ and
# the json files carry the tag: nothing overwrites the final outputs. Summary: sens_selkirk.py.
# Usage:  bash run_sens_selkirk.sh > run_sens_selkirk.log 2>&1
P="${PAPER_ROOT:?define PAPER_ROOT}"; A="$P/analysis"
NP=${NP:-2}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
MODELOS=$(sed -n '/^MODELOS="/,/"$/p' "$A/run_sens9.sh" | sed 's/^MODELOS="//; s/"$//')
uno() {  # folder tag excluded
  local d="$1" tag="$2" ex="$3" line rc prior cov
  line=$(echo "$MODELOS" | grep "^$d|"); IFS="|" read -r _ rc prior cov <<< "$line"
  mkdir -p "$P/models/$d/rot/sens9"
  ( cd "$P/models/$d" && env -u PAPER_HS_P0 -u PAPER_HS_FORCE PAPER_ROT_FILE="$P/models/$d/rot/$rc" PAPER_HS_PRIOR="$prior" \
      PAPER_HS_EXCLUDE="$ex" PAPER_HS_OUT_TAG="$tag" PAPER_HS_ROT_OUT="$P/models/$d/rot/sens9/$tag.rot" \
      py 03b_hotspot_fit.py > "hotspot_fit$tag.log" 2>&1 )
  echo "$(date +%T) $d $tag exit $?"
}
export -f uno; export P MODELOS
echo "== start $(date '+%F %T')"
J=""
for line in $MODELOS; do d="${line%%|*}"
  J+="$d _s9_sinSelkirk SO100-33DS-1,SO100-38DS-1"$'\n'
  J+="$d _s9_sinSelkirkB SO100-33DS-1,SO100-38DS-1,SO100-41DS-1"$'\n'
done
printf "%s" "$J" | xargs -P "$NP" -L 1 bash -c 'uno "$@"' _
echo "== end $(date '+%F %T')"
