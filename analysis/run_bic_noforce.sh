#!/bin/bash
# Suppl. S4 (Easter BIC): 03b fit WITHOUT forcing the Easter drift (no PAPER_HS_FORCE),
# for the fixed/drift BIC number. Final outputs are not overwritten: json and residuals carry the tag _s9_noFORCE, .rot in
# rot/sens9/. Read by manuscript_numbers.py (section bic).
# Base case without seed (default 03b), in the twelve models.
# Usage:  bash run_bic_noforce.sh
P="${PAPER_ROOT:?define PAPER_ROOT}"; A="$P/analysis"
MODELOS=$(sed -n '/^MODELOS="/,/"$/p' "$A/run_models.sh" | sed 's/^MODELOS="//; s/"$//')
uno() {
  local d="$1" rc="$2" prior="$3"
  mkdir -p "$P/models/$d/rot/sens9"; cd "$P/models/$d" || return
  env -u PAPER_HS_FORCE -u PAPER_HS_P0 PAPER_CIRCUIT_COV=b PAPER_TAG_EXTRA=_g13 PAPER_ROT_FILE="$P/models/$d/rot/$rc" \
    PAPER_HS_PRIOR="$prior" PAPER_HS_OUT_TAG=_s9_noFORCE PAPER_HS_ROT_OUT="$P/models/$d/rot/sens9/_s9_noFORCE.rot" \
    OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 py 03b_hotspot_fit.py > hotspot_fit_s9_noFORCE.log 2>&1
  echo "$d exit $?"
}
export -f uno; export P
for line in $MODELOS; do IFS="|" read -r d rc rp prior cov <<< "$line"; echo "$d $rc $prior"; done \
  | xargs -P "${NP:-2}" -L 1 bash -c 'uno "$@"' _
