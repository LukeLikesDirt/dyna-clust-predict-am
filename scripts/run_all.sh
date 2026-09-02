#!/bin/bash

# Script name:  run_all.sh
# Description:  Run the full pipeline for one primer set: prepare the
#               reference, predict cutoffs for each taxon group, then the
#               global Eukaryome cutoffs.
#
# Usage:        scripts/run_all.sh <primer_set> [taxon_group ...]
#               e.g. scripts/run_all.sh amv45nf_amdgr
#                    scripts/run_all.sh amv45nf_amdgr glomeromycota endogonomycetes fungi
#               (taxon_group defaults to: glomeromycota endogonomycetes)
#
# Note:         Must be run from the project root directory.
#               Does not run scripts/05_combine_cutoffs.sh -- run that
#               separately once every primer set you want is prepared.

set -eo pipefail

PRIMER_SET="${1:?Usage: $0 <primer_set> [taxon_group ...]}"
shift
TAXA=("$@")
if [[ ${#TAXA[@]} -eq 0 ]]; then
    TAXA=(glomeromycota endogonomycetes)
fi

echo "=============================================================================="
echo "RUN ALL: primer set = $PRIMER_SET, taxa = ${TAXA[*]}, plus eukaryome (global)"
echo "=============================================================================="

echo ""
echo ">>> [1/3] Preparing reference ($PRIMER_SET)..."
./scripts/01_prepare_reference.sh "$PRIMER_SET"

for taxon in "${TAXA[@]}"; do
    echo ""
    echo ">>> [2/3] Predicting cutoffs: $taxon ($PRIMER_SET)..."
    ./scripts/02_predict_cutoffs.sh "$PRIMER_SET" "$taxon"
done

echo ""
echo ">>> [3/3] Predicting global Eukaryome cutoffs ($PRIMER_SET)..."
./scripts/03_predict_euk_cutoffs.sh "$PRIMER_SET"

echo ""
echo "=============================================================================="
echo "RUN ALL COMPLETE: $PRIMER_SET"
echo "  Data: data/$PRIMER_SET/"
echo "  Run scripts/05_combine_cutoffs.sh to (re)write output/cutoffs_am*.txt"
echo "=============================================================================="
