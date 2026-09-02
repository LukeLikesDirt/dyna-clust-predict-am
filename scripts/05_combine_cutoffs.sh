#!/bin/bash

# Script name:  05_combine_cutoffs.sh
# Description:  Wrapper around R/combine_cutoffs.R that auto-discovers what
#               to combine, so adding a new primer pair needs no edit here:
#                 - primer sets: every config/primers/*.conf
#                 - taxa:        union of *.cutoffs.json.txt basenames found
#                                under data/<primer_set>/ for those primer
#                                sets (missing primer_set x taxon
#                                combinations are skipped by the R script,
#                                not an error)
#
# Usage:        scripts/05_combine_cutoffs.sh [--primer_sets a,b,...] [--taxa x,y,...]
#               e.g. scripts/05_combine_cutoffs.sh
#                    scripts/05_combine_cutoffs.sh --primer_sets wanda_aml2
#
# Note:         This script must be run from the project root directory,
#               after 02_predict_cutoffs.sh / 03_predict_euk_cutoffs.sh have
#               produced at least one *.cutoffs.json.txt.

set -eo pipefail

readonly COMBINE="./R/combine_cutoffs.R"
readonly PRIMER_CONF_DIR="./config/primers"
readonly DATA_DIR="./data"

PRIMER_SETS_OVERRIDE=""
TAXA_OVERRIDE=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --primer_sets) PRIMER_SETS_OVERRIDE="$2"; shift 2 ;;
        --taxa)        TAXA_OVERRIDE="$2"; shift 2 ;;
        *) echo "ERROR: Unknown argument: $1" >&2; exit 1 ;;
    esac
done

[[ -f "$COMBINE" ]] || { echo "ERROR: R script not found: $COMBINE" >&2; exit 1; }

# =============================================================================
# DISCOVER PRIMER SETS
# =============================================================================

# Portable across bash 3.2 (macOS default -- no mapfile/readarray/associative
# arrays, all bash 4.0+) and newer: build arrays with a plain read loop.

if [[ -n "$PRIMER_SETS_OVERRIDE" ]]; then
    IFS=',' read -ra PRIMER_SETS <<< "$PRIMER_SETS_OVERRIDE"
else
    PRIMER_SETS=()
    while IFS= read -r name; do
        PRIMER_SETS+=("$name")
    done < <(ls "$PRIMER_CONF_DIR"/*.conf 2>/dev/null | xargs -n1 basename | sed 's/\.conf$//' | sort)
fi

if [[ ${#PRIMER_SETS[@]} -eq 0 ]]; then
    echo "ERROR: No primer sets found in $PRIMER_CONF_DIR/*.conf" >&2
    exit 1
fi

# =============================================================================
# DISCOVER TAXA (union of *.cutoffs.json.txt across discovered primer sets)
# =============================================================================

if [[ -n "$TAXA_OVERRIDE" ]]; then
    IFS=',' read -ra TAXA <<< "$TAXA_OVERRIDE"
else
    TAXA=()
    while IFS= read -r name; do
        [[ -n "$name" ]] && TAXA+=("$name")
    done < <(
        for ps in "${PRIMER_SETS[@]}"; do
            for f in "$DATA_DIR/$ps"/*.cutoffs.json.txt; do
                [[ -f "$f" ]] || continue
                basename "$f" .cutoffs.json.txt
            done
        done | sort -u
    )
fi

if [[ ${#TAXA[@]} -eq 0 ]]; then
    echo "ERROR: No *.cutoffs.json.txt found under $DATA_DIR/<primer_set>/ for any of: ${PRIMER_SETS[*]}" >&2
    echo "Run 02_predict_cutoffs.sh / 03_predict_euk_cutoffs.sh first." >&2
    exit 1
fi

echo "=============================================================================="
echo "COMBINE CUTOFFS"
echo "  Primer sets (discovered): ${PRIMER_SETS[*]}"
echo "  Taxa (discovered):        ${TAXA[*]}"
echo "=============================================================================="

echo "Activating conda environment..."
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate dyna_clust_predict

IFS=',' PRIMER_SETS_CSV="${PRIMER_SETS[*]}"
IFS=',' TAXA_CSV="${TAXA[*]}"
unset IFS

Rscript "$COMBINE" --primer_sets "$PRIMER_SETS_CSV" --taxa "$TAXA_CSV"

conda deactivate
