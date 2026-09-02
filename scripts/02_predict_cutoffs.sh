#!/bin/bash

# Script name:  02_predict_cutoffs.sh
# Description:  Combined pipeline for taxon-group cutoff prediction, for any
#               primer set already prepared by 01_prepare_reference.sh:
#               1. Build per-rank subsets from eukaryome_full for the taxon
#                  group (optionally cascading from an identified base rank
#                  -- see CASCADE_BASE_RANK in config/taxa/<taxon>.conf).
#               2. Compute one master similarity matrix (reused across ranks).
#               3. Predict global cutoffs for each configured rank.
#
#               Replaces the old 02/03/04_predict_*_cutoffs_V4.sh, which were
#               ~95% identical copies differing only in the values now in
#               config/taxa/<taxon>.conf.
#
# Usage:        scripts/02_predict_cutoffs.sh <primer_set> <taxon_group>
#               e.g. scripts/02_predict_cutoffs.sh wanda_aml2     glomeromycota
#                    scripts/02_predict_cutoffs.sh amv45nf_amdgr  endogonomycetes
#
# Note:         This script must be run from the project root directory,
#               after 01_prepare_reference.sh <primer_set> has completed.

set -eo pipefail

# =============================================================================
# ARGUMENTS AND CONFIG
# =============================================================================

PRIMER_SET="${1:?Usage: $0 <primer_set> <taxon_group>}"
TAXON_GROUP="${2:?Usage: $0 <primer_set> <taxon_group>}"

PRIMER_CONF="./config/primers/${PRIMER_SET}.conf"
TAXON_CONF="./config/taxa/${TAXON_GROUP}.conf"

for f in "$PRIMER_CONF" "$TAXON_CONF"; do
    if [[ ! -f "$f" ]]; then
        echo "ERROR: Config not found: $f" >&2
        exit 1
    fi
done

# shellcheck source=/dev/null
source "$PRIMER_CONF"
# shellcheck source=/dev/null
source "$TAXON_CONF"

readonly V4_FASTA="./data/${PRIMER_SET}/eukaryome_full.fasta"
readonly V4_CLASS="./data/${PRIMER_SET}/eukaryome_full.classification"
readonly COMPUTE_SIM="./R/compute_sim.R"
readonly PREDICT="./R/predict.R"
readonly OUT_DIR="./data/${PRIMER_SET}"
readonly TMP_DIR="./tmp"
readonly LOG_DIR="./logs/${PRIMER_SET}"
readonly END_THRESH=1
readonly STEP=0.001
readonly N_CPUS="${SLURM_CPUS_PER_TASK:-$(sysctl -n hw.ncpu 2>/dev/null || nproc)}"
readonly RUN_PARALLEL="yes"

readonly SIM_FILE="$OUT_DIR/${TAXON_PREFIX}_${SIM_BASE_RANK}.sim"

mkdir -p "$OUT_DIR" "$TMP_DIR" "$LOG_DIR"
LOG_FILE="$LOG_DIR/02_predict_${TAXON_PREFIX}_cutoffs_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG_FILE") 2>&1

echo "=============================================================================="
echo "PREDICTING CUTOFFS: $TAXON_LABEL, primer set = $PRIMER_SET_NAME"
echo "  Ranks:       ${RANKS[*]}"
echo "  MIN_SIM:     $MIN_SIM"
echo "  Cascade base: ${CASCADE_BASE_RANK:-<none -- direct from taxon filter>}"
echo "=============================================================================="

# Universal rank -> classification column map (fixed schema, not
# taxon-specific: id, kingdom, phylum, class, order, family, genus, species).
get_rank_col() {
    case "$1" in
        kingdom) echo 2 ;;
        phylum)  echo 3 ;;
        class)   echo 4 ;;
        order)   echo 5 ;;
        family)  echo 6 ;;
        genus)   echo 7 ;;
        species) echo 8 ;;
        *) echo "ERROR: unknown rank '$1'" >&2; exit 1 ;;
    esac
}

# Mirrors is_identified() from R/utils.R.
filter_identified_classification() {
    local in_class="$1"
    local rank_col="$2"
    local out_class="$3"

    awk -F'\t' -v col="$rank_col" '
        NR == 1 { print; next }
        {
            val = $col
            gsub(/^[[:space:]]+|[[:space:]]+$/, "", val)
            l = tolower(val)

            if (val == "" || l == "na") next
            if (l == "unidentified" || l == "unclassified" || l ~ /^uncultured/) next
            if (l ~ /incertae[ _.-]sedis/) next
            if (l ~ /[ _.-]sp\./) next
            if (l ~ /unispike1|unispike2|unispike3/) next
            if (l ~ /archaea|bacteria/) next
            if (l ~ /mitochondrion|nucleomorph|plastid/) next

            print
        }
    ' "$in_class" > "$out_class"
}

downsample_dominant() {
    local class_in="$1"
    local fasta_in="$2"
    local rank_col="$3"
    local max_prop="$4"
    local seed="$5"
    local class_out="$6"
    local fasta_out="$7"

    local n_total n_dominant dominant_taxon n_other prop exceeds n_final

    n_total=$(tail -n+2 "$class_in" | wc -l | tr -d ' ')
    dominant_taxon=$(tail -n+2 "$class_in" | awk -F'\t' -v col="$rank_col" '{print $col}' \
        | sort | uniq -c | sort -rn | head -1 | awk '{$1=""; print}' | xargs)
    n_dominant=$(tail -n+2 "$class_in" | awk -F'\t' -v col="$rank_col" -v tax="$dominant_taxon" \
        '$col == tax' | wc -l | tr -d ' ')

    prop=$(echo "scale=4; $n_dominant / $n_total" | bc)
    exceeds=$(echo "$prop > $max_prop" | bc)

    if [[ "$exceeds" -eq 0 ]]; then
        echo "  No downsampling needed: '$dominant_taxon' = ${n_dominant}/${n_total} (${prop})"
        cp "$class_in" "$class_out"
        cp "$fasta_in" "$fasta_out"
        return 1
    fi

    n_other=$((n_total - n_dominant))
    echo "  Dominant '$dominant_taxon' = ${n_dominant}/${n_total} (${prop}); downsampling to ${n_other}"

    tail -n+2 "$class_in" | awk -F'\t' -v col="$rank_col" -v tax="$dominant_taxon" '$col == tax {print $1}' > "$TMP_DIR/dominant_ids_all.txt"
    tail -n+2 "$class_in" | awk -F'\t' -v col="$rank_col" -v tax="$dominant_taxon" '$col != tax {print $1}' > "$TMP_DIR/other_ids.txt"

    Rscript -e "
        set.seed(${seed})
        ids <- readLines('${TMP_DIR}/dominant_ids_all.txt')
        sampled <- sample(ids, ${n_other})
        writeLines(sampled, '${TMP_DIR}/dominant_ids_sampled.txt')
"

    cat "$TMP_DIR/other_ids.txt" "$TMP_DIR/dominant_ids_sampled.txt" > "$TMP_DIR/downsampled_ids.txt"

    awk 'NR == FNR { ids[$1] = 1; next } FNR == 1 { print; next } $1 in ids { print }' \
        "$TMP_DIR/downsampled_ids.txt" "$class_in" > "$class_out"
    seqkit grep -f "$TMP_DIR/downsampled_ids.txt" "$fasta_in" > "$fasta_out"

    n_final=$(tail -n+2 "$class_out" | wc -l | tr -d ' ')
    echo "  Downsampled dataset size: ${n_final}"

    rm -f "$TMP_DIR/dominant_ids_all.txt" "$TMP_DIR/other_ids.txt" \
          "$TMP_DIR/dominant_ids_sampled.txt" "$TMP_DIR/downsampled_ids.txt"

    return 0
}

is_downsample_rank() {
    local r="$1"
    for d in "${DOWNSAMPLE_RANKS[@]}"; do
        [[ "$r" == "$d" ]] && return 0
    done
    return 1
}

echo "Activating conda environment..."
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate dyna_clust_predict

for f in "$V4_FASTA" "$V4_CLASS" "$COMPUTE_SIM" "$PREDICT"; do
    [[ -f "$f" ]] || { echo "ERROR: Required file not found: $f (did you run 01_prepare_reference.sh $PRIMER_SET first?)" >&2; exit 1; }
done

echo ""
echo "=== STEP 1: PREPARE ${TAXON_LABEL} SUBSETS ==="
echo "$(date)"

all_class="$TMP_DIR/${TAXON_PREFIX}_all.classification"
all_fasta="$TMP_DIR/${TAXON_PREFIX}_all.fasta"
all_ids="$TMP_DIR/${TAXON_PREFIX}_all_ids.txt"

awk -F'\t' -v col="$TAXON_FILTER_COL" -v val="$TAXON_FILTER_VALUE" \
    'NR == 1 || $col == val' "$V4_CLASS" > "$all_class"
tail -n+2 "$all_class" | cut -f1 > "$all_ids"
seqkit grep -f "$all_ids" "$V4_FASTA" > "$all_fasta"

all_n=$(tail -n+2 "$all_class" | wc -l | tr -d ' ')
echo "${TAXON_LABEL} sequences: $all_n"

# Optional cascade: derive an identified-base subset at CASCADE_BASE_RANK
# first, and filter every configured rank from THAT instead of from the raw
# taxon-filtered "all" set. Used for broad groups (e.g. Fungi) where
# requiring at least phylum-level ID as a floor makes sense before going
# finer; narrow groups (Glomeromycota, Endogonomycetes) leave this empty and
# filter every rank directly from "all".
base_class="$all_class"
base_fasta="$all_fasta"
base_n="$all_n"

if [[ -n "${CASCADE_BASE_RANK:-}" ]]; then
    echo ""
    echo "Building ${CASCADE_BASE_RANK}-identified base dataset"
    base_col=$(get_rank_col "$CASCADE_BASE_RANK")
    base_class="$OUT_DIR/${TAXON_PREFIX}_${CASCADE_BASE_RANK}.classification"
    base_fasta="$OUT_DIR/${TAXON_PREFIX}_${CASCADE_BASE_RANK}.fasta"
    base_ids="$TMP_DIR/${CASCADE_BASE_RANK}_ids.txt"

    filter_identified_classification "$all_class" "$base_col" "$base_class"
    tail -n+2 "$base_class" | cut -f1 > "$base_ids"
    seqkit grep -f "$base_ids" "$all_fasta" > "$base_fasta"

    base_n=$(tail -n+2 "$base_class" | wc -l | tr -d ' ')
    echo "${TAXON_LABEL} sequences (${CASCADE_BASE_RANK}-identified): $base_n"
fi

echo ""
echo "Filtering with is_identified() logic from R/utils.R"
printf "%-10s %8s %8s  %-50s\n" "Rank" "Input" "Kept" "Output files"
echo "---------------------------------------------------------------------------------"

for rank in "${RANKS[@]}"; do
    if [[ -n "${CASCADE_BASE_RANK:-}" && "$rank" == "$CASCADE_BASE_RANK" ]]; then
        # Already built above as the cascade base -- just report it.
        printf "%-10s %8s %8s  %s, %s\n" "$rank" "$all_n" "$base_n" \
            "$(basename "$base_fasta")" "$(basename "$base_class")"
        continue
    fi

    col=$(get_rank_col "$rank")
    out_fasta="$OUT_DIR/${TAXON_PREFIX}_${rank}.fasta"
    out_class="$OUT_DIR/${TAXON_PREFIX}_${rank}.classification"
    ids_file="$TMP_DIR/${rank}_ids.txt"

    filter_identified_classification "$base_class" "$col" "$out_class"
    tail -n+2 "$out_class" | cut -f1 > "$ids_file"
    seqkit grep -f "$ids_file" "$base_fasta" > "$out_fasta"

    n_kept=$(tail -n+2 "$out_class" | wc -l | tr -d ' ')
    printf "%-10s %8s %8s  %s, %s\n" "$rank" "$base_n" "$n_kept" \
        "$(basename "$out_fasta")" "$(basename "$out_class")"
done

echo ""
echo "=== STEP 2: COMPUTE MASTER SIM MATRIX ==="
echo "$(date)"
sim_base_fasta="$OUT_DIR/${TAXON_PREFIX}_${SIM_BASE_RANK}.fasta"

if [[ -f "$SIM_FILE" ]]; then
    echo "Sim file already exists, reusing: $SIM_FILE"
else
    Rscript "$COMPUTE_SIM" --input "$sim_base_fasta" --out "$OUT_DIR" --min_sim "$MIN_SIM" --n_cpus "$N_CPUS" --tmp_dir "$TMP_DIR"
    echo "Sim file written: $SIM_FILE"
fi

echo ""
echo "=== STEP 3: PREDICT GLOBAL CUTOFFS ==="
echo "$(date)"
echo "Sim matrix : $SIM_FILE"
echo "Ranks      : ${RANKS[*]}"

FAILED_RANKS=()
for rank in "${RANKS[@]}"; do
    rank_fasta="$OUT_DIR/${TAXON_PREFIX}_${rank}.fasta"
    rank_class="$OUT_DIR/${TAXON_PREFIX}_${rank}.classification"

    pred_fasta="$rank_fasta"
    pred_class="$rank_class"

    echo ""
    echo "--- GLOBAL: ${rank} ---"
    echo "Start threshold: $START_THRESH"

    if is_downsample_rank "$rank"; then
        ds_fasta="$TMP_DIR/${TAXON_PREFIX}_${rank}_downsampled.fasta"
        ds_class="$TMP_DIR/${TAXON_PREFIX}_${rank}_downsampled.classification"
        downsample_dominant "$rank_class" "$rank_fasta" "$(get_rank_col "$rank")" "$MAX_PROPORTION" "$SEED" "$ds_class" "$ds_fasta" || true
        pred_fasta="$ds_fasta"
        pred_class="$ds_class"
    fi

    Rscript "$PREDICT" \
        --input "$pred_fasta" \
        --classification "$pred_class" \
        --rank "$rank" \
        --sim "$SIM_FILE" \
        --start_threshold "$START_THRESH" \
        --end_threshold "$END_THRESH" \
        --step "$STEP" \
        --prefix "$TAXON_PREFIX" \
        --id_col id \
        --run_parallel "$RUN_PARALLEL" \
        --n_cpus "$N_CPUS" \
        --tmp_dir "$TMP_DIR" \
        --out "$OUT_DIR" || FAILED_RANKS+=("$rank")

done

echo ""
echo "=== CLEANUP ==="
rm -f "$TMP_DIR/${TAXON_PREFIX}_all_ids.txt" "$TMP_DIR"/*_ids.txt \
      "$TMP_DIR/${TAXON_PREFIX}_all.classification" "$TMP_DIR/${TAXON_PREFIX}_all.fasta" \
      "$TMP_DIR/${TAXON_PREFIX}"_*_downsampled.fasta "$TMP_DIR/${TAXON_PREFIX}"_*_downsampled.classification

echo ""
if [[ ${#FAILED_RANKS[@]} -gt 0 ]]; then
    echo "=== PIPELINE COMPLETED WITH ${#FAILED_RANKS[@]} FAILURE(S) ==="
    printf '  - %s\n' "${FAILED_RANKS[@]}"
    exit 1
else
    echo "=== PIPELINE COMPLETED SUCCESSFULLY: ${TAXON_LABEL} / ${PRIMER_SET_NAME} ==="
fi

echo "$(date)"
echo "Log file: $LOG_FILE"
echo ""

conda deactivate
