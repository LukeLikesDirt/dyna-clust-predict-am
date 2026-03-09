#!/bin/bash

# Script name:  02_predict_glom_cutoffs_V4.sh
# Description:  Combined pipeline for Glomeromycota V4 cutoff prediction:
#               1. Build family/genus/species subsets from eukaryome_V4.
#               2. Compute family master similarity matrix.
#               3. Predict global cutoffs for family/genus/species.
# Note:         This script must be run from the project root directory.

set -euo pipefail

readonly V4_FASTA="./data/eukaryome_V4.fasta"
readonly V4_CLASS="./data/eukaryome_V4.classification"
readonly COMPUTE_SIM="./R/compute_sim.R"
readonly PREDICT="./R/predict.R"
readonly OUT_DIR="./data"
readonly TMP_DIR="./tmp"
readonly LOG_DIR="./logs"
readonly PREFIX="glomeromycota"
readonly STEP=0.001
readonly END_THRESH=1
readonly N_CPUS="${SLURM_CPUS_PER_TASK:-$(sysctl -n hw.ncpu 2>/dev/null || nproc)}"
readonly RUN_PARALLEL="yes"
readonly MAX_PROPORTION=0.5
readonly SEED=1986

readonly SIM_FILE="$OUT_DIR/glomeromycota_family_V4.sim"
RANKS=("family" "genus" "species")

mkdir -p "$OUT_DIR" "$TMP_DIR" "$LOG_DIR"
LOG_FILE="$LOG_DIR/02_predict_glom_cutoffs_V4_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG_FILE") 2>&1

get_start_thresh() {
    case "$1" in
        family)  echo "0.9" ;;
        genus)   echo "0.9" ;;
        species) echo "0.9" ;;
    esac
}

get_rank_col() {
    case "$1" in
        family)  echo 6 ;;
        genus)   echo 7 ;;
        species) echo 8 ;;
    esac
}

activate_conda() {
    set +u
    if command -v conda >/dev/null 2>&1; then
        eval "$(conda shell.bash hook)"
    elif [[ -f "$HOME/.bash_profile" ]]; then
        source "$HOME/.bash_profile"
    fi
    conda activate dyna_clust_predict
    set -u
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

echo "Activating conda environment..."
activate_conda

for f in "$V4_FASTA" "$V4_CLASS" "$COMPUTE_SIM" "$PREDICT"; do
    [[ -f "$f" ]] || { echo "ERROR: Required file not found: $f" >&2; exit 1; }
done

echo ""
echo "=== STEP 1: PREPARE GLOMEROMYCOTA SUBSETS ==="
echo "$(date)"

glom_class="$TMP_DIR/glomeromycota_all.classification"
glom_fasta="$TMP_DIR/glomeromycota_all.fasta"

aWK_tmp_ids="$TMP_DIR/glom_ids.txt"
awk -F'\t' 'NR == 1 || $3 == "Glomeromycota"' "$V4_CLASS" > "$glom_class"
tail -n+2 "$glom_class" | cut -f1 > "$aWK_tmp_ids"
seqkit grep -f "$aWK_tmp_ids" "$V4_FASTA" > "$glom_fasta"

glom_n=$(tail -n+2 "$glom_class" | wc -l | tr -d ' ')
echo "Glomeromycota sequences: $glom_n"

echo ""
echo "Filtering with is_identified() logic from R/utils.R"
printf "%-10s %8s %8s  %-50s\n" "Rank" "Input" "Kept" "Output files"
echo "---------------------------------------------------------------------------------"

for rank in "${RANKS[@]}"; do
    col=$(get_rank_col "$rank")
    out_fasta="$OUT_DIR/glomeromycota_${rank}_V4.fasta"
    out_class="$OUT_DIR/glomeromycota_${rank}_V4.classification"
    ids_file="$TMP_DIR/${rank}_ids.txt"

    filter_identified_classification "$glom_class" "$col" "$out_class"
    tail -n+2 "$out_class" | cut -f1 > "$ids_file"
    seqkit grep -f "$ids_file" "$glom_fasta" > "$out_fasta"

    n_kept=$(tail -n+2 "$out_class" | wc -l | tr -d ' ')
    printf "%-10s %8s %8s  %s, %s\n" "$rank" "$glom_n" "$n_kept" \
        "$(basename "$out_fasta")" "$(basename "$out_class")"
done

echo ""
echo "=== STEP 2: COMPUTE MASTER SIM MATRIX ==="
echo "$(date)"
family_fasta="$OUT_DIR/glomeromycota_family_V4.fasta"

if [[ -f "$SIM_FILE" ]]; then
    echo "Sim file already exists, reusing: $SIM_FILE"
else
    Rscript "$COMPUTE_SIM" --input "$family_fasta" --out "$OUT_DIR" --min_sim 0 --n_cpus "$N_CPUS" --tmp_dir "$TMP_DIR"
    echo "Sim file written: $SIM_FILE"
fi

echo ""
echo "=== STEP 3: PREDICT GLOBAL CUTOFFS ==="
echo "$(date)"
echo "Sim matrix : $SIM_FILE"
echo "Ranks      : ${RANKS[*]}"

FAILED_RANKS=()
for rank in "${RANKS[@]}"; do
    rank_fasta="$OUT_DIR/glomeromycota_${rank}_V4.fasta"
    rank_class="$OUT_DIR/glomeromycota_${rank}_V4.classification"
    st=$(get_start_thresh "$rank")

    pred_fasta="$rank_fasta"
    pred_class="$rank_class"

    echo ""
    echo "--- GLOBAL: ${rank} ---"
    echo "Start threshold: $st"

    if [[ "$rank" == "family" ]]; then
        ds_fasta="$TMP_DIR/glomeromycota_family_downsampled.fasta"
        ds_class="$TMP_DIR/glomeromycota_family_downsampled.classification"
        downsample_dominant "$rank_class" "$rank_fasta" "$(get_rank_col "$rank")" "$MAX_PROPORTION" "$SEED" "$ds_class" "$ds_fasta" || true
        pred_fasta="$ds_fasta"
        pred_class="$ds_class"
    fi

    Rscript "$PREDICT" \
        --input "$pred_fasta" \
        --classification "$pred_class" \
        --rank "$rank" \
        --sim "$SIM_FILE" \
        --start_threshold "$st" \
        --end_threshold "$END_THRESH" \
        --step "$STEP" \
        --prefix "$PREFIX" \
        --id_col id \
        --run_parallel "$RUN_PARALLEL" \
        --n_cpus "$N_CPUS" \
        --tmp_dir "$TMP_DIR" \
        --out "$OUT_DIR" || FAILED_RANKS+=("$rank")

done

echo ""
echo "=== CLEANUP ==="
rm -f "$TMP_DIR/glom_ids.txt" "$TMP_DIR"/*_ids.txt \
      "$TMP_DIR/glomeromycota_all.classification" "$TMP_DIR/glomeromycota_all.fasta" \
      "$TMP_DIR/glomeromycota_family_downsampled.fasta" "$TMP_DIR/glomeromycota_family_downsampled.classification"

echo ""
if [[ ${#FAILED_RANKS[@]} -gt 0 ]]; then
    echo "=== PIPELINE COMPLETED WITH ${#FAILED_RANKS[@]} FAILURE(S) ==="
    printf '  - %s\n' "${FAILED_RANKS[@]}"
    exit 1
else
    echo "=== PIPELINE COMPLETED SUCCESSFULLY ==="
fi

echo "$(date)"
echo "Log file: $LOG_FILE"

echo ""
set +u
conda deactivate
