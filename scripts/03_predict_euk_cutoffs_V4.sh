#!/bin/bash

# Script name:  03_predict_euk_cutoffs_V4.sh
# Description:  Combined pipeline for Eukaryome V4 global cutoff prediction:
#               1. Generate only global ID files for kingdom/phylum/class/order.
#               2. Subset FASTA/classification on the fly per rank.
#               3. Enforce max sequences and dominant-taxon cap.
#               4. Compute similarity matrices and predict global cutoffs.
# Note:         This script must be run from the project root directory.

set -euo pipefail

readonly SUBSET="./R/subset.R"
readonly COMPUTE_SIM="./R/compute_sim.R"
readonly PREDICT="./R/predict.R"
readonly DATA_DIR="./data"
readonly V4_FASTA="$DATA_DIR/eukaryome_V4.fasta"
readonly V4_CLASS="$DATA_DIR/eukaryome_V4.classification"
readonly TMP_DIR="./tmp"
readonly LOG_DIR="./logs"
readonly PREFIX="eukaryome"
readonly STEP=0.001
readonly END_THRESH=1
readonly N_CPUS="${SLURM_CPUS_PER_TASK:-$(sysctl -n hw.ncpu 2>/dev/null || nproc)}"
readonly RUN_PARALLEL="yes"
readonly MAX_PROPORTION=0.5
readonly MAX_SEQUENCES_IDS=25000
readonly MAX_SEQUENCES_PRED=10000
readonly MIN_SUBGROUPS=10
readonly MIN_SEQUENCES=30
readonly SEED=1986

RANKS=("kingdom" "phylum" "class" "order")

mkdir -p "$TMP_DIR" "$LOG_DIR" "$DATA_DIR"
LOG_FILE="$LOG_DIR/03_predict_euk_cutoffs_V4_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG_FILE") 2>&1

get_start_thresh() {
    case "$1" in
        kingdom) echo "0.8" ;;
        phylum)  echo "0.8" ;;
        class)   echo "0.8" ;;
        order)   echo "0.8" ;;
    esac
}

get_rank_col() {
    case "$1" in
        kingdom) echo 2 ;;
        phylum)  echo 3 ;;
        class)   echo 4 ;;
        order)   echo 5 ;;
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

subset_fasta_and_classification() {
    local id_file="$1"
    local global_fasta="$2"
    local global_class="$3"
    local out_fasta="$4"
    local out_class="$5"

    seqkit grep -f "$id_file" "$global_fasta" > "$out_fasta"
    awk 'NR == FNR { ids[$1] = 1; next } FNR == 1 { print; next } $1 in ids { print }' "$id_file" "$global_class" > "$out_class"
}

limit_to_max_sequences() {
    local class_in="$1"
    local fasta_in="$2"
    local seed="$3"
    local class_out="$4"
    local fasta_out="$5"

    local n_total
    n_total=$(tail -n+2 "$class_in" | wc -l | tr -d ' ')

    if (( n_total <= MAX_SEQUENCES_PRED )); then
        cp "$class_in" "$class_out"
        cp "$fasta_in" "$fasta_out"
        echo "  Dataset size within limit: ${n_total}/${MAX_SEQUENCES_PRED}"
        return 1
    fi

    echo "  Dataset exceeds limit: ${n_total}/${MAX_SEQUENCES_PRED}; random downsampling"
    tail -n+2 "$class_in" | cut -f1 > "$TMP_DIR/all_ids.txt"

    Rscript -e "
        set.seed(${seed})
        ids <- readLines('${TMP_DIR}/all_ids.txt')
        sel <- sample(ids, ${MAX_SEQUENCES_PRED})
        writeLines(sel, '${TMP_DIR}/maxseq_ids.txt')
"

    awk 'NR == FNR { ids[$1] = 1; next } FNR == 1 { print; next } $1 in ids { print }' \
        "$TMP_DIR/maxseq_ids.txt" "$class_in" > "$class_out"
    seqkit grep -f "$TMP_DIR/maxseq_ids.txt" "$fasta_in" > "$fasta_out"

    rm -f "$TMP_DIR/all_ids.txt" "$TMP_DIR/maxseq_ids.txt"
    return 0
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
        echo "  No dominant-taxon downsampling needed: '$dominant_taxon' = ${n_dominant}/${n_total} (${prop})"
        cp "$class_in" "$class_out"
        cp "$fasta_in" "$fasta_out"
        return 1
    fi

    n_other=$((n_total - n_dominant))
    echo "  Dominant taxon '$dominant_taxon' exceeds cap: ${n_dominant}/${n_total} (${prop})"
    echo "  Downsampling dominant taxon from ${n_dominant} to ${n_other}"

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

for f in "$SUBSET" "$COMPUTE_SIM" "$PREDICT" "$V4_FASTA" "$V4_CLASS"; do
    [[ -f "$f" ]] || { echo "ERROR: Required file not found: $f" >&2; exit 1; }
done

echo ""
echo "=== STEP 1: PREPARE REQUIRED GLOBAL ID FILES ==="
echo "$(date)"

tmp_subset_dir=$(mktemp -d "$TMP_DIR/euk_subset_ids.XXXXXX")
Rscript "$SUBSET" \
    --fasta_in "$V4_FASTA" \
    --classification_in "$V4_CLASS" \
    --output_dir "$tmp_subset_dir" \
    --min_subgroups "$MIN_SUBGROUPS" \
    --min_sequences "$MIN_SEQUENCES" \
    --max_sequences "$MAX_SEQUENCES_IDS" \
    --max_proportion "$MAX_PROPORTION"

echo "Copying required global ID files only..."
for rank in "${RANKS[@]}"; do
    src="$tmp_subset_dir/${rank}_pred_id_global.txt"
    dst="$DATA_DIR/${rank}_pred_id_global.txt"
    [[ -f "$src" ]] || { echo "ERROR: Missing staged ID file: $src" >&2; rm -rf "$tmp_subset_dir"; exit 1; }
    cp "$src" "$dst"
    echo "  $dst ($(wc -l < "$dst" | tr -d ' ') IDs)"
done
rm -rf "$tmp_subset_dir"

echo ""
echo "=== STEP 2: PREDICT GLOBAL CUTOFFS ==="
echo "$(date)"
echo "max_seq_no (ID generation) : $MAX_SEQUENCES_IDS"
echo "max_seq_no (prediction)    : $MAX_SEQUENCES_PRED"

FAILED_RANKS=()
for rank in "${RANKS[@]}"; do
    rank_col=$(get_rank_col "$rank")
    st=$(get_start_thresh "$rank")
    id_file="$DATA_DIR/${rank}_pred_id_global.txt"

    subset_fasta="$TMP_DIR/eukaryome_${rank}_global.fasta"
    subset_class="$TMP_DIR/eukaryome_${rank}_global.classification"
    maxseq_fasta="$TMP_DIR/eukaryome_${rank}_global_maxseq.fasta"
    maxseq_class="$TMP_DIR/eukaryome_${rank}_global_maxseq.classification"
    pred_fasta="$TMP_DIR/eukaryome_${rank}_global_pred.fasta"
    pred_class="$TMP_DIR/eukaryome_${rank}_global_pred.classification"
    sim_out="$TMP_DIR/eukaryome_${rank}_global_pred.sim"

    echo ""
    echo "--- GLOBAL: ${rank} ---"
    echo "  Start threshold: $st"

    subset_fasta_and_classification "$id_file" "$V4_FASTA" "$V4_CLASS" "$subset_fasta" "$subset_class"
    echo "  Initial subset size: $(( $(wc -l < "$subset_class") - 1 ))"

    limit_to_max_sequences "$subset_class" "$subset_fasta" "$SEED" "$maxseq_class" "$maxseq_fasta" || true
    downsample_dominant "$maxseq_class" "$maxseq_fasta" "$rank_col" "$MAX_PROPORTION" "$SEED" "$pred_class" "$pred_fasta" || true
    echo "  Final prediction size: $(( $(wc -l < "$pred_class") - 1 ))"

    Rscript "$COMPUTE_SIM" --input "$pred_fasta" --out "$TMP_DIR" --min_sim 0 --n_cpus "$N_CPUS" --tmp_dir "$TMP_DIR"

    generated_sim="$TMP_DIR/$(basename "${pred_fasta%.*}").sim"
    [[ -f "$generated_sim" ]] || { echo "  FAILED: Expected sim file missing: $generated_sim" >&2; FAILED_RANKS+=("$rank"); continue; }
    mv "$generated_sim" "$sim_out"

    Rscript "$PREDICT" \
        --input "$pred_fasta" \
        --classification "$pred_class" \
        --rank "$rank" \
        --sim "$sim_out" \
        --start_threshold "$st" \
        --end_threshold "$END_THRESH" \
        --step "$STEP" \
        --prefix "$PREFIX" \
        --id_col id \
        --run_parallel "$RUN_PARALLEL" \
        --n_cpus "$N_CPUS" \
        --tmp_dir "$TMP_DIR" \
        --out "$DATA_DIR" || FAILED_RANKS+=("$rank")

    rm -f "$subset_fasta" "$subset_class" "$maxseq_fasta" "$maxseq_class" "$pred_fasta" "$pred_class" "$sim_out"
done

echo ""
echo "=== CLEANUP ==="
rm -rf "$TMP_DIR"

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
