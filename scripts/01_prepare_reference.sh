#!/bin/bash

# Script name:  01_prepare_reference.sh
# Description:  Download the EUKARYOME General EUK SSU v2.0 database (once,
#               shared across all primer sets), trim reads with the given
#               primer pair to produce full-length reference files, then
#               generate a degenerate forward primer from the trimmed reads
#               and apply it to untrimmed reads (primer retained, truncated)
#               for downstream prediction.
#
# Usage:        scripts/01_prepare_reference.sh <primer_set>
#               e.g. scripts/01_prepare_reference.sh wanda_aml2
#                    scripts/01_prepare_reference.sh amv45nf_amdgr
#
# Note:         This script must be run from the project root directory.
#               Primer-specific parameters live in config/primers/<primer_set>.conf
#               -- this script itself has no primer-specific values.

set -eo pipefail

# =============================================================================
# ARGUMENTS AND CONFIG
# =============================================================================

PRIMER_SET="${1:?Usage: $0 <primer_set>  (e.g. wanda_aml2, amv45nf_amdgr)}"
PRIMER_CONF="./config/primers/${PRIMER_SET}.conf"

if [[ ! -f "$PRIMER_CONF" ]]; then
    echo "ERROR: Primer config not found: $PRIMER_CONF" >&2
    echo "Available primer sets:" >&2
    ls config/primers/*.conf 2>/dev/null | xargs -n1 basename | sed 's/\.conf$//' | sed 's/^/  - /' >&2
    exit 1
fi

# shellcheck source=/dev/null
source "$PRIMER_CONF"

readonly REF_SEQS_DIR="./data/ref_seqs"

# EUKARYOME download URL (EUKARYOME General EUK SSU v2.0) -- shared, not
# primer-specific, so it's cached once in data/ref_seqs/ regardless of which
# primer set is being prepared.
readonly EUKARYOME_URL="https://sisu.ut.ee/wp-content/uploads/sites/643/General_EUK_SSU_v2.0.zip"
readonly DOWNLOAD_FILE="$REF_SEQS_DIR/General_EUK_SSU_v2.0.zip"
readonly IN_FASTA="$REF_SEQS_DIR/General_EUK_SSU_v2.0.fasta"

# Output files (per primer set)
readonly OUT_DIR="./data/${PRIMER_SET}"
readonly TMP_DIR="./tmp"
readonly OUT_FASTA="$OUT_DIR/eukaryome_full.fasta"
readonly OUT_CLASS="$OUT_DIR/eukaryome_full.classification"
readonly OUT_FASTA_PARTIAL="$OUT_DIR/eukaryome_partial.fasta"
readonly OUT_CLASS_PARTIAL="$OUT_DIR/eukaryome_partial.classification"
readonly OUT_FASTA_COMBINED="$OUT_DIR/eukaryome_all.fasta"
readonly OUT_CLASS_COMBINED="$OUT_DIR/eukaryome_all.classification"

# Fixed cutadapt/pipeline parameters (not primer-specific)
readonly LOG_FILE="$TMP_DIR/logfile_cutadapt_eukaryome.txt"
readonly NUM_THREADS=10

# Helper R scripts (shared, not primer-specific)
readonly REFORMAT="./R/reformat.R"
readonly DEREP_LCA="./R/dereplicate_lca.R"
readonly CHECK_ANNOTATIONS="./R/check_annotations.R"
readonly COMBINE_LCA="./R/combine_lca.R"

echo "=============================================================================="
echo "PREPARING REFERENCE: primer set = $PRIMER_SET_NAME"
echo "  Forward primer: $PRIMER_FWD"
echo "  Reverse primer: $PRIMER_REV"
echo "  Length window:  ${MIN_LEN}-${MAX_LEN}bp"
echo "  Output dir:     $OUT_DIR"
echo "=============================================================================="

# =============================================================================
# DIRECTORY SETUP
# =============================================================================

mkdir -p "$REF_SEQS_DIR" "$OUT_DIR" "$TMP_DIR" "./logs/${PRIMER_SET}"

# =============================================================================
# ENVIRONMENT SETUP
# =============================================================================

echo "Activating conda environment..."
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate dyna_clust_predict

# =============================================================================
# INPUT VALIDATION
# =============================================================================

for f in "$REFORMAT" "$DEREP_LCA" "$CHECK_ANNOTATIONS" "$COMBINE_LCA"; do
    if [[ ! -f "$f" ]]; then
        echo "ERROR: R script not found: $f" >&2
        exit 1
    fi
done

# =============================================================================
# FILE DOWNLOAD: EUKARYOME (shared across primer sets)
# =============================================================================

echo ""
echo "=== DOWNLOADING EUKARYOME DATABASE (shared, once) ==="
echo "$(date)"

if [[ -f "$DOWNLOAD_FILE" && -f "$IN_FASTA" ]]; then
    echo "Files already exist in $REF_SEQS_DIR, skipping download and extraction."
else
    if [[ ! -f "$DOWNLOAD_FILE" ]]; then
        echo "URL: $EUKARYOME_URL"
        if ! curl -L -o "$DOWNLOAD_FILE" "$EUKARYOME_URL"; then
            echo "ERROR: Failed to download file from $EUKARYOME_URL" >&2
            exit 1
        fi
    else
        echo "Download file already exists, skipping download."
    fi

    if [[ ! -f "$IN_FASTA" ]]; then
        echo "Unzipping downloaded file..."
        if ! 7z x "$DOWNLOAD_FILE" -o"$REF_SEQS_DIR/" -y; then
            echo "ERROR: Failed to unzip $DOWNLOAD_FILE" >&2
            exit 1
        fi
    else
        echo "Extracted file already exists, skipping extraction."
    fi

    echo "Download and extraction completed successfully."
fi
echo ""

# =============================================================================
# PRIMER TRIMMING: EUKARYOME
# =============================================================================

if [[ ! -f "$IN_FASTA" ]]; then
    echo "ERROR: Input file $IN_FASTA not found!" >&2
    exit 1
fi

echo "=== PRIMER TRIMMING WITH $PRIMER_SET_NAME ==="
echo "Input file:     $IN_FASTA"
echo "Forward primer: $PRIMER_FWD (min overlap: $MIN_OVERLAP_FWD, max error rate: $MAX_ERROR_FWD)"
echo "Reverse primer: $PRIMER_REV (min overlap: $MIN_OVERLAP_REV, max error rate: $MAX_ERROR_REV)"
echo ""

TRIMMED_RAW="$TMP_DIR/eukaryome_trimmed_raw.fasta"
TRIMMED="$TMP_DIR/eukaryome_trimmed.fasta"
UNTRIMMED="$TMP_DIR/eukaryome_untrimmed.fasta"

cutadapt \
  -g "$PRIMER_FWD;min_overlap=$MIN_OVERLAP_FWD;max_error_rate=$MAX_ERROR_FWD;required"..."$PRIMER_REV;min_overlap=$MIN_OVERLAP_REV;max_error_rate=$MAX_ERROR_REV;required" \
  --cores "$NUM_THREADS" \
  --untrimmed-output "$UNTRIMMED" \
  -o "$TRIMMED_RAW" "$IN_FASTA" \
  >> "$LOG_FILE" 2>&1 < /dev/null

echo "Primer trimming completed."
if [[ -f "$LOG_FILE" ]]; then
    echo ""
    echo "=== TRIMMING SUMMARY ==="
    grep -E "(Total reads processed|Reads with adapters|Reads that were too short|Reads that were too long|Reads discarded as untrimmed|Reads written)" "$LOG_FILE" | tail -6
fi
echo ""

# =============================================================================
# LENGTH FILTER AND N REMOVAL
# =============================================================================

echo "=== LENGTH FILTER AND N REMOVAL ==="
echo "Length range: ${MIN_LEN}-${MAX_LEN}bp"
echo "Max N bases:  $MAX_N"
echo ""

LOG_FILTER="$TMP_DIR/logfile_cutadapt_filter.txt"
TOO_LONG_TRIMMED="$TMP_DIR/too_long_trimmed.fasta"

cutadapt \
  -m "$MIN_LEN" -M "$MAX_LEN" \
  --max-n "$MAX_N" \
  --cores "$NUM_THREADS" \
  --too-long-output "$TOO_LONG_TRIMMED" \
  -o "$TRIMMED" "$TRIMMED_RAW" \
  >> "$LOG_FILTER" 2>&1 < /dev/null

echo "Length/N filtering completed."
if [[ -f "$LOG_FILTER" ]]; then
    grep -E "(Total reads processed|Reads that were too short|Reads that were too long|Reads written)" "$LOG_FILTER" | tail -4
fi
n_too_long=$(grep -c "^>" "$TOO_LONG_TRIMMED" 2>/dev/null || true)
echo "Reads too long (rejoined to partial): ${n_too_long:-0}"
echo ""

# =============================================================================
# REFORMAT HEADERS
# =============================================================================

echo "=== REFORMATTING HEADERS ==="
echo "$(date)"

Rscript "$REFORMAT" \
    --fasta_in           "$TRIMMED" \
    --fasta_out          "$TMP_DIR/reformatted.fasta" \
    --classification_out "$TMP_DIR/reformatted.classification"

echo ""
echo "Reformatted FASTA written to     : $TMP_DIR/reformatted.fasta"
echo "Reformatted classification table : $TMP_DIR/reformatted.classification"
echo ""

# =============================================================================
# CHECK AND STANDARDISE ANNOTATIONS
# =============================================================================

echo "=== CHECKING AND STANDARDISING ANNOTATIONS ==="
echo "$(date)"

Rscript "$CHECK_ANNOTATIONS" \
    --classification_in  "$TMP_DIR/reformatted.classification" \
    --classification_out "$TMP_DIR/reformatted.classification"
echo ""

# =============================================================================
# DEREPLICATION WITH LCA ASSIGNMENT
# =============================================================================

echo "=== DEREPLICATE AND RESOLVE LCA TAXONOMY ==="
echo "$(date)"

Rscript "$DEREP_LCA" \
    --fasta_in             "$TMP_DIR/reformatted.fasta" \
    --fasta_out            "$OUT_FASTA" \
    --classification_in    "$TMP_DIR/reformatted.classification" \
    --classification_out   "$OUT_CLASS"

echo ""
echo "Full FASTA written to           : $OUT_FASTA"
echo "Full classification written to  : $OUT_CLASS"
echo ""

# =============================================================================
# GENERATE DEGENERATE FORWARD PRIMER FROM TRIMMED READS
# =============================================================================

echo "=== GENERATING DEGENERATE FORWARD PRIMER ==="
echo "$(date)"
echo ""

CONSENSUS_FWD=$(Rscript -e "
    library(Biostrings)
    seqs <- readDNAStringSet('${TRIMMED}')
    seqs <- seqs[width(seqs) >= ${FWD_CONSENSUS_LEN}]
    prefixes <- subseq(seqs, 1, ${FWD_CONSENSUS_LEN})
    cm <- consensusMatrix(prefixes)
    cons <- consensusString(cm, ambiguityMap = IUPAC_CODE_MAP, threshold = ${FWD_THRESHOLD})
    cat(toupper(cons))
" 2>/dev/null)

echo "Degenerate forward primer (${FWD_CONSENSUS_LEN}bp, threshold=${FWD_THRESHOLD}): $CONSENSUS_FWD"
echo ""

if [[ -z "$CONSENSUS_FWD" ]]; then
    echo "ERROR: Failed to generate degenerate forward primer." >&2
    exit 1
fi

# =============================================================================
# FILTER UNTRIMMED READS WITH DEGENERATE FORWARD PRIMER
# =============================================================================

echo "=== FILTERING UNTRIMMED READS WITH DEGENERATE FORWARD PRIMER ==="
echo "$(date)"
echo ""
echo "Input:          $UNTRIMMED"
echo "Forward primer: $CONSENSUS_FWD (${FWD_CONSENSUS_LEN}bp)"
echo "Action:         retain (primer kept in output)"
echo "Truncation:     ${TRUNCATION_LEN}bp"
echo ""

FWD_FILTERED_RAW="$TMP_DIR/eukaryome_fwd_filtered_raw.fasta"
FWD_FILTERED="$TMP_DIR/eukaryome_fwd_filtered.fasta"
LOG_FWD="$TMP_DIR/logfile_cutadapt_fwd_filter.txt"

cutadapt \
  -g "${CONSENSUS_FWD};min_overlap=${FWD_CONSENSUS_LEN};max_error_rate=1" \
  --action=retain \
  --discard-untrimmed \
  --length "$TRUNCATION_LEN" \
  --cores "$NUM_THREADS" \
  -o "$FWD_FILTERED_RAW" "$UNTRIMMED" \
  >> "$LOG_FWD" 2>&1 < /dev/null

echo "Forward primer filtering completed."
if [[ -f "$LOG_FWD" ]]; then
    echo ""
    echo "=== FORWARD PRIMER FILTER SUMMARY ==="
    grep -E "(Total reads processed|Reads with adapters|Reads discarded as untrimmed|Reads written)" "$LOG_FWD" | tail -4
fi

n_filtered=$(grep -c "^>" "$FWD_FILTERED_RAW" 2>/dev/null || true)
echo ""
echo "Reads passing forward primer filter: ${n_filtered:-0}"
echo ""

# =============================================================================
# JOIN TOO-LONG FULL READS TO PARTIAL DATASET
# =============================================================================

echo "=== JOINING TOO-LONG FULL READS TO PARTIAL DATASET ==="
echo "Truncating to ${TRUNCATION_LEN}bp and appending to partial reads..."
echo ""

TOO_LONG_TRUNCATED="$TMP_DIR/too_long_truncated.fasta"

cutadapt \
  --length "$TRUNCATION_LEN" \
  --cores "$NUM_THREADS" \
  -o "$TOO_LONG_TRUNCATED" "$TOO_LONG_TRIMMED" \
  > /dev/null 2>&1

cat "$TOO_LONG_TRUNCATED" >> "$FWD_FILTERED_RAW"

n_combined=$(grep -c "^>" "$FWD_FILTERED_RAW" 2>/dev/null || true)
echo "Partial reads after join: ${n_combined:-0}"
echo ""

# =============================================================================
# LENGTH FILTER AND N REMOVAL (partial)
# =============================================================================

echo "=== LENGTH FILTER AND N REMOVAL (partial) ==="
echo "Min length: ${MIN_LEN}bp"
echo "Max N bases: 0"
echo ""

vsearch \
  --fastx_filter "$FWD_FILTERED_RAW" \
  --fastaout "$FWD_FILTERED" \
  --fastq_maxns 0 \
  --fastq_minlen "$MIN_LEN" \
  --fasta_width 0 \
  2> "$TMP_DIR/logfile_vsearch_filter_partial.txt"

n_pass=$(grep -c "^>" "$FWD_FILTERED" 2>/dev/null || true)
n_removed=$(( ${n_filtered:-0} - ${n_pass:-0} ))
echo "Sequences passing filter: ${n_pass:-0}"
echo "Sequences removed:        ${n_removed}"
echo "Output written to: $FWD_FILTERED"
echo ""

# =============================================================================
# REFORMAT HEADERS (partial)
# =============================================================================

echo "=== REFORMATTING HEADERS (partial) ==="
echo "$(date)"

Rscript "$REFORMAT" \
    --fasta_in           "$FWD_FILTERED" \
    --fasta_out          "$TMP_DIR/reformatted_partial.fasta" \
    --classification_out "$TMP_DIR/reformatted_partial.classification"
echo ""

# =============================================================================
# CHECK AND STANDARDISE ANNOTATIONS (partial)
# =============================================================================

echo "=== CHECKING AND STANDARDISING ANNOTATIONS (partial) ==="
echo "$(date)"

Rscript "$CHECK_ANNOTATIONS" \
    --classification_in  "$TMP_DIR/reformatted_partial.classification" \
    --classification_out "$TMP_DIR/reformatted_partial.classification"
echo ""

# =============================================================================
# DEREPLICATION WITH LCA ASSIGNMENT (partial)
# =============================================================================

echo "=== DEREPLICATE AND RESOLVE LCA TAXONOMY (partial) ==="
echo "$(date)"

Rscript "$DEREP_LCA" \
    --fasta_in             "$TMP_DIR/reformatted_partial.fasta" \
    --fasta_out            "$OUT_FASTA_PARTIAL" \
    --classification_in    "$TMP_DIR/reformatted_partial.classification" \
    --classification_out   "$OUT_CLASS_PARTIAL"

echo ""
echo "Partial FASTA written to           : $OUT_FASTA_PARTIAL"
echo "Partial classification written to  : $OUT_CLASS_PARTIAL"
echo ""

# =============================================================================
# OVERLAP ANALYSIS: 100% IDENTITY BETWEEN FULL AND PARTIAL SEQUENCES
# =============================================================================

echo "=== OVERLAP ANALYSIS: full vs partial (100% identity) ==="
echo "$(date)"
echo "Query:    $OUT_FASTA_PARTIAL"
echo "Database: $OUT_FASTA"
echo ""

vsearch \
  --usearch_global "$OUT_FASTA_PARTIAL" \
  --db "$OUT_FASTA" \
  --id 1.0 \
  --strand plus \
  --maxaccepts 100 \
  --maxrejects 100 \
  --matched    "$TMP_DIR/overlap_partial_in_full.fasta" \
  --notmatched "$TMP_DIR/unique_partial.fasta" \
  --userout    "$TMP_DIR/usearch_overlap.txt" \
  --userfields query+target+id+ql+tl+alnlen+qcov+tcov+ids \
  --threads "$NUM_THREADS" \
  2> "$TMP_DIR/logfile_usearch_overlap.txt"

n_matched=$(grep -c "^>" "$TMP_DIR/overlap_partial_in_full.fasta" 2>/dev/null || true)
n_unique=$(grep -c  "^>" "$TMP_DIR/unique_partial.fasta"          2>/dev/null || true)
echo "Partial sequences with 100% identity match in full : ${n_matched:-0}"
echo "Partial sequences unique to the partial set          : ${n_unique:-0}"
echo ""

# =============================================================================
# COMBINE FULL AND PARTIAL WITH LCA FOR OVERLAPPING SEQUENCES
# =============================================================================

echo "=== COMBINING FULL AND PARTIAL SEQUENCES ==="
echo "$(date)"
echo ""

Rscript "$COMBINE_LCA" \
    --fasta_full    "$OUT_FASTA" \
    --class_full    "$OUT_CLASS" \
    --fasta_partial "$OUT_FASTA_PARTIAL" \
    --class_partial "$OUT_CLASS_PARTIAL" \
    --usearch_hits  "$TMP_DIR/usearch_overlap.txt" \
    --fasta_out     "$OUT_FASTA_COMBINED" \
    --class_out     "$OUT_CLASS_COMBINED"

echo ""
echo "Combined FASTA written to:          $OUT_FASTA_COMBINED"
echo "Combined classification written to: $OUT_CLASS_COMBINED"
echo ""

# =============================================================================
# CLEANUP
# =============================================================================

echo "=== CLEANUP ==="
echo "Removing $TMP_DIR contents..."
rm -rf "${TMP_DIR:?}"/*

echo ""
echo "=== REFERENCE PREPARATION COMPLETED SUCCESSFULLY: $PRIMER_SET_NAME ==="
echo "$(date)"

conda deactivate
