#!/bin/bash

# Script name:  01_new_2.sh
# Description:  Download the EUKARYOME General EUK SSU v2.0 database, trim
#               reads with WANDA/AML2 to produce eukaryome_V4_predict reference
#               files, then generate a degenerate forward primer from the
#               trimmed reads and apply it to untrimmed reads (primer retained,
#               truncated to 500bp) for downstream prediction.
# Note:         This script must be run from the project root directory.

set -eo pipefail

# =============================================================================
# PARAMETER SETUP
# =============================================================================

# EUKARYOME download URL (EUKARYOME General EUK SSU v2.0)
readonly EUKARYOME_URL="https://sisu.ut.ee/wp-content/uploads/sites/643/General_EUK_SSU_v2.0.zip"
readonly DOWNLOAD_FILE="./tmp/General_EUK_SSU_v2.0.zip"
readonly IN_FASTA="./tmp/General_EUK_SSU_v2.0.fasta"

# Output files
readonly OUT_FASTA="./data/eukaryome_V4_full.fasta"
readonly OUT_CLASS="./data/eukaryome_V4_full.classification"
readonly OUT_FASTA_PARTIAL="./data/eukaryome_V4_partial.fasta"
readonly OUT_CLASS_PARTIAL="./data/eukaryome_V4_partial.classification"
readonly OUT_FASTA_COMBINED="./data/eukaryome_V4_all.fasta"
readonly OUT_CLASS_COMBINED="./data/eukaryome_V4_all.classification"

# CUTADAPT PARAMETERS
readonly LOG_FILE="tmp/logfile_cutadapt_eukaryome.txt"
readonly NUM_THREADS=10
readonly MIN_LEN=440
readonly MAX_LEN=540
readonly MAX_N=0

# WANDA/AML2 trimming parameters
readonly MIN_OVERLAP_FWD=20
readonly MIN_OVERLAP_REV=22
readonly MAX_ERROR_FWD=2
readonly MAX_ERROR_REV=4

# WANDA 5' forward primer and AML2 3' reverse primer (reverse complement)
readonly PRIMER_FWD="CAGCCGCGGTAATTCCAGCT"
readonly PRIMER_REV="GGAAACCAAAGTGTTTGGGTTC"

# Degenerate forward primer parameters
readonly FWD_CONSENSUS_LEN=15
readonly FWD_THRESHOLD=0.25

# Truncation length for untrimmed read filtering
readonly TRUNCATION_LEN=540

# Helper R scripts
readonly REFORMAT="./R/reformat.R"
readonly DEREP_LCA="./R/dereplicate_lca.R"
readonly CHECK_ANNOTATIONS="./R/check_annotations.R"

# =============================================================================
# DIRECTORY SETUP
# =============================================================================

mkdir -p ./data/ tmp

# =============================================================================
# ENVIRONMENT SETUP
# =============================================================================

echo "Activating conda environment..."
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate dyna_clust_predict

# =============================================================================
# INPUT VALIDATION
# =============================================================================

if [[ ! -f "$REFORMAT" ]]; then
    echo "ERROR: R script not found: $REFORMAT" >&2
    exit 1
fi

if [[ ! -f "$DEREP_LCA" ]]; then
    echo "ERROR: R script not found: $DEREP_LCA" >&2
    exit 1
fi

# =============================================================================
# FILE DOWNLOAD: EUKARYOME
# =============================================================================

echo ""
echo "=== DOWNLOADING EUKARYOME DATABASE ==="
echo "$(date)"

if [[ -f "$DOWNLOAD_FILE" && -f "$IN_FASTA" ]]; then
    echo "Files already exist, skipping download and extraction."
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
        if ! 7z x "$DOWNLOAD_FILE" -o"./tmp/" -y; then
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
# PRIMER TRIMMING: EUKARYOME (WANDA / AML2)
# =============================================================================

if [[ ! -f "$IN_FASTA" ]]; then
    echo "ERROR: Input file $IN_FASTA not found!" >&2
    exit 1
fi

echo "=== PRIMER TRIMMING WITH WANDA / AML2 ==="
echo "Input file:     $IN_FASTA"
echo "Forward primer: $PRIMER_FWD (min overlap: $MIN_OVERLAP_FWD, max error rate: $MAX_ERROR_FWD)"
echo "Reverse primer: $PRIMER_REV (min overlap: $MIN_OVERLAP_REV, max error rate: $MAX_ERROR_REV)"
echo ""

TRIMMED_WANDA_RAW="./tmp/eukaryome_trimmed_wanda_aml2_raw.fasta"
TRIMMED_WANDA="./tmp/eukaryome_trimmed_wanda_aml2.fasta"
UNTRIMMED_WANDA="./tmp/eukaryome_untrimmed_wanda_aml2.fasta"

cutadapt \
  -g "$PRIMER_FWD;min_overlap=$MIN_OVERLAP_FWD;max_error_rate=$MAX_ERROR_FWD;required"..."$PRIMER_REV;min_overlap=$MIN_OVERLAP_REV;max_error_rate=$MAX_ERROR_REV;required" \
  --cores "$NUM_THREADS" \
  --untrimmed-output "$UNTRIMMED_WANDA" \
  -o "$TRIMMED_WANDA_RAW" "$IN_FASTA" \
  >> "$LOG_FILE" 2>&1 < /dev/null

echo "Primer trimming completed."
if [[ -f "$LOG_FILE" ]]; then
    echo ""
    echo "=== WANDA/AML2 TRIMMING SUMMARY ==="
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

LOG_FILTER="tmp/logfile_cutadapt_filter.txt"
TOO_LONG_TRIMMED="./tmp/too_long_trimmed.fasta"

cutadapt \
  -m "$MIN_LEN" -M "$MAX_LEN" \
  --max-n "$MAX_N" \
  --cores "$NUM_THREADS" \
  --too-long-output "$TOO_LONG_TRIMMED" \
  -o "$TRIMMED_WANDA" "$TRIMMED_WANDA_RAW" \
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

echo "Running reformat.R..."
Rscript "$REFORMAT" \
    --fasta_in           "$TRIMMED_WANDA" \
    --fasta_out          "tmp/reformatted.fasta" \
    --classification_out "tmp/reformatted.classification"

if [[ $? -ne 0 ]]; then
    echo "ERROR: reformat.R failed." >&2
    exit 1
fi

echo ""
echo "Reformatted FASTA written to     : tmp/reformatted.fasta"
echo "Reformatted classification table : tmp/reformatted.classification"
echo ""

# =============================================================================
# CHECK AND STANDARDISE ANNOTATIONS
# =============================================================================

echo "=== CHECKING AND STANDARDISING ANNOTATIONS ==="
echo "$(date)"

Rscript "$CHECK_ANNOTATIONS" \
    --classification_in  "tmp/reformatted.classification" \
    --classification_out "tmp/reformatted.classification"

if [[ $? -ne 0 ]]; then
    echo "ERROR: check_annotations.R failed." >&2
    exit 1
fi
echo ""

# =============================================================================
# DEREPLICATION WITH LCA ASSIGNMENT
# =============================================================================

echo "=== DEREPLICATE AND RESOLVE LCA TAXONOMY ==="
echo "$(date)"

Rscript "$DEREP_LCA" \
    --fasta_in             "tmp/reformatted.fasta" \
    --fasta_out            "$OUT_FASTA" \
    --classification_in    "tmp/reformatted.classification" \
    --classification_out   "$OUT_CLASS"

if [[ $? -ne 0 ]]; then
    echo "ERROR: dereplicate_lca.R failed." >&2
    exit 1
fi

echo ""
echo "Predict FASTA written to         : $OUT_FASTA"
echo "Predict classification written to: $OUT_CLASS"
echo ""

# =============================================================================
# GENERATE DEGENERATE FORWARD PRIMER FROM TRIMMED READS
# =============================================================================

echo "=== GENERATING DEGENERATE FORWARD PRIMER ==="
echo "$(date)"
echo ""

CONSENSUS_FWD=$(Rscript -e "
    library(Biostrings)
    seqs <- readDNAStringSet('${TRIMMED_WANDA}')
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
echo "Input:          $UNTRIMMED_WANDA"
echo "Forward primer: $CONSENSUS_FWD (${FWD_CONSENSUS_LEN}bp)"
echo "Action:         retain (primer kept in output)"
echo "Truncation:     ${TRUNCATION_LEN}bp"
echo ""

FWD_FILTERED_RAW="./tmp/eukaryome_fwd_filtered_raw.fasta"
FWD_FILTERED="./tmp/eukaryome_fwd_filtered.fasta"
LOG_FWD="tmp/logfile_cutadapt_fwd_filter.txt"

cutadapt \
  -g "${CONSENSUS_FWD};min_overlap=${FWD_CONSENSUS_LEN};max_error_rate=1" \
  --action=retain \
  --discard-untrimmed \
  --length "$TRUNCATION_LEN" \
  --cores "$NUM_THREADS" \
  -o "$FWD_FILTERED_RAW" "$UNTRIMMED_WANDA" \
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

TOO_LONG_TRUNCATED="./tmp/too_long_truncated.fasta"

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
# LENGTH FILTER AND N REMOVAL
# =============================================================================

echo "=== LENGTH FILTER AND N REMOVAL ==="
echo "Min length: ${MIN_LEN}bp"
echo "Max N bases: 0"
echo ""

vsearch \
  --fastx_filter "$FWD_FILTERED_RAW" \
  --fastaout "$FWD_FILTERED" \
  --fastq_maxns 0 \
  --fastq_minlen "$MIN_LEN" \
  --fasta_width 0 \
  2> tmp/logfile_vsearch_filter_partial.txt

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

echo "Running reformat.R..."
Rscript "$REFORMAT" \
    --fasta_in           "$FWD_FILTERED" \
    --fasta_out          "tmp/reformatted_partial.fasta" \
    --classification_out "tmp/reformatted_partial.classification"

if [[ $? -ne 0 ]]; then
    echo "ERROR: reformat.R failed." >&2
    exit 1
fi

echo ""
echo "Reformatted FASTA written to     : tmp/reformatted_partial.fasta"
echo "Reformatted classification table : tmp/reformatted_partial.classification"
echo ""

# =============================================================================
# CHECK AND STANDARDISE ANNOTATIONS (partial)
# =============================================================================

echo "=== CHECKING AND STANDARDISING ANNOTATIONS (partial) ==="
echo "$(date)"

Rscript "$CHECK_ANNOTATIONS" \
    --classification_in  "tmp/reformatted_partial.classification" \
    --classification_out "tmp/reformatted_partial.classification"

if [[ $? -ne 0 ]]; then
    echo "ERROR: check_annotations.R failed." >&2
    exit 1
fi
echo ""

# =============================================================================
# DEREPLICATION WITH LCA ASSIGNMENT (partial)
# =============================================================================

echo "=== DEREPLICATE AND RESOLVE LCA TAXONOMY (partial) ==="
echo "$(date)"

Rscript "$DEREP_LCA" \
    --fasta_in             "tmp/reformatted_partial.fasta" \
    --fasta_out            "$OUT_FASTA_PARTIAL" \
    --classification_in    "tmp/reformatted_partial.classification" \
    --classification_out   "$OUT_CLASS_PARTIAL"

if [[ $? -ne 0 ]]; then
    echo "ERROR: dereplicate_lca.R failed." >&2
    exit 1
fi

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
  --matched    tmp/overlap_partial_in_full.fasta \
  --notmatched tmp/unique_partial.fasta \
  --userout    tmp/usearch_overlap.txt \
  --userfields query+target+id+ql+tl+alnlen+qcov+tcov+ids \
  --threads "$NUM_THREADS" \
  2> tmp/logfile_usearch_overlap.txt

n_matched=$(grep -c "^>" tmp/overlap_partial_in_full.fasta 2>/dev/null || true)
n_unique=$(grep -c  "^>" tmp/unique_partial.fasta          2>/dev/null || true)
echo "Partial sequences with 100% identity match in full : ${n_matched:-0}"
echo "Partial sequences unique to the partial set          : ${n_unique:-0}"
echo "Overlap FASTA written to: tmp/overlap_partial_in_full.fasta"
echo "Unique  FASTA written to: tmp/unique_partial.fasta"
echo ""

# =============================================================================
# COMBINE FULL AND PARTIAL WITH LCA FOR OVERLAPPING SEQUENCES
# =============================================================================

echo "=== COMBINING FULL AND PARTIAL SEQUENCES ==="
echo "$(date)"
echo ""

COMBINE_LCA="R/combine_lca.R"

Rscript "$COMBINE_LCA" \
    --fasta_full    "$OUT_FASTA" \
    --class_full    "$OUT_CLASS" \
    --fasta_partial "$OUT_FASTA_PARTIAL" \
    --class_partial "$OUT_CLASS_PARTIAL" \
    --usearch_hits  "tmp/usearch_overlap.txt" \
    --fasta_out     "$OUT_FASTA_COMBINED" \
    --class_out     "$OUT_CLASS_COMBINED"

if [[ $? -ne 0 ]]; then
    echo "ERROR: combine_lca.R failed." >&2
    exit 1
fi

echo ""
echo "Combined FASTA written to:          $OUT_FASTA_COMBINED"
echo "Combined classification written to: $OUT_CLASS_COMBINED"
echo ""

# =============================================================================
# CLEANUP
# =============================================================================

echo "=== CLEANUP ==="
echo "Removing tmp/ contents..."
rm -rf tmp/*

echo ""
echo "=== PIPELINE COMPLETED SUCCESSFULLY ==="
echo "$(date)"

conda deactivate