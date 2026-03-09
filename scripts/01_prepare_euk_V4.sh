#!/bin/bash

# Script name:  01_prepare_V4_v2.sh
# Description:  Download the EUKARYOME SSU v2.0 and MAARJAM SSU v.2021 database,
#               trim EUKARYOME reads with WANDA/AML2, then generate a universal
#               15bp forward consensus and phylum-specific 12bp reverse consensus
#               primers from the trimmed reads. These are applied back to the
#               untrimmed EUKARYOME to maximise V4 recovery across all lineages.
#               MAARJAM sequences are filtered with Glomeromycota-specific primers.
#               Finally, headers are reformatted, databases merged, and duplicates
#               removed with reformat_V4.R.
# Note:         MAARJAM incorrectly states that V4 reads are forward-trimmed with
#               NS31, but they are actually trimmed with WANDA.
# Note:         This script must be run from the project root directory.

set -euo pipefail

# =============================================================================
# PARAMETER SETUP
# =============================================================================

# EUKARYOME download URL (EUKARYOME General EUK SSU v2.0)
readonly EUKARYOME_URL="https://sisu.ut.ee/wp-content/uploads/sites/643/General_EUK_SSU_v2.0.zip"
readonly DOWNLOAD_FILE="./tmp/General_EUK_SSU_v2.0.zip"
readonly IN_FASTA="./tmp/General_EUK_SSU_v2.0.fasta"

# Output files
readonly OUT_FASTA="./data/eukaryome_V4.fasta"
readonly OUT_CLASS="./data/eukaryome_V4.classification"

# Intermediate directories
readonly PRIMER_DIR="./tmp/primers"
readonly PHYLUM_DIR="./tmp/phyla"

# CUTADAPT PARAMETERS
readonly LOG_FILE="tmp/logfile_cutadapt_eukaryome.txt"
readonly NUM_THREADS=10
readonly MAX_ERROR_RATE=1
readonly MIN_LEN=440
readonly MAX_LEN=540

# WANDA/AML2 trimming parameters
readonly MIN_OVERLAP_FWD=20
readonly MIN_OVERLAP_REV=22

# WANDA 5' forward primer and AML2 3' reverse primer (reverse complement)
readonly PRIMER_FWD="CAGCCGCGGTAATTCCAGCT"
readonly PRIMER_REV="GGAAACCAAAGTGTTTGGGTTC"

# Consensus primer parameters
readonly FWD_CONSENSUS_LEN=15
readonly REV_CONSENSUS_LEN=12
readonly FWD_THRESHOLD=0.25
readonly REV_THRESHOLD=0.05
readonly MIN_SEQS_FOR_CONSENSUS=10

# Helper R script
readonly REFORMAT="./R/reformat.R"

# =============================================================================
# DIRECTORY SETUP
# =============================================================================

mkdir -p ./data/ tmp "$PRIMER_DIR" "$PHYLUM_DIR"

# =============================================================================
# ENVIRONMENT SETUP
# =============================================================================

echo "Activating conda environment..."
set +u
source ~/.bash_profile
conda activate dyna_clust_predict
set -u

# =============================================================================
# INPUT VALIDATION
# =============================================================================

if [[ ! -f "$REFORMAT" ]]; then
    echo "ERROR: R script not found: $REFORMAT" >&2
    exit 1
fi

# =============================================================================
# FILE DOWNLOAD: EUKARYOME
# =============================================================================

echo ""
echo "=== DOWNLOADING EUKARYOME DATABASE ==="
echo "$(date)"
echo "URL: $EUKARYOME_URL"

if ! curl -L -o "$DOWNLOAD_FILE" "$EUKARYOME_URL"; then
    echo "ERROR: Failed to download file from $EUKARYOME_URL" >&2
    exit 1
fi

echo "Unzipping downloaded file..."
if ! 7z x "$DOWNLOAD_FILE" -o"./tmp/" -y; then
    echo "ERROR: Failed to unzip $DOWNLOAD_FILE" >&2
    exit 1
fi

echo "Download and extraction completed successfully."
echo ""

# =============================================================================
# PRIMER TRIMMING: EUKARYOME (WANDA / AML2)
# =============================================================================

if [[ ! -f "$IN_FASTA" ]]; then
    echo 'ERROR: Input file $IN_FASTA not found!' >&2
    exit 1
fi

echo "=== PRIMER TRIMMING WITH WANDA / AML2 ==="
echo "Input file: $IN_FASTA"
echo "Output file: $OUT_FASTA"
echo "Forward primer: $PRIMER_FWD (min overlap: $MIN_OVERLAP_FWD)"
echo "Reverse primer: $PRIMER_REV (min overlap: $MIN_OVERLAP_REV)"
echo "Max error rate: $MAX_ERROR_RATE"
echo ""

TRIMMED_WANDA="./tmp/eukaryome_trimmed_wanda_aml2.fasta"

cutadapt \
  -g "$PRIMER_FWD;min_overlap=$MIN_OVERLAP_FWD;required"..."$PRIMER_REV;min_overlap=$MIN_OVERLAP_REV;required" \
  -e "$MAX_ERROR_RATE" \
  --discard-untrimmed \
  -M "$MAX_LEN" \
  -m "$MIN_LEN" \
  --cores "$NUM_THREADS" \
  -o "$TRIMMED_WANDA" "$IN_FASTA" \
  >> "$LOG_FILE" 2>&1

echo "Primer trimming completed."
if [[ -f "$LOG_FILE" ]]; then
    echo ""
    echo "=== WANDA/AML2 TRIMMING SUMMARY ==="
    grep -E "(Total reads processed|Reads with adapters|Reads that were too short|Reads that were too long|Reads discarded as untrimmed|Reads written)" "$LOG_FILE" | tail -6
fi
echo ""

# =============================================================================
# GENERATE CONSENSUS PRIMERS FROM TRIMMED READS
# =============================================================================

echo "=== GENERATING CONSENSUS PRIMERS ==="
echo "$(date)"
echo ""

# ── Universal 15bp forward primer ────────────────────────────────────────────

CONSENSUS_FWD=$(Rscript -e "
    library(Biostrings)
    seqs <- readDNAStringSet('${TRIMMED_WANDA}')
    seqs <- seqs[width(seqs) >= ${FWD_CONSENSUS_LEN}]
    prefixes <- subseq(seqs, 1, ${FWD_CONSENSUS_LEN})
    cm <- consensusMatrix(prefixes)
    cons <- consensusString(cm, ambiguityMap = IUPAC_CODE_MAP, threshold = ${FWD_THRESHOLD})
    cat(toupper(cons))
" 2>/dev/null)

echo "Universal forward consensus (${FWD_CONSENSUS_LEN}bp): $CONSENSUS_FWD"

# ── Per-phylum 12bp reverse primers ──────────────────────────────────────────

echo ""
echo "Generating per-phylum reverse primers (${REV_CONSENSUS_LEN}bp, threshold=${REV_THRESHOLD})..."
echo ""

# Extract phyla with >= MIN_SEQS_FOR_CONSENSUS sequences
grep "^>" "$TRIMMED_WANDA" \
  | sed 's/.*p__//;s/;.*//' \
  | sort | uniq -c | sort -rn \
  | awk -v min="$MIN_SEQS_FOR_CONSENSUS" '$1 >= min {print $2}' \
  > "$PRIMER_DIR/phyla_list.txt"

N_PHYLA=$(wc -l < "$PRIMER_DIR/phyla_list.txt" | tr -d ' ')
echo "Phyla with >= ${MIN_SEQS_FOR_CONSENSUS} trimmed sequences: ${N_PHYLA}"
echo ""

printf "%-25s %8s  %-${REV_CONSENSUS_LEN}s\n" "Phylum" "N_seqs" "REV_${REV_CONSENSUS_LEN}bp"
echo "-----------------------------------------------"

while read -r phylum; do
    
    phylum_fasta="${PHYLUM_DIR}/${phylum}_trimmed.fasta"
    seqkit grep -r -p "p__${phylum}" "$TRIMMED_WANDA" > "$phylum_fasta" 2>/dev/null
    n_seqs=$(grep -c "^>" "$phylum_fasta" 2>/dev/null || echo 0)

   
    rev_cons=$(Rscript -e "
        library(Biostrings)
        seqs <- readDNAStringSet('${phylum_fasta}')
        seqs <- seqs[width(seqs) >= ${REV_CONSENSUS_LEN}]
        suffixes <- subseq(seqs, width(seqs) - ${REV_CONSENSUS_LEN} + 1, width(seqs))
        cm <- consensusMatrix(suffixes)
        cons <- consensusString(cm, ambiguityMap = IUPAC_CODE_MAP, threshold = ${REV_THRESHOLD})
        cat(toupper(cons))
        " 2>/dev/null)

    
    echo "$rev_cons" > "${PRIMER_DIR}/${phylum}_rev.txt"
    printf "%-25s %8d  %s\n" "$phylum" "$n_seqs" "$rev_cons"

done < "$PRIMER_DIR/phyla_list.txt"

echo ""

# =============================================================================
# APPLY CONSENSUS PRIMERS TO UNTRIMMED EUKARYOME (PER PHYLUM)
# =============================================================================

echo "=== FILTERING UNTRIMMED EUKARYOME WITH CONSENSUS PRIMERS ==="
echo "$(date)"
echo ""
echo "Forward primer: $CONSENSUS_FWD (${FWD_CONSENSUS_LEN}bp, universal)"
echo "Reverse primers: per-phylum ${REV_CONSENSUS_LEN}bp consensus"
echo "Action: retain (keep primer + internal sequence)"
echo ""

# For each phylum, extract its untrimmed sequences and apply its own primers

# First empty the combined output file
> "$PHYLUM_DIR/all_filtered.fasta"

printf "%-25s %8s %8s %7s\n" "Phylum" "Untrim" "Recover" "%"
echo "-----------------------------------------------------"

while read -r phylum; do
    rev_cons=$(cat "${PRIMER_DIR}/${phylum}_rev.txt")

    untrimmed_phylum="${PHYLUM_DIR}/${phylum}_untrimmed.fasta"
    seqkit grep -r -p "p__${phylum}" "$IN_FASTA" > "$untrimmed_phylum" 2>/dev/null
    n_untrimmed=$(grep -c "^>" "$untrimmed_phylum" 2>/dev/null || echo 0)

    filtered_phylum="${PHYLUM_DIR}/${phylum}_filtered.fasta"
    cutadapt \
      -g "${CONSENSUS_FWD};min_overlap=${FWD_CONSENSUS_LEN};required...${rev_cons};min_overlap=${REV_CONSENSUS_LEN};required" \
      -e "$MAX_ERROR_RATE" \
      --action=retain \
      --discard-untrimmed \
      -m "$MIN_LEN" -M "$MAX_LEN" \
      --cores "$NUM_THREADS" \
      -o "$filtered_phylum" "$untrimmed_phylum" \
      > "${PHYLUM_DIR}/${phylum}_cutadapt.log" 2>&1

    n_recovered=$(grep -c "^>" "$filtered_phylum" 2>/dev/null || echo 0)
    if [[ "$n_untrimmed" -gt 0 ]]; then
        pct=$(echo "scale=1; $n_recovered * 100 / $n_untrimmed" | bc)
    else
        pct="0.0"
    fi
    printf "%-25s %8d %8d %6s%%\n" "$phylum" "$n_untrimmed" "$n_recovered" "$pct"

    cat "$filtered_phylum" >> "$PHYLUM_DIR/all_filtered.fasta"

done < "$PRIMER_DIR/phyla_list.txt"

# Also handle phyla below the MIN_SEQS_FOR_CONSENSUS threshold
# For these, fall back to the WANDA/AML2-trimmed version
echo ""
echo "Phyla below threshold (${MIN_SEQS_FOR_CONSENSUS} seqs): using WANDA/AML2-trimmed reads"
n_fallback=0

# Get the set of phyla already handled
handled_phyla=$(cat "$PRIMER_DIR/phyla_list.txt" | tr '\n' '|' | sed 's/|$//')

# Extract sequences from trimmed output that are NOT in any handled phylum
if [[ -n "$handled_phyla" ]]; then
    seqkit grep -r -v -p "p__(${handled_phyla})" "$TRIMMED_WANDA" > "$PHYLUM_DIR/fallback_trimmed.fasta" 2>/dev/null
else
    cp "$TRIMMED_WANDA" "$PHYLUM_DIR/fallback_trimmed.fasta"
fi

n_fallback=$(grep -c "^>" "$PHYLUM_DIR/fallback_trimmed.fasta" 2>/dev/null || echo 0)
echo "  Fallback sequences: $n_fallback"

cat "$PHYLUM_DIR/fallback_trimmed.fasta" >> "$PHYLUM_DIR/all_filtered.fasta"

# =============================================================================
# REFORMAT HEADERS
# =============================================================================

echo "=== REFORMATTING HEADERS ==="
echo "$(date)"

if [[ ! -f "$PHYLUM_DIR/all_filtered.fasta" ]]; then
    echo "ERROR: Input FASTA for reformatting headers not found: ${PHYLUM_DIR}/all_filtered.fasta" >&2
    exit 1
fi

echo "Running reformat.R..."
Rscript "$REFORMAT" \
    --fasta_in           "$PHYLUM_DIR/all_filtered.fasta" \
    --fasta_out          "$OUT_FASTA" \
    --classification_out "$OUT_CLASS"

if [[ $? -ne 0 ]]; then
    echo "ERROR: reformat.R failed." >&2
    exit 1
fi

echo ""
echo "Reformatted FASTA written to      : $OUT_FASTA"
echo "Classification table written to   : $OUT_CLASS"
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
