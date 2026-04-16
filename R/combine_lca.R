#!/usr/bin/env Rscript
# combine_lca.R — Combine full and partial V4 sequence sets with LCA for overlaps
#
# For partial sequences that match a full sequence at 100% identity with 100%
# query OR target coverage, the full sequence is retained as representative and
# its classification is updated via LCA using both classifications. All remaining
# partial sequences (no 100%-overlap match) are appended as-is.
#
# Usage:
#   Rscript R/combine_lca.R \
#     --fasta_full    data/eukaryome_V4_full.fasta \
#     --class_full    data/eukaryome_V4_full.classification \
#     --fasta_partial data/eukaryome_V4_partial.fasta \
#     --class_partial data/eukaryome_V4_partial.classification \
#     --usearch_hits  tmp/usearch_overlap.txt \
#     --fasta_out     data/eukaryome_V4_all.fasta \
#     --class_out     data/eukaryome_V4_all.classification
#
# Note: This script must be run from the project root directory.

suppressPackageStartupMessages({
  library(optparse)
  library(Biostrings)
  library(data.table)
  library(dplyr)
})

source("R/utils.R")

# ── Arguments ─────────────────────────────────────────────────────────────────

option_list <- list(
  make_option("--fasta_full",    type = "character", metavar = "FILE",
              help = "Full V4 FASTA [required]"),
  make_option("--class_full",    type = "character", metavar = "FILE",
              help = "Full V4 classification [required]"),
  make_option("--fasta_partial", type = "character", metavar = "FILE",
              help = "Partial V4 FASTA [required]"),
  make_option("--class_partial", type = "character", metavar = "FILE",
              help = "Partial V4 classification [required]"),
  make_option("--usearch_hits",  type = "character", metavar = "FILE",
              help = "vsearch --userout file (query+target+id+ql+tl+alnlen+qcov+tcov+ids) [required]"),
  make_option("--fasta_out",     type = "character", metavar = "FILE",
              help = "Output combined FASTA [required]"),
  make_option("--class_out",     type = "character", metavar = "FILE",
              help = "Output combined classification [required]")
)

opt <- parse_args(
  OptionParser(
    option_list = option_list,
    usage       = "%prog --fasta_full full.fa --class_full full.tsv ...",
    description = "Combine full and partial V4 sets with LCA for 100%-overlap matches."
  )
)

required_args <- c("fasta_full", "class_full", "fasta_partial", "class_partial",
                   "usearch_hits", "fasta_out", "class_out")
for (arg in required_args) {
  if (is.null(opt[[arg]])) stop("--", arg, " is required.")
}
for (arg in setdiff(required_args, c("fasta_out", "class_out"))) {
  if (!file.exists(opt[[arg]])) stop("File not found: ", opt[[arg]])
}

# ── Load inputs ────────────────────────────────────────────────────────────────

taxonomy_ranks <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")

cat("Loading full sequences from:", opt$fasta_full, "\n")
seqs_full <- readDNAStringSet(opt$fasta_full)
cat("Full sequences loaded:", length(seqs_full), "\n")

cat("Loading full classification from:", opt$class_full, "\n")
cls_full <- as.data.frame(fread(opt$class_full))

cat("Loading partial sequences from:", opt$fasta_partial, "\n")
seqs_partial <- readDNAStringSet(opt$fasta_partial)
cat("Partial sequences loaded:", length(seqs_partial), "\n")

cat("Loading partial classification from:", opt$class_partial, "\n")
cls_partial <- as.data.frame(fread(opt$class_partial))

cat("Loading usearch hits from:", opt$usearch_hits, "\n")
hits <- fread(opt$usearch_hits, header = FALSE,
              col.names = c("query", "target", "id", "ql", "tl", "alnlen",
                            "qcov", "tcov", "ids"))

# ── Filter for 100% identity with 100% query or target coverage ───────────────

hits_100 <- hits[id >= 100 & (qcov >= 100 | tcov >= 100)]
cat("Hits with 100% identity and 100% query/target coverage:", nrow(hits_100), "\n\n")

# ── LCA helper ────────────────────────────────────────────────────────────────

resolve_lca_rows <- function(rows) {
  result <- character(length(taxonomy_ranks))
  names(result) <- taxonomy_ranks
  stop_propagating <- FALSE
  for (rank in taxonomy_ranks) {
    if (stop_propagating) {
      result[rank] <- "unclassified"
      next
    }
    vals <- unique(rows[[rank]])
    vals <- vals[!vals %in% c("unclassified", "unidentified", NA, "")]
    if (length(vals) == 1) {
      result[rank] <- vals
    } else {
      result[rank] <- "unclassified"
      stop_propagating <- TRUE
    }
  }
  result
}

# ── Apply LCA updates to full classification ──────────────────────────────────

cls_full_updated <- cls_full

if (nrow(hits_100) > 0) {

  matched_partial_ids <- unique(hits_100$query)
  hits_by_target      <- split(hits_100$query, hits_100$target)

  for (target_id in names(hits_by_target)) {
    partial_ids  <- hits_by_target[[target_id]]
    full_row     <- cls_full_updated[cls_full_updated$id == target_id, , drop = FALSE]
    partial_rows <- cls_partial[cls_partial$id %in% partial_ids, , drop = FALSE]
    all_rows     <- rbind(full_row, partial_rows)

    lca <- resolve_lca_rows(all_rows)
    idx <- which(cls_full_updated$id == target_id)
    for (rank in taxonomy_ranks) {
      cls_full_updated[idx, rank] <- lca[rank]
    }
  }

  cat("Full sequences with LCA-updated classification:", length(hits_by_target), "\n")

} else {
  matched_partial_ids <- character(0)
  cat("No 100%-overlap matches found; no LCA updates applied.\n")
}

# ── Remove matched partial sequences ──────────────────────────────────────────

seqs_partial_unique <- seqs_partial[!names(seqs_partial) %in% matched_partial_ids]
cls_partial_unique  <- cls_partial[!cls_partial$id %in% matched_partial_ids, ]
cat("Partial sequences removed (matched to full):    ", length(matched_partial_ids), "\n")
cat("Partial sequences retained (unique to partial): ", length(seqs_partial_unique), "\n")

# ── Combine ───────────────────────────────────────────────────────────────────

seqs_all <- c(seqs_full, seqs_partial_unique)
cls_all  <- rbind(cls_full_updated, cls_partial_unique)

cat("\nCombined set summary\n")
cat("--------------------\n")
cat("Full sequences:           ", length(seqs_full), "\n")
cat("Unique partial sequences: ", length(seqs_partial_unique), "\n")
cat("Total combined:           ", length(seqs_all), "\n\n")

# ── Write outputs ─────────────────────────────────────────────────────────────

writeXStringSet(seqs_all, opt$fasta_out)
cat("Combined FASTA written to:          ", opt$fasta_out, "\n")

fwrite(cls_all, opt$class_out, sep = "\t")
cat("Combined classification written to: ", opt$class_out, "\n")
