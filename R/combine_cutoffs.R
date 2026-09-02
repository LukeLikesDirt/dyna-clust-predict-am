#!/usr/bin/env Rscript
# combine_cutoffs.R — Combine per-taxon predicted cutoffs into a
# per-primer-set table and a single master table across all primer sets.
#
# Helper script -- invoked by scripts/05_combine_cutoffs.sh, which
# auto-discovers --primer_sets and --taxa. Call directly only to override
# that discovery:
#   Rscript R/combine_cutoffs.R \
#     --primer_sets wanda_aml2,amv45nf_amdgr \
#     --taxa        glomeromycota,endogonomycetes,eukaryome
#
# For each <primer_set>, reads data/<primer_set>/<taxon>.cutoffs.json.txt
# (written by 02_predict_cutoffs.sh) for every taxon in --taxa, and:
#   1. Writes output/cutoffs_am_<primer_set>.txt -- the per-primer-set
#      combined table (family/genus/species rows only -- order is dropped,
#      it isn't used downstream by the classification pipeline).
#   2. Appends those same rows to output/cutoffs_am.txt -- the standalone
#      master table across every primer set passed in this invocation.
#
# Both tables carry a primers column (human-readable label, e.g.
# "WANDA–AML2", from PRIMER_LABEL in config/primers/<primer_set>.conf) and a
# primer_string column (each primer's own 5'->3' sequence -- PRIMER_REV in
# the config is stored as its reverse complement for cutadapt, so it is
# reverse-complemented back here for display). cut-off is rounded to the
# nearest 0.005.
#
# Each invocation fully rewrites both the per-primer-set files it touches
# and the master file (not an incremental append across separate runs) --
# pass every primer_set you want represented in one call.

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(readr)
  library(stringr)
  library(Biostrings)
})

option_list <- list(
  make_option("--primer_sets", type = "character",
              help = "Comma-separated primer set names (must match data/<name>/) [required]"),
  make_option("--taxa", type = "character", default = "glomeromycota,endogonomycetes,eukaryome",
              help = "Comma-separated taxon group names [default: %default]"),
  make_option("--data_dir", type = "character", default = "./data",
              help = "Directory containing data/<primer_set>/ [default: %default]"),
  make_option("--config_dir", type = "character", default = "./config/primers",
              help = "Directory containing <primer_set>.conf files [default: %default]"),
  make_option("--output_dir", type = "character", default = "./output",
              help = "Directory for output/cutoffs_am.txt and per-primer-set files [default: %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$primer_sets)) stop("--primer_sets is required (comma-separated).")

primer_sets <- str_split(opt$primer_sets, ",")[[1]] %>% str_trim()
taxa        <- str_split(opt$taxa, ",")[[1]] %>% str_trim()

if (length(primer_sets) < 1) stop("Need at least one primer set.")

message("Primer sets: ", paste(primer_sets, collapse = ", "))
message("Taxa:        ", paste(taxa, collapse = ", "))

dir.create(opt$output_dir, showWarnings = FALSE, recursive = TRUE)

# ── Helpers ───────────────────────────────────────────────────────────────────

round_nearest <- function(x, nearest = 0.005) round(round(x / nearest) * nearest, 3)

# Minimal parser for the primer .conf files: only pulls simple KEY="VALUE"
# lines (which is all these configs use -- no bash arrays, unlike taxon
# configs), so this doesn't need to shell out to actually source them.
read_primer_conf <- function(path) {
  if (!file.exists(path)) stop("Primer config not found: ", path)
  lines <- readLines(path)
  m <- regmatches(lines, regexec('^([A-Za-z_]+)="([^"]*)"$', lines))
  m <- m[lengths(m) == 3]
  vals <- vapply(m, `[[`, character(1), 3)
  keys <- vapply(m, `[[`, character(1), 2)
  as.list(setNames(vals, keys))
}

build_primer_string <- function(conf) {
  names_split <- str_split(conf$PRIMER_LABEL, "–")[[1]] %>% str_trim()  # en dash
  if (length(names_split) != 2) {
    stop("PRIMER_LABEL '", conf$PRIMER_LABEL, "' must be '<forward name>–<reverse name>'.")
  }
  fwd_name <- names_split[1]
  rev_name <- names_split[2]
  # PRIMER_REV is stored as its reverse complement (for cutadapt's 3' trim);
  # reverse-complement it back to the primer's own 5'->3' sequence.
  rev_seq_true <- as.character(reverseComplement(DNAString(conf$PRIMER_REV)))
  paste0(fwd_name, " (5'-", conf$PRIMER_FWD, "-3'); ",
         rev_name, " (5'-", rev_seq_true, "-3')")
}

# ── Read every primer set's per-taxon cutoffs, combine per primer set ────────

all_long <- list()

for (ps in primer_sets) {
  data_dir <- file.path(opt$data_dir, ps)

  conf          <- read_primer_conf(file.path(opt$config_dir, paste0(ps, ".conf")))
  primers       <- conf$PRIMER_LABEL
  primer_string <- build_primer_string(conf)

  ps_rows <- list()

  for (taxon in taxa) {
    f <- file.path(data_dir, paste0(taxon, ".cutoffs.json.txt"))
    if (!file.exists(f)) {
      warning("Missing (skipped): ", f)
      next
    }
    d <- read_tsv(f, show_col_types = FALSE) %>%
      dplyr::rename(cutoff = `cut-off`) %>%
      mutate(taxon = taxon)
    ps_rows[[taxon]] <- d
  }

  if (length(ps_rows) == 0) {
    warning("No cutoff files found for primer set '", ps, "' -- has 02_predict_cutoffs.sh been run for it?")
    next
  }

  ps_combined <- bind_rows(ps_rows) %>%
    filter(rank %in% c("family", "genus", "species")) %>%
    mutate(
      dataset       = str_to_sentence(taxon),
      cutoff        = round_nearest(cutoff),
      primers       = primers,
      primer_string = primer_string
    ) %>%
    select(rank, higher_rank, dataset, cutoff, confidence,
           `sequence number`, `group number`, `max proportion`,
           primers, primer_string) %>%
    dplyr::rename(`cut-off` = cutoff) %>%
    arrange(match(rank, c("family", "genus", "species")), dataset)

  all_long[[ps]] <- ps_combined

  out_file <- file.path(opt$output_dir, paste0("cutoffs_am_", ps, ".txt"))
  write_tsv(ps_combined, out_file)
  message("Wrote ", out_file, " (", nrow(ps_combined), " rows)")
}

if (length(all_long) == 0) stop("No cutoffs found for any primer set -- nothing to combine.")

# ── Master table: every primer set from this invocation, standalone ─────────

master_file <- file.path(opt$output_dir, "cutoffs_am.txt")
master <- bind_rows(all_long)
write_tsv(master, master_file)
message("\nWrote master table: ", master_file, " (", nrow(master), " rows)")
