#!/usr/bin/env Rscript
# test_script_params.R — Static checks that this repo's scripts pass the
# similarity settings the shared R/compute_sim.R and R/predict.R expect, and
# read predict.R's output under its current name. Not synced from
# dyna-clust-predict: these scripts are this repo's own.
#
# Usage: Rscript tests/test_script_params.R
# Note:  Must be run from the project root directory.

n_pass <- 0L
n_fail <- 0L
check <- function(desc, ok) {
  if (isTRUE(ok)) {
    n_pass <<- n_pass + 1L
    cat("PASS ", desc, "\n", sep = "")
  } else {
    n_fail <<- n_fail + 1L
    cat("FAIL ", desc, "\n", sep = "")
  }
}

# ── 02 and 03: --iddef 2 and --min_overlap from the primer config ─────────────

for (sc in c("scripts/02_predict_cutoffs.sh", "scripts/03_predict_euk_cutoffs.sh")) {
  lines <- readLines(sc)
  name  <- basename(sc)
  check(sprintf("%s: IDDEF=2", name), any(lines == "readonly IDDEF=2"))
  check(sprintf("%s: MIN_OVERLAP is the primer config's MIN_LEN", name),
        any(lines == 'readonly MIN_OVERLAP="$MIN_LEN"'))
  sim_calls <- grep('Rscript "\\$COMPUTE_SIM"', lines, value = TRUE)
  check(sprintf("%s: compute_sim.R gets --iddef and --min_overlap", name),
        length(sim_calls) == 1L &&
          grepl('--iddef "$IDDEF" --min_overlap "$MIN_OVERLAP"', sim_calls, fixed = TRUE))
  check(sprintf("%s: predict.R gets --iddef and --min_overlap", name),
        any(grepl('--iddef "$IDDEF"', lines, fixed = TRUE) & !grepl("COMPUTE_SIM", lines)) &&
          any(grepl('--min_overlap "$MIN_OVERLAP"', lines, fixed = TRUE) & !grepl("COMPUTE_SIM", lines)))
}

scripts <- list.files("scripts", pattern = "[.]sh$", full.names = TRUE)
check("no script passes the removed --min_cov",
      !any(vapply(scripts, function(f) any(grepl("min_cov|MIN_COV", readLines(f))), logical(1))))

# ── Every primer config gives a usable MIN_LEN ────────────────────────────────

for (conf in list.files("config/primers", pattern = "[.]conf$", full.names = TRUE)) {
  v <- sub("^MIN_LEN=", "", grep("^MIN_LEN=[0-9]+$", readLines(conf), value = TRUE))
  check(sprintf("%s: MIN_LEN is a positive integer", basename(conf)),
        length(v) == 1L && as.integer(v) > 0L)
}

# ── 05 reads predict.R's output under its current name ────────────────────────

check("predict.R writes <prefix>_raw_cutoffs.json (and .txt)",
      any(grepl('"_raw_cutoffs.json"', readLines("R/predict.R"), fixed = TRUE)))
for (f in c("scripts/05_combine_cutoffs.sh", "R/combine_cutoffs.R")) {
  lines <- readLines(f)
  check(sprintf("%s: reads _raw_cutoffs.txt, not .cutoffs.json.txt", basename(f)),
        any(grepl("_raw_cutoffs.txt", lines, fixed = TRUE)) &&
          !any(grepl("cutoffs.json", lines, fixed = TRUE)))
}

# ── Result ────────────────────────────────────────────────────────────────────

cat(sprintf("\nRESULT: %d passed, %d failed\n", n_pass, n_fail))
if (n_fail > 0) quit(status = 1)
