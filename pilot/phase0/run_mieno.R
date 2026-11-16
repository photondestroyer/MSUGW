# Task 0.3 runner: Mieno et al. (2024) analysis and figures on a WORKING COPY of their repository.
# Started by pilot/phase0/t0_3_benchmark.py as:  Rscript --vanilla run_mieno.R <work_dir> <n_boot> <out_dir>
# --vanilla: their .Rprofile and renv bootstrap (which install and download packages) never run.
# Their files are not edited. Two session-only adaptations, both recorded in the comparison note:
#   1. parallel::mclapply runs serially (Windows supports one core); the bootstrap draws are the same kind,
#      the random streams differ from their parallel run, so bootstrap intervals reproduce in distribution.
#   2. n_boot may be lower than their 1000 (argument), to bound the runtime; the point estimates do not use it.
args <- commandArgs(trailingOnly = TRUE)
work <- normalizePath(args[1], winslash = "/")
n_boot <- as.integer(args[2])
out_dir <- args[3]
setwd(work)

suppressPackageStartupMessages({
  library(here); library(dplyr); library(purrr); library(ggplot2); library(stringr); library(tidyr)
  library(data.table); library(sf); library(mgcv); library(fixest); library(cowplot); library(broom)
  library(lubridate); library(scales); library(patchwork)
})
here::i_am("DroughtRiskAquifer.Rproj")
assignInNamespace("mclapply", function(X, FUN, ..., mc.cores = 1L) lapply(X, FUN, ...), "parallel")
source(here("Codes/functions.R"))
dir.create("Results/Figures", recursive = TRUE, showWarnings = FALSE)

# ---- STEP 0.3-R1  their analysis (Codes/1_analysis.R), bootstrap count as given ----
t0 <- Sys.time()
src <- readLines(here("Codes/1_analysis.R"), warn = FALSE)
if (!is.na(n_boot) && n_boot != 1000L) src <- sub("1:1000,", sprintf("1:%d,", n_boot), src, fixed = TRUE)
eval(parse(text = src), envir = globalenv())
cat(sprintf("analysis done in %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))

# ---- STEP 0.3-R2  their figures (Codes/2_figures_tables.R) ----
fig_ok <- tryCatch({ source(here("Codes/2_figures_tables.R")); TRUE }, error = function(e) {
  cat("figures: ERROR ", conditionMessage(e), "\n"); FALSE })

# ---- STEP 0.3-R3  numbers for the comparison note ----
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
counts <- rbindlist(lapply(seq_len(nrow(all_results)), function(i) {
  d <- all_results$data[[i]]
  data.table(crop = all_results$crop[i], n_obs = nrow(d), n_counties = uniqueN(d$sc_code),
             year_min = min(d$year), year_max = max(d$year))
}))
fwrite(counts, file.path(out_dir, "mieno_sample_counts.csv"))
avg <- rbindlist(lapply(seq_len(nrow(all_results)), function(i) {
  x <- copy(all_results$avg_yield_sum[[i]]); x[, crop := all_results$crop[i]]; x
}), fill = TRUE)
fwrite(avg, file.path(out_dir, "mieno_avg_yield_sum.csv"))
dif <- rbindlist(lapply(seq_len(nrow(all_results)), function(i) {
  x <- copy(all_results$yield_dif_test[[i]]); x[, crop := all_results$crop[i]]; x
}), fill = TRUE)
fwrite(dif, file.path(out_dir, "mieno_yield_dif_test.csv"))
fwrite(data.table(n_boot = n_boot, figures_ok = fig_ok, r_version = R.version.string,
                  fixest = as.character(packageVersion("fixest")), mgcv = as.character(packageVersion("mgcv"))),
       file.path(out_dir, "mieno_run_info.csv"))
cat("done\n")
