# run_all.R — run steps 03 to 08 in order (after 01 + 02 and after you checked config/sample_metadata.csv)
# Rscript R/run_all.R            all steps
# Rscript R/run_all.R 05         start from step 05
args <- commandArgs(trailingOnly = TRUE); from <- if (length(args)) as.integer(args[1]) else 3
steps <- c("03_load_counts.R", "04_qc_doublets.R", "05_merge_integrate.R", "06_annotate_L1.R",
           "07_lineage_subclustering.R", "08_transfer_and_outputs.R")
for (s in steps[as.integer(substr(steps, 1, 2)) >= from]) {
  cat("\n================", s, "================\n"); t0 <- Sys.time()
  source(file.path("R", s), echo = FALSE, local = new.env())
  cat(s, "took", round(difftime(Sys.time(), t0, units = "mins"), 1), "min\n")
}
