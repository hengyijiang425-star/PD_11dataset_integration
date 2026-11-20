# 01_download_GEO.R — download GEO supplementary files for every dataset with loader = "auto"
# and a path under data/raw/ (the three old datasets are read from PD_3dataset_integration instead).
# download_mode = "series"    : the series' own supplementary files (usually <GSE>_RAW.tar)
# download_mode = "per_sample": each GSM's files separately (used when the RAW.tar is huge because of ATAC)
# Files matching geo_download_exclude_regex (ATAC, TCR, images, other species...) are skipped.
# Re-running skips complete files and retries failed ones.
source("R/utils.R")
suppressPackageStartupMessages(library(GEOquery))

ds <- cfg_datasets(); ds <- ds[ds$loader == "auto" & grepl("^data/raw", ds$path), ]
log <- list()
get_one <- function(g, acc, out, ex) {
  lst <- tryCatch(getGEOSuppFiles(acc, makeDirectory = FALSE, baseDir = out, fetch_files = FALSE),
                  error = function(e) { msg(acc, "cannot list files:", conditionMessage(e)); NULL })
  if (is.null(lst) || !nrow(lst)) return(NULL)
  skip <- if (is.na(ex)) rep(FALSE, nrow(lst)) else grepl(ex, lst$fname)
  status <- ifelse(skip, "skipped_by_regex", NA)
  for (j in which(!skip)) {
    dest <- file.path(out, lst$fname[j])
    if (file.exists(dest) && file.size(dest) > 0) { status[j] <- "already_there"; next }
    msg(g, "downloading", lst$fname[j])
    ok <- tryCatch({ download.file(lst$url[j], dest, mode = "wb", quiet = TRUE); TRUE },
                   error = function(e) { msg("  FAILED:", conditionMessage(e)); unlink(dest); FALSE })
    status[j] <- if (ok) "downloaded" else "FAILED"
  }
  data.frame(dataset = g, accession = acc, file = lst$fname, status = status)
}
for (i in seq_len(nrow(ds))) {
  g <- ds$dataset[i]; out <- ds$path[i]; ex <- ds$geo_download_exclude_regex[i]
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  if (!is.na(ds$download_mode[i]) && ds$download_mode[i] == "per_sample") {
    e <- tryCatch(getGEO(g, GSEMatrix = TRUE, getGPL = FALSE), error = function(e) NULL)
    if (is.null(e)) { msg(g, ": cannot read the sample list"); next }
    gsm <- unique(unlist(lapply(e, function(x) Biobase::pData(x)$geo_accession)))
    msg(g, ":", length(gsm), "samples, downloading per sample")
    log[[g]] <- do.call(rbind, lapply(gsm, function(a) get_one(g, a, out, ex)))
  } else log[[g]] <- get_one(g, g, out, ex)
}
log <- do.call(rbind, log)
write.csv(log, file.path(DIR$res, "01_download_log.csv"), row.names = FALSE)
print(table(log$dataset, log$status))
if (any(log$status == "FAILED", na.rm = TRUE)) cat("\nSome files FAILED (often a timeout on big files): just run this script again.\n")
cat("\nNext: Rscript R/02_inventory_and_metadata.R\n")
