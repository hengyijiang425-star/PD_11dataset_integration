# 04_qc_doublets.R — one QC rule for every sample + scDblFinder, per sample
# Rule (SET$qc in utils.R): nFeature >= 500, nCount >= 1000, percent.mt < 10 and <= sample median + 3 MAD (floor 1%),
#       nFeature <= 10^(median(log10 nFeature) + 3 * MAD)  (per-sample upper bound)
#       then remove scDblFinder doublets (required); drop samples with < 100 cells left.
# Output: data/processed/cells_qc.rds  (kept cells with QC values)
#         results/04_QC_per_sample.csv, plots/04_QC_*.png
source("R/utils.R")
man <- read.csv(file.path(DIR$proc, "manifest.csv")); cells <- readRDS(file.path(DIR$proc, "cells_raw.rds"))
has_dbl <- requireNamespace("scDblFinder", quietly = TRUE)
if (!has_dbl && isTRUE(SET$qc$require_doublets)) stop("scDblFinder does not load: run library(scDblFinder) to see why. ",
  "To run without doublet removal anyway, set require_doublets = FALSE in R/utils.R")
if (!has_dbl) warning("scDblFinder not installed: doublets will NOT be removed")
Q <- SET$qc; out <- list(); tab <- list()

for (j in seq_len(nrow(man))) {
  m <- load_counts(man$path[j]); g <- man$dataset[j]; k <- man$sample_key[j]
  mt <- grepl("^MT-", rownames(m))
  nc <- Matrix::colSums(m); nf <- n_detected(m)
  pmt <- if (any(mt)) 100 * Matrix::colSums(m[mt, , drop = FALSE]) / pmax(nc, 1) else rep(0, ncol(m))
  if (!any(mt)) msg(g, k, ": no MT- genes in this matrix (pre-filtered by authors?) -> percent.mt = 0")
  lf <- log10(pmax(nf, 1)); upper <- 10^(median(lf, na.rm = TRUE) + Q$mad_upper * mad(lf, na.rm = TRUE))
  nc <- unname(nc); nf <- unname(nf); pmt <- unname(pmt)
  if (anyNA(nc) || anyNA(nf)) msg(g, k, ":", sum(is.na(nc) | is.na(nf)), "cells with NA counts -> dropped")
  mt_lim <- if (any(mt)) max(Q$mt_floor, median(pmt, na.rm = TRUE) + Q$mt_mad * mad(pmt, na.rm = TRUE)) else Inf
  keep <- nf >= Q$min_features & nf <= upper & nc >= Q$min_counts & pmt < Q$max_mt & pmt <= mt_lim
  keep[is.na(keep)] <- FALSE
  cl <- rep("not_run", ncol(m)); scr <- rep(NA_real_, ncol(m))
  if (!sum(keep)) msg(g, k, ": no cell passed QC")
  if (has_dbl && sum(keep) >= 200) {
    x <- in_memory(m[, which(keep)])
    r <- tryCatch(scDblFinder::scDblFinder(x, returnType = "table"), error = function(e) { msg("  scDblFinder failed:", conditionMessage(e)); NULL })
    if (!is.null(r)) { r <- r[r$type == "real", ]; cl[which(keep)] <- as.character(r[colnames(x), "class"]); scr[which(keep)] <- r[colnames(x), "score"] }
  }
  final <- keep & cl != "doublet"; final[is.na(final)] <- FALSE
  tab[[j]] <- data.frame(row.names = NULL, dataset = g, sample_key = k, n_loaded = ncol(m), n_pass_qc = sum(keep), n_doublets = sum(cl == "doublet"),
    n_kept = sum(final), max_features_cut = round(upper), max_mt_cut = round(min(mt_lim, Q$max_mt), 2), median_features = median(nf[final], na.rm = TRUE), median_counts = median(nc[final], na.rm = TRUE),
    median_pct_mt = round(median(pmt[final], na.rm = TRUE), 2), dropped_sample = sum(final) < Q$min_cells_per_sample)
  if (sum(final) >= Q$min_cells_per_sample)
    out[[j]] <- data.frame(barcode = colnames(m)[final], dataset = g, sample_key = k, nCount = nc[final], nFeature = nf[final],
                           percent.mt = pmt[final], doublet_score = scr[final])
  if (j %% 10 == 0) msg(j, "/", nrow(man), "samples done")
}
tab <- do.call(rbind, tab); qc <- do.call(rbind, out)
qc <- merge(qc, cells, by = c("barcode", "dataset", "sample_key"), all.x = TRUE, sort = FALSE)
write.csv(tab, file.path(DIR$res, "04_QC_per_sample.csv"), row.names = FALSE)
saveRDS(qc, file.path(DIR$proc, "cells_qc.rds"))

s <- aggregate(cbind(n_loaded, n_pass_qc, n_doublets, n_kept) ~ dataset, tab, sum); print(s)
cat("Samples dropped (<", Q$min_cells_per_sample, "cells):", sum(tab$dropped_sample), "\n")
for (v in c("nFeature", "nCount", "percent.mt")) {
  p <- ggplot(qc, aes(x = dataset, y = .data[[v]], fill = dataset)) + geom_violin(scale = "width", show.legend = FALSE) +
    (if (v != "percent.mt") scale_y_log10() else NULL) + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle(paste("After QC:", v))
  ggsave(file.path(DIR$plot, paste0("04_QC_", v, ".png")), p, width = 9, height = 4.5, dpi = 200)
}
cat("\nNext: Rscript R/05_merge_integrate.R\n")
