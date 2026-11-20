# 02_inventory_and_metadata.R
# (1) unpacks the GEO files and reports what format each dataset is in
# (2) lists the sample keys the loader will create
# (3) pulls the GEO sample table and drafts config/sample_metadata.csv
#     (diagnosis / age / sex / PMI / region / enriched guessed from the GEO text)
# -> OPEN config/sample_metadata.csv AND FIX IT BY HAND before running 03. Everything downstream uses it.
source("R/utils.R")
suppressPackageStartupMessages(library(GEOquery))

ds <- cfg_datasets()
inv <- list(); keys <- list()

for (i in seq_len(nrow(ds))) {
  g <- ds$dataset[i]; p <- ds$path[i]; msg("==", g)
  if (ds$loader[i] == "seurat_rds") {
    if (!file.exists(p)) { msg("  missing", p); next }
    o <- readRDS(p); md <- o[[]]; sc <- first_col(colnames(md), ds$sample_col[i])
    if (is.null(sc)) stop(g, ": none of the sample columns '", ds$sample_col[i], "' exist; edit config/datasets.csv")
    k <- unique(as.character(md[[sc]]))
    d1 <- md[!duplicated(md[[sc]]), , drop = FALSE]
    pick <- function(c) { h <- first_col(colnames(d1), c); if (is.null(h)) NA else as.character(d1[[h]]) }
    keys[[g]] <- data.frame(dataset = g, sample_key = as.character(d1[[sc]]), from_object_diagnosis = pick("diagnosis|Diagnosis|condition|disease"),
                            from_object_age = pick("age|Age"), from_object_sex = pick("sex|Sex|gender"), from_object_PMI = pick("PMI_hours|PMI|pmi"))
    inv[[g]] <- data.frame(dataset = g, format = "seurat_rds", n_files = 1, n_samples = length(k),
                           example = paste(head(k, 4), collapse = "; "), n_cells = ncol(o))
    rm(o); gc(); next
  }
  if (!dir.exists(p) && !file.exists(p)) { msg("  missing", p, "- run 01_download_GEO.R"); next }
  ex <- ds$geo_download_exclude_regex[i]
  p <- dataset_dir(g, p, ex); f <- detect_files(p, ex)
  cf <- switch(f$format, h5 = f$h5, mtx = f$mtx, molinfo = f$molinfo, h5ad = f$h5ad, dense_text = f$txt, rds = f$rds, character(0))
  k <- unname(vapply(cf, sample_of_file, "")); bc_example <- NA; obs_cols <- NA
  if (f$format == "mtx" && length(cf) == 1 && !is.na(ds$barcode_sample_regex[i])) {     # one matrix, many samples
    bcf <- list.files(dirname(cf), pattern = "(barcodes|bcd).*\\.tsv(\\.gz)?$", full.names = TRUE)[1]
    bc <- fread(bcf, header = FALSE)[[1]]; bc_example <- paste(head(bc, 3), collapse = "; ")
    k <- unique(sub(paste0(ds$barcode_sample_regex[i], ".*$"), "\\1", bc[grepl(ds$barcode_sample_regex[i], bc)]))
    if (!length(k)) warning(g, ": barcode_sample_regex matched no barcodes; check the example barcodes in 02_inventory.csv")
  }
  if (f$format == "h5ad") {                    # samples and author labels live in obs: read them cheaply with rhdf5
    obs <- h5ad_obs_names(cf[1]); obs_cols <- paste(obs, collapse = "; ")
    sc <- first_col(obs, ds$sample_col[i]); ac <- first_col(obs, ds$author_label_col[i])
    msg("  h5ad obs columns:", obs_cols); msg("  sample column used:", sc, "| author label column:", ac)
    if (!is.null(sc)) k <- unique(h5ad_obs_col(cf[1], sc))
    if (is.null(sc) && !is.na(ds$barcode_sample_regex[i])) { bc <- h5ad_index(cf[1]); re <- ds$barcode_sample_regex[i]
      bc_example <- paste(head(bc, 3), collapse = "; "); k <- unique(sub(paste0(re, ".*$"), "\\1", bc[grepl(re, bc)])) }
    if (!is.null(ac)) { v <- h5ad_obs_col(cf[1], ac)
      write.csv(as.data.frame(table(author_label = v)), file.path(DIR$res, paste0("02_", g, "_author_labels.csv")), row.names = FALSE) }
  }
  if (!length(k)) { msg("  NO count files recognised for", g, "- see results/02_inventory.csv (files column)"); }
  if (length(cf) == 1 && f$format %in% c("rds", "dense_text"))
    warning(g, ": a single ", f$format, " file - samples must come from its cell metadata (set sample_col in datasets.csv)")
  if (length(k)) keys[[g]] <- data.frame(dataset = g, sample_key = unique(k), from_object_diagnosis = NA, from_object_age = NA,
                                         from_object_sex = NA, from_object_PMI = NA)
  inv[[g]] <- data.frame(dataset = g, format = f$format, n_files = length(cf), n_samples = length(unique(k)),
                         example = paste(head(basename(cf), 3), collapse = "; "), n_cells = NA, barcode_example = bc_example,
                         h5ad_obs_columns = obs_cols, files = paste(head(basename(f$all), 8), collapse = "; "))
}
inv <- data.table::rbindlist(inv, fill = TRUE); keys <- do.call(rbind, keys)
write.csv(inv, file.path(DIR$res, "02_inventory.csv"), row.names = FALSE); print(inv)

# ---------- GEO sample tables ----------
guess <- function(txt, pat) { m <- regmatches(txt, regexec(pat, txt, ignore.case = TRUE)); vapply(m, function(x) if (length(x) > 1) trimws(x[length(x)]) else NA_character_, "") }
diag_guess <- function(t) {
  r <- rep(NA_character_, length(t))
  rules <- list(PDD = "parkinson'?s? disease (with )?dementia|\\bPDD\\b", DLB = "dementia with lewy|\\bDLB\\b",
                ILBD = "incidental lewy|\\bILBD\\b", LBD = "lewy body disease|\\bLBD\\b", MSA = "multiple system atrophy|\\bMSA\\b",
                AD = "alzheimer|\\bAD\\b", PD = "parkinson|\\bI?PD\\b", Control = "control|healthy|normal|\\bCTRL?\\b|\\bHC\\b|\\bNC\\b|unaffected")
  for (nm in names(rules)) { hit <- is.na(r) & grepl(rules[[nm]], t, ignore.case = TRUE); r[hit] <- nm }
  r
}
geo <- list()
for (g in unique(keys$dataset)) {
  e <- tryCatch(getGEO(g, GSEMatrix = TRUE, getGPL = FALSE), error = function(e) NULL)
  if (is.null(e)) { msg(g, ": GEO sample table not available"); next }
  pd <- data.table::rbindlist(lapply(e, function(x) { d <- Biobase::pData(x); d[] <- lapply(d, as.character); d }), fill = TRUE)
  ch <- grep("^characteristics|source_name", names(pd), value = TRUE)
  pd$geo_text <- apply(pd[, c("title", ch), with = FALSE], 1, function(x) paste(x[!is.na(x) & x != ""], collapse = "; "))
  geo[[g]] <- data.frame(dataset = g, geo_accession = pd$geo_accession, title = pd$title, geo_text = pd$geo_text)
}
geo <- do.call(rbind, geo)

match_geo <- function(d, k) {
  s <- geo[geo$dataset == d, ]; if (is.null(s) || !nrow(s)) return(NA)
  h <- which(s$geo_accession == k); if (!length(h)) h <- which(s$title == k)
  if (!length(h)) h <- which(vapply(s$title, function(t) grepl(k, t, fixed = TRUE) || grepl(t, k, fixed = TRUE), TRUE))
  if (length(h) == 1) h else NA
}
ds_region <- setNames(ds$region, ds$dataset)
out <- do.call(rbind, lapply(seq_len(nrow(keys)), function(j) {
  d <- keys$dataset[j]; k <- keys$sample_key[j]; h <- match_geo(d, k)
  s <- if (!is.na(h)) geo[geo$dataset == d, ][h, ] else data.frame(geo_accession = NA, title = NA, geo_text = NA)
  t <- paste(s$geo_text, k)
  data.frame(dataset = d, sample_key = k, geo_accession = s$geo_accession, geo_title = s$title,
    diagnosis = ifelse(is.na(keys$from_object_diagnosis[j]), diag_guess(t), keys$from_object_diagnosis[j]),
    age = ifelse(is.na(keys$from_object_age[j]), guess(t, "\\bage[^:;]*:\\s*([0-9.]+)"), keys$from_object_age[j]),
    sex = ifelse(is.na(keys$from_object_sex[j]), guess(t, "(sex|gender)[^:;]*:\\s*(male|female|M|F)\\b"), keys$from_object_sex[j]),
    PMI_hours = ifelse(is.na(keys$from_object_PMI[j]), guess(t, "(pmi|post.?mortem interval)[^:;]*:\\s*([0-9.]+)"), keys$from_object_PMI[j]),
    region = { r <- guess(t, "(brain region|region|tissue)[^:;]*:\\s*([^;]+)"); ifelse(is.na(r), ds_region[[d]], r) },
    donor = NA,
    enriched = grepl("nurr|neun|sorted|FANS|Pos[A-H]?$|_pos\\b", t, ignore.case = TRUE) & !grepl("DAPI", t),
    include = TRUE, geo_text = s$geo_text)
}))
out$diagnosis[out$diagnosis %in% c("HC", "Ctrl", "CTRL", "control", "Control ")] <- "Control"
write.csv(out, file.path(DIR$cfg, "sample_metadata_DRAFT.csv"), row.names = FALSE)
if (!file.exists(file.path(DIR$cfg, "sample_metadata.csv"))) file.copy(file.path(DIR$cfg, "sample_metadata_DRAFT.csv"), file.path(DIR$cfg, "sample_metadata.csv"))

cat("\nSamples per dataset and diagnosis (guessed):\n"); print(table(out$dataset, out$diagnosis, useNA = "ifany"))
cat("\nNot matched to a GEO sample:", sum(is.na(out$geo_accession)), " | diagnosis NA:", sum(is.na(out$diagnosis)), "\n")
cat("\nNOW: open config/sample_metadata.csv, fix diagnosis/age/sex/PMI/region/donor/enriched, set include=FALSE to drop a sample.\n",
    "'donor' = same person sampled twice (e.g. GSE253462 SN + cingulate) -> give both rows the same donor id.\n",
    "Then: Rscript R/03_load_counts.R\n")
