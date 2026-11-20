# 03_load_counts.R — read every sample once, convert genes to symbols, save counts per sample
# Output: data/counts/<dataset>/<sample>     (on-disk BPCells, or .rds if BPCells is missing)
#         data/processed/cells_raw.rds       (one row per cell: barcode, dataset, sample_key, author_label)
#         data/processed/manifest.csv        (one row per sample: where its counts are, n cells)
# Re-running skips datasets already in the manifest (delete their rows to redo one).
source("R/utils.R")
ds <- cfg_datasets(); sm <- cfg_samples(); sm <- sm[as.logical(sm$include), ]
man_f <- file.path(DIR$proc, "manifest.csv"); cell_f <- file.path(DIR$proc, "cells_raw.rds")
man <- if (file.exists(man_f)) read.csv(man_f) else NULL
cells <- if (file.exists(cell_f)) readRDS(cell_f) else NULL
# datasets loaded with a different max_cells_per_sample are loaded again
if (!is.null(man)) {
  if (!"cap" %in% names(man)) man$cap <- NA
  redo <- unique(man$dataset[is.na(man$cap) | man$cap != SET$max_cells_per_sample])
  if (length(redo)) { msg("re-loading (cell cap changed):", paste(redo, collapse = ", "))
    man <- man[!man$dataset %in% redo, ]; cells <- cells[!cells$dataset %in% redo, ] }
}

cap <- function(bcs) if (length(bcs) > SET$max_cells_per_sample) sample(bcs, SET$max_cells_per_sample) else bcs
store <- function(m, g, k, author = NULL) {
  if (!k %in% sm$sample_key[sm$dataset == g]) { msg("  skip", k, "(not in sample_metadata or include=FALSE)"); return(invisible()) }
  keep <- cap(colnames(m)); if (length(keep) < ncol(m)) m <- m[, keep]
  path <- save_counts(m, g, make.names(k))
  man <<- rbind(man, data.frame(dataset = g, sample_key = k, path = path, n_cells = ncol(m), n_genes = nrow(m), cap = SET$max_cells_per_sample))
  al <- if (is.null(author)) rep(NA_character_, ncol(m)) else as.character(author[colnames(m)])
  cells <<- rbind(cells, data.frame(barcode = colnames(m), dataset = g, sample_key = k, author_label = al))
  msg("  ", k, ":", ncol(m), "cells,", nrow(m), "genes")
}
# optional author labels from a separate metadata file (e.g. Single Cell Portal export)
author_from_file <- function(i) {
  f <- ds$author_meta_file[i]; if (is.na(f) || !file.exists(f)) return(NULL)
  x <- fread(f, data.table = FALSE); x <- x[x[[1]] != "TYPE", ]          # SCP files have a TYPE row
  setNames(as.character(x[[ds$author_meta_label_col[i]]]), x[[ds$author_meta_barcode_col[i]]])
}

for (i in seq_len(nrow(ds))) {
  g <- ds$dataset[i]; if (!is.null(man) && g %in% man$dataset) { msg(g, "already loaded"); next }
  msg("==", g); p <- ds$path[i]; ex <- ds$geo_download_exclude_regex[i]

  if (ds$loader[i] == "seurat_rds") {
    o <- readRDS(p); DefaultAssay(o) <- "RNA"
    if (inherits(o[["RNA"]], "Assay5") && !"counts" %in% Layers(o[["RNA"]])) o[["RNA"]] <- JoinLayers(o[["RNA"]])
    cnt <- to_symbols(LayerData(o, assay = "RNA", layer = "counts")); md <- o[[]]
    sc <- first_col(colnames(md), ds$sample_col[i]); ac <- first_col(colnames(md), ds$author_label_col[i])
    if (is.null(ac)) msg("  no author label column found")
    au <- if (!is.null(ac)) setNames(as.character(md[[ac]]), rownames(md)) else NULL
    for (k in unique(as.character(md[[sc]]))) store(cnt[, rownames(md)[md[[sc]] == k], drop = FALSE], g, k, au)
    rm(o, cnt); gc()

  } else {
    p <- dataset_dir(g, p, ex); f <- detect_files(p, ex); au <- author_from_file(i)
    if (f$format == "h5") for (h in f$h5) store(read_h5_10x(h), g, sample_of_file(h), au)
    if (f$format == "molinfo") for (h in f$molinfo) store(read_molinfo(h), g, sample_of_file(h), au)
    if (f$format == "mtx") {
      if (length(f$mtx) == 1 && !is.na(ds$barcode_sample_regex[i])) {          # one matrix, many samples
        re <- ds$barcode_sample_regex[i]
        if (!SET$use_bpcells && file.size(f$mtx) > 1.5e9) {                      # too big for RAM: stream chosen cells
          if (!is.finite(SET$max_cells_per_sample)) stop(g, ": the full run needs BPCells (matrix is ", round(file.size(f$mtx) / 1e9, 1), " GB)")
          bc <- mtx_side_files(f$mtx)$bcs
          smp <- ifelse(grepl(re, bc), sub(paste0(re, ".*$"), "\\1", bc), NA)
          ok_s <- intersect(unique(na.omit(smp)), sm$sample_key[sm$dataset == g])
          keep <- unlist(lapply(ok_s, function(k) { x <- which(smp == k); if (length(x) > SET$max_cells_per_sample) sample(x, SET$max_cells_per_sample) else x }))
          m <- read_mtx_subset(f$mtx, keep)
        } else m <- read_mtx_triplet(f$mtx)
        smp <- ifelse(grepl(re, colnames(m)), sub(paste0(re, ".*$"), "\\1", colnames(m)), NA)
        for (k in unique(na.omit(smp))) store(m[, which(smp == k)], g, k, au)
      } else for (x in f$mtx) store(read_mtx_triplet(x), g, sample_of_file(x), au)
    }
    if (f$format == "dense_text") for (x in f$txt) store(read_dense_text(x), g, sample_of_file(x), au)
    if (f$format == "h5ad" && h5ad_is_dense(f$h5ad[1]) && !is.na(ds$barcode_sample_regex[i])) {   # big dense h5ad
      x <- f$h5ad[1]; bc <- h5ad_index(x); re <- ds$barcode_sample_regex[i]
      smp <- ifelse(grepl(re, bc), sub(paste0(re, ".*$"), "\\1", bc), NA)
      for (k in unique(na.omit(smp))) {
        if (!k %in% sm$sample_key[sm$dataset == g]) next
        idx <- which(smp == k); if (length(idx) > SET$max_cells_per_sample) idx <- sort(sample(idx, SET$max_cells_per_sample))
        msg("  reading", length(idx), "cells of", k, "from the dense h5ad"); store(read_h5ad_dense(x, idx), g, k, au)
      }
    } else if (f$format %in% c("h5ad", "rds")) for (x in c(f$h5ad, f$rds)) {
      r <- if (f$format == "h5ad") read_h5ad_counts(x) else { o <- readRDS(x)
             if (inherits(o, "Seurat")) list(counts = to_symbols(LayerData(JoinLayers(o), layer = "counts")), meta = o[[]])
             else list(counts = to_symbols(as(o, "dgCMatrix")), meta = NULL) }
      sc <- if (!is.null(r$meta)) first_col(colnames(r$meta), ds$sample_col[i]) else NULL
      ac <- if (!is.null(r$meta)) first_col(colnames(r$meta), ds$author_label_col[i]) else NULL
      a2 <- if (!is.null(ac)) setNames(as.character(r$meta[[ac]]), rownames(r$meta)) else au
      if (is.null(sc)) store(r$counts, g, sample_of_file(x), a2)
      else for (k in unique(as.character(r$meta[[sc]]))) store(r$counts[, rownames(r$meta)[r$meta[[sc]] == k], drop = FALSE], g, k, a2)
    }
    if (f$format == "none") msg("  no count files found in", p)
  }
  write.csv(man, man_f, row.names = FALSE); saveRDS(cells, cell_f)          # checkpoint after every dataset
}
cat("\nCells per dataset:\n"); print(aggregate(n_cells ~ dataset, man, sum))
miss <- setdiff(paste(sm$dataset, sm$sample_key), paste(man$dataset, man$sample_key))
if (length(miss)) cat("\nIn sample_metadata but not loaded (check names):\n", paste(miss, collapse = "\n "), "\n")
cat("\nNext: Rscript R/04_qc_doublets.R\n")
