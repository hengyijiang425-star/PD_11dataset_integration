# utils.R — shared settings and helpers for the 11-dataset PD integration
# Every numbered script starts with: source("R/utils.R")
# Run all scripts from the project root (the folder that holds PD_11dataset_integration.Rproj).

suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(ggplot2); library(data.table)
})
set.seed(20261110)
options(future.globals.maxSize = 200 * 1024^3, timeout = 36000, stringsAsFactors = FALSE)

# ============================ SETTINGS (edit here) ============================
SET <- list(
  # quick laptop test: e.g. 300 cells per sample. Inf = all cells (run on the server)
  max_cells_per_sample = 300,     # TEST RUN. Full run: Inf
  # "sketch": integrate a representative subset, then project every cell (needed for >1M nuclei)
  # "full"  : integrate all cells directly (only if RAM allows, or together with max_cells_per_sample)
  integration_mode = "full",      # TEST RUN. Full run: "sketch"
  integrate_cells_per_sample = 200,   # TEST RUN (laptop RAM): cells per DONOR used from step 05 on. Full run: Inf
  sketch_cells_per_sample = 800,   # cells kept per sample in the sketch
  use_bpcells = requireNamespace("BPCells", quietly = TRUE),   # on-disk counts; strongly recommended
  # one QC rule for every sample. Tightened 2026-11-17 (PI: "QC 拉紧"):
  #   genes >= 500, UMIs >= 1000, mito < 10% AND <= sample median + 3 MAD (floor 1%),
  #   genes <= sample median + 3 MAD (log scale), doublets removed (required)
  # mito stays 10% absolute: at 5% GSE184950 (Wang SN) loses most samples; the per-sample MAD rule tightens it where data are clean
  qc = list(min_features = 500, min_counts = 1000, max_mt = 10, mt_mad = 3, mt_floor = 1, mad_upper = 3,
            min_cells_per_sample = 100, require_doublets = TRUE),
  n_hvg = 3000, n_pcs = 30,
  res_main = 0.4, res_all = c(0.2, 0.4, 0.8),
  res_sub = 0.5,                   # resolution inside each lineage (neurons, myeloid, vascular)
  min_pct = 0.30,                  # a panel marker "supports" a call if detected in >=30% of the cluster
  min_support = 2,                 # ...and at least 2 panel markers must support it
  knn_k = 15,                      # label transfer from sketch to all cells
  min_shared_genes = 5000
)
DIR <- list(raw = "data/raw", proc = "data/processed", counts = "data/counts",
            res = "results", plot = "plots", cfg = "config")
for (d in DIR) dir.create(d, recursive = TRUE, showWarnings = FALSE)

cfg_datasets <- function() {
  x <- read.csv(file.path(DIR$cfg, "datasets.csv"), na.strings = c("", "NA"))
  x[as.logical(x$include), ]
}
cfg_map    <- function() read.csv(file.path(DIR$cfg, "celltype_harmonization_map_11datasets.csv"), na.strings = c("", "NA"))
cfg_panels <- function() read.csv(file.path(DIR$cfg, "marker_panels_v2.csv"))
cfg_samples <- function() {
  f <- file.path(DIR$cfg, "sample_metadata.csv")
  if (!file.exists(f)) stop("config/sample_metadata.csv missing: run R/02_inventory_and_metadata.R, review the draft, save it as config/sample_metadata.csv")
  read.csv(f, na.strings = c("", "NA"))
}
first_col <- function(cols, candidates) {           # candidates like "sample|patient|orig.ident"
  if (is.null(candidates) || is.na(candidates)) return(NULL)
  h <- intersect(strsplit(candidates, "\\|")[[1]], cols); if (length(h)) h[1] else NULL
}
msg <- function(...) cat(format(Sys.time(), "%H:%M:%S"), "|", ..., "\n")

# ============================ count storage ============================
# Each sample is saved once as data/counts/<dataset>/<sample> (BPCells dir) or <sample>.rds (dgCMatrix)
save_counts <- function(m, dataset, sample) {
  d <- file.path(DIR$counts, dataset); dir.create(d, recursive = TRUE, showWarnings = FALSE)
  if (SET$use_bpcells) {
    out <- file.path(d, sample); if (dir.exists(out)) unlink(out, recursive = TRUE)
    if (!inherits(m, "IterableMatrix")) m <- as(m, "dgCMatrix")
    BPCells::write_matrix_dir(BPCells::convert_matrix_type(m, "uint32_t"), out, overwrite = TRUE)
    out
  } else { out <- file.path(d, paste0(sample, ".rds")); saveRDS(as(m, "dgCMatrix"), out); out }
}
load_counts <- function(path) if (dir.exists(path)) BPCells::open_matrix_dir(path) else readRDS(path)
n_detected <- function(m) if (inherits(m, "IterableMatrix")) Matrix::colSums(BPCells::binarize(m, threshold = 0)) else Matrix::colSums(m > 0)
in_memory <- function(m) if (inherits(m, "IterableMatrix")) as(m, "dgCMatrix") else m

# ============================ gene IDs -> symbols ============================
# Ensembl IDs (with or without version) are mapped with org.Hs.eg.db; duplicated symbols are summed
# (in-memory) or the first copy kept (on-disk). Mitochondrial genes keep the "MT-" prefix.
to_symbols <- function(m, feature_symbols = NULL) {
  rn <- rownames(m)
  if (!is.null(feature_symbols)) rn <- ifelse(is.na(feature_symbols) | feature_symbols == "", rn, feature_symbols)
  ens <- grepl("^ENSG[0-9]+", rn)
  if (mean(ens) > 0.5) {
    if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) stop("Install org.Hs.eg.db to map Ensembl IDs")
    key <- sub("\\..*$", "", rn[ens])
    sym <- suppressMessages(AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, key, "SYMBOL", "ENSEMBL"))
    rn[ens] <- ifelse(is.na(sym), key, sym)
  }
  rn <- sub("^MT\\.", "MT-", rn)                      # some files write MT.ND1
  rownames(m) <- rn
  if (anyDuplicated(rn)) {
    if (inherits(m, "IterableMatrix")) m <- m[!duplicated(rn), ]
    else m <- as(Matrix::fac2sparse(factor(rn, levels = unique(rn))) %*% m, "dgCMatrix")
  }
  m
}

# ============================ readers for GEO supplementary files ============================
# Unpack *_RAW.tar and nested .tar.gz / .tgz once (only files not matching the exclude regex)
unpack_all <- function(dir, exclude = NA) {
  for (i in 1:3) {
    tars <- list.files(dir, pattern = "\\.(tar|tar\\.gz|tgz)$", recursive = TRUE, full.names = TRUE)
    if (!length(tars)) break
    tars <- tars[!file.exists(paste0(tars, ".unpacked"))]
    if (!length(tars)) break
    for (t in tars) {
      inside <- untar(t, list = TRUE)
      keep <- if (is.na(exclude)) inside else inside[!grepl(exclude, inside)]
      ex <- file.path(dirname(t), sub("\\.(tar|tar\\.gz|tgz)$", "", basename(t)))
      dir.create(ex, showWarnings = FALSE); msg("untar", basename(t), "->", length(keep), "files")
      if (length(keep)) untar(t, files = keep, exdir = ex)
      file.create(paste0(t, ".unpacked"))
    }
  }
  for (z in list.files(dir, pattern = "\\.h5ad\\.gz$", recursive = TRUE, full.names = TRUE)) {   # h5ad must be unzipped
    if (!is.na(exclude) && grepl(exclude, z)) next
    out <- sub("\\.gz$", "", z); done <- paste0(out, ".gunzip_done")
    if (!file.exists(done)) {                      # write to a temp name, rename only when complete
      msg("gunzip", basename(z), "(can take a few minutes)")
      tmp <- paste0(out, ".part"); st <- system2("gzip", c("-dc", shQuote(z)), stdout = tmp)
      if (st != 0) stop("gunzip failed for ", z, " - re-download it with 01_download_GEO.R")
      file.rename(tmp, out); file.create(done)
    }
  }
}

# Folder that holds a dataset's raw files. If config path is a .tar (e.g. the GSE184950_RAW.tar of the
# 3-dataset project), it is unpacked into data/raw/<dataset>/ so the old project folder is never touched.
dataset_dir <- function(dataset, path, exclude = NA) {
  if (file.exists(path) && !dir.exists(path) && grepl("\\.tar$", path)) {
    d <- file.path(DIR$raw, dataset); dir.create(d, recursive = TRUE, showWarnings = FALSE)
    flag <- file.path(d, ".from_tar_done")
    if (!file.exists(flag)) {
      inside <- untar(path, list = TRUE); keep <- if (is.na(exclude)) inside else inside[!grepl(exclude, inside)]
      msg("untar", basename(path), "->", d); untar(path, files = keep, exdir = d); file.create(flag)
    }
    path <- d
  }
  unpack_all(path, exclude); path
}

# Find count files and decide the format of a dataset folder
detect_files <- function(dir, exclude = NA) {
  f <- list.files(dir, recursive = TRUE, full.names = TRUE)
  f <- f[!grepl("\\.unpacked$|\\.from_tar_done$|\\.gunzip_done$|\\.part$|\\.(tar|tar\\.gz|tgz)$|\\.h5ad\\.gz$", f)]
  if (!is.na(exclude)) f <- f[!grepl(exclude, f)]
  molinfo <- f[grepl("molecule_info\\.h5$", f)]
  h5   <- f[grepl("\\.h5$", f) & !grepl("raw_feature_bc_matrix|molecule_info", f)]
  mtx  <- f[grepl("(^|_|/)matrix\\.mtx(\\.gz)?$|\\.mtx(\\.gz)?$", f)]
  h5ad <- f[grepl("\\.h5ad$", f)]
  rds  <- f[grepl("\\.(rds|RDS)$", f)]
  txt  <- f[grepl("\\.(csv|tsv|txt)(\\.gz)?$", f) &
            !grepl("barcode|bcd|feature|genes\\.tsv|meta|annot|cluster", basename(f), ignore.case = TRUE)]
  named <- txt[grepl("count|umi|matrix|expr", basename(txt), ignore.case = TRUE)]
  txt <- if (length(named)) named else txt[grepl("^GSM[0-9]+", basename(txt))]   # e.g. GSM8021038_OJ7_nigra.txt.gz
  fmt <- if (length(h5)) "h5" else if (length(mtx)) "mtx" else if (length(molinfo)) "molinfo" else if (length(h5ad)) "h5ad" else
         if (length(txt)) "dense_text" else if (length(rds)) "rds" else "none"
  list(format = fmt, h5 = h5, mtx = mtx, molinfo = molinfo, h5ad = h5ad, txt = txt, rds = rds, all = f)
}

# sample id from a file name: GSM accession if present, else the file prefix
sample_of_file <- function(path) {
  b <- basename(path); g <- regmatches(b, regexpr("GSM[0-9]+", b))
  if (length(g)) return(g)
  g <- regmatches(path, regexpr("GSM[0-9]+", path)); if (length(g)) return(g)
  k <- sub("[_.-]?(filtered_feature_bc_matrix|matrix|counts?|umi|UMI|expression).*$", "", b)
  if (!nzchar(k)) k <- basename(dirname(path))
  if (k %in% c("filtered_feature_bc_matrix", "outs")) k <- basename(dirname(dirname(path)))
  k
}

# 10x triplet (matrix.mtx + barcodes + features/genes) in the same folder, or as <prefix>_matrix.mtx.gz
# Without BPCells a 5 GB matrix (Kamath) cannot be read into RAM. For a capped run we stream the file once
# through gzip + awk and keep only the chosen cells (columns), so memory stays small. keep = column numbers (1-based).
mtx_side_files <- function(mtx) {
  d <- dirname(mtx); pre <- sub("matrix\\.mtx(\\.gz)?$", "", basename(mtx))
  pick <- function(pats) { for (p in pats) { x <- list.files(d, pattern = paste0("^", gsub("([.|()^${}+*?\\[\\]])", "\\\\\\1", pre), p), full.names = TRUE); if (length(x)) return(x[1]) }; NA }
  bc <- pick(c("barcodes\\.tsv(\\.gz)?$", "bcd\\.tsv(\\.gz)?$", "barcodes.*"))
  ft <- pick(c("features\\.tsv(\\.gz)?$", "genes\\.tsv(\\.gz)?$", "features.*", "genes.*"))
  if (is.na(bc) || is.na(ft)) stop("Cannot find barcodes/features next to ", mtx)
  list(feats = fread(ft, header = FALSE, sep = "\t"), bcs = fread(bc, header = FALSE, sep = "\t")[[1]])
}
read_mtx_subset <- function(mtx, keep) {
  sf <- mtx_side_files(mtx); keep <- sort(unique(keep))
  kf <- tempfile(fileext = ".tsv"); out <- tempfile(fileext = ".tsv"); on.exit(unlink(c(kf, out)))
  writeLines(paste(keep, seq_along(keep), sep = "\t"), kf)
  awk <- paste0("BEGIN{while((getline l < \"", kf, "\")>0){split(l,a,\"\\t\"); k[a[1]]=a[2]}} ",
                "/^%/{next} !h{h=1; next} ($2 in k){print $1\"\\t\"k[$2]\"\\t\"$3}")
  cmd <- paste(if (grepl("\\.gz$", mtx)) "gzip -dc" else "cat", shQuote(mtx), "| awk", shQuote(awk), ">", shQuote(out))
  msg("  streaming", basename(mtx), "for", length(keep), "cells (one pass over the file; can take 15-40 min)")
  st <- system(cmd); if (st != 0) stop("streaming the matrix failed")
  t <- fread(out, header = FALSE, col.names = c("i", "j", "x"))
  m <- sparseMatrix(i = t$i, j = t$j, x = as.numeric(t$x), dims = c(nrow(sf$feats), length(keep)))
  rownames(m) <- sf$feats[[1]]; colnames(m) <- sf$bcs[keep]
  sym <- if (ncol(sf$feats) >= 2) sf$feats[[2]] else NULL
  if (ncol(sf$feats) >= 3) { ge <- sf$feats[[3]] == "Gene Expression"; if (any(!ge)) { m <- m[which(ge), ]; sym <- sym[ge] } }
  to_symbols(as(m, "CsparseMatrix"), sym)
}
read_mtx_triplet <- function(mtx) {
  d <- dirname(mtx); pre <- sub("matrix\\.mtx(\\.gz)?$", "", basename(mtx))
  pick <- function(pats) { for (p in pats) { x <- list.files(d, pattern = paste0("^", gsub("([.|()^${}+*?\\[\\]])", "\\\\\\1", pre), p), full.names = TRUE); if (length(x)) return(x[1]) }; NA }
  bc <- pick(c("barcodes\\.tsv(\\.gz)?$", "bcd\\.tsv(\\.gz)?$", "barcodes.*"))
  ft <- pick(c("features\\.tsv(\\.gz)?$", "genes\\.tsv(\\.gz)?$", "features.*", "genes.*"))
  if (is.na(bc) || is.na(ft)) stop("Cannot find barcodes/features next to ", mtx)
  feats <- fread(ft, header = FALSE, sep = "\t"); bcs <- fread(bc, header = FALSE, sep = "\t")[[1]]
  sym <- if (ncol(feats) >= 2) feats[[2]] else NULL
  big <- file.size(mtx) > 1.5e9
  if (big && SET$use_bpcells) {
    tmp <- file.path(DIR$counts, "_import", basename(d)); unlink(tmp, recursive = TRUE)
    m <- BPCells::import_matrix_market(mtx, outdir = tmp, row_names = feats[[1]], col_names = bcs)
  } else {
    if (big) warning(basename(mtx), " is large; install BPCells to avoid running out of memory")
    m <- as(Matrix::readMM(mtx), "dgCMatrix"); rownames(m) <- feats[[1]]; colnames(m) <- bcs
  }
  if (ncol(feats) >= 3) {                                     # multiome: keep Gene Expression rows only
    ge <- feats[[3]] == "Gene Expression"; if (any(!ge)) { m <- m[which(ge), ]; sym <- sym[ge] }
  }
  if (grepl("raw_feature_bc_matrix", mtx)) {                  # unfiltered Cell Ranger output: drop empty droplets
    nf <- n_detected(m); m <- m[, which(nf >= 200)]; msg("  raw matrix: kept", ncol(m), "barcodes with >= 200 genes")
  }
  to_symbols(m, sym)
}
read_h5_10x <- function(h5) {
  if (SET$use_bpcells) {
    m <- tryCatch(BPCells::open_matrix_10x_hdf5(h5, feature_type = "Gene Expression"), error = function(e) NULL)
    if (!is.null(m)) {                                           # rownames are gene IDs -> use names
      nm <- tryCatch(rhdf5::h5read(h5, "matrix/features/name"), error = function(e) NULL)
      ft <- tryCatch(rhdf5::h5read(h5, "matrix/features/feature_type"), error = function(e) NULL)
      if (!is.null(nm) && !is.null(ft)) nm <- nm[ft == "Gene Expression"]
      if (!is.null(nm) && length(nm) == nrow(m)) return(to_symbols(m, nm))
      return(to_symbols(m))
    }
  }
  x <- Read10X_h5(h5, use.names = TRUE); if (is.list(x)) x <- x[["Gene Expression"]]
  to_symbols(x)
}
# Dense gene x cell table (one per sample). Read in blocks of genes so a 36k x 15k table never sits in RAM as dense.
read_dense_text <- function(f, block = 3000) {
  tmp <- f
  if (grepl("\\.gz$", f)) { tmp <- tempfile(fileext = ".txt"); R.utils::gunzip(f, destname = tmp, remove = FALSE) }
  on.exit(if (tmp != f) unlink(tmp))
  hdr <- scan(tmp, what = "", nlines = 1, quiet = TRUE, sep = if (grepl("\\.csv$", sub("\\.gz$", "", f))) "," else "")
  n_lines <- R.utils::countLines(tmp) - 1
  parts <- list(); genes <- character(0); start <- 0
  while (start < n_lines) {
    x <- fread(tmp, skip = start + 1, nrows = min(block, n_lines - start), header = FALSE, data.table = FALSE)
    if (start == 0 && ncol(x) == length(hdr)) hdr <- hdr[-1]     # header also names the gene column
    g <- x[[1]]; x <- x[, -1, drop = FALSE]
    x <- as.matrix(x)
    if (start == 0 && any(abs(x - round(x)) > 1e-6, na.rm = TRUE)) warning(basename(f), ": values are not integer counts (normalized?)")
    parts[[length(parts) + 1]] <- as(Matrix(x, sparse = TRUE), "dgCMatrix"); genes <- c(genes, as.character(g))
    start <- start + nrow(x)
  }
  m <- do.call(rbind, parts); rownames(m) <- genes; colnames(m) <- make.unique(hdr)
  to_symbols(m)
}

# Cell Ranger (v3+) molecule_info.h5 -> UMI counts for the barcodes Cell Ranger called as cells (pass_filter)
read_molinfo <- function(f) {
  rd <- function(x) as.vector(rhdf5::h5read(f, x, bit64conversion = "double"))   # plain vectors (h5read gives 1-d arrays)
  ls <- rhdf5::h5ls(f)$name
  if (!"barcode_info" %in% ls) stop(basename(f), ": not a Cell Ranger v3 molecule_info.h5 (no barcode_info/pass_filter)")
  bcs <- as.character(rd("barcodes")); pf <- matrix(as.numeric(rd("barcode_info/pass_filter")), nrow = 3)
  cell_idx <- sort(unique(pf[1, ]))                       # 0-based barcode indices called as cells
  bi <- as.numeric(rd("barcode_idx")); keep <- bi %in% cell_idx
  fi <- as.numeric(rd("feature_idx"))[keep]
  ok <- if ("umi_type" %in% ls) as.numeric(rd("umi_type"))[keep] == 1 else rep(TRUE, sum(keep))
  j <- match(bi[keep], cell_idx); rm(bi)
  ft <- as.character(rd("features/feature_type")); nm <- as.character(rd("features/name")); id <- as.character(rd("features/id"))
  rhdf5::h5closeAll()
  m <- Matrix::sparseMatrix(i = fi[ok] + 1, j = j[ok], x = rep(1, sum(ok)), dims = c(length(nm), length(cell_idx)))  # one row per UMI -> summed
  rownames(m) <- id; colnames(m) <- paste0(bcs[cell_idx + 1], "-1")
  ge <- which(ft == "Gene Expression"); m <- m[ge, , drop = FALSE]
  msg("  molecule_info:", ncol(m), "cells called by Cell Ranger")
  to_symbols(as(m, "CsparseMatrix"), nm[ge])
}

# h5ad obs columns without loading the matrix
h5ad_obs_names <- function(f) { x <- rhdf5::h5ls(f, recursive = 2); n <- x$name[x$group == "/obs"]; rhdf5::h5closeAll(); setdiff(n, c("_index", "__categories")) }
h5ad_obs_col <- function(f, col) {
  x <- rhdf5::h5ls(f, recursive = 3); kids <- x$name[x$group == paste0("/obs/", col)]
  v <- if (all(c("categories", "codes") %in% kids)) {
         cat <- rhdf5::h5read(f, paste0("obs/", col, "/categories")); cod <- rhdf5::h5read(f, paste0("obs/", col, "/codes")); cat[cod + 1]
       } else if (paste0("/obs/__categories") %in% x$group && col %in% x$name[x$group == "/obs/__categories"]) {
         cat <- rhdf5::h5read(f, paste0("obs/__categories/", col)); cat[rhdf5::h5read(f, paste0("obs/", col)) + 1]
       } else rhdf5::h5read(f, paste0("obs/", col))
  rhdf5::h5closeAll(); as.character(v)
}

# h5ad whose X is a plain dense matrix (GSE235330: 151k cells x 58k genes, 35 GB): read only the cells we need,
# in blocks, and keep them sparse. idx = row numbers (cells) in obs order.
h5ad_is_dense <- function(f) { x <- rhdf5::h5ls(f, recursive = 1); rhdf5::h5closeAll(); x$otype[x$name == "X"] == "H5I_DATASET" }
h5ad_index <- function(f, what = "obs") { v <- as.character(rhdf5::h5read(f, paste0(what, "/_index"))); rhdf5::h5closeAll(); v }
read_h5ad_dense <- function(f, idx, block = 2000) {
  genes <- h5ad_index(f, "var"); cells <- h5ad_index(f, "obs"); idx <- sort(idx); parts <- list()
  for (b in split(idx, ceiling(seq_along(idx) / block))) {
    rng <- min(b):max(b)
    x <- if (length(rng) > 4 * length(b)) rhdf5::h5read(f, "X", index = list(NULL, b))
         else rhdf5::h5read(f, "X", index = list(NULL, rng))[, match(b, rng), drop = FALSE]
    parts[[length(parts) + 1]] <- as(Matrix(x, sparse = TRUE), "dgCMatrix"); rm(x); gc(verbose = FALSE)
  }
  rhdf5::h5closeAll()
  m <- do.call(cbind, parts); rownames(m) <- genes; colnames(m) <- cells[idx]; to_symbols(m)
}

# h5ad: use raw integer counts (layers['counts'], else .raw, else X if integer)
read_h5ad_counts <- function(f) {
  if (!requireNamespace("zellkonverter", quietly = TRUE)) stop("Install zellkonverter (Bioconductor) to read ", f)
  sce <- zellkonverter::readH5AD(f, reader = "R", raw = TRUE)
  an <- SummarizedExperiment::assayNames(sce); isint <- function(m) { v <- m[, seq_len(min(200, ncol(m)))]; v <- v@x[seq_len(min(1e5, length(v@x)))]; all(abs(v - round(v)) < 1e-6) }
  m <- NULL
  if ("counts" %in% an) m <- SummarizedExperiment::assay(sce, "counts")
  if (is.null(m) && "raw" %in% SingleCellExperiment::altExpNames(sce)) {
    r <- SingleCellExperiment::altExp(sce, "raw"); r <- as(SummarizedExperiment::assay(r, 1), "dgCMatrix"); if (isint(r)) m <- r }
  if (is.null(m)) { m <- SummarizedExperiment::assay(sce, an[1]); if (!isint(as(m, "dgCMatrix"))) warning(basename(f), ": X is not integer and no raw counts found") }
  list(counts = to_symbols(as(m, "dgCMatrix")), meta = as.data.frame(SummarizedExperiment::colData(sce)))
}

# ============================ multi-marker annotation ============================
# panels: named list of gene vectors. Returns a per-cluster table with the call, support and confidence.
top_de <- function(o, col, n_per = 2000) {
  cl <- as.character(o[[col]][, 1])
  cells <- unlist(lapply(split(colnames(o), cl), function(x) if (length(x) > n_per) sample(x, n_per) else x))
  s <- subset(o, cells = cells); X <- in_memory(LayerData(s, layer = "data")); g <- as.character(s[[col]][cells, 1])
  if (requireNamespace("presto", quietly = TRUE)) {
    de <- presto::wilcoxauc(X, g)
    de <- de[de$logFC > 0.25 & de$padj < 0.05 & de$pct_in > 10, c("group", "feature", "logFC", "auc", "pct_in", "pct_out", "padj")]
    names(de)[1:2] <- c("cluster", "gene"); de[order(de$cluster, -de$logFC), ]
  } else {
    Idents(s) <- col
    de <- FindAllMarkers(s, only.pos = TRUE, logfc.threshold = 0.25, min.pct = 0.1, verbose = FALSE)
    de <- de[de$p_val_adj < 0.05, ]; de[order(de$cluster, -de$avg_log2FC), ]
  }
}
# markers that are up in the cluster within >= min_ds datasets (computed per dataset, then intersected)
conserved_markers <- function(o, col, min_ds = 2, n_per = 1000) {
  if (!requireNamespace("presto", quietly = TRUE)) return(NULL)
  out <- list()
  for (ds in unique(o$dataset)) {
    cells <- colnames(o)[o$dataset == ds]
    if (length(unique(o[[col]][cells, 1])) < 2) next
    cl <- as.character(o[[col]][cells, 1])
    cells <- unlist(lapply(split(cells, cl), function(x) if (length(x) > n_per) sample(x, n_per) else x))
    if (length(cells) < 50) next
    de <- presto::wilcoxauc(in_memory(LayerData(o, layer = "data")[, cells]), as.character(o[[col]][cells, 1]))
    de <- de[de$logFC > 0.25 & de$padj < 0.05 & de$pct_in > 10, ]; de$dataset <- ds
    out[[ds]] <- de[, c("group", "feature", "logFC", "dataset")]
  }
  if (!length(out)) return(NULL)
  de <- do.call(rbind, out)
  agg <- aggregate(cbind(n_datasets = dataset) ~ group + feature, de, function(x) length(unique(x)))
  lf  <- aggregate(logFC ~ group + feature, de, mean)
  agg <- merge(agg, lf); agg <- agg[agg$n_datasets >= min_ds, ]
  agg <- agg[order(agg$group, -agg$n_datasets, -agg$logFC), ]; names(agg)[1:2] <- c("cluster", "gene"); agg
}

annotate_clusters <- function(o, col, panels, tag, author_col = NULL, call_map = NULL) {
  # author_col: metadata column with author labels already mapped to unified names (e.g. "author_L1")
  # call_map  : named vector panel -> unified name at the same level, used to compare with authors
  panels <- lapply(panels, function(g) intersect(g, rownames(o))); panels <- panels[lengths(panels) >= 2]
  genes <- unique(unlist(panels)); cl <- as.character(o[[col]][, 1])
  X <- in_memory(LayerData(o, layer = "data")[genes, , drop = FALSE])
  pct <- sapply(split(seq_along(cl), cl), function(i) Matrix::rowMeans(X[, i, drop = FALSE] > 0))   # genes x clusters
  o <- AddModuleScore(o, features = panels, name = "pnl_", seed = 1)
  sc <- aggregate(o[[]][, paste0("pnl_", seq_along(panels)), drop = FALSE], by = list(cluster = cl), FUN = mean)
  colnames(sc)[-1] <- names(panels)
  de <- top_de(o, col); cm <- conserved_markers(o, col)
  au <- if (!is.null(author_col) && author_col %in% colnames(o[[]])) as.character(o[[author_col]][, 1]) else rep(NA, ncol(o))
  dsl <- sort(unique(o$dataset))
  rows <- lapply(sc$cluster, function(k) {
    v <- unlist(sc[sc$cluster == k, names(panels)]); ord <- order(v, decreasing = TRUE)
    best <- names(panels)[ord[1]]; second <- if (length(ord) > 1) names(panels)[ord[2]] else NA
    p <- pct[panels[[best]], k]; support <- sum(p >= SET$min_pct)
    call <- if (support >= SET$min_support) best else "Unresolved"
    a <- au[cl == k]; a <- a[!is.na(a) & !grepl("mixed|unsplit|TBD", a)]
    amaj <- if (length(a)) names(sort(table(a), decreasing = TRUE))[1] else NA
    afr  <- if (length(a)) round(mean(a == amaj), 2) else NA
    call_u <- if (!is.null(call_map) && call %in% names(call_map)) call_map[[call]] else call
    agree <- if (is.na(amaj)) NA else amaj == call_u
    conf <- if (call == "Unresolved") "Low" else if (support >= 3 && !isFALSE(agree) && (v[ord[1]] - v[ord[min(2, length(ord))]]) > 0.1) "High" else
            if (!isFALSE(agree)) "Medium" else "Low"
    sm <- table(o$sample_uid[cl == k]); dsf <- round(prop.table(table(factor(o$dataset[cl == k], levels = dsl))), 2)
    cmk <- if (!is.null(cm)) head(cm$gene[cm$cluster == k], 10) else character(0)
    data.frame(cluster = k, n_cells = sum(cl == k), call = call, confidence = conf,
      score_best = round(v[ord[1]], 3), second = second, margin = round(v[ord[1]] - v[ord[min(2, length(ord))]], 3),
      panel_markers = paste0(names(p), "(", round(100 * p), "%)", collapse = ", "), n_panel_markers_supporting = support,
      top10_DE_markers = paste(head(de$gene[de$cluster == k], 10), collapse = ", "),
      conserved_markers_2plus_datasets = paste(cmk, collapse = ", "),
      author_majority = amaj, author_fraction = afr, agrees_with_author = agree,
      n_datasets_ge5pct = sum(dsf >= 0.05), datasets = paste0(names(dsf)[dsf > 0], "=", dsf[dsf > 0], collapse = "; "),
      top_donor = names(which.max(sm)), top_donor_fraction = round(max(sm) / sum(sm), 2),
      mean_doublet_score = round(mean(o$doublet_score[cl == k], na.rm = TRUE), 3),
      median_nFeature = round(median(o$nFeature[cl == k], na.rm = TRUE)), median_pct_mt = round(median(o$percent.mt[cl == k], na.rm = TRUE), 2),
      frac_midbrain = if ("region_group" %in% colnames(o[[]])) round(mean(o$region_group[cl == k] == "Midbrain"), 2) else NA)
  })
  tab <- do.call(rbind, rows); tab <- tab[order(suppressWarnings(as.numeric(tab$cluster))), ]
  write.csv(tab, file.path(DIR$res, paste0(tag, "_cluster_annotation.csv")), row.names = FALSE)
  write.csv(de, file.path(DIR$res, paste0(tag, "_cluster_DE_markers.csv")), row.names = FALSE)
  if (!is.null(cm)) write.csv(cm, file.path(DIR$res, paste0(tag, "_cluster_conserved_markers.csv")), row.names = FALSE)
  list(table = tab, panels = panels)
}

dot_two_ways <- function(o, col_cluster, col_annot, panels, tag) {
  feat <- panels[order(names(panels))]
  thm <- theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8), strip.text.x = element_text(size = 7))
  w <- max(14, 0.32 * length(unlist(feat)) + 4)
  sv <- function(p, f, n) ggsave(file.path(DIR$plot, f), p, width = w, height = 0.32 * n + 4, dpi = 200, bg = "white", limitsize = FALSE)
  sv(DotPlot(o, features = feat, group.by = col_cluster) + thm + ggtitle(paste(tag, "| markers by cluster")),
     paste0(tag, "_dotplot_by_cluster.png"), length(unique(o[[col_cluster]][, 1])))
  sv(DotPlot(o, features = feat, group.by = col_annot) + thm + ggtitle(paste(tag, "| markers by annotation")),
     paste0(tag, "_dotplot_by_annotation.png"), length(unique(o[[col_annot]][, 1])))
  o$annot_ds <- paste(o[[col_annot]][, 1], o$dataset, sep = " | ")
  sv(DotPlot(o, features = feat, group.by = "annot_ds") + thm + ggtitle(paste(tag, "| by annotation and dataset")),
     paste0(tag, "_dotplot_by_annotation_dataset.png"), length(unique(o$annot_ds)))
}
# coarse region of a sample, used for plausibility checks (DA neurons only exist in midbrain samples)
region_group <- function(r) ifelse(grepl("nigra|midbrain|SNpc|^SN$|mesencephal", r, ignore.case = TRUE), "Midbrain",
                          ifelse(grepl("cortex|cingulate|frontal|BA9", r, ignore.case = TRUE), "Cortex",
                          ifelse(grepl("putamen|caudate|accumbens|NBM|striat|basal", r, ignore.case = TRUE), "Striatum/basal forebrain", "Other")))
panels_for <- function(step) { p <- cfg_panels(); p <- p[p$step == step, ]; split(p$gene, p$panel) }
panel_to <- function(step, level) { p <- unique(cfg_panels()[cfg_panels()$step == step, c("panel", level)]); setNames(p[[level]], p$panel) }
