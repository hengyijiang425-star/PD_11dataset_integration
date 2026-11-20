# 05_merge_integrate.R — shared genes, one Seurat v5 object (one layer per sample), Harmony
# integration_mode = "sketch" (default, for ~1.5M nuclei):
#   normalize every sample -> HVGs -> LeverageScore sketch (SET$sketch_cells_per_sample per sample)
#   -> PCA + Harmony (batch = sample, so dataset is corrected too) on the sketch -> clusters + UMAP model
#   -> ProjectIntegration / ProjectData give every cell a Harmony embedding ("harmony.full") and UMAP ("full.umap")
# integration_mode = "full": PCA + Harmony(dataset + sample) on all cells (use with max_cells_per_sample)
# Output: results/05_full.rds      all cells, counts stay on disk (do not move data/counts/)
#         results/05_analysis.rds  in-memory object used for annotation (the sketch, or all cells in "full" mode)
source("R/utils.R"); suppressPackageStartupMessages(library(harmony))
man <- read.csv(file.path(DIR$proc, "manifest.csv")); qc <- readRDS(file.path(DIR$proc, "cells_qc.rds"))
sm <- cfg_samples(); map <- cfg_map()
if (is.finite(SET$integrate_cells_per_sample)) {           # laptop test: cap per DONOR (Kamath has 97 libraries from 22 donors)
  dn <- sm$donor[match(paste(qc$dataset, qc$sample_key), paste(sm$dataset, sm$sample_key))]
  grp <- paste(qc$dataset, ifelse(is.na(dn), qc$sample_key, dn))
  qc <- do.call(rbind, lapply(split(qc, grp), function(x)
    if (nrow(x) > SET$integrate_cells_per_sample) x[sample(nrow(x), SET$integrate_cells_per_sample), ] else x))
  msg("using", nrow(qc), "cells (", SET$integrate_cells_per_sample, "per donor )")
}
man <- man[paste(man$dataset, man$sample_key) %in% unique(paste(qc$dataset, qc$sample_key)), ]

# ---------- shared genes ----------
genes <- lapply(man$path, function(p) rownames(load_counts(p)))
by_ds <- lapply(split(genes, man$dataset), function(l) Reduce(intersect, l))
shared <- Reduce(intersect, by_ds)
cat("Genes per dataset:\n"); print(sapply(by_ds, length)); cat("Shared genes:", length(shared), "\n")
if (length(shared) < SET$min_shared_genes) stop("Too few shared genes - check gene IDs (see Genes per dataset above)")
write.csv(data.frame(gene = shared), file.path(DIR$res, "05_shared_genes.csv"), row.names = FALSE)

# ---------- one layer per sample ----------
man$sample_uid <- paste(man$dataset, make.names(man$sample_key), sep = "_")
qc$sample_uid <- paste(qc$dataset, make.names(qc$sample_key), sep = "_")
qc$cell <- paste(qc$sample_uid, qc$barcode, sep = "_")
mats <- setNames(lapply(seq_len(nrow(man)), function(j) {
  m <- load_counts(man$path[j]); b <- qc$barcode[qc$sample_uid == man$sample_uid[j]]
  m <- m[shared, b]; colnames(m) <- paste(man$sample_uid[j], b, sep = "_"); m }), man$sample_uid)
msg("combining", length(mats), "samples into one matrix")
cnt_all <- do.call(cbind, mats); rm(mats); gc()                 # one matrix: much faster than ~300 Seurat layers
obj <- CreateSeuratObject(counts = cnt_all); rm(cnt_all); gc()

# ---------- cell metadata ----------
md <- qc[match(colnames(obj), qc$cell), ]
smk <- match(paste(md$dataset, md$sample_key), paste(sm$dataset, sm$sample_key))
for (v in c("diagnosis", "age", "sex", "PMI_hours", "region", "donor", "enriched")) md[[v]] <- sm[[v]][smk]
md$donor <- ifelse(is.na(md$donor), md$sample_uid, paste(md$dataset, md$donor, sep = "_"))
md$region_group <- region_group(md$region)
key <- paste(md$dataset, tolower(trimws(md$author_label))); mk <- paste(map$dataset, tolower(trimws(map$author_label)))
md$author_L1 <- map$unified_L1[match(key, mk)]; md$author_L2 <- map$unified_L2[match(key, mk)]
un <- unique(md[!is.na(md$author_label) & is.na(md$author_L1), c("dataset", "author_label")])
if (nrow(un)) { cat("\nAuthor labels not in the crosswalk (add them to config/celltype_harmonization_map_11datasets.csv):\n"); print(un) }
rownames(md) <- md$cell
obj <- AddMetaData(obj, md[, setdiff(colnames(md), c("cell", "barcode"))])
cat("\nCells per dataset:\n"); print(table(obj$dataset))

npc <- 1:SET$n_pcs
if (SET$integration_mode == "sketch") {
  obj[["RNA"]] <- split(obj[["RNA"]], f = obj$sample_uid)      # sketch + Harmony work per sample layer
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, nfeatures = SET$n_hvg, verbose = FALSE)
  obj <- SketchData(obj, ncells = SET$sketch_cells_per_sample, method = "LeverageScore", sketched.assay = "sketch", seed = 1)
  DefaultAssay(obj) <- "sketch"; msg("sketch cells:", length(Cells(obj[["sketch"]])))
  obj <- FindVariableFeatures(obj, nfeatures = SET$n_hvg, verbose = FALSE)
  obj <- ScaleData(obj, verbose = FALSE); obj <- RunPCA(obj, npcs = SET$n_pcs, verbose = FALSE)
  obj <- RunUMAP(obj, reduction = "pca", dims = npc, reduction.name = "umap.pre", verbose = FALSE)
  obj <- IntegrateLayers(obj, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony", verbose = FALSE)
  obj <- FindNeighbors(obj, reduction = "harmony", dims = npc, verbose = FALSE)
  for (r in SET$res_all) obj <- FindClusters(obj, resolution = r, verbose = FALSE)
  obj$cluster_L1 <- obj[[paste0("sketch_snn_res.", SET$res_main)]][, 1]
  obj <- RunUMAP(obj, reduction = "harmony", dims = npc, return.model = TRUE, reduction.name = "umap", verbose = FALSE)
  msg("projecting all cells")
  obj <- ProjectIntegration(obj, sketched.assay = "sketch", assay = "RNA", reduction = "harmony")
  obj <- ProjectData(obj, sketched.assay = "sketch", assay = "RNA", sketched.reduction = "harmony.full",
                     full.reduction = "harmony.full", dims = npc, refdata = list(cluster_L1_full = "cluster_L1"), umap.model = "umap")
  obj[["sketch"]] <- JoinLayers(obj[["sketch"]])
  cells <- Cells(obj[["sketch"]]); cnt <- LayerData(obj, assay = "sketch", layer = "counts")
  saveRDS(obj, file.path(DIR$res, "05_full.rds"))
} else {
  # one layer per DATASET (11) instead of per sample (~300): same HVG logic as the 3-dataset project, much faster
  msg("splitting by dataset"); obj[["RNA"]] <- split(obj[["RNA"]], f = obj$dataset)
  msg("normalize + HVGs per dataset"); obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, nfeatures = SET$n_hvg, verbose = FALSE)
  msg("scale + PCA"); obj <- ScaleData(obj, verbose = FALSE); obj <- RunPCA(obj, npcs = SET$n_pcs, verbose = FALSE)
  obj <- JoinLayers(obj)
  msg("UMAP before integration"); obj <- RunUMAP(obj, reduction = "pca", dims = npc, reduction.name = "umap.pre", verbose = FALSE)
  msg("Harmony (dataset + sample)")
  obj <- RunHarmony(obj, group.by.vars = c("dataset", "sample_uid"), reduction.use = "pca", dims.use = npc, reduction.save = "harmony", verbose = FALSE)
  msg("clustering"); obj <- FindNeighbors(obj, reduction = "harmony", dims = npc, verbose = FALSE)
  for (r in SET$res_all) obj <- FindClusters(obj, resolution = r, verbose = FALSE)
  obj$cluster_L1 <- obj[[paste0("RNA_snn_res.", SET$res_main)]][, 1]
  msg("UMAP after integration"); obj <- RunUMAP(obj, reduction = "harmony", dims = npc, reduction.name = "umap", verbose = FALSE)
  obj[["RNA"]]$scale.data <- NULL; an <- obj; rm(obj); gc()   # full mode: one object only (saves RAM)
}

# ---------- slim in-memory analysis object (sketch mode) ----------
if (SET$integration_mode == "sketch") {
  an <- CreateSeuratObject(counts = in_memory(cnt[, cells]), meta.data = obj[[]][cells, ])
  an <- NormalizeData(an, verbose = FALSE)
  for (r in c("pca", "harmony", "umap.pre", "umap")) {
    e <- Embeddings(obj, r)[cells, , drop = FALSE]
    an[[r]] <- CreateDimReducObject(e, key = paste0(gsub("[^A-Za-z]", "", r), "_"), assay = "RNA")
  }
  rm(obj, cnt); gc()
}
an@misc$mode <- SET$integration_mode
saveRDS(an, file.path(DIR$res, "05_analysis.rds"))

for (g in c("dataset", "region", "diagnosis", "cluster_L1"))
  ggsave(file.path(DIR$plot, paste0("05_umap_", g, ".png")), DimPlot(an, reduction = "umap", group.by = g, raster = TRUE, label = g == "cluster_L1"),
         width = 9, height = 7, dpi = 200, bg = "white")
ggsave(file.path(DIR$plot, "05_umap_BEFORE_harmony_dataset.png"), DimPlot(an, reduction = "umap.pre", group.by = "dataset", raster = TRUE), width = 9, height = 7, dpi = 200, bg = "white")
write.csv(as.data.frame.matrix(round(prop.table(table(an$cluster_L1, an$dataset), 1), 3)), file.path(DIR$res, "05_cluster_dataset_fraction.csv"))
cat("\nNext: Rscript R/06_annotate_L1.R\n")
