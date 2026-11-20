# 07_lineage_subclustering.R — re-run Harmony INSIDE each lineage, then level-2 annotation
# Global Harmony is not strong enough for fine types, so Neuron, Microglia/Myeloid and Vascular cells
# are re-normalized, re-integrated (Harmony on dataset + sample) and re-clustered on their own.
# Panels: config/marker_panels_v2.csv with step = lineage name.
# Neuron calls: DA / CADPS2-RIT2 low-TH / Excitatory / Inhibitory interneuron / MSN / Cholinergic
#   (MSN is stored as L2 = Inhibitory (GABAergic) neuron, L3 = MSN; Smooth muscle -> L2 Pericyte, L3 SMC)
# Other lineages keep L2 = L1 (Astrocyte, Oligodendrocyte, OPC, T cell, Ependymal).
source("R/utils.R"); suppressPackageStartupMessages(library(harmony))
USE_MANUAL <- FALSE      # TRUE = take celltype_L1 from an edited results/06_L1_cluster_annotation.csv
an <- readRDS(file.path(DIR$res, "06_analysis_L1.rds"))
if (USE_MANUAL) { t1 <- read.csv(file.path(DIR$res, "06_L1_cluster_annotation.csv"))
  an$celltype_L1 <- t1$celltype_L1[match(as.character(an$cluster_L1), as.character(t1$cluster))] }
an$cluster_fine <- paste0(an$celltype_L1, "_c", an$cluster_L1)
L3_of <- c("MSN" = "MSN", "Smooth muscle" = "SMC", "Macrophage/BAM" = "Macrophage/BAM", "Cholinergic neuron" = NA)

subcluster <- function(lin) {
  cells <- colnames(an)[an$celltype_L1 == lin]
  if (length(cells) < 200) { msg(lin, ": only", length(cells), "cells - skipped"); return(NULL) }
  msg(lin, ":", length(cells), "cells")
  s <- CreateSeuratObject(counts = LayerData(an, layer = "counts")[, cells], meta.data = an[[]][cells, ])
  n_ds <- table(s$dataset); s$split_ds <- ifelse(s$dataset %in% names(n_ds)[n_ds >= 50], s$dataset, "small_datasets")
  s[["RNA"]] <- split(s[["RNA"]], f = s$split_ds)                # HVGs per dataset, then combined
  s <- NormalizeData(s, verbose = FALSE); s <- FindVariableFeatures(s, nfeatures = 2000, verbose = FALSE)
  s <- ScaleData(s, verbose = FALSE); s <- RunPCA(s, npcs = SET$n_pcs, verbose = FALSE); s <- JoinLayers(s)
  s <- RunHarmony(s, group.by.vars = c("dataset", "sample_uid"), reduction.use = "pca", dims.use = 1:SET$n_pcs, verbose = FALSE)
  s <- RunUMAP(s, reduction = "harmony", dims = 1:SET$n_pcs, verbose = FALSE)
  s <- FindNeighbors(s, reduction = "harmony", dims = 1:SET$n_pcs, verbose = FALSE)
  s <- FindClusters(s, resolution = SET$res_sub, verbose = FALSE); s$cluster_sub <- as.character(s$seurat_clusters)
  tag <- paste0("07_", gsub("[^A-Za-z]", "", lin))
  pnl <- panels_for(lin); to2 <- panel_to(lin, "unified_L2"); to1 <- panel_to(lin, "unified_L1")
  # contamination check: other lineages' level-1 panels. A cluster that scores best for one of them is a likely
  # doublet / ambient-RNA cluster, not a new cell type.
  l1 <- cfg_panels(); l1 <- l1[l1$step == "L1" & l1$unified_L1 != lin, ]
  for (q in setdiff(unique(l1$panel), names(pnl))) { nm <- paste0("Contam: ", q); pnl[[nm]] <- l1$gene[l1$panel == q]
    to2[nm] <- paste(lin, "(doublet/ambient?)"); to1[nm] <- lin }
  a <- annotate_clusters(s, "cluster_sub", pnl, tag, author_col = "author_L2", call_map = to2)
  tab <- a$table; tab$celltype_L2 <- ifelse(tab$call == "Unresolved", paste(lin, "(unresolved)"), to2[tab$call])
  tab$celltype_L1 <- ifelse(tab$call == "Unresolved", lin, to1[tab$call])
  # plausibility: DA-lineage neurons only exist in midbrain samples
  imp <- tab$celltype_L2 %in% c("DA neuron", "CADPS2/RIT2 low-TH neuron") & !is.na(tab$frac_midbrain) & tab$frac_midbrain < 0.5
  tab$flag <- ifelse(imp, "DA-lineage call but < 50% midbrain cells", ifelse(grepl("^Contam", tab$call), "likely doublet/ambient", ""))
  # cortical CADPS2/RIT2-high clusters are deep-layer excitatory neurons (THEMIS, FEZF2): use the next-best panel instead
  dal <- c("DA neuron", "CADPS2/RIT2 low-TH neuron")
  alt <- ifelse(!is.na(tab$second) & tab$second %in% names(to2) & !grepl("^Contam", tab$second) & !(to2[tab$second] %in% dal), to2[tab$second], NA)
  tab$celltype_L2[imp] <- ifelse(is.na(alt[imp]), paste(lin, "(unresolved)"), alt[imp])
  tab$flag[imp & !is.na(alt)] <- paste0("DA-lineage call in non-midbrain cluster -> next-best panel: ", tab$second[imp & !is.na(alt)])
  tab$celltype_L3 <- unname(L3_of[tab$call]); write.csv(tab, file.path(DIR$res, paste0(tag, "_cluster_annotation.csv")), row.names = FALSE)
  i <- match(s$cluster_sub, tab$cluster)
  s$celltype_L1 <- tab$celltype_L1[i]; s$celltype_L2 <- tab$celltype_L2[i]
  bad <- s$celltype_L2 %in% c("DA neuron", "CADPS2/RIT2 low-TH neuron") & s$region_group != "Midbrain"
  if (any(bad)) { msg("  ", sum(bad), "cortex/striatum cells in DA-lineage clusters -> unresolved"); s$celltype_L2[bad] <- paste(lin, "(unresolved)") }
  s$celltype_L3 <- tab$celltype_L3[i]; s$panel_call_L2 <- tab$call[i]; s$confidence_L2 <- tab$confidence[i]
  dot_two_ways(s, "cluster_sub", "panel_call_L2", a$panels, tag)
  for (g in c("dataset", "region", "cluster_sub", "panel_call_L2", "diagnosis"))
    ggsave(file.path(DIR$plot, paste0(tag, "_umap_", g, ".png")), DimPlot(s, group.by = g, label = g %in% c("cluster_sub", "panel_call_L2"), repel = TRUE, raster = TRUE),
           width = 8, height = 6, dpi = 200, bg = "white")
  if (any(!is.na(s$author_L2)))
    write.csv(as.data.frame.matrix(table(author = s$author_L2, ours = s$celltype_L2)), file.path(DIR$res, paste0(tag, "_author_vs_ours_L2.csv")))
  saveRDS(s, file.path(DIR$res, paste0(tag, ".rds")))
  data.frame(cell = colnames(s), L1 = s$celltype_L1, L2 = s$celltype_L2, L3 = s$celltype_L3, fine = paste0(lin, "_s", s$cluster_sub), conf = s$confidence_L2)
}
res <- do.call(rbind, lapply(c("Neuron", "Microglia/Myeloid", "Vascular", "Astrocyte"), subcluster))
i <- match(res$cell, colnames(an))
an$celltype_L1[i] <- res$L1; an$celltype_L2[i] <- res$L2; an$celltype_L3[i] <- res$L3; an$cluster_fine[i] <- res$fine; an$confidence_L2 <- an$confidence_L1
an$confidence_L2[i] <- res$conf
saveRDS(an, file.path(DIR$res, "07_analysis_L2.rds"))
print(table(an$celltype_L1, an$celltype_L2))
cat("\nNext: Rscript R/08_transfer_and_outputs.R\n")
