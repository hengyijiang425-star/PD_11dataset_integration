# 06_annotate_L1.R — level-1 (lineage) annotation of every cluster with MARKER SETS, not single genes
# For each cluster: module score for every panel in config/marker_panels_v2.csv (step = "L1"),
# best panel is called only if >= 2 of its genes are detected in >= 30% of the cluster's cells.
# Also reported per cluster: top-10 DE markers, markers conserved in >= 2 datasets, author agreement,
# dataset mix and top-donor share (to spot study-specific or one-donor clusters).
source("R/utils.R")
an <- readRDS(file.path(DIR$res, "05_analysis.rds"))
pnl <- panels_for("L1"); to_L1 <- panel_to("L1", "unified_L1"); to_L2 <- panel_to("L1", "unified_L2")

a1 <- annotate_clusters(an, "cluster_L1", pnl, "06_L1", author_col = "author_L1", call_map = to_L1)
tab <- a1$table
tab$celltype_L1 <- ifelse(tab$call == "Unresolved", "Unresolved", to_L1[tab$call])
tab$celltype_L2_from_L1 <- ifelse(tab$call == "Unresolved", "Unresolved", to_L2[tab$call])
write.csv(tab, file.path(DIR$res, "06_L1_cluster_annotation.csv"), row.names = FALSE)

i <- match(as.character(an$cluster_L1), tab$cluster)
an$L1_panel_call <- tab$call[i]; an$celltype_L1 <- tab$celltype_L1[i]; an$confidence_L1 <- tab$confidence[i]
an$celltype_L2 <- tab$celltype_L2_from_L1[i]          # refined for Neuron / Myeloid / Vascular in step 07
an$celltype_L3 <- NA_character_

dot_two_ways(an, "cluster_L1", "celltype_L1", a1$panels, "06_L1")
for (g in c("celltype_L1", "L1_panel_call", "confidence_L1"))
  ggsave(file.path(DIR$plot, paste0("06_umap_", g, ".png")), DimPlot(an, reduction = "umap", group.by = g, label = g != "confidence_L1", repel = TRUE, raster = TRUE),
         width = 9, height = 7, dpi = 200, bg = "white")
if (any(!is.na(an$author_L1))) {
  x <- table(author = an$author_L1, ours = an$celltype_L1)
  write.csv(as.data.frame.matrix(x), file.path(DIR$res, "06_author_vs_ours_L1.csv"))
  cat("Agreement with authors (cells that have an author label):",
      round(100 * mean(an$author_L1 == an$celltype_L1, na.rm = TRUE), 1), "%\n")
}
saveRDS(an, file.path(DIR$res, "06_analysis_L1.rds"))
print(tab[, c("cluster", "n_cells", "call", "confidence", "n_datasets_ge5pct", "top10_DE_markers")], row.names = FALSE)
cat("\nCheck results/06_L1_cluster_annotation.csv. Low-confidence or one-dataset/one-donor clusters need a look before step 07.\n",
    "To override a call by hand, edit celltype_L1 in that CSV and re-run 07 with USE_MANUAL <- TRUE.\n",
    "Next: Rscript R/07_lineage_subclustering.R\n")
