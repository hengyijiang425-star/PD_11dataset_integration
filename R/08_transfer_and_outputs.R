# 08_transfer_and_outputs.R — give EVERY cell its labels, then write the final tables and plots
# Sketch mode: each non-sketch cell takes the majority label of its k = 15 nearest sketch cells in the
# projected Harmony space ("harmony.full"); transfer_score = share of neighbours that agree.
# Outputs (results/):
#   PD11_final_cell_annotations.csv.gz   one row per cell: dataset, sample, donor, diagnosis, region, L1/L2/L3, UMAP
#   PD11_marker_sets_per_cluster.csv     per cluster: call, confidence, panel markers, top-10 DE, conserved markers
#   PD11_donor_composition_L2.csv        donor x cell type proportions (enriched/sorted samples excluded)
#   PD11_cluster_dataset_mixing.csv      dataset entropy per cluster (1 = evenly mixed)
#   PD11_full_annotated.rds              Seurat object, all cells (counts on disk in data/counts)
source("R/utils.R")
an <- readRDS(file.path(DIR$res, "07_analysis_L2.rds"))
full <- if (an@misc$mode == "sketch") readRDS(file.path(DIR$res, "05_full.rds")) else an
labs <- c("celltype_L1", "celltype_L2", "celltype_L3", "cluster_fine")
key <- do.call(paste, c(an[[]][, labs], sep = "||"))

if (an@misc$mode == "sketch") {
  ref <- Embeddings(full, "harmony.full")[colnames(an), 1:SET$n_pcs]
  qry_cells <- setdiff(colnames(full), colnames(an))
  msg("transferring labels to", length(qry_cells), "cells")
  lab <- setNames(key, colnames(an)); score <- setNames(rep(1, ncol(an)), colnames(an))
  use_bn <- requireNamespace("BiocNeighbors", quietly = TRUE)
  for (ch in split(qry_cells, ceiling(seq_along(qry_cells) / 2e5))) {
    q <- Embeddings(full, "harmony.full")[ch, 1:SET$n_pcs, drop = FALSE]
    nn <- if (use_bn) BiocNeighbors::queryKNN(ref, q, k = SET$knn_k, BNPARAM = BiocNeighbors::AnnoyParam())$index
          else RANN::nn2(ref, q, k = SET$knn_k)$nn.idx
    v <- matrix(key[nn], nrow = nrow(nn))
    best <- apply(v, 1, function(x) { t <- table(x); names(t)[which.max(t)] })
    lab[ch] <- best; score[ch] <- rowMeans(v == best); msg("  ", length(ch), "done")
  }
  lab <- lab[colnames(full)]; score <- score[colnames(full)]
  umap <- Embeddings(full, "full.umap")
} else {
  lab <- setNames(key, colnames(an)); score <- setNames(rep(1, ncol(an)), colnames(an)); umap <- Embeddings(an, "umap")
}
sp <- do.call(rbind, strsplit(lab, "||", fixed = TRUE))
for (j in seq_along(labs)) full[[labs[j]]] <- ifelse(sp[, j] == "NA", NA, sp[, j])
full$transfer_score <- score; full$in_sketch <- colnames(full) %in% colnames(an)
# short display names in the PI's style (config/display_names.csv): EN, IN, Astro, Oligo, Micro_PVM, Adaptive, Mural, Endo ...
dn <- read.csv(file.path(DIR$cfg, "display_names.csv"))
ct <- dn$cell_type[match(full$celltype_L2, dn$unified_L2)]
ct[!is.na(full$celltype_L3) & full$celltype_L3 == "MSN"] <- "MSN"
ct[is.na(ct)] <- "Unresolved"                       # "(unresolved)", "(doublet/ambient?)"
ct_levels <- unique(dn$cell_type); full$cell_type <- factor(ct, levels = ct_levels[ct_levels %in% ct])
ct_cols <- setNames(dn$color[!duplicated(dn$cell_type)], dn$cell_type[!duplicated(dn$cell_type)])
saveRDS(full, file.path(DIR$res, "PD11_full_annotated.rds"))

md <- full[[]]
fin <- data.frame(cell = colnames(full), md[, intersect(c("dataset", "sample_key", "sample_uid", "donor", "diagnosis", "age", "sex", "PMI_hours", "region",
                  "enriched", "nCount", "nFeature", "percent.mt", "doublet_score", "author_label", "author_L1", "author_L2",
                  labs, "cell_type", "transfer_score", "in_sketch"), colnames(md))], UMAP_1 = umap[colnames(full), 1], UMAP_2 = umap[colnames(full), 2])
fwrite(fin, file.path(DIR$res, "PD11_final_cell_annotations.csv.gz"))

# ---------- marker sets per cluster (L1 + every lineage) ----------
mk <- list(read.csv(file.path(DIR$res, "06_L1_cluster_annotation.csv")))
mk[[1]]$level <- "L1 (all cells)"
for (f in list.files(DIR$res, pattern = "^07_.*_cluster_annotation\\.csv$", full.names = TRUE)) {
  x <- read.csv(f); x$level <- sub("^07_(.*)_cluster_annotation\\.csv$", "L2 within \\1", basename(f)); mk[[length(mk) + 1]] <- x }
mk <- data.table::rbindlist(mk, fill = TRUE)
keep <- c("level", "cluster", "n_cells", "call", "celltype_L1", "celltype_L2", "celltype_L3", "confidence", "panel_markers",
          "n_panel_markers_supporting", "top10_DE_markers", "conserved_markers_2plus_datasets", "author_majority", "agrees_with_author",
          "n_datasets_ge5pct", "datasets", "top_donor_fraction")
fwrite(mk[, intersect(keep, names(mk)), with = FALSE], file.path(DIR$res, "PD11_marker_sets_per_cluster.csv"))

# ---------- mixing, composition, plots ----------
ent <- function(x) { p <- prop.table(table(x)); p <- p[p > 0]; if (length(unique(md$dataset)) < 2) NA else -sum(p * log(p)) / log(length(unique(md$dataset))) }
mix <- aggregate(dataset ~ cluster_fine, md, ent); names(mix)[2] <- "dataset_entropy"
mix$n_cells <- as.vector(table(md$cluster_fine)[mix$cluster_fine])
mix$n_datasets <- as.vector(tapply(md$dataset, md$cluster_fine, function(x) length(unique(x)))[mix$cluster_fine])
fwrite(mix[order(mix$dataset_entropy), ], file.path(DIR$res, "PD11_cluster_dataset_mixing.csv"))

cm <- md[!(md$enriched %in% TRUE) & !is.na(md$celltype_L2), ]
comp <- as.data.frame(prop.table(table(cm$donor, cm$celltype_L2), 1)); names(comp) <- c("donor", "celltype_L2", "proportion")
dm <- unique(cm[, intersect(c("donor", "dataset", "region", "diagnosis", "age", "sex", "PMI_hours"), names(cm))])
dm <- dm[!duplicated(dm$donor), ]
fwrite(merge(comp, dm, by = "donor"), file.path(DIR$res, "PD11_donor_composition_L2.csv"))

pl <- function(g, f, lab = FALSE) {
  p <- ggplot(fin[sample(nrow(fin)), ], aes(UMAP_1, UMAP_2, colour = .data[[g]])) + geom_point(size = 0.05, alpha = 0.4) +
    theme_void() + guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) + ggtitle(g)
  if (g == "cell_type") p <- p + scale_colour_manual(values = ct_cols)
  if (lab) { cen <- aggregate(cbind(UMAP_1, UMAP_2) ~ get(g), fin, median); names(cen)[1] <- g
    p <- p + geom_text(data = cen, aes(label = .data[[g]]), colour = "black", size = 3) }
  ggsave(file.path(DIR$plot, f), p, width = 10, height = 7.5, dpi = 200, bg = "white")
}
pl("dataset", "08_umap_all_cells_dataset.png"); pl("region", "08_umap_all_cells_region.png")
pl("diagnosis", "08_umap_all_cells_diagnosis.png"); pl("celltype_L1", "08_umap_all_cells_L1.png", TRUE)
pl("celltype_L2", "08_umap_all_cells_L2.png", TRUE); pl("cell_type", "08_umap_all_cells_cell_type.png", TRUE)
p <- ggplot(fin, aes(x = dataset, fill = cell_type)) + geom_bar(position = "fill", colour = "white", linewidth = 0.2) + theme_bw() +
  scale_fill_manual(values = ct_cols) + ylab("fraction of nuclei") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(DIR$plot, "08_cell_type_fraction_by_dataset.png"), p, width = 10, height = 6, dpi = 200, bg = "white")
for (smp in head(names(sort(table(fin$sample_uid), decreasing = TRUE)), 3)) {     # one-sample UMAPs, like the PI's example
  x <- fin[fin$sample_uid == smp, ]
  p <- ggplot(fin, aes(UMAP_1, UMAP_2)) + geom_point(colour = "#ececec", size = 0.05) +
    geom_point(data = x, aes(colour = cell_type), size = 0.4) + scale_colour_manual(values = ct_cols) + theme_void() +
    guides(colour = guide_legend(override.aes = list(size = 3))) + ggtitle(smp)
  ggsave(file.path(DIR$plot, paste0("08_umap_sample_", make.names(smp), ".png")), p, width = 9, height = 7, dpi = 200, bg = "white")
}
p <- ggplot(fin, aes(x = dataset, fill = celltype_L2)) + geom_bar(position = "fill") + theme_bw() + ylab("fraction of nuclei") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(DIR$plot, "08_celltype_L2_fraction_by_dataset.png"), p, width = 10, height = 6, dpi = 200, bg = "white")

cat("\nCells per L2 type:\n"); print(sort(table(fin$celltype_L2), decreasing = TRUE))
cat("\nLabel transfer: median score", round(median(score), 2), "; cells with score < 0.5:", sum(score < 0.5), "\n")
cat("\nDONE. Key files: results/PD11_marker_sets_per_cluster.csv, results/PD11_final_cell_annotations.csv.gz, plots/08_*.png\n")
