# 00_install_packages.R — install only what is missing, binaries only (no compile questions)
# Run this file ON ITS OWN, then run the next script. Do not paste several commands at once:
# an install question ("Do you want to install from sources?") would swallow the next line.
options(install.packages.compile.from.source = "never")
Sys.unsetenv(c("GITHUB_PAT", "GITHUB_TOKEN"))
cran <- c("Seurat", "harmony", "Matrix", "data.table", "ggplot2", "RANN", "remotes", "BiocManager", "patchwork", "R.utils", "hexbin")
bioc <- c("GEOquery", "scDblFinder", "SingleCellExperiment", "org.Hs.eg.db", "AnnotationDbi", "rhdf5", "zellkonverter")
ok <- function(p) suppressWarnings(suppressMessages(requireNamespace(p, quietly = TRUE)))
miss <- cran[!vapply(cran, ok, logical(1))]
if (length(miss)) install.packages(miss, type = "binary")
miss <- bioc[!vapply(bioc, ok, logical(1))]
if (length(miss)) BiocManager::install(miss, update = FALSE, ask = FALSE, type = "binary", force = TRUE)

# Optional, only if a ready-made binary exists for your R version (no compiling):
#   BPCells = counts on disk (needed for the full 1.5M-nuclei run, not for the laptop test)
#   presto  = fast marker tests (FindAllMarkers is used without it)
for (p in c("BPCells", "presto")) if (!ok(p))
  tryCatch(install.packages(p, repos = c(sprintf("https://%s.r-universe.dev", if (p == "BPCells") "bnprks" else "immunogenomics"),
                                         "https://cloud.r-project.org"), type = "binary"),
           error = function(e) NULL, warning = function(w) NULL)

for (p in c(cran, bioc, "presto", "BPCells")) {
  v <- if (ok(p)) as.character(packageVersion(p)) else "MISSING"
  cat(sprintf("%-22s %s\n", p, v))
}
if (!ok("scDblFinder")) {
  cat("\nscDblFinder is installed but does not load. The reason:\n")
  print(tryCatch({ library(scDblFinder); "loads now" }, error = function(e) conditionMessage(e)))
}
cat("\nR", R.version$major, ".", R.version$minor, "| platform:", R.version$platform, "\n")
cat("BPCells / presto MISSING is OK for the laptop test (SET$max_cells_per_sample = 300, integration_mode = 'full').\n")
