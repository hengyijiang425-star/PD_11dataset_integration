# PD 11-dataset snRNA-seq atlas

Re-processes 11 public Parkinson's disease single-nucleus RNA-seq datasets (midbrain/SN, cortex, striatum/basal forebrain) with **one** pipeline: the same QC, one Harmony integration, and one cell-type naming scheme. Every cluster is labelled from a **marker set**, never from a single gene. Two downstream analyses then compare PD and control donors: cell-type composition and pseudobulk differential expression.

> **Status (19 Nov 2026): test run.** Every step runs end to end on a laptop, but on a subsample: at most 300 nuclei loaded per sample, and 200 per donor from integration onward. All numbers below come from this test run. The full run (~1.5 M nuclei) is planned for the Cayuga cluster.
>
> **No data are in this repository**: no raw or processed counts, Seurat objects or result tables. `data/`, `results/` and `plots/` are git-ignored. The notebooks hold the outputs as they were when the notebooks ran.

## Datasets

| GEO | Region | Paper |
|---|---|---|
| GSE157783 | Midbrain | Smajić 2022 Brain |
| GSE184950 | SN | Wang 2024 |
| GSE243639 | SNpc | Martirosyan 2024 Mol Neurodegener |
| GSE148434 | SN | Lee 2023 Sci Adv |
| GSE253462 | SN + cingulate cortex | Jakubiak 2025 Nat Commun |
| GSE235330 | NBM + NAc + putamen | Pech 2025 eLife |
| GSE178146 | Anterior cingulate cortex | Feleke 2021 Acta Neuropathol |
| GSE161045 | Putamen | Xu 2023 Nat Aging |
| GSE202210 | Prefrontal cortex (BA9) | Zhu 2024 Sci Transl Med |
| GSE193688 | Midbrain/SN (Multiome, RNA part) | Adams 2024 Nat Aging |
| GSE178265 | SNpc (+ caudate) | Kamath 2022 Nat Neurosci |

How to read each dataset is set in `config/datasets.csv`. GSE157783, GSE184950 and GSE243639 are read from the earlier `PD_3dataset_integration` project, which must sit next to this folder.

## Repository layout

```
R/            the pipeline, steps 00–08 (Seurat v5); settings in R/utils.R → SET
analysis/     Python steps after the pipeline (pseudobulk, DE, composition test, figures)
config/       dataset list, marker panels, display names/colours, author-label crosswalk
notebooks/    executed notebooks with outputs: QC, atlas overview, composition, DE
```

## Run order

### 1. Pipeline in R (open `PD_11dataset_integration.Rproj` in RStudio)

```r
source("R/00_install_packages.R")        # installs missing packages (binary only)
source("R/01_download_GEO.R")            # GEO supplementary files -> data/raw/  (Kamath alone is ~5 GB)
source("R/02_inventory_and_metadata.R")  # detects file formats, drafts config/sample_metadata.csv
# >>> check config/sample_metadata.csv by hand: diagnosis, donor, region, age, sex, PMI, enriched, include
source("R/run_all.R")                    # steps 03–08; from the shell: Rscript R/run_all.R 05 to resume at 05
```

| Step | What it does | Main output |
|---|---|---|
| 03 `load_counts` | Reads every sample (h5 / mtx / h5ad / rds / text), converts genes to symbols, saves counts per sample | `data/counts/`, manifest |
| 04 `qc_doublets` | One QC rule for every sample (below) + scDblFinder | `results/04_QC_per_sample.csv` |
| 05 `merge_integrate` | Caps nuclei per **donor**, keeps the genes shared by all datasets (16,345), Harmony on dataset + sample | `results/05_*.rds` |
| 06 `annotate_L1` | Major lineages, each called from a marker panel | `results/06_L1_cluster_annotation.csv` |
| 07 `lineage_subclustering` | Neuron, Microglia/Myeloid, Vascular, Astrocyte: re-integration, sub-clustering, L2 labels; contamination panels; region-plausibility check | `results/07_*` |
| 08 `transfer_and_outputs` | kNN label transfer to every nucleus, `cell_type` from `config/display_names.csv`, final tables and plots | `results/PD11_*`, `plots/08_*` |

**QC (step 04):** genes ≥ 500, UMIs ≥ 1,000, mito < 10 % **and** ≤ sample median + 3 MAD (floor 1 %), genes ≤ sample median + 3 MAD, doublets removed; samples with fewer than 100 nuclei left are dropped.

**Annotation rule (06/07):** a cluster gets a label only if ≥ 2 markers from that panel are detected in ≥ 30 % of its nuclei. Contamination panels flag doublet or ambient-RNA clusters, which are labelled "Unresolved". DA-lineage labels are allowed only in midbrain samples; elsewhere the next-best panel is used. Panels and their literature sources are in `config/marker_panels_v2.csv`.

**Test vs full run** (`SET` in `R/utils.R`):

| Setting | Test (laptop) | Full (cluster) |
|---|---|---|
| `max_cells_per_sample` | 300 | `Inf` |
| `integrate_cells_per_sample` (per donor) | 200 | `Inf` |
| `integration_mode` | `"full"` | `"sketch"` (needs BPCells) |

### 2. Downstream analyses in Python (run from the project root)

```bash
python analysis/01_pseudobulk.py         # counts summed per donor x cell type (sorted libraries + Unresolved excluded)
python analysis/02_pseudobulk_DE.py      # pydeseq2, ~ dataset + condition, per region x cell type
python analysis/03_DE_figures.py         # DE count bars, volcanoes, shared-gene heatmap
python analysis/04_composition_test.py   # PD vs control proportions, empirical logit ~ PD + dataset, BH within region
```

Python packages: pandas, numpy, scipy, statsmodels, matplotlib, adjustText, pydeseq2, rdata (reads the per-sample `.rds` counts).

### 3. Notebooks (outputs saved)

| Notebook | Content |
|---|---|
| `01_QC.ipynb` | Nuclei kept per dataset, dropped samples |
| `02_atlas_overview.ipynb` | Cell-type counts, UMAPs by cell type and dataset, cell-type mix per dataset, agreement with author labels, dataset mixing per cluster, marker dot plots |
| `03_composition_test.ipynb` | PD vs control composition test; midbrain DA neurons per donor |
| `04_pseudobulk_DE.ipynb` | DE summary, volcano plots, top genes, stress-gene caveat |

## Results so far (test run)

- 96,019 nuclei loaded → **67,979 kept** after QC (71 %), 3,023 doublets; 37 small samples dropped.
- **38,895 annotated nuclei**, 209 donors, 14 cell types (DA_Neu, CADPS2_Neu, EN, IN, MSN, Astro, OPC, Oligo, Micro_PVM, Adaptive, Mural, Endo, Fibro, Ependymal), plus 6.7 % Unresolved.
- **99.5 % agreement** with author L1 labels (on the 2,200 nuclei that have author labels).
- Glia mix across datasets (94–100 % of glial nuclei sit in clusters shared by ≥ 8 datasets). Only 50 % of neurons do, which is expected because each region has different neurons.
- **Midbrain DA neurons are lower in PD**: 2.27 % → 0.90 %, lower in 6 of 7 datasets, p = 0.004, FDR = 0.048 (`04_composition_test.py`).
- **Pseudobulk DE** is concentrated in midbrain glia: Oligo 2,179 genes, Micro_PVM 123, Astro 106, OPC 76. Cortex and striatum: 0–3 genes.

## Open issues / questions

1. **Test scale only.** With ≤ 200 nuclei per donor, many donors have 0–5 DA neurons. The DA_Neu DE result (9 donors) is underpowered. Everything has to be repeated in the full run.
2. **Full run not set up yet.** Plan: Cayuga (`scu-cpu`, 768 GB nodes), `integration_mode = "sketch"`, which needs BPCells (not installed on the laptop). `presto` is also missing, so marker DE falls back to `FindAllMarkers` (slow).
3. **DE confounders.** The model is `~ dataset + condition` only, with no PMI, age or sex. Heat-shock genes (DNAJA1, DNAJB6, CHORDC1, PTGES3) go up and CIRBP goes down in every glial type: a shared stress / agonal / PMI signature, not necessarily PD. To do: add covariates where the metadata exist, and focus on cell-type-specific genes.
4. **Sorted libraries.** Neuron-sorted libraries inflate neuron fractions. For now they are excluded only in the composition analysis (metadata `enriched = TRUE`, plus midbrain libraries with > 40 % neurons). The 4 GSE253462 sorted libraries are not yet flagged in `config/sample_metadata.csv`, and the > 40 % rule is not yet part of the R pipeline.
5. **GSE235330** has no donor or region per cell yet, so it is used for annotation only, not for PD vs control.
6. **Striatum** has one dataset with 4 PD vs 4 control donors, which is too few for DE.
7. **Composition numbers vs slides.** The 18 Nov lab-meeting slide gave DA_Neu FDR = 0.034. That came from an earlier version of the test that tested fewer cell types. The script here tests every cell type present in ≥ 10 donors, which gives FDR = 0.048. The direction and p-value are essentially the same.
8. **OPC vs fibroblast (question from lab meeting).** PDGFRA is also expressed in fibroblasts. OPC calls require other panel markers (VCAN, OLIG2) and fibroblasts have their own panel (COL1A2, DCN), but this still needs an explicit check.
9. **Shallow datasets lose most nuclei at QC** (GSE184950 keeps 29 %, GSE178146 44 %). Should these datasets get their own thresholds?
10. **Other diagnoses** (PDD, DLB, LBD, AD in GSE178146, GSE178265 and GSE161045) are integrated but not yet analysed. A Lewy-body spectrum vs control comparison is a possible next step.
11. **Author labels** are available for only a few datasets. More can be added through the `author_meta_*` columns in `config/datasets.csv`.
