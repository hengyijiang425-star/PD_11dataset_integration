"""PDGFRA question: pseudobulk CPM of OPC and fibroblast marker genes in every cell type (all donors pooled).
Input results/pseudobulk/pb_counts.npy, pb_meta.csv, pb_genes.csv. Output results/followup/pdgfra_check_cpm.csv. Run from the project root."""
import os, numpy as np, pandas as pd
os.makedirs("results/followup", exist_ok=True)
X = np.load("results/pseudobulk/pb_counts.npy"); M = pd.read_csv("results/pseudobulk/pb_meta.csv")
G = pd.read_csv("results/pseudobulk/pb_genes.csv").gene.tolist(); gi = {g: i for i, g in enumerate(G)}
genes = [g for g in ["PDGFRA","VCAN","OLIG2","OLIG1","SOX10","PTPRZ1","BCAN","COL1A2","COL1A1","DCN","LAMA2","SLC6A13","CEMIP","PDGFRB","RGS5","CLDN5","MBP","GFAP"] if g in gi]
S = pd.DataFrame(X).groupby(M.cell_type.values).sum(); cpm = S.div(S.sum(axis=1), axis=0) * 1e6
T = cpm[[gi[g] for g in genes]]; T.columns = genes
out = T.round(1).join(M.groupby("cell_type").agg(nuclei=("n_cells", "sum"), donors=("donor", "nunique")))
out.to_csv("results/followup/pdgfra_check_cpm.csv"); print(out.to_string())
