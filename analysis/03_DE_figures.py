"""Figures for pseudobulk DE: DE-count bar chart, midbrain volcano plots, genes shared across midbrain glia.
Input: results/pseudobulk/DE_*.csv from 02_pseudobulk_DE.py."""
import pandas as pd, numpy as np, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from adjustText import adjust_text
import os
# Run from the project root: python analysis/03_DE_figures.py
D="results/pseudobulk/"
OUT="results/pseudobulk/figures/"; os.makedirs(OUT,exist_ok=True)
plt.rcParams.update({"font.family":"DejaVu Sans","font.size":11,"axes.spines.top":False,"axes.spines.right":False})
UP,DN,GR="#C2410C","#1D4ED8","#B8B8B8"
s=pd.read_csv(D+"DE_summary.csv")
s["lab"]=s.region+" · "+s.cell_type+"  ("+s.donors_PD.astype(str)+" vs "+s.donors_Ctrl.astype(str)+")"
s=s.sort_values(["region","DE_padj05"],ascending=[True,True])
order=["Striatum","Cortex","Midbrain"]; s["r"]=s.region.map({k:i for i,k in enumerate(order)})
s=s.sort_values(["r","DE_padj05"])
fig,ax=plt.subplots(figsize=(9,6.2))
y=np.arange(len(s))
ax.barh(y,s.up_in_PD,color=UP,height=0.7,label="Up in PD")
ax.barh(y,-s.down_in_PD,color=DN,height=0.7,label="Down in PD")
ax.set_yticks(y); ax.set_yticklabels(s.lab,fontsize=10)
ax.set_xscale("symlog",linthresh=10)
ax.set_xticks([-1000,-100,-10,0,10,100,1000]); ax.set_xticklabels(["1000","100","10","0","10","100","1000"])
ax.axvline(0,color="#444",lw=0.8)
for yi,(u,d) in enumerate(zip(s.up_in_PD,s.down_in_PD)):
    t=u+d
    ax.text(max(u,0.5)*1.25 if u>0 else 1.2,yi,str(t),va="center",fontsize=9,color="#333")
ax.set_xlabel("Number of DE genes (padj < 0.05), symlog scale")
ax.set_title("Pseudobulk DE, PD vs control (DESeq2, ~dataset + condition)",loc="left",fontsize=12)
ax.legend(frameon=False,loc="lower right"); fig.tight_layout(); fig.savefig(OUT+"de_counts.png",dpi=200)

fig,axs=plt.subplots(1,3,figsize=(15,5.2))
for ax,ct in zip(axs,["Micro_PVM","Astro","Oligo"]):
    d=pd.read_csv(D+f"DE_Midbrain_{ct}.csv").dropna(subset=["padj"])
    d["nl"]=-np.log10(d.padj.clip(lower=1e-300))
    sig=d.padj<0.05; up=sig&(d.log2FoldChange>0); dn=sig&(d.log2FoldChange<0)
    ax.scatter(d.log2FoldChange[~sig],d.nl[~sig],s=4,c=GR,alpha=.5,lw=0)
    ax.scatter(d.log2FoldChange[up],d.nl[up],s=8,c=UP,lw=0)
    ax.scatter(d.log2FoldChange[dn],d.nl[dn],s=8,c=DN,lw=0)
    ax.axhline(-np.log10(0.05),ls="--",c="#888",lw=.8)
    top=pd.concat([d[up].nsmallest(6,"padj"),d[dn].nsmallest(6,"padj")])
    ts=[ax.text(r.log2FoldChange,r.nl,r.gene,fontsize=8.5) for _,r in top.iterrows()]
    adjust_text(ts,ax=ax,x=d.log2FoldChange[sig].values,y=d.nl[sig].values,arrowprops=dict(arrowstyle='-',color='#666',lw=.5),expand=(1.3,1.6))
    ax.set_xlim(-3.5,3.5)
    ax.set_title(f"Midbrain {ct}: {up.sum()} up / {dn.sum()} down",loc="left",fontsize=11.5)
    ax.set_xlabel("log2 fold change (PD / control)")
axs[0].set_ylabel("-log10 adjusted p")
fig.tight_layout(); fig.savefig(OUT+"volcano_midbrain.png",dpi=180)

# shared genes across midbrain glia
sets={}
for ct in ["Micro_PVM","Astro","Oligo","OPC"]:
    d=pd.read_csv(D+f"DE_Midbrain_{ct}.csv").dropna(subset=["padj"])
    sets[ct]=d[d.padj<0.05].set_index("gene").log2FoldChange
genes=set.union(*[set(v.index) for v in sets.values()])
cnt={g:sum(g in v.index for v in sets.values()) for g in genes}
shared=sorted([g for g,c in cnt.items() if c>=3],key=lambda g:-cnt[g])
rows=[]
for g in shared: rows.append({"gene":g,"n_celltypes":cnt[g],**{ct:round(v.get(g,np.nan),2) for ct,v in sets.items()}})
sh=pd.DataFrame(rows).sort_values(["n_celltypes","gene"],ascending=[False,True]); sh.to_csv(OUT+"shared_DE_midbrain_glia.csv",index=False)
print(sh.to_string())
# direction concordance for genes in >=2
c2=[g for g,c in cnt.items() if c>=2]
conc=sum(len(set(np.sign([v[g] for v in sets.values() if g in v.index])))==1 for g in c2)
print("genes in >=2 celltypes",len(c2),"same direction",conc)

# heatmap of genes DE in >=3 midbrain glial types
cts=["Micro_PVM","Astro","Oligo","OPC"]; M=sh[cts].values
fig,ax=plt.subplots(figsize=(6.2,6.4))
im=ax.imshow(np.ma.masked_invalid(M),cmap="RdBu_r",vmin=-2,vmax=2,aspect="auto"); ax.set_facecolor("#EEEEEE")
ax.set_xticks(range(len(cts))); ax.set_xticklabels(cts); ax.set_yticks(range(len(sh))); ax.set_yticklabels(sh.gene)
for i in range(M.shape[0]):
    for j in range(M.shape[1]):
        v=M[i,j]
        if not np.isnan(v): ax.text(j,i,f"{v:.1f}",ha="center",va="center",fontsize=9,color="white" if abs(v)>1.2 else "#222")
for s_ in ax.spines.values(): s_.set_visible(False)
ax.set_xticks(np.arange(-.5,len(cts)),minor=True); ax.set_yticks(np.arange(-.5,len(sh)),minor=True)
ax.grid(which="minor",color="white",lw=2); ax.tick_params(which="minor",length=0)
cb=fig.colorbar(im,ax=ax,shrink=.6); cb.set_label("log2FC (PD/control); grey = not DE")
ax.set_title("Genes DE in ≥3 midbrain glial types",loc="left",fontsize=12)
fig.tight_layout(); fig.savefig(OUT+"shared_genes.png",dpi=200)
