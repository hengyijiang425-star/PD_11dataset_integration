# Follow-up figures: PDGFRA check, DE with covariates, DA subtypes. Run from the project root after 05 and 06.
import pandas as pd,numpy as np,matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import os
U="results/followup/"; os.makedirs(U+"figures",exist_ok=True)
plt.rcParams.update({"font.family":"DejaVu Sans","font.size":11,"axes.spines.top":False,"axes.spines.right":False})
UP,DN="#C2410C","#1D4ED8"
# A: PDGFRA check
c=pd.read_csv(U+"pdgfra_check_cpm.csv",index_col=0)
cts=["OPC","Fibro","Oligo","Astro","Mural","Endo","Micro_PVM"]; gs=["PDGFRA","VCAN","OLIG2","SOX10","PTPRZ1","COL1A2","COL1A1","DCN","LAMA2","SLC6A13"]
L=np.log10(c.loc[cts,gs]+1)
fig,ax=plt.subplots(figsize=(10,4.6))
im=ax.imshow(L.values,cmap="Greys",aspect="auto",vmin=0,vmax=4)
ax.set_xticks(range(len(gs))); ax.set_xticklabels(gs,rotation=35,ha="right"); ax.set_yticks(range(len(cts)))
ax.set_yticklabels([f"{t}  (n={int(c.loc[t,'nuclei']):,})" for t in cts])
for i in range(len(cts)):
    for j in range(len(gs)):
        v=c.loc[cts[i],gs[j]]; ax.text(j,i,f"{v:,.0f}",ha="center",va="center",fontsize=8.5,color="white" if L.values[i,j]>2.4 else "#222")
ax.axvline(4.5,color="#C2410C",lw=2)
ax.text(2,-0.85,"OPC panel",ha="center",fontsize=10.5,color="#444"); ax.text(7.5,-0.85,"Fibroblast panel",ha="center",fontsize=10.5,color="#444")
for s in ax.spines.values(): s.set_visible(False)
cb=fig.colorbar(im,ax=ax,shrink=.8); cb.set_label("log10(CPM + 1)")
ax.set_title("Pseudobulk expression (CPM), all donors pooled",loc="left",fontsize=11.5,pad=22)
fig.tight_layout(); fig.savefig(U+"figures/pdgfra_check.png",dpi=200)

# B: DE with covariates
cts=["Micro_PVM","Astro","Oligo","OPC"]; a=[];b=[];k=[]
for ct in cts:
    x=pd.read_csv(U+f"DE_cov_Midbrain_{ct}_nocov.csv"); y=pd.read_csv(U+f"DE_cov_Midbrain_{ct}_cov.csv")
    sa=set(x.gene[x.padj<0.05]); sb=set(y.gene[y.padj<0.05]); a.append(len(sa)); b.append(len(sb)); k.append(len(sa&sb))
    print(ct,x.n_PD[0],x.n_Ctrl[0],x.datasets[0],len(sa),len(sb),len(sa&sb))
fig,axs=plt.subplots(1,2,figsize=(14,5.2),gridspec_kw={"width_ratios":[1,1.35]})
ax=axs[0]; xx=np.arange(len(cts)); w=.38
ax.bar(xx-w/2,a,w,color="#9CA3AF",label="~ dataset + PD"); ax.bar(xx+w/2,b,w,color="#2f6fb0",label="+ age, sex, PMI")
for i in range(len(cts)): ax.text(i-w/2,a[i]*1.08,str(a[i]),ha="center",fontsize=9); ax.text(i+w/2,b[i]*1.08 if b[i] else 1.2,str(b[i]),ha="center",fontsize=9)
ax.set_yscale("symlog",linthresh=10); ax.set_yticks([0,10,100,1000]); ax.set_yticklabels(["0","10","100","1000"])
ax.set_xticks(xx); ax.set_xticklabels(cts); ax.set_ylabel("DE genes (padj < 0.05)"); ax.legend(frameon=False,loc="upper left")
ax.set_title("Same donors, two models (midbrain)",loc="left",fontsize=11.5)
w_=pd.read_csv(U+"watch_genes_cov_vs_nocov.csv")
keep=[("Oligo","HSPA1A"),("Oligo","DNAJA1"),("Oligo","CIRBP"),("Micro_PVM","HSPH1"),("Micro_PVM","HSPB1"),("Astro","HSPH1"),
      ("Micro_PVM","CALCOCO2"),("Astro","PSAT1"),("Astro","SLC7A11"),("Oligo","SLC38A2"),("Oligo","IRS2"),("Oligo","ABCC8")]
ax=axs[1]
for i,(ct,g) in enumerate(keep):
    r=w_[(w_.cell_type==ct)&(w_.gene==g)].iloc[0]
    stress=i<6; col="#B45309" if stress else "#2f6fb0"
    ax.plot([r.lfc_nocov,r.lfc_cov],[i,i],color="#bbb",lw=2,zorder=1)
    ax.scatter(r.lfc_nocov,i,s=40,facecolor="white",edgecolor=col,lw=1.6,zorder=2)
    ax.scatter(r.lfc_cov,i,s=48,color=col,zorder=3)
    ax.text(2.25,i,("padj %.2g" % r.padj_cov),va="center",fontsize=9,color="#111" if r.padj_cov<0.05 else "#999")
ax.set_yticks(range(len(keep))); ax.set_yticklabels([f"{g} · {c}" for c,g in keep]); ax.invert_yaxis()
ax.axvline(0,color="#444",lw=.8); ax.set_xlim(-2.2,2.9); ax.set_xlabel("log2 fold change PD / control   (open = before, filled = with age, sex, PMI)")
ax.axhline(5.5,color="#ddd",lw=1)
ax.text(-2.15,0.35,"Stress genes",color="#B45309",fontsize=10); ax.text(-2.15,6.1,"Cell-type genes",color="#2f6fb0",fontsize=10)
ax.set_title("Effect sizes before and after covariates",loc="left",fontsize=11.5,pad=10)
fig.tight_layout(); fig.savefig(U+"figures/de_covariates.png",dpi=180)

# C: DA subtypes
d=pd.read_csv(U+"DA_subtype_per_donor.csv")
s=d[(d.n>=10)&d.diagnosis.isin(["PD","Control"])].copy()
s=s[s.dataset.isin([x for x,g in s.groupby("dataset") if g.diagnosis.nunique()==2])]
fig,axs=plt.subplots(1,3,figsize=(13,4.6),sharey=True)
dsc=dict(zip(sorted(s.dataset.unique()),["#0F766E","#7C3AED","#B45309"]))
rng=np.random.default_rng(1)
for ax,(v,lab) in zip(axs,[("AGTR1p","AGTR1+"),("SOX6p","SOX6+ (CALB1−)"),("CALB1p","CALB1+ (SOX6−)")]):
    for j,dx in enumerate(["Control","PD"]):
        g=s[s.diagnosis==dx]; f=100*g[v]/g.n
        ax.scatter(j+rng.uniform(-.15,.15,len(g)),f,c=[dsc[x] for x in g.dataset],s=36,alpha=.85,lw=0)
        ax.hlines(f.mean(),j-.28,j+.28,color="#111",lw=2)
    ax.set_xticks([0,1]); ax.set_xticklabels(["Control","PD"]); ax.set_xlim(-.6,1.6); ax.set_title(lab+" among DA nuclei",loc="left",fontsize=11.5)
axs[0].set_ylabel("% of the donor's DA nuclei")
for k,c_ in dsc.items(): axs[2].scatter([],[],c=c_,label=k)
axs[2].legend(frameon=False,fontsize=9,loc="upper right")
fig.tight_layout(); fig.savefig(U+"figures/da_subtypes.png",dpi=180)
print(s.groupby(["dataset","diagnosis"]).apply(lambda g:(100*g.AGTR1p/g.n).mean().round(1)).to_string())
