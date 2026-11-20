# Detection of DA-subtype genes in every DA_Neu nucleus (test run); reads the per-sample count .rds files.
import os, re, time, numpy as np, pandas as pd, rdata, scipy.sparse as sp
# Run from the project root.
GENES=['TH','SLC6A3','SLC18A2','KCNJ6','ALDH1A1','SOX6','AGTR1','CALB1','CALCR','CRYM','GFRA1','RET','NR4A2']
ann=pd.read_csv('results/PD11_final_cell_annotations.csv.gz',low_memory=False)
da=ann[ann.cell_type=='DA_Neu']
man=pd.read_csv('data/processed/manifest.csv')
man['sample_uid']=man.dataset+'_'+man.sample_key.map(lambda k: re.sub(r'[^A-Za-z0-9._]','.',k if re.match(r'^[A-Za-z.]',k) else 'X'+k))
man=man[man.sample_uid.isin(set(da.sample_uid))]
print('samples',len(man),'DA nuclei',len(da),flush=True)
conv={**rdata.conversion.DEFAULT_CLASS_MAP,'dgCMatrix':lambda obj,attrs: attrs}
out=[]; t0=time.time()
for k,(_,r) in enumerate(man.iterrows()):
    a=rdata.conversion.convert(rdata.parser.parse_file(r.path),constructor_dict=conv)
    nr,nc=[int(v) for v in a['Dim']]; m=sp.csc_matrix((np.asarray(a['x'],float),np.asarray(a['i'],int),np.asarray(a['p'],int)),shape=(nr,nc))
    g=[str(x) for x in a['Dimnames'][0]]; b=[str(x) for x in a['Dimnames'][1]]; gpos={x:i for i,x in enumerate(g)}; pos={x:i for i,x in enumerate(b)}
    sub=da[da.sample_uid==r.sample_uid]; bc=sub.cell.str.slice(len(r.sample_uid)+1)
    keep=[(c,pos[x]) for c,x in zip(sub.cell,bc) if x in pos]
    if not keep: continue
    cols=[j for _,j in keep]; tot=np.asarray(m[:,cols].sum(axis=0)).ravel()
    d=pd.DataFrame({'cell':[c for c,_ in keep],'nUMI':tot})
    for G in GENES: d[G]=np.asarray(m[gpos[G],cols].todense()).ravel() if G in gpos else np.nan
    out.append(d)
    if k%10==0: print(k,round(time.time()-t0),'s',flush=True)
D=pd.concat(out).merge(da[['cell','sample_uid','dataset','donor','diagnosis','region','enriched']],on='cell')
D.to_csv('results/followup/DA_subtype_genes_per_nucleus.csv',index=False); print('done',D.shape,flush=True)
