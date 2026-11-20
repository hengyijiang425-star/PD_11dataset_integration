# Pseudobulk counts per donor x cell type from the test-run annotations (python; reads the per-sample count .rds files)
import os, re, csv, numpy as np, pandas as pd, rdata, scipy.sparse as sp, time
root=os.path.expanduser('~/mnt/PD_11dataset_integration'); os.chdir(root)
ann=pd.read_csv('results/PD11_final_cell_annotations.csv.gz',low_memory=False)
enr=set(pd.read_csv('results/composition_test/neuron_enriched_libraries.csv').sample_uid)
ann=ann[(ann.enriched!=True)&(~ann.sample_uid.isin(enr))&(ann.cell_type!='Unresolved')]
genes=pd.read_csv('results/05_shared_genes.csv').gene.tolist(); gi={g:i for i,g in enumerate(genes)}
man=pd.read_csv('data/processed/manifest.csv')
man['sample_uid']=man.dataset+'_'+man.sample_key.map(lambda k: re.sub(r'[^A-Za-z0-9._]','.',k if re.match(r'^[A-Za-z.]',k) else 'X'+k))
need=set(ann.sample_uid); man=man[man.sample_uid.isin(need)]
print('samples to read',len(man),flush=True)
conv={**rdata.conversion.DEFAULT_CLASS_MAP,'dgCMatrix':lambda obj,attrs: attrs}
rows=[]; meta=[]; t0=time.time()
for k,(_,r) in enumerate(man.iterrows()):
    a=rdata.conversion.convert(rdata.parser.parse_file(r.path),constructor_dict=conv)
    nr,nc=[int(v) for v in a['Dim']]; m=sp.csc_matrix((np.asarray(a['x'],float),np.asarray(a['i'],int),np.asarray(a['p'],int)),shape=(nr,nc))
    g=[str(x) for x in a['Dimnames'][0]]; b=[str(x) for x in a['Dimnames'][1]]
    sel=[i for i,x in enumerate(g) if x in gi]; tgt=[gi[g[i]] for i in sel]
    m=m[sel,:]
    sub=ann[ann.sample_uid==r.sample_uid]; bc=sub.cell.str.slice(len(r.sample_uid)+1)
    pos={x:i for i,x in enumerate(b)}
    for ct,grp in sub.groupby('cell_type'):
        cols=[pos[x] for x in bc[grp.index] if x in pos]
        if not cols: continue
        v=np.zeros(len(genes)); v[tgt]=np.asarray(m[:,cols].sum(axis=1)).ravel()
        rows.append(v); x=grp.iloc[0]
        meta.append(dict(sample_uid=r.sample_uid,dataset=x.dataset,donor=x.donor,diagnosis=x.diagnosis,region=x.region,cell_type=ct,n_cells=len(cols)))
    if k%25==0: print(k,round(time.time()-t0),'s',flush=True)
M=pd.DataFrame(meta); X=np.vstack(rows).astype(np.int64)
# sum libraries of the same donor
M['key']=M.donor+'|'+M.cell_type; keys=M.key.unique()
agg=pd.DataFrame(X).groupby(M.key.values).sum().loc[keys]
Md=M.groupby('key').agg(dataset=('dataset','first'),donor=('donor','first'),diagnosis=('diagnosis','first'),region=('region','first'),cell_type=('cell_type','first'),n_cells=('n_cells','sum')).loc[keys]
agg.columns=genes; agg.to_parquet('results/pseudobulk/pseudobulk_counts.parquet') if False else None
np.save('results/pseudobulk/pb_counts.npy',agg.values); Md.to_csv('results/pseudobulk/pb_meta.csv'); pd.Series(genes).to_csv('results/pseudobulk/pb_genes.csv',index=False,header=['gene'])
print('done',agg.shape,flush=True)
