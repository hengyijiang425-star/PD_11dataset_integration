# PD vs Control pseudobulk DESeq2 (pydeseq2) per region x cell type; design ~ dataset + condition
import os, numpy as np, pandas as pd, warnings, time; warnings.filterwarnings('ignore')
from pydeseq2.dds import DeseqDataSet; from pydeseq2.ds import DeseqStats
from pydeseq2.default_inference import DefaultInference
os.chdir(os.path.expanduser('~/mnt/PD_11dataset_integration/results/pseudobulk'))
X=np.load('pb_counts.npy'); M=pd.read_csv('pb_meta.csv',index_col=0); genes=pd.read_csv('pb_genes.csv').gene.tolist()
def rg(r):
    r=str(r).lower()
    if any(k in r for k in ['nigra','midbrain','snpc']): return 'Midbrain'
    if any(k in r for k in ['cortex','cingulate','frontal']): return 'Cortex'
    return 'Striatum'
M['region_group']=M.region.map(rg); M['condition']=M.diagnosis
keep=M.condition.isin(['PD','Control'])&(M.n_cells>=10); M=M[keep]; X=X[keep.values]
import sys, glob
summ=[]
for (reg,ct),s in M.groupby(['region_group','cell_type']):
    if os.path.exists(f'DE_{reg}_{ct}.csv'): continue
    if len(sys.argv)>1 and f'{reg}|{ct}'!=sys.argv[1]: continue
    n=s.condition.value_counts()
    if n.get('PD',0)<4 or n.get('Control',0)<4: continue
    ds_ok=[d for d,x in s.groupby('dataset') if x.condition.nunique()==2]; s=s[s.dataset.isin(ds_ok)]
    n=s.condition.value_counts()
    if n.get('PD',0)<4 or n.get('Control',0)<4: continue
    C=pd.DataFrame(X[[M.index.get_loc(i) for i in s.index]],index=s.index,columns=genes)
    C=C.loc[:,(C>=5).sum(0)>=max(3,int(0.25*len(C)))]
    meta=s[['dataset','condition']].copy()
    t0=time.time()
    try:
        dds=DeseqDataSet(counts=C,metadata=meta,design='~dataset + condition',refit_cooks=True,inference=DefaultInference(n_cpus=4),quiet=True)
        dds.deseq2(); st=DeseqStats(dds,contrast=['condition','PD','Control'],inference=DefaultInference(n_cpus=4),quiet=True); st.summary()
        r=st.results_df.sort_values('pvalue'); r.insert(0,'gene',r.index); r['region']=reg; r['cell_type']=ct
        r.to_csv(f'DE_{reg}_{ct}.csv',index=False)
        sig=r[(r.padj<0.05)]
        summ.append(dict(region=reg,cell_type=ct,donors_PD=int(n.get('PD',0)),donors_Ctrl=int(n.get('Control',0)),datasets=len(ds_ok),genes_tested=C.shape[1],
            DE_padj05=len(sig),up_in_PD=int((sig.log2FoldChange>0).sum()),down_in_PD=int((sig.log2FoldChange<0).sum()),
            top_up=', '.join(sig[sig.log2FoldChange>0].gene.head(8)),top_down=', '.join(sig[sig.log2FoldChange<0].gene.head(8))))
        print(reg,ct,C.shape,len(sig),round(time.time()-t0),'s',flush=True)
    except Exception as e: print('FAILED',reg,ct,e,flush=True)
    old=pd.read_csv('DE_summary.csv') if os.path.exists('DE_summary.csv') else pd.DataFrame()
    pd.concat([old,pd.DataFrame(summ)]).drop_duplicates(['region','cell_type'],keep='last').to_csv('DE_summary.csv',index=False); summ=[]
print('ALL DONE',flush=True)
