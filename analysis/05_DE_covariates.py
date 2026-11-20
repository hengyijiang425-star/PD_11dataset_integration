# Midbrain glia pseudobulk DE with donor covariates (age, sex, PMI), and the same donors without covariates.
# usage: python analysis/05_DE_covariates.py <cell_type> <cov|nocov>
import os, sys, time, warnings, numpy as np, pandas as pd; warnings.filterwarnings('ignore')
from pydeseq2.dds import DeseqDataSet; from pydeseq2.ds import DeseqStats; from pydeseq2.default_inference import DefaultInference
# Run from the project root.
ct, mode = sys.argv[1], sys.argv[2]
out=f'results/followup/DE_cov_Midbrain_{ct}_{mode}.csv'
X=np.load('results/pseudobulk/pb_counts.npy'); M=pd.read_csv('results/pseudobulk/pb_meta.csv',index_col=0); genes=pd.read_csv('results/pseudobulk/pb_genes.csv').gene.tolist()
cov=pd.read_csv('results/followup/donor_covariates.csv').set_index('donor')
M=M.join(cov,on='donor')
M['sex']=M.sex.astype(str).str.upper().str[0].map({'M':'M','F':'F'})
mid=M.region.astype(str).str.lower().str.contains('nigra|midbrain|snpc')
s=M[mid&M.diagnosis.isin(['PD','Control'])&(M.n_cells>=10)&(M.cell_type==ct)].dropna(subset=['age','sex','PMI'])
s=s[s.dataset.isin([d for d,x in s.groupby('dataset') if x.diagnosis.nunique()==2])].copy()
s['condition']=s.diagnosis
for v in ['age','PMI']: s[v+'_z']=(s[v]-s[v].mean())/s[v].std()
C=pd.DataFrame(X[[M.index.get_loc(i) for i in s.index]],index=s.index,columns=genes)
C=C.loc[:,(C>=5).sum(0)>=max(3,int(0.25*len(C)))]
design='~dataset + age_z + sex + PMI_z + condition' if mode=='cov' else '~dataset + condition'
t0=time.time()
dds=DeseqDataSet(counts=C,metadata=s[['dataset','condition','sex','age_z','PMI_z']],design=design,refit_cooks=True,inference=DefaultInference(n_cpus=4),quiet=True)
dds.deseq2(); st=DeseqStats(dds,contrast=['condition','PD','Control'],inference=DefaultInference(n_cpus=4),quiet=True); st.summary()
r=st.results_df.sort_values('pvalue'); r.insert(0,'gene',r.index); r['n_PD']=(s.condition=='PD').sum(); r['n_Ctrl']=(s.condition=='Control').sum(); r['datasets']=s.dataset.nunique()
r.to_csv(out,index=False)
print(ct,mode,s.condition.value_counts().to_dict(),s.dataset.nunique(),'datasets',C.shape,'DE',(r.padj<0.05).sum(),round(time.time()-t0),'s')
