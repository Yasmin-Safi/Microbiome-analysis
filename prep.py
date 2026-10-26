import pandas as pd, numpy as np, warnings
warnings.filterwarnings('ignore')
U='/mnt/user-data/uploads/'
cnt=pd.read_excel(U+'Taxonomy_abundance_count.xlsx',sheet_name='Species')
tax=['Kingdom','Phylum','Class','Order','Family','Genus','Species']
samp=[c for c in cnt.columns if c[:3].isdigit()]
sp=cnt['Species'].str.replace('__','',regex=False).values
gen=cnt['Genus'].str.replace('__','',regex=False).values
CT=cnt[samp].T.astype(float); CT.columns=sp
meta=pd.DataFrame({'sample':samp})
meta['pid']=meta['sample'].str[:3].astype(int); meta['tp']=meta['sample'].str[-2:]
meta=meta.set_index('sample')
depth=CT.sum(axis=1)
R=CT.div(depth,axis=0)
complete=meta.groupby('pid')['tp'].nunique(); complete=complete[complete==3].index
meta['complete']=meta.pid.isin(complete)
# patient data
p=pd.read_excel(U+'patient_data.xlsx',sheet_name='Patient data'); p.columns=[c.strip() for c in p.columns]
sn=pd.read_excel(U+'patient_data.xlsx',sheet_name='Sample Numb').iloc[:62,:4]
