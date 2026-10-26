exec(open('prep.py').read())
from scipy import stats
import json
from statsmodels.stats.multitest import multipletests
pd.set_option('display.width',250,'display.max_columns',30)
out={}
# ---- QC ----
q=pd.read_excel(U+'Qulaity_control_stats.xlsx',header=3).iloc[:,:7].dropna(subset=['SampleName'])
q=q[q.SampleName.astype(str).str.contains(r'^\d{3}\.S\d$')].set_index('SampleName')
for c in q.columns: q[c]=pd.to_numeric(q[c])
q['tp']=[s[-2:] for s in q.index]; q['len']=q['Total Bases']/q['Read Count']
q=q.join(depth.rename('table_reads'),how='left'); q['frac']=q.table_reads/q['Read Count']
print(len(q), q.table_reads.isna().sum())
g=q.groupby('tp').agg(n=('Read Count','size'),reads_med=('Read Count','median'),reads_q1=('Read Count',lambda x:x.quantile(.25)),reads_q3=('Read Count',lambda x:x.quantile(.75)),len_mean=('len','mean'),Q20=('Q20 (%)','mean'),Q30=('Q30 (%)','mean'),GC=('GC (%)','mean'),tab_med=('table_reads','median'),frac_med=('frac','median'),frac_q1=('frac',lambda x:x.quantile(.25)),frac_q3=('frac',lambda x:x.quantile(.75)))
print(g.round(3)); g.to_csv('qc_by_tp.csv')
print('KW merged reads',stats.kruskal(*[q['Read Count'][q.tp==t] for t in ['S1','S2','S3']]).pvalue, 'KW frac',stats.kruskal(*[q.frac[q.tp==t] for t in ['S1','S2','S3']]).pvalue)
print('overall merged reads',q['Read Count'].median(),q['Read Count'].min(),q['Read Count'].max(),' frac',q.frac.median(),q.frac.min(),q.frac.max())
print('Q30 range',q['Q30 (%)'].min(),q['Q30 (%)'].max(), 'len', q['len'].mean())
# ---- attrition ----
p2=p.copy(); p2['pid']=p2['Case ID']
p2['group']=np.where(p2.pid.isin(complete),'complete','incomplete')
sn.columns=['pid','S1','S2','S3']; sn['S3']=sn.S3.str.strip(); sn['S2']=sn.S2.str.strip()
reason=sn.loc[~sn.pid.isin(complete),['pid','S2','S3']]
print(reason.assign(S3=reason.S3.str.lower()).S3.value_counts())
# sample presence
pres=meta.reset_index().groupby('pid').tp.apply(set)
miss=p2[~p2.pid.isin(complete)].pid.tolist()
print('incomplete ids',miss, {i:sorted(pres[i]) for i in miss})
rows=[]
cont=['Age','BMI','TUMOR SIZE']
for c in cont:
    a=p2[c][p2.group=='complete']; b=p2[c][p2.group=='incomplete']
    rows.append((c,f'{a.mean():.1f} ± {a.std():.1f}',f'{b.mean():.1f} ± {b.std():.1f}',stats.mannwhitneyu(a,b).pvalue))
def fisher(col,pos,label):
    t=pd.crosstab(p2.group,p2[col].astype(str).str.strip().str.lower()==pos); 
    a=int(t.loc['complete',True]) if True in t.columns else 0; b=int(t.loc['incomplete',True]) if True in t.columns else 0
    rows.append((label,f'{a}/46 ({100*a/46:.0f}%)',f'{b}/16 ({100*b/16:.0f}%)',stats.fisher_exact([[a,46-a],[b,16-b]])[1]))
fisher('Gender','male','Male sex'); fisher('Smoking','yes','Current/ex-smoker'); fisher('COMOR(Y/N)','yes','Any comorbidity')
p2['cN2']=p2.Cn=='N2'; p2['T34']=p2.Ct.isin(['T3','T4'])
fisher('Cn','n2','cN2'); fisher('Ct','t4','cT4')
# chemo regimen, cycles
p2['cyc']=p2.cycle.astype(str).str.extract(r'(\d+)').astype(float)
a=p2.cyc[p2.group=='complete']; b=p2.cyc[p2.group=='incomplete']; rows.append(('Chemotherapy cycles, median',f'{a.median():.0f} ({a.min():.0f}-{a.max():.0f})',f'{b.median():.0f} ({b.min():.0f}-{b.max():.0f})',stats.mannwhitneyu(a,b).pvalue))
# baseline microbiome: S1 alpha & depth
A=pd.read_excel(U+'Microbiome_Analysis_Results.xlsx',sheet_name='Alpha_Diversity').set_index('Sample')
s1=A[A.Timepoint=='S1'].copy(); s1['pid']=[int(s[:3]) for s in s1.index]; s1['grp']=np.where(s1.pid.isin(complete),'complete','incomplete')
for m,l in [('Shannon','Baseline Shannon'),('Observed_Species','Baseline observed richness')]:
    a=s1[m][s1.grp=='complete']; b=s1[m][s1.grp=='incomplete']; rows.append((l,f'{a.mean():.2f} ± {a.std():.2f}',f'{b.mean():.2f} ± {b.std():.2f}',stats.mannwhitneyu(a,b).pvalue))
A2=A.copy(); A2['pid']=[int(s[:3]) for s in A2.index]
# S2 values in those lacking S3
s2=A2[A2.Timepoint=='S2']; s2['grp']=np.where(s2.pid.isin(complete),'complete','incomplete')
for m,l in [('Shannon','Post-RT (S2) Shannon'),('Observed_Species','Post-RT (S2) observed richness')]:
    a=s2[m][s2.grp=='complete']; b=s2[m][s2.grp=='incomplete']; rows.append((l,f'{a.mean():.2f} ± {a.std():.2f} (n={len(a)})',f'{b.mean():.2f} ± {b.std():.2f} (n={len(b)})',stats.mannwhitneyu(a,b).pvalue))
AT=pd.DataFrame(rows,columns=['Variable','Complete (n=46)','Incomplete (n=16)','p']); print(AT.round(3).to_string()); AT.to_csv('attrition_compare.csv',index=False)
# ---- delta-adjusted tipping point (pattern mixture) for S2->S3 change ----
Ar=pd.read_csv('alpha_rarefied.csv',index_col=0)
def tipping(metric,Adf):
    s=Adf[metric]; W=pd.DataFrame({'S1':s[meta.tp=='S1'].values,},index=meta.pid[meta.tp=='S1'])
    W['S2']=pd.Series(s[meta.tp=='S2'].values,index=meta.pid[meta.tp=='S2']); W['S3']=pd.Series(s[meta.tp=='S3'].values,index=meta.pid[meta.tp=='S3'])
    obs=W.loc[complete]; dchg=(obs.S3-obs.S2); sd=dchg.std()
    res=[]
    for k in [0,0.5,1.0,1.5,2.0]:
        w=W.copy(); mS3=w.S3.isna()&w.S2.notna(); mx=w.S3.isna()&w.S2.isna()
        w.loc[mS3,'S3']=w.loc[mS3,'S2']+dchg.mean()-k*sd   # impute: observed mean S2->S3 change shifted down by k SD
        # patients with no S2 and no S3: not included in S2vS3 test (need S2)
        ww=w.dropna(subset=['S2','S3'])
        pw=stats.wilcoxon(ww.S2,ww.S3).pvalue; md=float((ww.S3-ww.S2).median())
        res.append((metric,k,len(ww),round(md,3),pw))
    return res
tp_rows=[]
for m in ['Observed','Shannon']: tp_rows+=tipping(m,Ar)
TP=pd.DataFrame(tp_rows,columns=['metric','delta_SD','n_pairs','median_change','p_S2vsS3']); print(TP.round(4)); TP.to_csv('tipping.csv',index=False)
