exec(open('prep.py').read())
import json
from scipy import stats
rng=np.random.default_rng(7)
pd.set_option('display.width',250,'display.max_columns',30,'display.max_rows',60)
tps=['S1','S2','S3']; tpv=meta.tp.values; cm=meta.complete.values
# original definition (all samples, detected >0)
orig={t:int(((CT[tpv==t]>0).mean(axis=0)>0.5).sum()) for t in tps}
orig20={t:int(((CT[tpv==t]>0).mean(axis=0)>0.2).sum()) for t in tps}
print('orig >50%',orig,'>20%',orig20)
# complete cohort only, raw
cc={t:int(((CT[(tpv==t)&cm]>0).mean(axis=0)>0.5).sum()) for t in tps}; print('complete raw >50%',cc)
# rarefied sensitivity: 100 iterations at 7500; thresholds: >0 reads, >=3 reads (~0.04%)
def rarefy(M,d):
    return np.array([rng.multivariate_hypergeometric(r.astype(int),d) for r in M])
res={}
NIT=50
for lab,mask in [('all samples',np.ones(len(tpv),bool)),('complete cohort (n=46/TP)',cm)]:
    for thr in [1,3]:
        counts={t:[] for t in tps}; counts20={t:[] for t in tps}
        for it in range(NIT):
            M=rarefy(CT.values,7500)
            for t in tps:
                sel=M[(tpv==t)&mask]
                prev=(sel>=thr).mean(axis=0); counts[t].append(int((prev>0.5).sum())); counts20[t].append(int((prev>0.2).sum()))
        res[(lab,thr)]={t:(np.mean(counts[t]),np.percentile(counts[t],2.5),np.percentile(counts[t],97.5)) for t in tps}
        res[(lab,thr,'20')]={t:(np.mean(counts20[t]),) for t in tps}
core=[]
for k,v in res.items():
    if len(k)==2: core.append({'set':k[0],'min_reads':k[1],**{f'{t}_core50':round(v[t][0],1) for t in tps},**{f'{t}_core20':round(res[(k[0],k[1],'20')][t][0],1) for t in tps}})
CD=pd.DataFrame(core); print(CD.to_string()); CD.to_csv('core_sensitivity.csv',index=False)
# overlap: of S1 core species (original), fraction still >50% at S3 (complete, rarefied one iteration)
M=rarefy(CT.values,7500)
c1=set(CT.columns[((M[(tpv=='S1')&cm]>=1).mean(axis=0)>0.5)]); c3=set(CT.columns[((M[(tpv=='S3')&cm]>=1).mean(axis=0)>0.5)])
print(len(c1),len(c3),len(c1&c3))
json.dump({'orig':orig,'orig20':orig20,'complete_raw':cc,'overlap':[len(c1),len(c3),len(c1&c3)]},open('core.json','w'))
# genus table
g=pd.read_csv('da_genus.csv').set_index('feature')
kg=['Fusobacterium','Bacteroides','Streptococcus','Lactobacillus','Limosilactobacillus','Lactiplantibacillus','Parvimonas','Peptostreptococcus','Gemella','Porphyromonas','Turicibacter','Ruminococcus','Faecalibacterium','Akkermansia','Veillonella','Lacticaseibacillus','Lacticaseibacillus']
kg=[k for k in kg if k in g.index]
print(g.loc[kg,['prev','mean_S1','mean_S2','mean_S3','KW_p','KW_q_filt','W_S1S2_q','W_S2S3_q','W_S1S3_q','LMM_S2vS1_q','LMM_S3vS2_q','LMM_S3vS1_q','LMMd_S3vS1_q']].round(4).to_string())
print(g[g.KW_q_filt<0.05][['KW_q_filt','mean_S1','mean_S2','mean_S3']].round(5).to_string())
# alpha SD of S2->S3 changes for tipping
Ar=pd.read_csv('alpha_rarefied.csv',index_col=0)
for m in ['Observed','Shannon']:
    s=Ar[m]; a=s[(tpv=='S2')&cm].values; b=s[(tpv=='S3')&cm].values
    pid2=meta.pid[(tpv=='S2')&cm].values; pid3=meta.pid[(tpv=='S3')&cm].values
    dd=pd.Series(b,index=pid3)-pd.Series(a,index=pid2); print(m,'SD S2->S3 change',round(dd.std(),3),'mean',round(dd.mean(),3))
