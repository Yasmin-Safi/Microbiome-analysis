exec(open('prep.py').read())
from scipy import stats
from statsmodels.stats.multitest import multipletests
import statsmodels.formula.api as smf, json, warnings
warnings.filterwarnings('ignore')
rng=np.random.default_rng(5)
tps=['S1','S2','S3']; tpv=meta.tp.values; pid=meta.pid.values; cm=meta.complete.values
def bh(p):
    p=np.array(p,float); q=np.full(len(p),np.nan); m=~np.isnan(p); q[m]=multipletests(p[m],method='fdr_bh')[1]; return q
prev=(CT>0).mean(axis=0)
def clr(M,pc): 
    X=np.log(M+pc); return X-X.mean(1,keepdims=True)
def lmm_counts(pc,pmin):
    keep=np.where(prev.values>=pmin)[0]; X=clr(CT.values,pc)[:,keep]
    d=pd.DataFrame({'tp':pd.Categorical(tpv,categories=tps),'pid':pid}); r2=[];r3=[];r32=[]
    for j in range(len(keep)):
        d['y']=X[:,j]
        try:
            r=smf.mixedlm('y ~ C(tp)',d,groups=d['pid']).fit(reml=True)
            r2.append(r.pvalues['C(tp)[T.S2]']); r3.append(r.pvalues['C(tp)[T.S3]']); r32.append(float(r.t_test(np.array([0,-1,1.]).reshape(1,-1)).pvalue))
        except: r2.append(np.nan);r3.append(np.nan);r32.append(np.nan)
    q2,q3,q32=bh(r2),bh(r3),bh(r32)
    return {'n_tested':len(keep),'S2vS1':int((q2<.05).sum()),'S3vS2':int((q32<.05).sum()),'S3vS1':int((q3<.05).sum())}, pd.Series(q3,index=CT.columns[keep])
rows=[];qs={}
for pc in [0.1,0.5,1.0]:
    for pmin in [0.05,0.10,0.20]:
        if pc!=0.5 and pmin!=0.10: continue
        r,q3=lmm_counts(pc,pmin); rows.append({'pseudocount':pc,'min_prevalence':pmin,**r}); qs[(pc,pmin)]=q3; print(rows[-1],flush=True)
SENS=pd.DataFrame(rows); SENS.to_csv('da_sensitivity.csv',index=False)
# overlap of S3vS1 significant sets across settings
base=set(qs[(0.5,0.10)][qs[(0.5,0.10)]<.05].index)
ov=[{'setting':f'pc={k[0]}, prev>={k[1]}','n_sig':int((v<.05).sum()),'overlap_with_primary':len(base & set(v[v<.05].index))} for k,v in qs.items()]
pd.DataFrame(ov).to_csv('da_sensitivity_overlap.csv',index=False); print(pd.DataFrame(ov))
# ---------- ALDEx2-like: Dirichlet MC -> CLR -> paired Wilcoxon (complete cohort) ----------
keep=np.where(prev.values>=0.10)[0]; names=CT.columns[keep]
ci=np.where(cm)[0]; Cc=CT.values[ci]; tpc=tpv[ci]; pidc=pid[ci]
order={t:[np.where((tpc==t)&(pidc==p_))[0][0] for p_ in complete] for t in tps}
NMC=64; P={k:np.zeros((NMC,len(keep))) for k in ['S1S2','S2S3','S1S3']}; E={k:np.zeros((NMC,len(keep))) for k in P}
for m in range(NMC):
    G=rng.gamma(Cc+0.5); Pr=G/G.sum(1,keepdims=True); L=np.log(Pr); Z=(L-L.mean(1,keepdims=True))[:,keep]
    for (a,b),k in [(('S1','S2'),'S1S2'),(('S2','S3'),'S2S3'),(('S1','S3'),'S1S3')]:
        x=Z[order[a]]; y=Z[order[b]]; dlt=y-x
        for j in range(len(keep)):
            P[k][m,j]=stats.wilcoxon(x[:,j],y[:,j]).pvalue if np.any(dlt[:,j]!=0) else 1.0
        E[k][m]=dlt.mean(0)
al=pd.DataFrame({'feature':names})
for k in P: 
    al[f'ALDEx_{k}_p']=P[k].mean(0); al[f'ALDEx_{k}_q']=bh(al[f'ALDEx_{k}_p']); al[f'ALDEx_{k}_effect']=E[k].mean(0)
al.to_csv('aldex_like.csv',index=False)
print({k:int((al[f'ALDEx_{k}_q']<.05).sum()) for k in P})
# concordance with LMM
s=pd.read_csv('da_species.csv').set_index('feature'); lm=set(s.index[s.LMM_S3vS1_q<.05]); ad=set(al.feature[al.ALDEx_S1S3_q<.05])
print('S3vS1 LMM',len(lm),'ALDEx-like',len(ad),'both',len(lm&ad))
json.dump({'lmm':len(lm),'aldex':len(ad),'both':len(lm&ad),'aldex_counts':{k:int((al[f"ALDEx_{k}_q"]<.05).sum()) for k in P}},open('concordance.json','w'))
