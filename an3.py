exec(open('prep.py').read())
from scipy import stats
from statsmodels.stats.multitest import multipletests
import statsmodels.formula.api as smf, json, warnings
warnings.filterwarnings('ignore')
tps=['S1','S2','S3']; tpv=meta.tp.values; pid=meta.pid.values
def bh(p): 
    p=np.array(p,float); q=np.full(len(p),np.nan); m=~np.isnan(p); q[m]=multipletests(p[m],method='fdr_bh')[1]; return q
def clr(M,pc=0.5):
    X=np.log(M+pc); return X-X.mean(1,keepdims=True)
def analyse(Cnt,label,prev_min=0.10):
    names=list(Cnt.columns); Rl=Cnt.div(Cnt.sum(axis=1),axis=0)
    prev=(Cnt>0).mean(axis=0); keep=prev>=prev_min
    # KW all features
    kw=[]
    for s in names:
        g=[Rl[s].values[tpv==t] for t in tps]
        try: kw.append(stats.kruskal(*g).pvalue)
        except: kw.append(np.nan)
    df=pd.DataFrame({'feature':names,'prev':prev.values,'KW_p':kw,'mean_S1':[Rl[s].values[tpv=='S1'].mean() for s in names],'mean_S2':[Rl[s].values[tpv=='S2'].mean() for s in names],'mean_S3':[Rl[s].values[tpv=='S3'].mean() for s in names]})
    df['KW_q_all']=bh(df.KW_p)
    df['tested_filtered']=keep.values
    df['KW_q_filt']=np.nan; df.loc[df.tested_filtered,'KW_q_filt']=bh(df.loc[df.tested_filtered,'KW_p'])
    # filtered: paired (complete) & LMM
    F=[s for s,k in zip(names,keep) if k]
    X=pd.DataFrame(clr(Cnt[F].values),index=Cnt.index,columns=F)
    ld=np.log10(Cnt.sum(axis=1).values)
    cm=meta.complete.values
    fr=[];wl=[]
    Xc=X[cm]; mc=meta[cm]
    wide={t:Xc[mc.tp==t].set_index(mc.pid[mc.tp==t]).loc[complete] for t in tps}
    lm=[]
    for s in F:
        a,b,c=[wide[t][s].values for t in tps]
        try: fr_p=stats.friedmanchisquare(a,b,c).pvalue
        except: fr_p=np.nan
        row={'feature':s,'Friedman_p':fr_p}
        for (n,x,y) in [('S1S2',a,b),('S2S3',b,c),('S1S3',a,c)]:
            try: row['W_'+n+'_p']=stats.wilcoxon(x,y).pvalue
            except: row['W_'+n+'_p']=np.nan
            row['clrdiff_'+n]=float(np.mean(y-x))
        # LMM
        d=pd.DataFrame({'y':X[s].values,'tp':pd.Categorical(tpv,categories=['S1','S2','S3']),'pid':pid,'ld':ld})
        try:
            r=smf.mixedlm('y ~ C(tp)',d,groups=d['pid']).fit(reml=True)
            row['LMM_S2vS1_est']=r.params['C(tp)[T.S2]']; row['LMM_S2vS1_p']=r.pvalues['C(tp)[T.S2]']
            row['LMM_S3vS1_est']=r.params['C(tp)[T.S3]']; row['LMM_S3vS1_p']=r.pvalues['C(tp)[T.S3]']
            # params order: Intercept, S2, S3, Group Var
            tt=r.t_test(np.array([0,-1,1.]).reshape(1,-1))
            row['LMM_S3vS2_est']=float(tt.effect[0]); row['LMM_S3vS2_p']=float(tt.pvalue)
            # depth-adjusted
            r2=smf.mixedlm('y ~ C(tp) + ld',d,groups=d['pid']).fit(reml=True)
            row['LMMd_S2vS1_p']=r2.pvalues['C(tp)[T.S2]']; row['LMMd_S3vS1_p']=r2.pvalues['C(tp)[T.S3]']
            tt2=r2.t_test(np.array([0,-1,1,0.]).reshape(1,-1)); row['LMMd_S3vS2_p']=float(tt2.pvalue)
        except Exception as e:
            pass
        lm.append(row)
    LM=pd.DataFrame(lm)
    for c in [c for c in LM.columns if c.endswith('_p')]: LM[c.replace('_p','_q')]=bh(LM[c])
    out=df.merge(LM,on='feature',how='left')
    return out
sp_res=analyse(CT,'species'); sp_res.to_csv('da_species.csv',index=False)
# genus
gcnt=CT.copy(); gcnt.columns=[g.replace('__','') for g in cnt['Genus']]
gcnt=gcnt.T.groupby(level=0).sum().T
g_res=analyse(gcnt,'genus'); g_res.to_csv('da_genus.csv',index=False)
def summ(r,label):
    s={'n_features':len(r),'n_tested_filtered':int(r.tested_filtered.sum()),
       'KW_p05':int((r.KW_p<0.05).sum()),'KW_q05_all':int((r.KW_q_all<0.05).sum()),'KW_q05_filt':int((r.KW_q_filt<0.05).sum()),'KW_q10_filt':int((r.KW_q_filt<0.10).sum()),
       'Friedman_q05':int((r.Friedman_q<0.05).sum()),'LMM_S2vS1_q05':int((r.LMM_S2vS1_q<0.05).sum()),'LMM_S3vS1_q05':int((r.LMM_S3vS1_q<0.05).sum()),'LMM_S3vS2_q05':int((r.LMM_S3vS2_q<0.05).sum()),
       'LMMd_S2vS1_q05':int((r.LMMd_S2vS1_q<0.05).sum()),'LMMd_S3vS1_q05':int((r.LMMd_S3vS1_q<0.05).sum()),'LMMd_S3vS2_q05':int((r.LMMd_S3vS2_q<0.05).sum()),
       'W_S1S2_q05':int((r.W_S1S2_q<0.05).sum()),'W_S2S3_q05':int((r.W_S2S3_q<0.05).sum()),'W_S1S3_q05':int((r.W_S1S3_q<0.05).sum())}
    print(label,json.dumps(s)); return s
S={'species':summ(sp_res,'species'),'genus':summ(g_res,'genus')}
json.dump(S,open('da_summary.json','w'),indent=1)
