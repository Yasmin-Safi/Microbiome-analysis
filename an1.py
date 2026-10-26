exec(open('prep.py').read())
from scipy import stats
from statsmodels.stats.multitest import multipletests
import statsmodels.formula.api as smf, json, itertools
rng=np.random.default_rng(42)
out={}
tps=['S1','S2','S3']
# ---------- rarefaction ----------
RD=7500
def rarefy(counts,d,rng):
    arr=np.zeros_like(counts)
    for i in range(counts.shape[0]):
        c=counts[i].astype(int)
        arr[i]=rng.multivariate_hypergeometric(c,d)
    return arr
NIT=20
Cv=CT.values
def alpha_from_counts(M):
    tot=M.sum(1,keepdims=True); pr=M/tot
    with np.errstate(divide='ignore',invalid='ignore'):
        sh=-np.nansum(np.where(pr>0,pr*np.log(pr),0),1)
    sob=(M>0).sum(1); si=1-(pr**2).sum(1)
    n1=(M==1).sum(1); n2=(M==2).sum(1)
    ch=np.where(n2>0,sob+n1**2/(2*np.maximum(n2,1)),sob+n1*(n1-1)/2)
    pe=sh/np.log(sob)
    return pd.DataFrame({'Shannon':sh,'Simpson':si,'Observed':sob,'Chao1':ch,'Evenness':pe})
def bray(M):
    P=M/M.sum(1,keepdims=True)
    n=P.shape[0]; D=np.zeros((n,n))
    for i in range(n): D[i]=np.abs(P[i]-P).sum(1)/(P[i]+P).sum(1)
    return D
def jac(M):
    B=(M>0).astype(float); inter=B@B.T; s=B.sum(1); un=s[:,None]+s[None,:]-inter
    return 1-inter/un
alphas=[];DB=0;DJ=0;RAR=[]
for it in range(NIT):
    M=rarefy(Cv,RD,rng); RAR.append(M)
    alphas.append(alpha_from_counts(M).values)
    DB=DB+bray(M)/NIT; DJ=DJ+jac(M)/NIT
A_r=pd.DataFrame(np.mean(alphas,0),index=CT.index,columns=['Shannon','Simpson','Observed','Chao1','Evenness'])
A_raw=alpha_from_counts(Cv); A_raw.index=CT.index
# raw Chao1 per Rmd; use original file for raw alpha to match manuscript
orig=pd.read_excel(U+'Microbiome_Analysis_Results.xlsx',sheet_name='Alpha_Diversity').set_index('Sample')
A_raw=orig[['Shannon','Simpson','Observed_Species','Chao1','Evenness']].rename(columns={'Observed_Species':'Observed'}).loc[CT.index]
np.save('DB.npy',DB); np.save('DJ.npy',DJ)
A_r.to_csv('alpha_rarefied.csv'); 
import pickle; pickle.dump(RAR[0],open('rar0.pkl','wb'))
# ---------- depth ----------
d={}
d['depth_by_tp']={t:[float(depth[meta.tp==t].median()),float(depth[meta.tp==t].quantile(.25)),float(depth[meta.tp==t].quantile(.75)),float(depth[meta.tp==t].mean())] for t in tps}
d['depth_KW_p']=float(stats.kruskal(*[depth[meta.tp==t] for t in tps]).pvalue)
wide=depth.to_frame('d').join(meta); W=wide.pivot(index='pid',columns='tp',values='d')
W3=W.loc[complete]
d['depth_paired_S1S3_p']=float(stats.wilcoxon(W3.S1,W3.S3).pvalue); d['depth_paired_S1S2_p']=float(stats.wilcoxon(W3.S1,W3.S2).pvalue); d['depth_paired_S2S3_p']=float(stats.wilcoxon(W3.S2,W3.S3).pvalue)
d['spearman_depth_vs_obs_raw']={t:[float(x) for x in stats.spearmanr(depth[meta.tp==t],A_raw.Observed[meta.tp==t])] for t in tps}
d['spearman_depth_vs_obs_all']=[float(x) for x in stats.spearmanr(depth,A_raw.Observed)]
d['spearman_depth_vs_shannon_all']=[float(x) for x in stats.spearmanr(depth,A_raw.Shannon)]
d['spearman_depth_vs_obs_rar']=[float(x) for x in stats.spearmanr(depth,A_r.Observed)]
d['n_samples_below_rd']=int((depth<RD).sum()); d['min_depth']=float(depth.min())
out['depth']=d
# ---------- alpha: effect sizes, tests ----------
def rb_unpaired(a,b):
    u=stats.mannwhitneyu(a,b,alternative='two-sided').statistic
    return 2*u/(len(a)*len(b))-1  # P(a>b)-P(a<b)
def rb_paired(a,b):
    dlt=np.array(a)-np.array(b); dlt=dlt[dlt!=0]
    r=stats.rankdata(np.abs(dlt)); Wp=r[dlt>0].sum(); Wm=r[dlt<0].sum()
    return (Wp-Wm)/(Wp+Wm)
def alpha_tests(A,label):
    rows=[]
    for m in ['Shannon','Simpson','Observed','Chao1','Evenness']:
        s=A[m]; g=[s[meta.tp==t] for t in tps]
        kw=stats.kruskal(*g)
        row={'set':label,'metric':m,'KW_H':kw.statistic,'KW_p':kw.pvalue}
        for a,b in [('S1','S2'),('S2','S3'),('S1','S3')]:
            x=s[meta.tp==a]; y=s[meta.tp==b]
            row[f'MW_{a}{b}_p']=stats.mannwhitneyu(x,y).pvalue
            row[f'MW_{a}{b}_rb']=rb_unpaired(y,x)  # later vs earlier; negative = decline
            xs=s[(meta.tp==a)&meta.complete].groupby(meta.pid).first().loc[complete]; ys=s[(meta.tp==b)&meta.complete].groupby(meta.pid).first().loc[complete]
            row[f'W_{a}{b}_p']=stats.wilcoxon(xs,ys).pvalue
            row[f'W_{a}{b}_rb']=rb_paired(ys,xs)
            row[f'W_{a}{b}_medchange']=float(np.median(ys.values-xs.values))
        rows.append(row)
    df=pd.DataFrame(rows)
    # BH across 5 metrics within each test family
    for c in [c for c in df.columns if c.endswith('_p')]:
        df[c.replace('_p','_q')]=multipletests(df[c],method='fdr_bh')[1]
    return df
at_raw=alpha_tests(A_raw,'raw'); at_rar=alpha_tests(A_r,'rarefied')
at_raw.to_csv('alpha_tests_raw.csv',index=False); at_rar.to_csv('alpha_tests_rar.csv',index=False)
# ---------- LMM alpha (all samples, random intercept) ----------
def lmm(A,label,covar=False):
    rows=[]
    for m in ['Shannon','Observed','Chao1']:
        df=pd.DataFrame({'y':A[m].values,'tp':meta.tp.values,'pid':meta.pid.values,'ld':np.log10(depth.values)})
        f='y ~ C(tp)'+(' + ld' if covar else '')
        r=smf.mixedlm(f,df,groups=df['pid']).fit(reml=True)
        for t in ['S2','S3']:
            k=f'C(tp)[T.{t}]'
            rows.append({'set':label,'metric':m,'contrast':f'{t} vs S1','est':r.params[k],'se':r.bse[k],'p':r.pvalues[k],'lo':r.conf_int().loc[k,0],'hi':r.conf_int().loc[k,1]})
        # S3 vs S2 via re-leveling
        df['tp2']=pd.Categorical(df.tp,categories=['S2','S1','S3'])
        r2=smf.mixedlm('y ~ C(tp2)'+(' + ld' if covar else ''),df,groups=df['pid']).fit(reml=True)
        k='C(tp2)[T.S3]'
        rows.append({'set':label,'metric':m,'contrast':'S3 vs S2','est':r2.params[k],'se':r2.bse[k],'p':r2.pvalues[k],'lo':r2.conf_int().loc[k,0],'hi':r2.conf_int().loc[k,1]})
    return pd.DataFrame(rows)
L=pd.concat([lmm(A_raw,'raw'),lmm(A_raw,'raw+logdepth',True),lmm(A_r,'rarefied')])
L.to_csv('alpha_lmm.csv',index=False)
out['alpha_raw_vs_rar']=at_rar.round(4).to_dict('records')
json.dump(out,open('out1.json','w'),indent=1,default=float)
print(json.dumps(d,indent=1,default=float))
pd.set_option('display.width',250,'display.max_columns',50)
print(at_raw[['metric','KW_p','MW_S1S2_p','MW_S1S3_p','W_S1S2_p','W_S2S3_p','W_S1S3_p','W_S1S3_rb','W_S1S3_q']].round(4))
print(at_rar[['metric','KW_p','MW_S1S2_p','MW_S2S3_p','MW_S1S3_p','W_S1S2_p','W_S2S3_p','W_S1S3_p','W_S1S2_rb','W_S2S3_rb','W_S1S3_rb','W_S1S2_q','W_S2S3_q','W_S1S3_q']].round(4))
print(L.round(4).to_string())
