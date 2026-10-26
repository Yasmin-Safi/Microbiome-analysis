exec(open('prep.py').read())
from scipy import stats
from statsmodels.stats.multitest import multipletests
import statsmodels.formula.api as smf, statsmodels.api as sm, json, warnings
warnings.filterwarnings('ignore')
rng=np.random.default_rng(11)
tps=['S1','S2','S3']; tpv=meta.tp.values; pid=meta.pid.values; cm=meta.complete.values
OUT={}
def bh(p):
    p=np.array(p,float); q=np.full(len(p),np.nan); m=~np.isnan(p); q[m]=multipletests(p[m],method='fdr_bh')[1]; return q
# ---------- A. within-patient dissimilarity ----------
DB=np.load('DB.npy'); DJ=np.load('DJ.npy')
from scipy.spatial.distance import pdist,squareform
DBraw=squareform(pdist(R.values,'braycurtis')); DJraw=squareform(pdist(R.values>0,'jaccard'))
ix={(p_,t):i for i,(p_,t) in enumerate(zip(pid,tpv))}
rows=[]
for name,D in [('Bray-Curtis (raw)',DBraw),('Jaccard (raw)',DJraw),('Bray-Curtis (rarefied)',DB),('Jaccard (rarefied)',DJ)]:
    d12=np.array([D[ix[(p_,'S1')],ix[(p_,'S2')]] for p_ in complete]); d23=np.array([D[ix[(p_,'S2')],ix[(p_,'S3')]] for p_ in complete]); d13=np.array([D[ix[(p_,'S1')],ix[(p_,'S3')]] for p_ in complete])
    # between-patient reference at S1 (median distance to other patients at S1)
    s1=[ix[(p_,'S1')] for p_ in complete]; btw=np.median(D[np.ix_(s1,s1)][np.triu_indices(len(s1),1)])
    rows.append({'distance':name,'median_S1_S2':np.median(d12),'median_S2_S3':np.median(d23),'median_S1_S3':np.median(d13),'median_between_patients_S1':btw,
     'p_S1S2_vs_S2S3':stats.wilcoxon(d12,d23).pvalue,'p_S1S2_vs_between':stats.mannwhitneyu(d12,D[np.ix_(s1,s1)][np.triu_indices(len(s1),1)]).pvalue,
     'frac_patients_S1S2_gt_S2S3':float((d12>d23).mean())})
WP=pd.DataFrame(rows); WP.to_csv('within_patient_dist.csv',index=False); print(WP.round(4).to_string())
# store per-patient for figure
pd.DataFrame({'pid':complete,'BC_S1S2':[DB[ix[(p_,'S1')],ix[(p_,'S2')]] for p_ in complete],'BC_S2S3':[DB[ix[(p_,'S2')],ix[(p_,'S3')]] for p_ in complete],'JAC_S1S2':[DJ[ix[(p_,'S1')],ix[(p_,'S2')]] for p_ in complete],'JAC_S2S3':[DJ[ix[(p_,'S2')],ix[(p_,'S3')]] for p_ in complete]}).to_csv('within_patient_pairs.csv',index=False)
# ---------- B. inverse-probability weighting ----------
pp=p.copy(); pp['pid']=pp['Case ID']; pp['complete']=pp.pid.isin(complete).astype(int)
A=pd.read_csv('alpha_rarefied.csv',index_col=0)
base=A.loc[meta.index[tpv=='S1']].copy(); base.index=meta.pid[tpv=='S1'].values
pp=pp.set_index('pid').join(base[['Observed','Shannon']].rename(columns={'Observed':'b_obs','Shannon':'b_sh'}))
X=pd.DataFrame({'age':pp.Age,'male':(pp.Gender.str.strip().str.lower()=='male').astype(int),'bmi':pp.BMI,'cN2':(pp.Cn=='N2').astype(int),'cT4':(pp.Ct=='T4').astype(int),'smoker':(pp.Smoking.astype(str).str.strip().str.lower()=='yes').astype(int),'comor':(pp['COMOR(Y/N)'].astype(str).str.strip().str.lower()=='yes').astype(int),'b_sh':pp.b_sh})
Xz=(X-X.mean())/X.std(); Xz=sm.add_constant(Xz)
lg=sm.Logit(pp.complete.values,Xz.values).fit_disp=False if False else sm.Logit(pp.complete.values,Xz.values).fit(disp=0,maxiter=200,method='bfgs')
ps=lg.predict(Xz.values); w=1/ps
from sklearn.metrics import roc_auc_score
auc=roc_auc_score(pp.complete.values,ps)
print('IPW AUC',round(auc,3),'weights range',w[pp.complete.values==1].min().round(2),w[pp.complete.values==1].max().round(2))
def wmean(x,wt): return np.sum(x*wt)/np.sum(wt)
res=[]
for m in ['Observed','Shannon']:
    W=pd.DataFrame({t:A[m][tpv==t].values for t in tps},index=[meta.pid[tpv==t].values for t in tps][0] if False else None) if False else None
    s=A[m]; Wd={t:pd.Series(s[tpv==t].values,index=meta.pid[tpv==t].values) for t in tps}
    for a,b in [('S1','S2'),('S2','S3'),('S1','S3')]:
        dch=(Wd[b]-Wd[a]).loc[complete]; wt=pd.Series(w,index=pp.index).loc[complete]
        est=wmean(dch.values,wt.values); un=dch.mean()
        bs=[]
        for _ in range(2000):
            k=rng.integers(0,len(dch),len(dch)); bs.append(wmean(dch.values[k],wt.values[k]))
        res.append({'metric':m,'contrast':f'{a}->{b}','unweighted_mean_change':un,'IPW_mean_change':est,'IPW_CI_lo':np.percentile(bs,2.5),'IPW_CI_hi':np.percentile(bs,97.5)})
IP=pd.DataFrame(res); IP.to_csv('ipw.csv',index=False); print(IP.round(3).to_string())
OUT['ipw_auc']=auc
# ---------- C. covariate-adjusted alpha mixed models ----------
cov=pp[['Age','BMI']].copy(); cov['male']=X.male
d=pd.DataFrame({'y':A['Observed'].values,'sh':A['Shannon'].values,'tp':pd.Categorical(tpv,categories=tps),'pid':pid})
d=d.join(cov,on='pid'); d['Age_c']=d.Age-d.Age.mean(); d['BMI_c']=d.BMI-d.BMI.mean()
rows=[]
for resp in ['y','sh']:
    for lab,f in [('timepoint only',f'{resp} ~ C(tp)'),('+ age, sex, BMI',f'{resp} ~ C(tp)+Age_c+male+BMI_c')]:
        r=smf.mixedlm(f,d,groups=d['pid']).fit(reml=True)
        for t in ['S2','S3']:
            k=f'C(tp)[T.{t}]'; rows.append({'metric':'Observed (rarefied)' if resp=='y' else 'Shannon (rarefied)','model':lab,'contrast':f'{t} vs S1','est':r.params[k],'lo':r.conf_int().loc[k,0],'hi':r.conf_int().loc[k,1],'p':r.pvalues[k]})
CV=pd.DataFrame(rows); CV.to_csv('alpha_covariate.csv',index=False); print(CV.round(4).to_string())
# ---------- D. rarefaction depth sensitivity ----------
def rarefy_mat(M,dp):
    keep=M.sum(1)>=dp; out=np.zeros((keep.sum(),M.shape[1])); 
    for j,i in enumerate(np.where(keep)[0]): out[j]=rng.multivariate_hypergeometric(M[i].astype(int),dp)
    return out,keep
rows=[]
for dp in [3000,5000,7500,10000,15000]:
    ob=[];sh=[];kp=None
    for it in range(10):
        M,keep=rarefy_mat(CT.values,dp); kp=keep
        pr=M/M.sum(1,keepdims=True); ob.append((M>0).sum(1)); sh.append(-np.sum(np.where(pr>0,pr*np.log(np.where(pr>0,pr,1)),0),1))
    ob=np.mean(ob,0); sh=np.mean(sh,0)
    dd=pd.DataFrame({'y':ob,'sh':sh,'tp':pd.Categorical(tpv[kp],categories=tps),'pid':pid[kp]})
    for resp,lab in [('y','Observed richness'),('sh','Shannon')]:
        r=smf.mixedlm(f'{resp} ~ C(tp)',dd,groups=dd['pid']).fit(reml=True)
        e2=r.params['C(tp)[T.S2]']; e3=r.params['C(tp)[T.S3]']
        t=r.t_test(np.array([0,-1,1.]).reshape(1,-1))
        rows.append({'depth':dp,'metric':lab,'n_samples':int(kp.sum()),'S1':int((tpv[kp]=='S1').sum()),'S2':int((tpv[kp]=='S2').sum()),'S3':int((tpv[kp]=='S3').sum()),'S2-S1':e2,'S2-S1_p':r.pvalues['C(tp)[T.S2]'],'S3-S1':e3,'S3-S1_p':r.pvalues['C(tp)[T.S3]'],'S3-S2':float(t.effect[0]),'S3-S2_p':float(t.pvalue)})
RS=pd.DataFrame(rows); RS.to_csv('rarefaction_depth_sens.csv',index=False); print(RS.round(4).to_string())
