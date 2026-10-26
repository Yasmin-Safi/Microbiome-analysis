exec(open('prep.py').read())
from scipy import stats
from statsmodels.stats.multitest import multipletests
import json
from skbio import DistanceMatrix
from skbio.stats.distance import permdisp
rng=np.random.default_rng(2026)
def bray(P):
    n=P.shape[0]; D=np.zeros((n,n))
    for i in range(n): D[i]=np.abs(P[i]-P).sum(1)/(P[i]+P).sum(1)
    return D
def jacc(M):
    B=(M>0).astype(float); inter=B@B.T; s=B.sum(1); un=s[:,None]+s[None,:]-inter
    return 1-inter/un
dists={'BC_raw':bray(R.values),'Jac_raw':jacc(CT.values),'BC_rar':np.load('DB.npy'),'Jac_rar':np.load('DJ.npy')}
idx=meta.index
def gower(D):
    n=D.shape[0]; J=np.eye(n)-1/n; return -0.5*J@(D**2)@J
def hat(X):
    X=X-X.mean(0); 
    Q,_=np.linalg.qr(X); 
    # rank-safe
    u,s,vt=np.linalg.svd(X,full_matrices=False); r=(s>1e-8).sum(); U_=u[:,:r]; return U_@U_.T
def dummies(v):
    return pd.get_dummies(pd.Series(v),drop_first=True).values.astype(float)
def permanova_blocked(D,tp,pid,nperm=9999,restrict=True):
    n=len(tp); G=gower(D); trG=np.trace(G)
    Hp=hat(dummies(pid)); trHp=np.trace(Hp@G)
    def stat(t):
        Ht=hat(np.column_stack([dummies(pid),dummies(t)]))
        ss_t=np.trace(Ht@G)-trHp; ss_r=trG-np.trace(Ht@G)
        dft=len(set(t))-1; dfr=n-np.linalg.matrix_rank(np.column_stack([dummies(pid),dummies(t)]))-1
        return (ss_t/dft)/(ss_r/dfr), ss_t/trG
    F0,R20=stat(tp); 
    pid=np.asarray(pid); tp=np.asarray(tp); groups=[np.where(pid==u)[0] for u in np.unique(pid)]
    ge=1
    # speed: precompute Hp; recompute Ht each perm (n<=167: fine)
    for _ in range(nperm):
        t=tp.copy()
        for g in groups: t[g]=rng.permutation(tp[g])
        if stat(t)[0]>=F0-1e-12: ge+=1
    return F0,R20,ge/(nperm+1)
def permanova_simple(D,tp,nperm=999):
    tp=np.asarray(tp); n=len(tp); G=gower(D); trG=np.trace(G)
    def stat(t):
        Ht=hat(dummies(t)); ss_t=np.trace(Ht@G); ss_r=trG-ss_t; dft=len(set(t))-1
        return (ss_t/dft)/(ss_r/(n-dft-1)), ss_t/trG
    F0,R20=stat(tp); ge=1
    for _ in range(nperm):
        if stat(rng.permutation(tp))[0]>=F0-1e-12: ge+=1
    return F0,R20,ge/(nperm+1)
res=[]
cm=meta.complete.values
NP=2999
for name,D in dists.items():
    # (a) original-style: unrestricted, all 167
    F,R2,pv=permanova_simple(D,meta.tp.values,999); res.append(dict(dist=name,design='All samples (n=167), unrestricted permutations (original analysis)',F=F,R2=R2,p=pv))
    # (b) all 167 with patient as block, permutation within patient
    F,R2,pv=permanova_blocked(D,meta.tp.values,meta.pid.values,NP); res.append(dict(dist=name,design='All samples (n=167), patient-restricted permutations, pid + timepoint',F=F,R2=R2,p=pv))
    # (c) complete cohort n=138
    ii=np.where(cm)[0]; Dc=D[np.ix_(ii,ii)]
    F,R2,pv=permanova_blocked(Dc,meta.tp.values[ii],meta.pid.values[ii],NP); res.append(dict(dist=name,design='Complete cohort (n=46x3), patient-restricted permutations, pid + timepoint',F=F,R2=R2,p=pv))
    print(name,'done',flush=True)
PR=pd.DataFrame(res); PR.to_csv('permanova.csv',index=False); print(PR.round(4).to_string())
# pairwise (complete cohort, restricted)
pw=[]
for name in ['BC_rar','Jac_rar','BC_raw','Jac_raw']:
    D=dists[name]
    for a,b in [('S1','S2'),('S2','S3'),('S1','S3')]:
        ii=np.where(cm&meta.tp.isin([a,b]).values)[0]; Dc=D[np.ix_(ii,ii)]
        F,R2,pv=permanova_blocked(Dc,meta.tp.values[ii],meta.pid.values[ii],NP)
        pw.append(dict(dist=name,pair=f'{a} vs {b}',F=F,R2=R2,p=pv))
PW=pd.DataFrame(pw)
for name in PW.dist.unique():
    m=PW.dist==name; PW.loc[m,'q']=multipletests(PW.loc[m,'p'],method='fdr_bh')[1]
PW.to_csv('permanova_pairwise.csv',index=False); print(PW.round(4).to_string())
# PERMDISP both metrics
pdisp=[]
for name,D in dists.items():
    for lab,ii in [('All (n=167)',np.arange(len(idx))),('Complete (n=138)',np.where(cm)[0])]:
        dm=DistanceMatrix(D[np.ix_(ii,ii)],ids=[str(i) for i in range(len(ii))])
        r=permdisp(dm,meta.tp.values[ii],test='centroid',permutations=9999)
        # group mean distances to centroid via PCoA
        pdisp.append(dict(dist=name,set=lab,F=float(r['test statistic']),p=float(r['p-value'])))
PD=pd.DataFrame(pdisp); PD.to_csv('permdisp.csv',index=False); print(PD.round(4).to_string())
# mean distance to centroid per tp (PCoA-based) for BC_raw & Jaccard_raw
def disp_by_tp(D,ii):
    from skbio.stats.ordination import pcoa
    dm=DistanceMatrix(D[np.ix_(ii,ii)],ids=[str(i) for i in range(len(ii))])
    o=pcoa(dm,number_of_dimensions=0).samples.values  # all positive axes
    # include negative eigenvalues ignored (centroid method on positive axes)
    t=meta.tp.values[ii]; res={}
    for g in ['S1','S2','S3']:
        X=o[t==g]; c=np.median(X,axis=0) if False else X.mean(0); res[g]=float(np.sqrt(((X-c)**2).sum(1)).mean())
    return res
dd={}
for name in ['BC_raw','Jac_raw','BC_rar','Jac_rar']:
    dd[name]={'all':disp_by_tp(dists[name],np.arange(len(idx))),'complete':disp_by_tp(dists[name],np.where(cm)[0])}
json.dump(dd,open('disp_by_tp.json','w'),indent=1); print(json.dumps(dd,indent=1))
