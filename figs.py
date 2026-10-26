exec(open('prep.py').read())
import json, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from scipy import stats
from scipy.spatial.distance import pdist,squareform
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,'axes.spines.right':False,'axes.titlesize':10,'axes.titleweight':'bold','figure.dpi':100})
COL={'S1':'#3B7DD8','S2':'#E8742C','S3':'#4C9A5B'}; LAB={'S1':'S1 baseline','S2':'S2 post-RT','S3':'S3 post-CT'}
tps=['S1','S2','S3']; tpv=meta.tp.values; pid=meta.pid.values; cm=meta.complete.values
def pcoa(D,k=2):
    n=len(D); J=np.eye(n)-1/n; G=-0.5*J@(D**2)@J; w,v=np.linalg.eigh(G); o=np.argsort(w)[::-1]; w=w[o]; v=v[:,o]
    pos=w>0; return v[:,:k]*np.sqrt(w[:k]), w[:k]/w[pos].sum()
def ell(ax,X,c):
    cov=np.cov(X.T); vals,vecs=np.linalg.eigh(cov); a=np.degrees(np.arctan2(*vecs[:,1][::-1])); 
    ax.add_patch(Ellipse(X.mean(0),2*2*np.sqrt(vals[1]*1.0),2*2*np.sqrt(vals[0]*1.0),angle=a,fc=c,ec=c,alpha=0.12,lw=1))
def panel(l,ax): ax.text(-0.12,1.06,l,transform=ax.transAxes,fontsize=13,fontweight='bold')
PM=pd.read_csv('permanova.csv'); PW=pd.read_csv('permanova_pairwise.csv'); PD=pd.read_csv('permdisp.csv')
def getp(dist,design_sub):
    r=PM[(PM.dist==dist)&PM.design.str.contains(design_sub,regex=False)].iloc[0]; return r.F,r.R2,r.p
def fmtp(p): return 'p < 0.001' if p<0.001 else f'p = {p:.3f}'
# ================= Revised Figure 3 (beta) =================
DBraw=squareform(pdist(R.values,'braycurtis')); DJraw=squareform(pdist(R.values>0,'jaccard'))
DBr=np.load('DB.npy'); DJr=np.load('DJ.npy')
fig,axs=plt.subplots(2,3,figsize=(14,8.6))
for ax,D,nm,key,l in [(axs[0,0],DBraw,'Bray-Curtis','BC_raw','A'),(axs[0,1],DJraw,'Jaccard','Jac_raw','B')]:
    X,ev=pcoa(D)
    for t in tps:
        m=tpv==t; ax.scatter(X[m,0],X[m,1],s=28,c=COL[t],alpha=0.75,label=LAB[t],edgecolor='white',lw=0.4); ell(ax,X[m],COL[t]); ax.scatter(*X[m].mean(0),marker='*',s=170,c=COL[t],edgecolor='k',zorder=5)
    F,R2,p_=getp(key,'All samples (n=167), patient-restricted')
    F2,R22,p2=getp(key,'Complete cohort')
    dp=PD[(PD.dist==key)&(PD.set=='All (n=167)')].iloc[0]
    ax.set_title(f'PCoA: {nm} dissimilarity'); ax.set_xlabel(f'PCo1 ({100*ev[0]:.1f}%)'); ax.set_ylabel(f'PCo2 ({100*ev[1]:.1f}%)')
    ax.text(0.02,0.02,f'PERMANOVA (within-patient permutations)\nR² = {R2:.3f}, {fmtp(p_)}\nPERMDISP: {fmtp(dp.p)}',transform=ax.transAxes,fontsize=7.5,va='bottom',bbox=dict(fc='white',ec='0.8',alpha=.9)); panel(l,ax)
axs[0,0].legend(frameon=False,fontsize=8,loc='upper right')
# C dispersion (distance to centroid, Jaccard & BC) 
ax=axs[0,2]; pos=0; ticks=[];labs=[]
for nm,D,key in [('Bray-Curtis',DBraw,'BC_raw'),('Jaccard',DJraw,'Jac_raw')]:
    X,_=pcoa(D,k=min(60,len(D)-1)); 
    for t in tps:
        m=tpv==t; c=X[m].mean(0); dd=np.sqrt(((X[m]-c)**2).sum(1)); bp=ax.boxplot(dd,positions=[pos],widths=0.65,patch_artist=True,showfliers=False)
        bp['boxes'][0].set_facecolor(COL[t]); bp['boxes'][0].set_alpha(.8)
        for e in bp['medians']: e.set_color('k')
        pos+=1
    ticks.append(pos-2); labs.append(nm); pos+=1
ax.set_xticks(ticks); ax.set_xticklabels(labs); ax.set_ylabel('Distance to timepoint centroid'); ax.set_ylim(0.2,0.85)
dpb=PD[(PD.dist=='BC_raw')&(PD.set=='All (n=167)')].p.iloc[0]; dpj=PD[(PD.dist=='Jac_raw')&(PD.set=='All (n=167)')].p.iloc[0]
ax.set_title('Multivariate dispersion (PERMDISP)'); ax.text(0.5,0.97,f'Bray-Curtis {fmtp(dpb)}   Jaccard {fmtp(dpj)}',transform=ax.transAxes,ha='center',va='top',fontsize=8); panel('C',ax)
# D variance partition (complete cohort, rarefied)
def gower(D): n=len(D); J=np.eye(n)-1/n; return -0.5*J@(D**2)@J
ii=np.where(cm)[0]; res={}
for nm,D in [('Bray-Curtis',DBr),('Jaccard',DJr)]:
    Dc=D[np.ix_(ii,ii)]; G=gower(Dc); X=pd.get_dummies(pid[ii],drop_first=True).values.astype(float); X-=X.mean(0)
    Q,_=np.linalg.qr(X); rp=np.trace(Q@Q.T@G)/np.trace(G)
    T=pd.get_dummies(tpv[ii],drop_first=True).values.astype(float); XT=np.column_stack([X,T-T.mean(0)]); Q2,_=np.linalg.qr(XT); rfull=np.trace(Q2@Q2.T@G)/np.trace(G)
    res[nm]=(rp,rfull-rp,1-rfull)
ax=axs[1,0]; x=np.arange(2); bot=np.zeros(2)
for k,(lab,c) in enumerate([('Patient identity','#7F7F7F'),('Timepoint','#C0392B'),('Residual','#D9D9D9')]):
    v=np.array([res[n][k] for n in res]); ax.bar(x,v,bottom=bot,color=c,label=lab,width=0.55)
    for xi,vi,bi in zip(x,v,bot): ax.text(xi,bi+vi/2,f'{100*vi:.1f}%',ha='center',va='center',fontsize=8,color='white' if k<2 else 'k')
    bot+=v
ax.set_xticks(x); ax.set_xticklabels(list(res)); ax.set_ylabel('Share of total dissimilarity'); ax.legend(frameon=False,fontsize=8,loc='upper center',bbox_to_anchor=(.5,-0.12),ncol=3); ax.set_title('Variance partition (46 patients, rarefied)'); panel('D',ax)
json.dump({k:list(v) for k,v in res.items()},open('varpart.json','w'))
# E within-patient distances
wp=pd.read_csv('within_patient_pairs.csv'); WPd=pd.read_csv('within_patient_dist.csv').set_index('distance')
ax=axs[1,1]
for k,(a,b,nm,key) in enumerate([('BC_S1S2','BC_S2S3','Bray-Curtis','Bray-Curtis (rarefied)'),('JAC_S1S2','JAC_S2S3','Jaccard','Jaccard (rarefied)')]):
    x0=k*3
    for i in range(len(wp)): ax.plot([x0,x0+1],[wp[a][i],wp[b][i]],c='0.75',lw=0.6,zorder=1)
    ax.scatter(np.full(len(wp),x0),wp[a],s=14,c=COL['S2'],zorder=2); ax.scatter(np.full(len(wp),x0+1),wp[b],s=14,c=COL['S3'],zorder=2)
    ax.plot([x0,x0+1],[wp[a].median(),wp[b].median()],c='k',lw=2.2,zorder=3)
    ax.text(x0+0.5,1.03,fmtp(WPd.loc[key,'p_S1S2_vs_S2S3']),ha='center',fontsize=8)
ax.set_xticks([0,1,3,4]); ax.set_xticklabels(['S1→S2','S2→S3','S1→S2','S2→S3']); ax.set_ylim(0.2,1.1); ax.set_ylabel('Within-patient dissimilarity (rarefied)'); ax.set_title('Change within patients, by treatment phase')
ax.text(0.5,-0.2,'Bray-Curtis                                      Jaccard',transform=ax.transAxes,ha='center',fontsize=8.5); panel('E',ax)
# F pairwise R2
ax=axs[1,2]; sub=PW[PW.dist.isin(['BC_rar','Jac_rar'])]; pairs=['S1 vs S2','S2 vs S3','S1 vs S3']; w=0.36
for k,(d_,c,nm) in enumerate([('BC_rar','#5B6C8F','Bray-Curtis'),('Jac_rar','#B8860B','Jaccard')]):
    s=sub[sub.dist==d_].set_index('pair').loc[pairs]
    ax.bar(np.arange(3)+(k-.5)*w,s.R2*100,w,color=c,label=nm)
    for i,(r2,q) in enumerate(zip(s.R2,s.q)): ax.text(i+(k-.5)*w,r2*100+0.03,('q<0.001' if q<0.001 else f'q={q:.3f}'),ha='center',fontsize=7,rotation=90,va='bottom')
ax.set_xticks(range(3)); ax.set_xticklabels(pairs); ax.set_ylabel('PERMANOVA R² (%)'); ax.set_ylim(0,3.6); ax.legend(frameon=False,fontsize=8); ax.set_title('Pairwise contrasts (restricted, BH q)'); panel('F',ax)
fig.suptitle('Figure 3 (revised). Beta diversity across treatment timepoints, accounting for repeated measures',fontsize=12,fontweight='bold',y=0.995)
fig.tight_layout(rect=[0,0.02,1,0.97]); fig.savefig('figs/Figure3_revised_beta_diversity.png',dpi=250); plt.close(fig); print('fig3')
# ================= Fig S1 depth =================
Ar=pd.read_csv('alpha_rarefied.csv',index_col=0)
orig=pd.read_excel(U+'Microbiome_Analysis_Results.xlsx',sheet_name='Alpha_Diversity').set_index('Sample').loc[CT.index]
fig,axs=plt.subplots(1,4,figsize=(17,4.4))
ax=axs[0]; bp=ax.boxplot([depth[tpv==t] for t in tps],patch_artist=True,showfliers=True,flierprops=dict(markersize=3)); 
for b_,t in zip(bp['boxes'],tps): b_.set_facecolor(COL[t]); b_.set_alpha(.8)
for e in bp['medians']: e.set_color('k')
ax.set_yscale('log'); ax.set_xticklabels(tps); ax.set_ylabel('Reads per sample (log scale)'); ax.axhline(7500,ls='--',c='0.4',lw=1); ax.text(3.45,7800,'7,500',fontsize=7,ha='right'); o1=json.load(open('out1.json'))['depth']; ax.set_title('Sequencing depth by timepoint'); ax.text(0.03,0.04,'Kruskal-Wallis p < 0.001',transform=ax.transAxes,fontsize=8); panel('A',ax)
ax=axs[1]
for t in tps: m=tpv==t; ax.scatter(depth[m],orig.Observed_Species[m],s=14,c=COL[t],alpha=.7,label=t)
rho=stats.spearmanr(depth,orig.Observed_Species); ax.set_xscale('log'); ax.set_xlabel('Reads per sample'); ax.set_ylabel('Observed species (unrarefied)'); ax.set_title('Richness depends on depth'); ax.text(0.03,0.93,f'Spearman ρ = {rho[0]:.2f}, p < 0.001',transform=ax.transAxes,fontsize=8); ax.legend(frameon=False,fontsize=8,loc='lower right'); panel('B',ax)
ax=axs[2]; L=pd.read_csv('alpha_lmm.csv'); L=L[L.metric=='Observed']; ypos=0; yt=[];yl=[]
for con in ['S2 vs S1','S3 vs S2','S3 vs S1']:
    for st,c,nm in [('raw','#7F7F7F','unrarefied'),('raw+logdepth','#5B6C8F','depth-adjusted'),('rarefied','#C0392B','rarefied')]:
        r=L[(L.contrast==con)&(L.set==st)].iloc[0]; ax.errorbar(r.est,ypos,xerr=[[r.est-r.lo],[r.hi-r.est]],fmt='o',c=c,capsize=2,ms=5); ypos-=1
    yt.append(ypos+2); yl.append(con); ypos-=0.8
ax.axvline(0,c='k',lw=.8); ax.set_yticks(yt); ax.set_yticklabels(yl); ax.set_xlabel('Mixed-model difference in observed richness (95% CI)')
for st,c,nm in [('raw','#7F7F7F','unrarefied'),('raw+logdepth','#5B6C8F','depth-adjusted'),('rarefied','#C0392B','rarefied')]: ax.scatter([],[],c=c,label=nm)
ax.legend(frameon=False,fontsize=8,loc='lower left'); ax.set_title('Effect of depth handling'); panel('C',ax)
ax=axs[3]; RS=pd.read_csv('rarefaction_depth_sens.csv'); RS=RS[RS.metric=='Observed richness']
for col,c,lab in [('S2-S1','#E8742C','S2 − S1'),('S3-S2','#4C9A5B','S3 − S2'),('S3-S1','#3B7DD8','S3 − S1')]:
    ax.plot(RS.depth,RS[col],'-o',c=c,label=lab)
for _,r in RS.iterrows(): ax.text(r.depth,r['S3-S2']-3.2,f"p={r['S3-S2_p']:.2f}",ha='center',fontsize=6.5,color='#4C9A5B')
ax.axhline(0,c='k',lw=.8); ax.set_xlabel('Rarefaction depth (reads)'); ax.set_ylabel('Mixed-model difference (species)'); ax.set_title('Sensitivity to rarefaction depth'); ax.legend(frameon=False,fontsize=8); panel('D',ax)
fig.suptitle('Supplementary Figure S1. Sequencing depth and its influence on richness estimates',fontsize=11,fontweight='bold'); fig.tight_layout(rect=[0,0,1,0.94]); fig.savefig('figs/FigureS1_depth_sensitivity.png',dpi=220); plt.close(fig); print('s1')
# ================= Fig S2 paired rarefied alpha =================
at=pd.read_csv('alpha_tests_rar.csv').set_index('metric')
fig,axs=plt.subplots(1,3,figsize=(13,4.4))
for ax,m,lab,l in [(axs[0],'Observed','Observed richness (rarefied, 7,500 reads)','A'),(axs[1],'Shannon','Shannon diversity (rarefied)','B')]:
    W=pd.DataFrame({t:pd.Series(Ar[m][tpv==t].values,index=pid[tpv==t]) for t in tps}).loc[complete]
    for p_ in complete: ax.plot([0,1,2],W.loc[p_],c='0.8',lw=.7,zorder=1)
    for k,t in enumerate(tps): ax.scatter(np.full(len(W),k)+np.random.default_rng(1).uniform(-.06,.06,len(W)),W[t],s=14,c=COL[t],zorder=2)
    ax.plot([0,1,2],W.median(),'k-o',lw=2.2,zorder=3)
    r=at.loc[m]
    ax.text(0.5,1.01,f"q = {r.W_S1S2_q:.3f}" if r.W_S1S2_q>=0.001 else 'q < 0.001',transform=ax.get_xaxis_transform(),ha='center',fontsize=8)
    ax.text(1.5,1.01,f"q = {r.W_S2S3_q:.2f}",transform=ax.get_xaxis_transform(),ha='center',fontsize=8)
    ax.text(1,1.07,('q < 0.001' if r.W_S1S3_q<0.001 else f"q = {r.W_S1S3_q:.3f}"),transform=ax.get_xaxis_transform(),ha='center',fontsize=8)
    ax.set_xticks([0,1,2]); ax.set_xticklabels([LAB[t] for t in tps]); ax.set_ylabel(lab); ax.set_title(lab.split(' (')[0],pad=22); panel(l,ax)
ax=axs[2]; rr=pd.read_csv('ipw.csv'); rr=rr[rr.metric=='Observed']; x=np.arange(3); w=.35
ax.bar(x-w/2,rr.unweighted_mean_change,w,color='#9AA5B1',label='Complete cases'); ax.bar(x+w/2,rr.IPW_mean_change,w,color='#5B6C8F',label='Inverse-probability weighted')
ax.errorbar(x+w/2,rr.IPW_mean_change,yerr=[rr.IPW_mean_change-rr.IPW_CI_lo,rr.IPW_CI_hi-rr.IPW_mean_change],fmt='none',c='k',capsize=3)
ax.axhline(0,c='k',lw=.8); ax.set_xticks(x); ax.set_xticklabels(rr.contrast); ax.set_ylabel('Mean change in observed richness (rarefied)'); ax.legend(frameon=False,fontsize=8); ax.set_title('Weighting for attrition',pad=22); panel('C',ax)
fig.suptitle('Supplementary Figure S2. Within-patient change in rarefied alpha diversity (46 patients)',fontsize=11,fontweight='bold'); fig.tight_layout(rect=[0,0,1,0.94]); fig.savefig('figs/FigureS2_paired_alpha_rarefied.png',dpi=220); plt.close(fig); print('s2')
