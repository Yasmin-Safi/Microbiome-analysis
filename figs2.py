exec(open('prep.py').read())
import json, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,'axes.spines.right':False,'axes.titlesize':10,'axes.titleweight':'bold'})
COL={'S1':'#3B7DD8','S2':'#E8742C','S3':'#4C9A5B'}; tps=['S1','S2','S3']
def panel(l,ax): ax.text(-0.12,1.06,l,transform=ax.transAxes,fontsize=13,fontweight='bold')
s=pd.read_csv('da_species.csv').set_index('feature'); ge=pd.read_csv('da_genus.csv').set_index('feature')
SEN=pd.read_csv('da_sensitivity.csv'); _da=json.load(open('da_summary.json'))['species']; SEN=pd.concat([pd.DataFrame([{'pseudocount':0.5,'min_prevalence':0.10,'n_tested':323,'S2vS1':_da['LMM_S2vS1_q05'],'S3vS2':_da['LMM_S3vS2_q05'],'S3vS1':_da['LMM_S3vS1_q05'],'basis':'primary'}]),SEN.assign(basis='full')],ignore_index=True); conc=json.load(open('concordance.json')); da=json.load(open('da_summary.json'))['species']
# ============ S3: differential abundance =============
fig=plt.figure(figsize=(16,6.4)); gs=fig.add_gridspec(1,3,width_ratios=[1.05,1.35,1.0])
ax=fig.add_subplot(gs[0]); labs=['Kruskal-Wallis\nunadjusted p<0.05','KW, BH across\n805 species','KW, BH across\n323 filtered','Mixed models\nS3 vs S1 (FDR)','Mixed models\n+ depth adj.','Dirichlet-CLR\npaired Wilcoxon\nS3 vs S1']
vals=[da['KW_p05'],da['KW_q05_all'],da['KW_q05_filt'],da['LMM_S3vS1_q05'],da['LMMd_S3vS1_q05'],conc['aldex']]
cols=['#C9CED6','#5B6C8F','#5B6C8F','#C0392B','#C0392B','#B8860B']
ax.barh(range(6)[::-1],vals,color=cols); 
for i,v in zip(range(6)[::-1],vals): ax.text(v+1,i,str(v),va='center',fontsize=9)
ax.set_yticks(range(6)[::-1]); ax.set_yticklabels(labs,fontsize=8); ax.set_xlabel('Number of species called significant'); ax.set_title('Effect of statistical approach'); panel('A',ax)
ax=fig.add_subplot(gs[1]); rob=s[(s.LMM_S3vS1_q<.05)&(s.LMMd_S3vS1_q<.05)].copy(); rob['abs']=rob.LMM_S3vS1_est.abs(); rob=rob.sort_values('abs',ascending=False).head(30).sort_values('LMM_S3vS1_est')
ax.barh(range(len(rob)),rob.LMM_S3vS1_est,color=np.where(rob.LMM_S3vS1_est>0,'#4C9A5B','#C0392B'))
ax.set_yticks(range(len(rob))); ax.set_yticklabels(rob.index,fontsize=7,style='italic'); ax.axvline(0,c='k',lw=.8)
ax.set_xlabel('CLR difference, S3 vs S1 (mixed model)'); ax.set_title(f'Species robust to depth adjustment (top 30 of {int(((s.LMM_S3vS1_q<.05)&(s.LMMd_S3vS1_q<.05)).sum())})'); panel('B',ax)
ax=fig.add_subplot(gs[2]); x=np.arange(len(SEN)); w=.27
for k,(c,col,lab) in enumerate([('S2vS1','#E8742C','S2 vs S1'),('S3vS2','#4C9A5B','S3 vs S2'),('S3vS1','#3B7DD8','S3 vs S1')]):
    ax.bar(x+(k-1)*w,SEN[c],w,color=col,label=lab)
    for xi,v in zip(x+(k-1)*w,SEN[c]): ax.text(xi,v+1,str(v),ha='center',fontsize=7)
ax.set_xticks(x); ax.set_xticklabels([('PRIMARY\n' if bs=='primary' else 'full-comp.\n')+f'pc {a}, ≥{int(b*100)}%' for a,b,bs in zip(SEN.pseudocount,SEN.min_prevalence,SEN.basis)],fontsize=6.8); ax.set_ylabel('Species with FDR q < 0.05'); ax.legend(frameon=False,fontsize=8); ax.set_title('Sensitivity to CLR specification'); ax.text(0.5,-0.2,'PRIMARY = CLR of prevalence-filtered table; others = CLR on full composition, then filtered',transform=ax.transAxes,ha='center',fontsize=6.5); panel('C',ax)
fig.suptitle('Supplementary Figure S3. Differential abundance after multiplicity correction and under alternative compositional specifications',fontsize=11,fontweight='bold'); fig.tight_layout(rect=[0,0,1,0.94]); fig.savefig('figs/FigureS3_differential_abundance.png',dpi=220); plt.close(fig); print('s3')
# ============ S4: attrition =============
fig,axs=plt.subplots(1,3,figsize=(14,4.6))
ax=axs[0]; reasons=[('Death',5),('Metastasis / progression',6),('Withdrew',2),('Stopped treatment, travelled',1),('Insufficient sample',2)]
c=['#7F7F7F','#C0392B','#E8A33D','#5B6C8F','#9AA5B1']; left=0
for (l,n),cc in zip(reasons,c): ax.barh(0,n,left=left,color=cc,label=f'{l} (n={n})'); ax.text(left+n/2,0,str(n),ha='center',va='center',color='white',fontsize=9,fontweight='bold'); left+=n
ax.barh(1,3,color='#B8860B',label='S2 not received (n=3)'); ax.text(1.5,1,'3',ha='center',va='center',color='white',fontweight='bold')
ax.set_yticks([0,1]); ax.set_yticklabels(['S3 missing\n(16)','S2 missing\n(3)']); ax.set_xlabel('Patients'); ax.legend(frameon=False,fontsize=7,loc='upper center',bbox_to_anchor=(.5,-0.18),ncol=2); ax.set_title('Reasons for missing samples'); ax.text(0.98,0.55,'11 of 16 S3 losses:\ndeath or progression',transform=ax.transAxes,ha='right',fontsize=8); panel('A',ax)
ax=axs[1]; T=pd.read_csv('tipping.csv')
for m,c_ in [('Observed','#5B6C8F'),('Shannon','#C0392B')]:
    t=T[T.metric==m]; ax.plot(t.delta_SD,t.p_S2vsS3,'-o',c=c_,label='Observed richness' if m=='Observed' else 'Shannon')
ax.axhline(0.05,ls='--',c='k',lw=1); ax.text(2.0,0.06,'p = 0.05',ha='right',fontsize=8); ax.set_yscale('log'); ax.set_xlabel('Additional decline assumed for the 13 patients\nwith S2 but no S3 (SD of observed S2→S3 change)'); ax.set_ylabel('Paired Wilcoxon p, S2 vs S3 (n = 59)'); ax.legend(frameon=False,fontsize=8); ax.set_title('Tipping-point analysis'); panel('B',ax)
ax=axs[2]; at=pd.read_csv('attrition_compare.csv')
v=at[at.Variable.isin(['Age','BMI','Male sex','Baseline Shannon','Baseline observed richness','Post-RT (S2) Shannon','Post-RT (S2) observed richness','Any comorbidity','Current/ex-smoker','cN2','cT4','TUMOR SIZE','Chemotherapy cycles, median'])].copy()
v['lab']=v.Variable.replace({'TUMOR SIZE':'Tumor size','Chemotherapy cycles, median':'Chemo cycles'}); v=v.iloc[::-1]
ax.scatter(-np.log10(v.p),range(len(v)),c=np.where(v.p<.05,'#C0392B','#5B6C8F'),s=36); ax.axvline(-np.log10(.05),ls='--',c='k',lw=1); ax.set_yticks(range(len(v))); ax.set_yticklabels(v.lab,fontsize=7.5); ax.set_xlabel('−log10 p (complete vs incomplete patients)'); ax.set_title('Completers vs non-completers'); panel('C',ax)
fig.suptitle('Supplementary Figure S4. Sample attrition and its potential influence on the S2→S3 inference',fontsize=11,fontweight='bold'); fig.tight_layout(rect=[0,0,1,0.93]); fig.savefig('figs/FigureS4_attrition.png',dpi=220); plt.close(fig); print('s4')
# ============ S5: core + genus =============
fig,axs=plt.subplots(1,2,figsize=(13,4.8))
ax=axs[0]; co=pd.read_csv('core_sensitivity.csv'); jc=json.load(open('core.json'))
rows=[('Unrarefied, all samples\n(n = 62/59/46)',[jc['orig'][t] for t in tps],'#9AA5B1'),('Rarefied, 46 patients\n(≥1 read)',[float(co[(co.set=='complete cohort (n=46/TP)')&(co.min_reads==1)][f'{t}_core50'].iloc[0]) for t in tps],'#5B6C8F'),('Rarefied, 46 patients\n(≥3 reads)',[float(co[(co.set=='complete cohort (n=46/TP)')&(co.min_reads==3)][f'{t}_core50'].iloc[0]) for t in tps],'#C0392B')]
x=np.arange(3); w=.26
for k,(lab,v,c) in enumerate(rows):
    ax.bar(x+(k-1)*w,v,w,color=c,label=lab)
    for xi,vi in zip(x+(k-1)*w,v): ax.text(xi,vi+1.5,f'{vi:.0f}',ha='center',fontsize=8)
ax.set_xticks(x); ax.set_xticklabels(tps); ax.set_ylabel('Species detected in >50% of samples'); ax.legend(frameon=False,fontsize=8); ax.set_title('Operational "core" membership'); panel('A',ax)
ax=axs[1]; kg=['Fusobacterium','Parvimonas','Peptostreptococcus','Porphyromonas','Gemella','Streptococcus','Bacteroides','Limosilactobacillus']
cm_=plt.cm.tab10(np.arange(len(kg)))
for g,c in zip(kg,cm_):
    r=ge.loc[g]; y=[r.mean_S1*100,r.mean_S2*100,r.mean_S3*100]; q=r.LMM_S3vS1_q
    ax.plot([0,1,2],y,'-o',c=c,label=f'{g} (S3 vs S1 q {"<0.001" if q<0.001 else f"={q:.3f}"})')
ax.set_yscale('log'); ax.set_xticks([0,1,2]); ax.set_xticklabels(['S1','S2','S3']); ax.set_ylabel('Mean relative abundance (%, log scale)'); ax.legend(frameon=False,fontsize=7.2,loc='lower left'); ax.set_title('Genus-level trajectories (mixed-model q)'); panel('B',ax)
fig.suptitle('Supplementary Figure S5. Core-membership sensitivity and genus-level robustness',fontsize=11,fontweight='bold'); fig.tight_layout(rect=[0,0,1,0.93]); fig.savefig('figs/FigureS5_core_genus.png',dpi=220); plt.close(fig); print('s5')
