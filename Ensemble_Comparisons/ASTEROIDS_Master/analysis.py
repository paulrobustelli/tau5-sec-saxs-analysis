from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
DATA=Path('data')
META=json.loads((DATA/'metadata.json').read_text())
Z={p:dict(np.load(DATA/f'{p}_structural.npz')) for p in ['WT','AA']}
MD=dict(np.load(DATA/'MD_structural.npz'))
COLORS=['#777777','#2876b8','#df762f','#24946d','#ab3bb2']
plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False,'figure.dpi':100})
def curves(key,stem):return [pd.read_csv(DATA/f'{key}_ens{i}_{stem}.csv') for i in range(1,4)]
def members(key,view='full',with_md=True):
 c=META['conditions'][key];z=Z[c['protein']];sl=slice(None) if view=='full' else slice(63,119);cols=[0,1] if view=='full' else [2,3];pairs=z['pairs'];ans=[]
 for label,ids in [('Raw pool',c['pool'])]+[(f'Ens {i}',v) for i,v in enumerate(c['selected'],1)]:
  mat=np.zeros((120,120));mat[pairs[:,0],pairs[:,1]]=z['contacts'][ids].mean(0);mat+=mat.T
  ans.append({'name':label,'v':z['values'][ids][:,cols],'h':z['helix'][ids,sl],'c':mat[sl,sl],'n':len(ids),'ids':ids})
 if view=='R2R3' and with_md:ans.append({'name':'2022 apo WT MD','v':MD['values'],'h':MD['helix'],'c':MD['contacts'],'n':len(MD['values'])})
 return ans

def run_counts(h):
 hist=np.zeros(h.shape[1]+1);longest=[]
 for row in h:
  edges=np.diff(np.r_[False,row,False].astype(int));lens=np.where(edges==-1)[0]-np.where(edges==1)[0]
  hist+=np.bincount(lens,minlength=len(hist));longest.append(max(lens,default=0))
 return hist/len(h),np.array(longest)

def summaries():
 out=[];qout=[]
 for key,c in META['conditions'].items():
  shifts=curves(key,'shifts');saxs=curves(key,'saxs')
  for view in ['full','R2R3']:
   for i,g in enumerate(members(key,view,False)):
    d={'ensemble':key+('_raw' if i==0 else f'_ens{i}'),'view':view,'n':g['n'],'Rg_mean_nm':g['v'][:,1].mean(),'Rg_SD_nm':g['v'][:,1].std(ddof=0),'Sa_mean':g['v'][:,0].mean(),'Sa_SD':g['v'][:,0].std(ddof=0)}
    for sep in [5,10,20]:
     u,v=np.triu_indices(len(g['c']),sep+1);d[f'contact_gt{sep}_prob']=g['c'][u,v].mean();d[f'contact_gt{sep}_pairs']=len(u);d[f'contact_gt{sep}_count']=g['c'][u,v].sum()
    sh=shifts[max(0,i-1)];sc=saxs[max(0,i-1)];shiftcol='Exp_minus_pool_ppm' if i==0 else 'Exp_minus_selected_ppm';calc='I_pool' if i==0 else 'I_selected'
    if view=='R2R3':sh=sh[sh.AR_residue.between(391,446)]
    for atom,a in sh.groupby('atom'):d[atom+'_RMSD_ppm']=np.sqrt(np.mean(a[shiftcol]**2))
    # SAXS restraints are full-chain measurements; never assign them to the cropped segment.
    if view=='full':
     e=sc.I_exp-sc[calc];r=e/sc.sigma_exp;d['SAXS_RMSD_intensity']=np.sqrt(np.mean(e**2));d['SAXS_chi2_per_point']=np.mean(r**2)
     for lo,hi in [(0,.04),(.04,.1),(.1,.25),(.25,np.inf)]:
      mask=(sc.q_Ainv>=lo)&(sc.q_Ainv<hi);qout.append({'ensemble':d['ensemble'],'q_start_Ainv':lo,'q_end_Ainv':hi,'n_points':int(mask.sum()),'chi2_sum':float(np.sum(r[mask]**2)),'chi2_per_point':float(np.mean(r[mask]**2)),'fraction_total_chi2':float(np.sum(r[mask]**2)/np.sum(r**2))})
    out.append(d)
 g={'v':MD['values'],'c':MD['contacts']};d={'ensemble':'2022 apo WT MD','view':'R2R3','n':len(g['v']),'Rg_mean_nm':g['v'][:,1].mean(),'Rg_SD_nm':g['v'][:,1].std(),'Sa_mean':g['v'][:,0].mean(),'Sa_SD':g['v'][:,0].std()}
 for sep in [5,10,20]:
  u,v=np.triu_indices(56,sep+1);d[f'contact_gt{sep}_prob']=g['c'][u,v].mean();d[f'contact_gt{sep}_pairs']=len(u);d[f'contact_gt{sep}_count']=g['c'][u,v].sum()
 out.append(d);return pd.DataFrame(out),pd.DataFrame(qout)
SUMMARY,QSUMMARY=summaries()
def summary_table(view):
 d=SUMMARY[SUMMARY.view==view].copy();o=d[['ensemble','n']].copy()
 for atom in ['CA','C','CB','N','NH','HA']:
  col=atom+'_RMSD_ppm'
  if col in d:o[atom+' RMSD (ppm)']=d[col].map(lambda x:'—' if pd.isna(x) else f'{x:.3f}')
 if view=='full':o['SAXS χ²/point']=d.SAXS_chi2_per_point.map(lambda x:f'{x:.3f}')
 o['Rg ± SD (nm)']=[f'{a:.3f} ± {b:.3f}' for a,b in zip(d.Rg_mean_nm,d.Rg_SD_nm)];o['Sα ± SD']=[f'{a:.3f} ± {b:.3f}' for a,b in zip(d.Sa_mean,d.Sa_SD)]
 for sep in [5,10,20]:o[f'P(contact), |i−j|>{sep}']=d[f'contact_gt{sep}_prob'].map(lambda x:f'{x:.4f}')
 return o.reset_index(drop=True)

def sa_rg(key,view):
 groups=members(key,view);fig,axs=plt.subplots(1,len(groups),figsize=(4*len(groups),4.2),sharex=True,sharey=True,layout='constrained');extent=[.9,2.5] if view=='R2R3' else [1.5,6.5]
 for ax,g in zip(axs,groups):
  a,xe,ye=np.histogram2d(g['v'][:,1],g['v'][:,0],30,[extent,[0,25]],density=True)
  a=-.001987*300*np.log(np.flipud(a)+.000001)
  im=ax.imshow(a,interpolation='gaussian',extent=[ye[0],ye[-1],xe[0],xe[-1]],cmap='jet',aspect='auto',vmin=.1,vmax=3)
  inside=((g['v'][:,1]>=extent[0])&(g['v'][:,1]<=extent[1])&(g['v'][:,0]<=25)).sum()
  ax.set(xlabel='Sα',title=f"{g['name']} (n={g['n']:,})\n{inside/g['n']:.1%} inside display range",xlim=(0,24.9));fig.colorbar(im,ax=ax,shrink=.8,ticks=[1,2,3])
 axs[0].set_ylabel('Cα Rg (nm)');fig.suptitle(key+' | '+view+' | original AR histogram/log transform; not pool free energy');plt.show()

def contacts(key,view):
 groups=members(key,view);fig,axs=plt.subplots(1,len(groups),figsize=(4.3*len(groups),4.5),layout='constrained')
 for ax,g in zip(axs,groups):
  sns.heatmap(g['c'],ax=ax,cmap='jet',vmin=0,vmax=1,cbar=True,cbar_kws={'shrink':.75,'label':'Contact probability'},rasterized=True);ax.invert_yaxis();ax.grid(alpha=.5)
  if view=='full':ticks=np.array([0,12,32,52,72,92,112,119]);labels=['GP','340','360','380','400','420','440','447']
  else:ticks=np.arange(0,56,10);labels=[str(391+i) for i in ticks]
  ax.set_xticks(ticks+.5,labels,rotation=45);ax.set_yticks(ticks+.5,labels);ax.set(title=f"{g['name']} (n={g['n']:,})",xlabel='AR residue',ylabel='AR residue')
 fig.suptitle(key+' | '+view+' | closest-heavy <1.2 nm; adjacent pairs retained');plt.show()

def distributions(key,view):
 groups=members(key,view);fig,axs=plt.subplots(1,3,figsize=(17,4.6),layout='constrained');maxsa=max(g['v'][:,0].max() for g in groups);bins=np.linspace(0,np.ceil(maxsa)+1,31);maxlen=0
 for g,col in zip(groups,COLORS):
  axs[0].hist(g['v'][:,0],bins=bins,weights=np.full(g['n'],1/g['n']),histtype='step',lw=1.8,label=g['name'],color=col)
  hist,long=run_counts(g['h']);xs=np.flatnonzero(hist);maxlen=max(maxlen,long.max());axs[1].step(np.arange(len(hist)),hist,where='mid',label=g['name'],color=col)
  h=np.bincount(long,minlength=g['h'].shape[1]+1)/g['n'];axs[2].step(np.arange(len(h)),h,where='mid',label=g['name'],color=col)
 axs[0].set(xlabel='Sα',ylabel='Fraction of conformers / bin');axs[1].set(xlabel='DSSP-H run length (residues)',ylabel='Mean number of helices / conformer',xlim=(0,maxlen+1));axs[2].set(xlabel='Longest DSSP-H run (0 = none)',ylabel='Fraction of conformers',xlim=(0,maxlen+1))
 for ax in axs:ax.legend(fontsize=9)
 fig.suptitle(key+' | '+view+' | realized DSSP helix lengths, not imposed torsion-window lengths');plt.show()

def helicity(key):
 c=META['conditions'][key];d2d=pd.read_csv(DATA/f"{c['protein']}_d2d.csv");fig,axs=plt.subplots(2,1,figsize=(16,7),layout='constrained')
 for ax,view in zip(axs,['full','R2R3']):
  xx=np.arange(328,448) if view=='full' else np.arange(391,447)
  for g,col in zip(members(key,view),COLORS):ax.plot(xx,g['h'].mean(0),label=g['name'],color=col)
  d=d2d if view=='full' else d2d[d2d.author_residue.between(391,446)];ax.plot(d.author_residue,d.helix,'k--',lw=2,label='δ2D '+c['protein']);ax.set(ylim=(0,1),xlim=(330,446) if view=='full' else (391,446),ylabel='DSSP-H population',xlabel='AR residue',title=view);ax.legend(ncol=6,fontsize=9)
 fig.suptitle(key+' | before/after equal-weight subset selection');plt.show()

def chemical_shifts(key,secondary=False):
 ds=curves(key,'shifts');atoms=[a for a in ['CA','C','CB','N','NH','HA'] if a in set(ds[0].atom)];fig,axs=plt.subplots(len(atoms),1 if secondary else 2,figsize=(17,3*len(atoms)),squeeze=False,layout='constrained')
 for row,atom in enumerate(atoms):
  s=[d[d.atom==atom].set_index('AR_residue') for d in ds];x=s[0].index.to_numpy()
  if secondary:
   ax=axs[row,0];ax.plot(x,s[0].experimental_secondary_ppm,'k.-',label='Experiment');ax.plot(x,s[0].pool_secondary_ppm,color=COLORS[0],label='Raw pool')
   for i,d in enumerate(s,1):ax.plot(d.index,d.selected_secondary_ppm,color=COLORS[i],label=f'Ens {i}')
   ax.set(ylabel=f'Δδ {atom} (ppm)',xlim=(330,446));ax.legend(ncol=5,fontsize=9)
  else:
   vals=[s[0].Exp_minus_pool_ppm]+[v.Exp_minus_selected_ppm for v in s]
   for ax,absolute in zip(axs[row],[False,True]):
    for i,y in enumerate(vals):ax.bar(x+(i-1.5)*.2,np.abs(y) if absolute else y,width=.2,color=COLORS[i],label='Raw pool' if i==0 else f'Ens {i}')
    ax.axhline(0,color='k',lw=.5);ax.set(ylabel=f"{atom}: {'|Exp − Calc|' if absolute else 'Exp − Calc'} (ppm)",xlim=(329,447));ax.legend(ncol=4,fontsize=8)
 for ax in axs[-1]:ax.set_xlabel('AR residue')
 fig.suptitle(key+' | '+('secondary shifts' if secondary else 'signed and absolute residuals'));plt.show()

def guinier_estimate(q,I,sigma):
 mask=(q<=.03)&(I>0);x=np.asarray(q[mask])**2;y=np.log(np.asarray(I[mask]));w=np.asarray(I[mask]/sigma[mask]);slope,intercept=np.polyfit(x,y,1,w=w)
 return np.sqrt(-3*slope) if slope<0 else np.nan,mask

def saxs(key):
 ds=curves(key,'saxs');d=ds[0];q=d.q_Ainv;exp=d.I_exp;err=d.sigma_exp;vals=[d.I_pool]+[s.I_selected for s in ds];fig,axs=plt.subplots(2,2,figsize=(16,9),layout='constrained')
 axs[0,0].errorbar(q,exp,yerr=err,fmt='k.',ms=2,alpha=.5,label='Experiment');mask=q<=.06;axs[1,0].plot(q[mask]**2,np.log(exp[mask]),'k.',label='Experiment')
 er,gm=guinier_estimate(q,exp,err)
 for i,v in enumerate(vals):
  label='Raw pool' if i==0 else f'Ens {i}';r=(exp-v)/err;rg,_=guinier_estimate(q,v,err)
  axs[0,0].plot(q,v,color=COLORS[i],label=label);axs[0,1].plot(q,r,color=COLORS[i],label=label);axs[1,0].plot(q[mask]**2,np.log(v[mask]),color=COLORS[i],label=f'{label}: Rg={rg:.1f} Å');axs[1,1].plot(q,np.cumsum(r*r)/np.sum(r*r),color=COLORS[i],label=label)
 axs[0,0].set(yscale='log',xlabel='q (Å⁻¹)',ylabel='I(q), experimental scale');axs[0,1].set(xlabel='q (Å⁻¹)',ylabel='(Exp − Calc)/σ');axs[0,1].axhline(0,color='k',lw=.5)
 axs[1,0].axvline(.03**2,color='k',ls=':',label='Slope fit qmax=.03 Å⁻¹');axs[1,0].set(xlabel='q² (Å⁻²)',ylabel='ln I(q)',title=f'Low-q diagnostic: exp Rg={er:.1f} Å; qmax·Rg={.03*er:.2f}')
 axs[1,1].axvline(.04,color='k',ls=':');axs[1,1].set(xlabel='q (Å⁻¹)',ylabel='Cumulative fraction of SAXS χ²',ylim=(0,1));
 for ax in axs.flat:ax.legend(fontsize=9)
 fig.suptitle(key+' | raw and selected curves, each with its own analytical intensity scale');plt.show()
