from pathlib import Path
import json
import numpy as np,pandas as pd
from scipy.stats import spearmanr,rankdata
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
O=Path(__file__).resolve().parent;R=O.parents[1];C=R/'outputs/Ensemble_Comparison_121_2026-10-15';OLD=Path('/Users/paulrobustelli/Documents/Codex/2026-09-28/is-th')
NAMES=['REST_WT','IDPForge_WT','bAIes_WT','ASTEROIDS_WT','IDPForge_AA','bAIes_AA','ASTEROIDS_AA']
FIELDS=['buried_HSA_nonlocal_A2','hydrophobic_pairs_sep6_50','aromatic_hydrophobic_pairs_sep6_50']
def load(name):
 d=pd.read_csv(O/name/'metrics.csv');w=pd.read_csv(O/name/'weights.csv');variant=name.split('_')[-1]
 q,y,sy=np.loadtxt(R/f'outputs/IDPForge_D2D_MaxEnt_2026-10-15/experimental/{variant}_SAXS.dat',unpack=True)
 if name.startswith('ASTEROIDS'):
  case='WT' if variant=='WT' else 'W397A_W433A';p=pd.read_csv(OLD/f'work/refined-{case}-saxs/predicted_dro5.csv');assert np.array_equal(p.conformer_id.astype(str),d.conformer_id.astype(str));nativeq=np.array([float(c.split('=')[1]) for c in p.columns[1:]]);raw=np.array([np.interp(q,nativeq,row) for row in p.iloc[:,1:].to_numpy()]).T
 else:
  z=np.load(C/name/'observable_matrix.npz');raw=z['raw_saxs'];assert np.allclose(z['q'],q)
 assert raw.shape==(len(q),len(d))
 return d,w,raw,q,y,sy

def score(curve,q,y,sy):
 scale=np.sum(curve*y/sy**2)/np.sum((curve/sy)**2);pred=scale*curve;z=(pred-y)/sy
 return scale,{'chi2_all':np.mean(z*z),'chi2_mid':np.mean(z[(q>.08)&(q<.25)]**2),'chi2_low':np.mean(z[q<=.08]**2),'chi2_high':np.mean(z[q>=.25]**2)},pred

def main():
 summary=[];corr=[];groupstats=[];curverows=[];matched=[]
 fig,axs=plt.subplots(7,4,figsize=(18,23),layout='constrained')
 fig2,ax2=plt.subplots(2,3,figsize=(15,8),layout='constrained')
 colors={'REST':'#cc79a7','IDPForge':'#0072b2','bAIes':'#d55e00','ASTEROIDS':'#009e73'}
 for row,name in enumerate(NAMES):
  d,w,raw,q,y,sy=load(name);mid=(q>.08)&(q<.25);v=name.split('_')[-1];method=name.rsplit('_',1)[0]
  for stage in ['uniform','combined']:
   ww=w[stage].to_numpy();ss={'ensemble':name,'stage':stage,'N':len(d),'fraction_Rg_lt19':float(ww@(d.Rg_CA_A<19)),'Neff':1/(ww@ww)}
   for f in ['Rg_CA_A']+FIELDS+['hydrophobic_mean_coordination','hydrophobic_residues_ge3_contacts']:ss[f]=float(ww@d[f])
   summary.append(ss)
  # Correlations are structural descriptors, not independent per-frame SAXS fits.
  shape=raw[np.argmin(abs(q-.15))]/raw[np.argmin(abs(q-.03))]
  for f in FIELDS:
   X=np.c_[np.ones(len(d)),rankdata(d.Rg_CA_A),rankdata(d.Rg_CA_A)**2/len(d)];aa=rankdata(d[f]);bb=rankdata(shape);aa-=X@np.linalg.lstsq(X,aa,rcond=None)[0];bb-=X@np.linalg.lstsq(X,bb,rcond=None)[0]
   corr.append({'ensemble':name,'metric':f,'rho_with_Rg':spearmanr(d[f],d.Rg_CA_A).statistic,'rho_with_I015_I003':spearmanr(d[f],shape).statistic,'rank_residual_correlation_after_Rg_and_Rg2':np.corrcoef(aa,bb)[0,1]})
  splits=[('Rg_CA_A',[19,25,35],['<19 Å','19–25 Å','25–35 Å','≥35 Å'])]
  for f in FIELDS:
   cuts=np.unique(np.quantile(d[f],[1/3,2/3]));splits.append((f,cuts,[f'group {i+1}' for i in range(len(cuts)+1)]))
  for col,(f,cuts,labels) in enumerate(splits):
   group=np.digitize(d[f],cuts)
   for g,label in enumerate(labels):
    use=group==g;n=use.sum()
    if not n:continue
    curve=raw[:,use].mean(1);scale,met,pred=score(curve,q,y,sy)
    axs[row,col].plot(q,pred/y,label=f'{label}: n={n}, mid χ²={met["chi2_mid"]:.1f}',lw=1)
    groupstats.append({'ensemble':name,'split':f,'group':g,'label':label,'N':int(n),'lower':-np.inf if g==0 else cuts[g-1],'upper':np.inf if g==len(cuts) else cuts[g],'Rg_mean':d.loc[use,'Rg_CA_A'].mean(),'burial_mean':d.loc[use,FIELDS[0]].mean(),'combined_weight_fraction':float(w.loc[use,'combined'].sum()),'global_scale':scale,**met})
    for qi,ri,pi in zip(q,curve,pred):curverows.append({'ensemble':name,'split':f,'group':g,'q':qi,'raw_intensity':ri,'scaled_intensity':pi})
   ax=axs[row,col];ax.axhline(1,color='k',lw=.7);ax.axvspan(.08,.25,color='grey',alpha=.1);ax.set_ylim(.2,1.8);ax.set_title(name+' • '+f.replace('_',' '),fontsize=9);ax.legend(fontsize=6);ax.set(xlabel='q (Å⁻¹)',ylabel='Prediction / experiment')
  # Match high/low burial using same counts per 2-A Rg bin, no curve-fitting in selection.
  groups=[[],[]];rgids=np.floor(d.Rg_CA_A/2).astype(int)
  for rgid in sorted(set(rgids)):
   ii=np.flatnonzero(rgids==rgid);ii=ii[np.argsort(d.loc[ii,FIELDS[0]].to_numpy(),kind='stable')];half=len(ii)//2
   if half>=3:groups[0].extend(ii[:half]);groups[1].extend(ii[-half:])
  for lab,ii in zip(['low_burial_Rg_matched','high_burial_Rg_matched'],groups):
   if not ii:continue
   scale,met,pred=score(raw[:,ii].mean(1),q,y,sy);matched.append({'ensemble':name,'group':lab,'N':len(ii),'Rg_mean':d.loc[ii,'Rg_CA_A'].mean(),'burial_mean':d.loc[ii,FIELDS[0]].mean(),**met})
  for col,f in enumerate(FIELDS):
   ax=ax2[0 if v=='WT' else 1,col];bins=np.arange(14,62,2);cent=[];means=[]
   for lo,hi in zip(bins[:-1],bins[1:]):
    mask=(d.Rg_CA_A>=lo)&(d.Rg_CA_A<hi)
    if mask.sum()>=5:cent.append(d.loc[mask,'Rg_CA_A'].mean());means.append(d.loc[mask,f].mean())
   ax.plot(cent,means,'o-',ms=3,label=method,color=colors[method]);ax.axvline(19,color='grey',ls=':');ax.set(title=v+' • '+f.replace('_',' '),xlabel='Cα Rg (Å)',ylabel=f);ax.legend(fontsize=8)
 fig.suptitle('Conditional subensemble SAXS curves: one fitted full-q scale per group; no per-conformer scaling\nHydrophobic/contact tertiles are within each pool; shaded region q=0.08–0.25 Å⁻¹. These are uniform group averages, not MaxEnt fits.',fontsize=12)
 for ext in ['png','pdf']:fig.savefig(O/f'conditional_SAXS.{ext}',dpi=160);fig2.savefig(O/f'packing_at_matched_Rg.{ext}',dpi=160)
 plt.close('all')
 pd.DataFrame(summary).to_csv(O/'ensemble_packing_before_after.csv',index=False);pd.DataFrame(corr).to_csv(O/'structure_SAXS_correlations.csv',index=False);pd.DataFrame(groupstats).to_csv(O/'conditional_group_scores.csv',index=False);pd.DataFrame(curverows).to_csv(O/'conditional_SAXS_curves.csv',index=False);pd.DataFrame(matched).to_csv(O/'matched_Rg_burial_group_scores.csv',index=False)
 print(pd.DataFrame(summary).to_string(index=False));print(pd.DataFrame(matched).to_string(index=False))
if __name__=='__main__':main()
