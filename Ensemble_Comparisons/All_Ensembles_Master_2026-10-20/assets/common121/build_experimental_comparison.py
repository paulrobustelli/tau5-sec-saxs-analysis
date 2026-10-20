"""Common experimental observations; archived ASTEROIDS predictions and validated new fits."""
from pathlib import Path
import json
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
O=Path(__file__).resolve().parent;R=O.parents[1];P=O/'master_figures';P.mkdir(exist_ok=True)
exp=pd.read_csv(R/'outputs/IDPForge_D2D_MaxEnt_2026-10-15/experimental/chemical_shifts.csv');metrics=[];curves=[];shifts=[]
for name in ['ASTEROIDS_WT','ASTEROIDS_AA','REST_WT','bAIes_WT','bAIes_AA','IDPForge_WT','IDPForge_AA']:
 d=O/name;variant=name.rsplit('_',1)[1]
 if name.startswith('ASTEROIDS'):
  sf=pd.read_csv(d/'matched_saxs_fit.csv');sf=sf.rename(columns={'I_prior_dro5':'I_uniform','I_selected_dro5':'I_MaxEnt'});after='selected'
  records=[]
  for stage in ['prior','selected']:
   a=pd.read_csv(d/(stage+'_archived_predictions.csv'));a=a[a.dataset!='SAXS'].copy();a['AR_residue']=a.id.str.split(':').str[0].astype(int);a['atom']=a.id.str.split(':').str[1].replace({'HN':'H','NH':'H'});a=a[['AR_residue','atom','predicted']];assert not a.duplicated(['AR_residue','atom']).any()
   e=exp[exp.variant==variant].merge(a,on=['AR_residue','atom'],how='inner',validate='one_to_one');e['stage']=stage;e['ensemble']=name;records.append(e)
  cf=pd.concat(records,ignore_index=True);cf.to_csv(d/'matched_chemical_shift_fit.csv',index=False)
 else:
  if not (d/'fit_complete.json').exists():continue
  assert (d/'fit_validation.json').exists();sf=pd.read_csv(d/'saxs_fit.csv');after='MaxEnt';c=pd.read_csv(d/'chemical_shift_fit.csv');records=[]
  for stage,col in [('prior','unweighted_ppm'),('MaxEnt','reweighted_ppm')]:
   e=c.rename(columns={col:'predicted'}).copy();e['stage']=stage;e['ensemble']=name;records.append(e)
  cf=pd.concat(records,ignore_index=True)
 shifts.append(cf)
 for stage,col in [('prior','I_uniform'),(after,'I_MaxEnt')]:
  metrics.append({'ensemble':name,'stage':stage,'datatype':'SAXS','N_observations':len(sf),'metric':'chi2_mean_experimental_sigma','value':float(np.mean(((sf[col]-sf.I_exp)/sf.sigma_exp)**2))})
  for atom,g in cf[cf.stage==stage].groupby('atom'):metrics.append({'ensemble':name,'stage':stage,'datatype':atom,'N_observations':len(g),'metric':'RMSE_ppm','value':float(np.sqrt(np.mean((g.predicted-g.value)**2)))})
 fig,ax=plt.subplots(2,1,figsize=(7,7),sharex=True,gridspec_kw={'height_ratios':[3,1]});ax[0].errorbar(sf.q_Ainv,sf.I_exp,yerr=sf.sigma_exp,fmt='.',color='k',ms=3,label='experiment')
 for stage,col,color in [('prior','I_uniform','#777777'),(after,'I_MaxEnt','#d95f02')]:
  ax[0].plot(sf.q_Ainv,sf[col],color=color,label=stage);ax[1].plot(sf.q_Ainv,(sf[col]-sf.I_exp)/sf.sigma_exp,color=color)
 ax[0].set(yscale='log',ylabel='Intensity (experimental scale)',title=name.replace('_',' '));ax[0].legend(frameon=False);ax[1].axhline(0,color='k',lw=.6);ax[1].set(xlabel='q (Å⁻¹)',ylabel='Residual / σ');fig.tight_layout();fig.savefig(P/(name+'_SAXS_before_after.png'),dpi=180);fig.savefig(P/(name+'_SAXS_before_after.pdf'));plt.close(fig)
 sf['ensemble']=name;curves.append(sf)
pd.DataFrame(metrics).to_csv(O/'experimental_agreement.csv',index=False)
if shifts:pd.concat(shifts,ignore_index=True).to_csv(O/'matched_chemical_shift_predictions.csv',index=False)
if curves:pd.concat(curves,ignore_index=True).to_csv(O/'matched_saxs_curves.csv',index=False)
print('Comparison includes',len(curves),'ensembles',flush=True)
if shifts:
 allcs=pd.concat(shifts,ignore_index=True)
 for name,g in allcs.groupby('ensemble',sort=False):
  atoms=sorted(g.atom.unique());fig,axes=plt.subplots(len(atoms),1,figsize=(11,2*len(atoms)),squeeze=False)
  for ax,atom in zip(axes.ravel(),atoms):
   for stage,h in g[g.atom==atom].groupby('stage',sort=False):
    h=h.sort_values('AR_residue');ax.plot(h.AR_residue,h.predicted-h.value,'.-',ms=3,lw=.8,label=stage,color='#777777' if stage=='prior' else '#d95f02')
   ax.axhline(0,color='k',lw=.5);ax.set(ylabel=atom+' residual (ppm)',xlim=(330,446));ax.legend(frameon=False,loc='upper right',fontsize=8)
  axes[-1,0].set_xlabel('AR residue');fig.suptitle(name.replace('_',' ')+' • predicted minus experimental chemical shifts');fig.tight_layout();fig.savefig(P/(name+'_CS_before_after.png'),dpi=160);fig.savefig(P/(name+'_CS_before_after.pdf'));plt.close(fig)
 # Explicit coverage matrix prevents missing predictions from being treated as zero errors.
 coverage=allcs.groupby(['ensemble','stage','atom']).size().reset_index(name='N_observations');coverage.to_csv(O/'chemical_shift_comparison_coverage.csv',index=False)
