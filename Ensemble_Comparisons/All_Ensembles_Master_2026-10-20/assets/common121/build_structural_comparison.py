"""Curated common-construct structural comparison; include MaxEnt only after validated completion."""
from pathlib import Path
import json
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
O=Path(__file__).resolve().parent;P=O/'master_figures';P.mkdir(exist_ok=True)
rows=[];sets=[]
for name in ['ASTEROIDS_WT','ASTEROIDS_AA','REST_WT','bAIes_WT','bAIes_AA','IDPForge_WT','IDPForge_AA']:
 d=O/name
 if name.startswith('ASTEROIDS'):
  stages=[('prior',d/'prior_structural_cache.npz'),('selected',d/'structural_cache.npz')]
 else:stages=[('prior',d/'structural_cache.npz')]
 for stage,path in stages:
  if not path.exists():continue
  z=np.load(path);w=z['weights'] if 'weights' in z else np.full(len(z['values']),1/len(z['values']));sets.append((name,stage,z,w))
 if not name.startswith('ASTEROIDS') and (d/'fit_complete.json').exists():
  assert (d/'fit_validation.json').exists();w=pd.read_csv(d/'weights.csv').combined.to_numpy();assert np.isfinite(w).all() and np.isclose(w.sum(),1) and 1/(len(w)*np.sum(w*w))>=.1;sets.append((name,'MaxEnt',np.load(d/'structural_cache.npz'),w))
for name,stage,z,w in sets:
 row={'ensemble':name,'stage':stage,'N':len(w),'Neff':float(1/np.sum(w*w))}
 for j,col in enumerate(z['columns']):row[str(col)]=float(np.einsum('i,i->',w,z['values'][:,j]))
 for reg,lo,hi in [('R2_main',390,404),('R2_minor',404,410),('R3',430,440)]:row[reg+'_H_percent']=float(100*np.einsum('i,i->',w,z['helix'][:,lo-328:hi-327].mean(1)))
 rows.append(row)
pd.DataFrame(rows).to_csv(O/'structural_comparison.csv',index=False)
fig,axes=plt.subplots(4,2,figsize=(12,13));axes=axes.ravel();colors={'prior':'#777777','selected':'#1b9e77','MaxEnt':'#d95f02'}
for ax,name in zip(axes,dict.fromkeys(x[0] for x in sets)):
 for nn,stage,z,w in sets:
  if nn!=name:continue
  ax.hist(z['values'][:,0],bins=np.arange(10,65,1),weights=w,density=True,histtype='step',lw=1.7,label=stage,color=colors[stage])
 ax.set(title=name.replace('_',' '),xlabel='Cα Rg (Å)',ylabel='Probability density');ax.legend(frameon=False)
for ax in axes[len(set(x[0] for x in sets)):]:ax.set_visible(False)
fig.suptitle('121-residue constructs • ASTEROIDS selection or validated Kish MaxEnt',fontsize=13);fig.tight_layout();fig.savefig(P/'Rg_before_after.png',dpi=180);fig.savefig(P/'Rg_before_after.pdf');plt.close(fig)
# Export all contact matrices with explicit AR coordinate mapping, including tag residues.
contacts={};profiles=[]
for name,stage,z,w in sets:
 v=np.einsum('i,ij->j',w,z['contacts']);m=np.zeros((121,121));pairs=z['contact_pairs'];m[pairs[:,0],pairs[:,1]]=v;m+=m.T;contacts[name+'__'+stage]=m
 h=np.einsum('i,ij->j',w,z['helix'])*100
 for idx,value in enumerate(h):profiles.append({'ensemble':name,'stage':stage,'model_residue':idx+1,'AR_residue':idx+328 if idx>=2 else None,'DSSP_H_percent':value})
np.savez_compressed(O/'contact_matrices_before_after.npz',**contacts);pd.DataFrame(profiles).to_csv(O/'helicity_profiles_before_after.csv',index=False)
(O/'structural_comparison_status.json').write_text(json.dumps({'datasets':[{ 'ensemble':n,'stage':s} for n,s,z,w in sets],'contact_definition':'closest heavy atom distance <12 Å; diagonal zero; 121 residues GP+AR330448','Rg':'geometric C-alpha; not SAXS forward-model Rg','incomplete_MaxEnt':[n for n in ['REST_WT','bAIes_WT','bAIes_AA','IDPForge_WT','IDPForge_AA'] if not (O/n/'fit_complete.json').exists()]},indent=2))
print('Saved',len(rows),'ensemble stages',flush=True)
# Separate panels retain before/after distributions without obscuring changes.
from matplotlib.colors import LogNorm
for name in dict.fromkeys(x[0] for x in sets):
 subset=[x for x in sets if x[0]==name];fig,axes=plt.subplots(2,len(subset),figsize=(5*len(subset),8),squeeze=False)
 for j,(nn,stage,z,w) in enumerate(subset):
  hist,xe,ye=np.histogram2d(z['values'][:,0],z['values'][:,2],bins=[np.linspace(10,70,61),np.linspace(0,35,71)],weights=w)
  outside=1-hist.sum();assert outside<1e-8,(name,stage,'Sa/Rg plot excludes probability',outside)
  im=axes[0,j].pcolormesh(xe,ye,np.ma.masked_equal(hist.T,0),norm=LogNorm(1e-4,.1),cmap='viridis');axes[0,j].set(title=stage,xlabel='Cα Rg (Å)',ylabel='Sα (115 windows)');fig.colorbar(im,ax=axes[0,j],label='Probability per bin')
  m=contacts[nn+'__'+stage];im=axes[1,j].imshow(np.maximum(m,5e-5),origin='lower',extent=[327.5,448.5,327.5,448.5],norm=LogNorm(5e-5,1),cmap='magma');axes[1,j].set(xlabel='AR residue (328–329 represent GP tag)',ylabel='AR residue',title='Closest heavy atom <12 Å');fig.colorbar(im,ax=axes[1,j],label='Contact population')
 fig.suptitle(name.replace('_',' ')+' • 121 residues');fig.tight_layout();fig.savefig(P/(name+'_SaRg_contacts.png'),dpi=170);fig.savefig(P/(name+'_SaRg_contacts.pdf'));plt.close(fig)
