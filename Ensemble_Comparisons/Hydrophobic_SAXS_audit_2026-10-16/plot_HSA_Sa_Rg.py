from pathlib import Path
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter
from matplotlib.lines import Line2D
O=Path(__file__).resolve().parent;C=O.parent/'Ensemble_Comparison_121_2026-10-15'
colors={'ASTEROIDS':'#009e73','REST':'#cc79a7','bAIes':'#d55e00','IDPForge':'#0072b2'}
frames=[]
for name in ['ASTEROIDS_WT','REST_WT','IDPForge_WT','bAIes_WT','ASTEROIDS_AA','IDPForge_AA','bAIes_AA']:
 t=pd.read_csv(O/name/'metrics.csv');w=pd.read_csv(O/name/'weights.csv');z=np.load(C/name/('prior_structural_cache.npz' if name.startswith('ASTEROIDS') else 'structural_cache.npz'),allow_pickle=True)
 ids=z['conformer_ids'].astype(str) if name.startswith('ASTEROIDS') else [Path(x).stem for x in z['files'].astype(str)]
 assert list(t.conformer_id)==list(ids)
 assert np.allclose(t.Rg_CA_A,z['values'][:,0],atol=.002)
 t['Sa_full']=z['values'][:,2];t['ensemble']=name;t['prior_weight']=w.uniform;t['after_weight']=w.combined
 frames.append(t[['ensemble','conformer_id','Rg_CA_A','Sa_full','buried_HSA_nonlocal_A2','prior_weight','after_weight']])
all=pd.concat(frames,ignore_index=True);all.to_csv(O/'HSA_Sa_Rg_plot_data.csv',index=False)
for xcol,xlabel,tag in [('Sa_full','Sα (full chain)','Sa'),('Rg_CA_A','Cα Rg (Å)','Rg')]:
 fig,axs=plt.subplots(2,2,figsize=(14,10),layout='constrained')
 xb=np.linspace(0,np.ceil(all[xcol].max()/5)*5,61);yb=np.linspace(0,np.ceil(all.buried_HSA_nonlocal_A2.max()/100)*100,61)
 for i,v in enumerate(['WT','AA']):
  for j,stage in enumerate(['prior','after']):
   ax=axs[i,j];handles=[]
   for name,t in all.groupby('ensemble'):
    if not name.endswith('_'+v):continue
    method=name.rsplit('_',1)[0];w=t[stage+'_weight'].to_numpy();w=w/w.sum()
    h,_,_=np.histogram2d(t[xcol],t.buried_HSA_nonlocal_A2,bins=[xb,yb],weights=w);assert np.isclose(h.sum(),1)
    h=gaussian_filter(h,1.2);h/=h.sum();rank=np.sort(h.ravel())[::-1];cum=rank.cumsum()
    for mass,ls in [(.75,'-'),(.5,':')]:
     level=rank[np.searchsorted(cum,mass)];ax.contour((xb[1:]+xb[:-1])/2,(yb[1:]+yb[:-1])/2,h.T,levels=[level],colors=[colors[method]],linestyles=[ls],linewidths=1.7)
    handles.append(Line2D([],[],color=colors[method],label=method))
   ax.set(title=v+' — '+('Unweighted' if stage=='prior' else 'Reweighted / ASTEROIDS selected'),xlabel=xlabel,ylabel='Nonlocal buried HSA (Å²)',xlim=(xb[0],xb[-1]),ylim=(yb[0],yb[-1]));ax.legend(handles=handles,fontsize=9)
 fig.suptitle('Buried hydrophobic surface area versus '+xlabel+'\nDotted: 50%; solid: 75% probability contours • common bins, 1.2-bin Gaussian smoothing\nNonlocal burial relative to local ±6-residue environment; common121 ensembles')
 for ext in ['png','pdf']:fig.savefig(O/f'buried_HSA_vs_{tag}.{ext}',dpi=170)
 plt.close(fig)
print('Saved matched data and two before/after contour figures')
