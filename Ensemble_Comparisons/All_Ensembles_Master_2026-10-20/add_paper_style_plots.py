from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
O=Path(__file__).resolve().parent;R=O.parents[1];S=R/'outputs/Ensemble_Comparison_121_2026-10-15';D=O/'assets/paper_style';D.mkdir(exist_ok=True)
names=['ASTEROIDS_WT','ASTEROIDS_AA','REST_WT','bAIes_WT','bAIes_AA','IDPForge_WT','IDPForge_AA'];coverage=[];compact={}
cm=np.load(O/'assets/common121/contact_matrices_before_after.npz')
for name in names:
 stages=['prior','selected' if name.startswith('ASTEROIDS') else 'MaxEnt']
 for stage in stages:
  key=name+'__'+stage
  if (D/'plot_inputs.npz').exists():
   cache=np.load(D/'plot_inputs.npz');values=cache[key+'__values'];w=cache[key+'__weights']
  else:
   p=S/name/('prior_structural_cache.npz' if name.startswith('ASTEROIDS') and stage=='prior' else 'structural_cache.npz');z=np.load(p);values=z['values'];w=(pd.read_csv(S/name/'weights.csv').combined.to_numpy() if stage=='MaxEnt' else z['weights'] if 'weights' in z else np.ones(len(values))/len(values))
  assert np.isclose(w.sum(),1) and np.isfinite(values).all();compact[key+'__values']=values;compact[key+'__weights']=w
 for view in ['full','R2R3']:
  extent=[1.5,6.5] if view=='full' else [.9,2.5];rgcol,sacol=(0,2) if view=='full' else (1,3)
  fig,axs=plt.subplots(1,2,figsize=(10,4.8),sharex=True,sharey=True,layout='constrained')
  for ax,stage in zip(axs,stages):
   key=name+'__'+stage;values=compact[key+'__values'];w=compact[key+'__weights'];rg=values[:,rgcol]/10;sa=values[:,sacol]
   inside=(rg>=extent[0])&(rg<=extent[1])&(sa>=0)&(sa<=25);retained=float(w[inside].sum());assert retained>0
   hist,xe,ye=np.histogram2d(rg,sa,bins=30,range=[extent,[0,25]],weights=w,density=True)
   field=-.001987*300*np.log(np.flipud(hist)+1e-6)
   im=ax.imshow(field,interpolation='gaussian',extent=[ye[0],ye[-1],xe[0],xe[-1]],cmap='jet',aspect='auto',vmin=.1,vmax=3)
   ax.set(xlabel='Sα',xlim=(0,24.9),title=f"{'Before: prior' if stage=='prior' else 'After: '+stage}\n{retained:.1%} weight inside display range")
   fig.colorbar(im,ax=ax,shrink=.8,ticks=[1,2,3],label='−RT ln(density + 10⁻⁶)')
   coverage.append(dict(ensemble=name,stage=stage,view=view,retained_weight=retained,excluded_weight=1-retained,N=len(w)))
  axs[0].set_ylabel('Cα Rg (nm)');fig.suptitle(f'{name} | '+('full GP + AR330–448' if view=='full' else 'AR391–446')+'\nOriginal AR histogram/log display; not thermodynamic free energy',fontsize=11)
  for ext in ['png','pdf']:fig.savefig(D/f'{name}_{view}_SaRg.{ext}',dpi=170)
  plt.close(fig)
  fig,axs=plt.subplots(1,2,figsize=(10,4.8),layout='constrained')
  for ax,stage in zip(axs,stages):
   mat=cm[name+'__'+stage];mat=mat if view=='full' else mat[63:119,63:119];assert mat.shape==((121,121) if view=='full' else (56,56))
   sns.heatmap(mat,ax=ax,cmap='jet',vmin=0,vmax=1,cbar=True,cbar_kws={'shrink':.75,'label':'Contact probability'},rasterized=True);ax.invert_yaxis();ax.grid(alpha=.5)
   if view=='full':ticks=np.array([0,12,32,52,72,92,112,120]);labels=['GP','340','360','380','400','420','440','448']
   else:ticks=np.array([0,9,19,29,39,49,55]);labels=[str(391+i) for i in ticks]
   ax.set_xticks(ticks+.5,labels,rotation=45);ax.set_yticks(ticks+.5,labels);ax.set(xlabel='AR residue',ylabel='AR residue',title='Before: prior' if stage=='prior' else 'After: '+stage)
  fig.suptitle(f'{name} | {view} | closest-heavy <1.2 nm; adjacent pairs retained',fontsize=11)
  for ext in ['png','pdf']:fig.savefig(D/f'{name}_{view}_contacts.{ext}',dpi=170)
  plt.close(fig)
np.savez_compressed(D/'plot_inputs.npz',**compact);pd.DataFrame(coverage).to_csv(D/'display_range_coverage.csv',index=False)
(D/'definitions.json').write_text(json.dumps({'reference':'work/asteroids_reference/analysis.py: sa_rg and contacts','SaRg':'30x30 weighted density histogram; -0.001987*300*ln(density+1e-6); gaussian interpolation; jet; color limits 0.1-3','ranges':{'Sa':[0,25],'full_Rg_nm':[1.5,6.5],'R2R3_Rg_nm':[.9,2.5]},'normalization':'Density normalized within display window; retained/excluded full-ensemble weights explicitly reported. Not thermodynamic free energy.','contact':'Closest heavy atom <1.2 nm; diagonal zero; adjacent pairs retained; jet 0-1','region':'R2R3 AR391-446 inclusive; full current121 construct; cropped Sa includes 50 windows fully inside fragment','after':'Combined Kish MaxEnt for 5 new datasets; historical subset selection for ASTEROIDS'},indent=2))
p=O/'Tau5_All_Ensembles_Master.ipynb';nb=json.loads(p.read_text());nb['cells']=[c for c in nb['cells'] if not c.get('id','').startswith('paperstyle-')];extra=[]
def cell(kind,s):
 d=dict(cell_type=kind,metadata={},source=s,id='paperstyle-'+str(len(extra)))
 if kind=='code':d.update(execution_count=None,outputs=[])
 extra.append(d)
cell('markdown','## Paper-style Sα–Rg and contact maps: full chain and R2–R3\nThese reproduce the plotting conventions from the earlier AR Nature Communications analysis: 30×30 density histograms, −RT ln(density+10⁻⁶), Gaussian display interpolation, jet with 0.1–3 limits; contacts use closest-heavy <1.2 nm, linear jet 0–1, adjacent pairs retained.\n\nFull-chain views use the current 121-residue construct. R2–R3 crops are AR391–446, matching the previous paper comparison. This crop is distinct from the main/minor helix summary regions. Each panel reports the **ensemble weight retained within the original display bounds**; histograms normalize within that window. The transform is a display convention, not a thermodynamic free-energy estimate.\n\nGray/historical methods elsewhere keep their original conventions. Here every current ensemble receives separate before/after full and cropped panels; ASTEROIDS “after” is subset selection, the other five are combined Kish MaxEnt.')
cell('code',"display(pd.read_csv(assets/'paper_style/display_range_coverage.csv'))")
for name in names:
 cell('markdown','### '+name.replace('_',' '))
 cell('code',f"figures('paper_style',[{name+'_full_SaRg.png'!r},{name+'_R2R3_SaRg.png'!r},{name+'_full_contacts.png'!r},{name+'_R2R3_contacts.png'!r}])")
idx=next(i for i,c in enumerate(nb['cells']) if '## Full Jeffreys bAIes production' in ''.join(c['source']));nb['cells'][idx:idx]=extra;p.write_text(json.dumps(nb,indent=1))
print('Created 28 paired figures in PNG/PDF and inserted 16 notebook cells',flush=True)
