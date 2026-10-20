from pathlib import Path
import json,shutil,hashlib
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
O=Path(__file__).resolve().parent;R=O.parents[1];S=R/'outputs/bAIes_harmonic_2026-10-16';D=O/'assets/harmonic_production';D.mkdir(exist_ok=True)
for folder in ['comparison_analysis','final_region_analysis']:
 for f in (S/folder).iterdir():
  if f.is_file() and f.suffix in ['.csv','.json','.png','.pdf']:shutil.copy2(f,D/f.name)
(D/'status.json').write_text(json.dumps({'production':'WT and AA completed 2 us each; 20000 frames per sequence','structural_analysis':'complete','forward_calculations':'started 2026-10-19;1000 chronological frames per sequence','reweighting':'pending; not launched','native_construct':'120 residues GP+AR330-447','SAXS_construct':'121 with terminal Cys extension; same prior procedure','plot_comparison':'Current prior pools use121; harmonic full-chain metrics use native120. Crop AR391-446 matches exactly.'},indent=2))
z=np.load(O/'assets/paper_style/plot_inputs.npz');hel=pd.read_csv(O/'assets/common121/helicity_profiles_before_after.csv');colors={'ASTEROIDS':'#009e73','REST':'#cc79a7','bAIes':'#d55e00','IDPForge':'#0072b2','Harmonic':'#222222'};refs=json.loads((O/'assets/method_overlays/experimental_Rg_reference.json').read_text());names=sorted(set(k.split('__')[0] for k in z.files))
fig,axs=plt.subplots(2,2,figsize=(15,10),layout='constrained')
for row,v in enumerate(['WT','AA']):
 for name in names:
  if not name.endswith('_'+v):continue
  vals=z[name+'__prior__values'];w=z[name+'__prior__weights'];method=name.rsplit('_',1)[0];label=method+' (121)';color=colors[method]
  axs[row,0].hist(vals[:,0],bins=np.linspace(10,75,66),weights=w,density=True,histtype='step',color=color,label=label);axs[row,0].axvline(w@vals[:,0],color=color,ls=':',lw=.8)
  h=hel[(hel.ensemble==name)&(hel.stage=='prior')];axs[row,1].plot(h.AR_residue,h.DSSP_H_percent,color=color,label=method)
 t=pd.read_csv(D/f'{v}_timeseries.csv');h=pd.read_csv(D/f'{v}_residue_helicity.csv');axs[row,0].hist(t.Rg_CA_A,bins=np.linspace(10,75,66),density=True,histtype='step',color='k',lw=2,label='bAIes harmonic (native120)');axs[row,0].axvline(t.Rg_CA_A.mean(),ls=':',color='k');axs[row,0].axvline(refs[v]['Rg_A'],ls=':',color='grey',lw=2,label='Experimental SAXS Rg');axs[row,1].plot(h.AR_residue,h.DSSP_H*100,color='k',lw=2,label='bAIes harmonic')
 for lo,hi in [(390,404),(404,410),(430,440)]:axs[row,1].axvspan(lo,hi,color='grey',alpha=.08)
 axs[row,0].set(title=v+' • unweighted Rg',xlabel='CA Rg (Å)',ylabel='Probability density');axs[row,1].set(title=v+' • unweighted helicity',xlabel='AR residue',ylabel='DSSP H (%)',ylim=(0,100))
 for ax in axs[row]:ax.legend(fontsize=8)
fig.suptitle('All methods including completed harmonic bAIes: structural priors\nDotted lines: ensemble mean CA Rg and experimental SAXS Rg reference (distinct metrics). No harmonic reweighting yet.')
for ext in ['png','pdf']:fig.savefig(D/f'all_methods_with_harmonic.{ext}',dpi=170)
plt.close(fig)
# Prior Sa/Rg overlays include full and exact matching crop.
for view in ['full','R2R3']:
 fig,axs=plt.subplots(1,2,figsize=(14,6),layout='constrained')
 for ax,v in zip(axs,['WT','AA']):
  t=pd.read_csv(D/f'{v}_timeseries.csv');s=pd.read_csv(D/f'{v}_Sa_timeseries.csv');datasets=[]
  for name in names:
   if name.endswith('_'+v):
    vals=z[name+'__prior__values'];datasets.append((name.rsplit('_',1)[0],vals[:,0 if view=='full' else 1],vals[:,2 if view=='full' else 3],z[name+'__prior__weights']))
  datasets.append(('Harmonic',t.Rg_CA_A.to_numpy() if view=='full' else s.Rg_R2R3_CA_A.to_numpy(),s.Sa_full.to_numpy() if view=='full' else s.Sa_AR391_446.to_numpy(),np.ones(len(t))/len(t)))
  from scipy.ndimage import gaussian_filter
  for label,rg,sa,w in datasets:
   hist,xe,ye=np.histogram2d(sa,rg,bins=(np.linspace(0,60,101),np.linspace(0,90 if view=='full' else 60,101)),weights=w);assert np.isclose(hist.sum(),1);hist=gaussian_filter(hist,1.2);rank=np.sort(hist.ravel())[::-1];cut=rank[np.searchsorted(np.cumsum(rank),.5*rank.sum())];outer=rank[np.searchsorted(np.cumsum(rank),.75*rank.sum())];ax.contour((xe[:-1]+xe[1:])/2,(ye[:-1]+ye[1:])/2,hist.T,levels=[outer,cut],colors=[colors[label]],linestyles=['-',':'],linewidths=1.8);ax.plot([],[],color=colors[label],label=label)
  ax.set(title=v,xlabel='Sα',ylabel='CA Rg (Å)',xlim=(0,30),ylim=(10,70) if view=='full' else (5,35));ax.legend()
 fig.suptitle('Unweighted '+view+' • dotted 50% / solid 75% probability contours, common bins and smoothing\nHarmonic full chain native120; other pools121. R2R3 crop AR391–446 identical.')
 for ext in ['png','pdf']:fig.savefig(D/f'all_methods_SaRg_{view}.{ext}',dpi=170)
 plt.close(fig)
# Harmonic full/cropped contacts, same closest-heavy12A linear convention.
for view in ['full','R2R3']:
 fig,axs=plt.subplots(1,2,figsize=(12,5),layout='constrained')
 for ax,v in zip(axs,['WT','AA']):
  mat=np.loadtxt(D/f'{v}_contact_map.csv',delimiter=',');assert mat.shape==(120,120)
  if view=='R2R3':mat=mat[63:119,63:119];ticks=[0,9,19,29,39,49,55];labels=[391+i for i in ticks]
  else:ticks=[2,22,42,62,82,102,119];labels=[330,350,370,390,410,430,447]
  im=ax.imshow(mat,origin='lower',vmin=0,vmax=1,cmap='jet');ax.set_xticks(ticks,labels,rotation=45);ax.set_yticks(ticks,labels);ax.set(title=v+' harmonic • prior',xlabel='AR residue',ylabel='AR residue');fig.colorbar(im,ax=ax,label='Contact probability')
 fig.suptitle('Completed harmonic bAIes '+view+' • closest heavy <12 Å; adjacent pairs retained')
 for ext in ['png','pdf']:fig.savefig(D/f'harmonic_contacts_{view}.{ext}',dpi=170)
 plt.close(fig)
# Update notebook in place with a reproducible idempotent insertion, preserving existing chapters.
p=O/'Tau5_All_Ensembles_Master.ipynb';nb=json.loads(p.read_text());backup=O/'Tau5_All_Ensembles_Master_before_harmonic_2026-10-19.ipynb'
if not backup.exists():shutil.copy2(p,backup)
nb['cells']=[c for c in nb['cells'] if not c.get('id','').startswith('harmonic-production-')]
for c in nb['cells']:
 if c['cell_type']!='markdown':continue
 s=''.join(c['source']);s=s.replace('New harmonic production is ongoing.','New harmonic production and structural analysis are complete; native CS/SAXS calculations started October 19; harmonic reweighting is pending.').replace('Harmonic equilibrations and ongoing production','Harmonic equilibrations and completed production').replace('equilibrations complete; production ongoing','equilibrations and production complete').replace('2 μs additional production is ongoing.','2 μs additional production is complete. See the new structural comparison chapter above.').replace('No CS/SAXS reweighting exists for these new runs yet.','CS/SAXS backcalculations started October 19; no reweighting exists yet.')
 c['source']=s
cells=[]
def cell(kind,s):
 c={'cell_type':kind,'metadata':{},'source':s,'id':'harmonic-production-'+str(len(cells))}
 if kind=='code':c.update(execution_count=None,outputs=[])
 cells.append(c)
cell('markdown','## Completed harmonic bAIes — structural comparison added October 19\nBoth 2 μs productions are complete (20,000 frames each, following 20 ns equilibration). This chapter adds the full trajectories to the structural priors and compares them with Jeffreys bAIes. **Chemical shifts and SAXS are computing; harmonic reweighted curves and structures are unavailable.** Existing before/after experimental panels below still refer to the earlier Jeffreys bAIes.\n\nHarmonic native topology is GP+AR330–447 (120 residues); existing current forward-model pools have121. This difference is explicit for full-chain metrics. The AR391–446 crop is identical. Region summaries retain main R2 AR390–404, minor R2 AR404–410, R3 AR430–440. Persistent cis states and time-block results are included. No transfer of prior weights to the new trajectories.')
cell('code',"table('harmonic_production','likelihood_comparison.csv')\nfigures('harmonic_production',['all_methods_with_harmonic.png','all_methods_SaRg_full.png','all_methods_SaRg_R2R3.png'])")
cell('code',"figures('harmonic_production',['helicity_Rg_comparison.png','Sa_Rg_full.png','Sa_Rg_R2R3.png','harmonic_contacts_full.png','harmonic_contacts_R2R3.png','contact_comparison_linear.png','restraint_comparison.png','time_block_comparison.png','packing_comparison.png','cis_comparison.png'])")
cell('code',"display(json.loads((assets/'harmonic_production/status.json').read_text()))\ntable('harmonic_production','packing_comparison.csv')\ntable('harmonic_production','cis_well_change_summary.csv')")
nb['cells'][2:2]=cells;p.write_text(json.dumps(nb,indent=1))
# Replace stale old status file as well.
(O/'assets/harmonic/status.json').write_text((D/'status.json').read_text())
(D/'source_manifest.json').write_text(json.dumps({f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in D.iterdir() if f.is_file() and f.name!='source_manifest.json'},indent=2))
print('Harmonic production structure chapter and all-method prior overlays added',flush=True)
