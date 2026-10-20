from pathlib import Path
import json
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
def gaussian_filter(a,sigma,mode="reflect"):
    x=np.arange(-5,6);k=np.exp(-0.5*(x/sigma)**2);k/=k.sum()
    for axis in [0,1]:a=np.apply_along_axis(lambda v:np.convolve(np.pad(v,5,mode="symmetric"),k,mode="valid"),axis,a)
    return a
O=Path(__file__).resolve().parent;A=O/'assets';D=A/'method_overlays';D.mkdir(exist_ok=True)
saxs=pd.read_csv(A/'common121/matched_saxs_curves.csv');cs=pd.read_csv(A/'common121/matched_chemical_shift_predictions.csv');hel=pd.read_csv(A/'common121/helicity_profiles_before_after.csv');z=np.load(A/'paper_style/plot_inputs.npz')
colors={'ASTEROIDS':'#009e73','REST':'#cc79a7','bAIes':'#d55e00','IDPForge':'#0072b2'};names=list(saxs.ensemble.unique());refs={};stats=[]
def method(n):return n.rsplit('_',1)[0]
def members(v):return [n for n in names if n.endswith('_'+v)]
def arrays(n,after):
 stage=('selected' if n.startswith('ASTEROIDS') else 'MaxEnt') if after else 'prior';key=n+'__'+stage;return z[key+'__values'],z[key+'__weights']
def setup(title):
 fig,ax=plt.subplots(2,2,figsize=(14,9),layout='constrained')
 for i,v in enumerate(['WT','AA']):
  for j in range(2):ax[i,j].set_title(v+' — '+('Reweighted / selected' if j else 'Unweighted prior'))
 fig.suptitle(title);return fig,ax
def save(fig,name):
 for ext in ['png','pdf']:fig.savefig(D/f'{name}.{ext}',dpi=180)
 plt.close(fig)
# Experimental Rg from exactly the curves used for fitting, with explicit common window.
for v in ['WT','AA']:
 g=saxs[saxs.ensemble==members(v)[0]].sort_values('q_Ainv')
 for n in members(v):assert np.allclose(saxs[saxs.ensemble==n].I_exp,g.I_exp) and np.allclose(saxs[saxs.ensemble==n].q_Ainv,g.q_Ainv)
 m=(g.q_Ainv>=.012)&(g.q_Ainv<=.035)&(g.I_exp>0);q=g.q_Ainv[m].to_numpy();y=np.log(g.I_exp[m]);sig=(g.sigma_exp/g.I_exp)[m].to_numpy();X=np.c_[np.ones(len(q)),q*q];cov=np.linalg.inv((X/sig[:,None]).T@(X/sig[:,None]));beta=cov@((X/sig[:,None]).T@(y/sig));rg=float(np.sqrt(-3*beta[1]));assert q[-1]*rg<1.3
 refs[v]=dict(Rg_A=rg,conditional_SE_A=float(1.5*np.sqrt(cov[1,1])/rg),qmin=float(q[0]),qmax=float(q[-1]),n=len(q),qmaxRg=float(q[-1]*rg),method='Weighted log-I Guinier fit to same experimental curve, common q window .012-.035 A^-1; reference is SAXS Rg, not geometric CA Rg')
(D/'experimental_Rg_reference.json').write_text(json.dumps(refs,indent=2))
fig,axs=setup('All methods: SAXS versus experiment — same forward protocol')
for i,v in enumerate(['WT','AA']):
 for j in range(2):
  ax=axs[i,j];g=saxs[saxs.ensemble==members(v)[0]];ax.errorbar(g.q_Ainv,g.I_exp,yerr=g.sigma_exp,fmt='.',ms=2,color='black',alpha=.6,label='Experiment',zorder=5)
  for n in members(v):
   h=saxs[saxs.ensemble==n];ax.plot(h.q_Ainv,h['I_MaxEnt' if j else 'I_uniform'],color=colors[method(n)],lw=1.5,label=method(n))
  ax.set(yscale='log',xlabel='q (Å⁻¹)',ylabel='Intensity (experimental scale)');ax.legend(fontsize=8)
save(fig,'SAXS_all_methods')
fig,axs=setup('All methods: Rg distributions — dotted lines show means and experimental SAXS reference')
for i,v in enumerate(['WT','AA']):
 for j in range(2):
  ax=axs[i,j]
  for n in members(v):
   vals,w=arrays(n,j);mean=float(w@vals[:,0]);ax.hist(vals[:,0],bins=np.arange(10,71,1),weights=w,density=True,histtype='step',lw=1.7,color=colors[method(n)],label=f'{method(n)}: mean {mean:.2f} Å');ax.axvline(mean,color=colors[method(n)],ls=':',lw=1.6);stats.append(dict(ensemble=n,after=bool(j),mean_CA_Rg_A=mean))
  ref=refs[v];ax.axvline(ref['Rg_A'],color='black',ls=':',lw=2.5,label=f"Exp. Guinier: {ref['Rg_A']:.2f} Å");ax.set(xlabel='Geometric Cα Rg (Å); black line = SAXS Rg reference',ylabel='Probability density');ax.legend(fontsize=8)
save(fig,'Rg_all_methods');pd.DataFrame(stats).to_csv(D/'Rg_means.csv',index=False)
fig,axs=setup('All methods: residue helicity (DSSP H)')
for i,v in enumerate(['WT','AA']):
 for j in range(2):
  ax=axs[i,j]
  for n in members(v):
   h=hel[(hel.ensemble==n)&((hel.stage!='prior') if j else (hel.stage=='prior'))].dropna(subset=['AR_residue']);ax.plot(h.AR_residue,h.DSSP_H_percent,color=colors[method(n)],label=method(n),lw=1.5)
  for a,b in [(390,404),(404,410),(430,440)]:ax.axvspan(a,b,color='gray',alpha=.08)
  ax.set(xlabel='AR residue',ylabel='DSSP H (%)',xlim=(330,448),ylim=(0,100));ax.legend(fontsize=8)
save(fig,'helicity_all_methods')
# Match CA observations across methods within each variant, no zero filling.
fig,axs=setup('All methods: Cα chemical-shift residuals (prediction − experiment)');coverage=[]
for i,v in enumerate(['WT','AA']):
 common=set.intersection(*[set(cs[(cs.ensemble==n)&(cs.atom=='CA')].AR_residue) for n in members(v)])
 for j in range(2):
  ax=axs[i,j];reference=None
  for n in members(v):
   h=cs[(cs.ensemble==n)&(cs.atom=='CA')&cs.AR_residue.isin(common)&((cs.stage!='prior') if j else (cs.stage=='prior'))].sort_values('AR_residue');assert not h.AR_residue.duplicated().any()
   if reference is None:reference=h.value.to_numpy()
   else:assert np.allclose(h.value,reference)
   r=h.predicted-h.value;ax.plot(h.AR_residue,r,'.-',ms=2.5,lw=1,color=colors[method(n)],label=f'{method(n)} (RMSE {np.sqrt(np.mean(r*r)):.3f})');coverage.append(dict(variant=v,ensemble=n,after=bool(j),n_common_CA=len(h)))
  ax.axhline(0,color='black',lw=.8);ax.set(xlabel='AR residue',ylabel='Cα residual (ppm)',xlim=(330,446));ax.legend(fontsize=8)
save(fig,'CA_residuals_all_methods');pd.DataFrame(coverage).to_csv(D/'CA_common_coverage.csv',index=False)
# Weighted density contours let all methods be visible on the same axes.
for view,rgcol,sacol in [('full',0,2),('R2R3',1,3)]:
 fig,axs=setup('All methods: Sα–Rg '+view+' — solid 50%, dashed 80% contours (1.2-bin Gaussian smoothing)')
 allv=np.concatenate([arrays(n,j)[0][:,[rgcol,sacol]] for n in names for j in [0,1]]);rb=np.linspace(0,np.ceil(allv[:,0].max()/5)*5,51);sb=np.linspace(0,np.ceil(allv[:,1].max()/5)*5,51)
 for i,v in enumerate(['WT','AA']):
  for j in range(2):
   ax=axs[i,j]
   for n in members(v):
    vals,w=arrays(n,j);hist,_,_=np.histogram2d(vals[:,sacol],vals[:,rgcol],bins=[sb,rb],weights=w);assert np.isclose(hist.sum(),1)
    hist=gaussian_filter(hist,1.2,mode="reflect");hist/=hist.sum()
    ordered=np.sort(hist.ravel())[::-1];cum=ordered.cumsum()
    for mass,ls in [(.8,'--'),(.5,'-')]:
     level=ordered[min(np.searchsorted(cum,mass),len(ordered)-1)]
     if level>0:ax.contour((sb[1:]+sb[:-1])/2,(rb[1:]+rb[:-1])/2,hist.T,levels=[level],colors=[colors[method(n)]],linestyles=[ls],linewidths=1.2)
    ax.scatter(w@vals[:,sacol],w@vals[:,rgcol],color=colors[method(n)],s=22,marker='x')
   ax.set(xlabel='Sα',ylabel='Cα Rg (Å)',xlim=(sb[0],sb[-1]),ylim=(rb[0],rb[-1]));ax.legend(handles=[Line2D([0],[0],color=colors[method(n)],label=method(n)) for n in members(v)],fontsize=8)
 save(fig,'SaRg_all_methods_'+view)
nbp=O/'Tau5_All_Ensembles_Master.ipynb';nb=json.loads(nbp.read_text());nb['cells']=[c for c in nb['cells'] if not c.get('id','').startswith('overlay-')];extra=[]
def add(kind,text):
 c=dict(cell_type=kind,source=text,metadata={},id='overlay-'+str(len(extra)))
 if kind=='code':c.update(execution_count=None,outputs=[])
 extra.append(c)
add('markdown','## All-method overlays: before versus after\nWT and AA are separated to retain their correct experiments. Each panel overlays all completed methods under the common121 protocol. “Unweighted” means uniform weights in the supplied pool, not unbiased molecular dynamics; bAIes pools already contain AF2 restraints. SAXS prior curves retain their fitted ensemble intensity scale. ASTEROIDS after is selected, other after is combined Kish MaxEnt. Historical incompatible protocols remain in their separate sections.\n\nRg means and the experimental reference use dotted lines. Experimental Rg is a weighted Guinier fit to the exact experimental SAXS curve used here over q=0.012–0.035 Å⁻¹, with qmaxRg<1.3; this is explicitly distinct from geometric Cα Rg and from earlier P(r)/other-window estimates. Cα residuals use the common observed residue intersection within each sequence. Sα–Rg overlays show 50%/80% weighted histogram probability contours (solid/dashed; 1.2-bin Gaussian smoothing for display), x marks means, and include the full data range; individual paper-style heatmaps remain below.')
for name in ['SAXS_all_methods','Rg_all_methods','helicity_all_methods','SaRg_all_methods_full','SaRg_all_methods_R2R3','CA_residuals_all_methods']:add('code',f"figures('method_overlays',['{name}.png'])")
add('code',"display(json.loads((assets/'method_overlays/experimental_Rg_reference.json').read_text()))\ntable('method_overlays','CA_common_coverage.csv')")
idx=next(i for i,c in enumerate(nb['cells']) if '## Full Jeffreys bAIes production' in ''.join(c['source']));# place earlier, immediately after setup
idx=2
nb['cells'][idx:idx]=extra;nbp.write_text(json.dumps(nb,indent=1));print('Six overlay figures and notebook section complete')
