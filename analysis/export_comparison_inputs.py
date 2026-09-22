from pathlib import Path
import json,sys,numpy as np,nbformat
R=Path(__file__).resolve().parents[1]; sys.path.insert(0,str(R/'analysis'))
from inspect_data import load
from component_recovery import SPECS,transform
O=R/'results_comparison';O.mkdir(exist_ok=True);(O/'curves').mkdir(exist_ok=True)
registry=[]
def add(path,group,label,sample,error='Exported pointwise uncertainties; unknown q covariance',**kw):
 d=np.loadtxt(R/path);assert d.ndim==2 and d.shape[1]>=3
 registry.append(dict(id=str(path).replace('/','__').replace('.dat',''),path=str(path),group=group,label=label,sample=sample,error_model=error,**kw))
# Export the original averages for a pure I(q)-reading notebook.
for sample,lo,hi in [('WT_092025',175,185),('AA_052026',172,182)]:
 a=np.array([np.loadtxt(R/'BNL_SAXS/dat_files'/sample/f'profile_{j}.dat') for j in range(lo,hi+1)])
 p=f'results_comparison/curves/{sample}_Natalie_average.dat';np.savetxt(R/p,np.c_[a[0,:,0],a[:,:,1].mean(0),np.sqrt((a[:,:,2]**2).sum(0))/len(a)],header='q I sigma_original_independent_frame_assumption')
 add(p,'Original averages',f'{sample[:2]} Natalie frames {lo}–{hi}',sample,'Original independent-frame error propagation; shared baseline uncertainty unavailable')
for folder,pattern,group in [('poster_results','*_selected_average.dat','Selected non-EFA'),('independent_baseline','*_average.dat','Baseline variants'),('component_results','*.dat','Audited EFA and matched average'),('results','*rank1_profile.dat','Legacy rank-one'),('results','*rotation_component*.dat','Legacy EFA')]:
 for p in sorted((R/folder).glob(pattern)):
  if 'Pr' in p.name:continue
  sample='WT_092025' if p.name.startswith('WT') else 'AA_052026'
  error='Exported conditional errors; model and baseline choice excluded'
  if 'candidate_main' in p.name:error='Pointwise EFA bootstrap SD; paired bootstrap curves retained separately'
  add(p.relative_to(R),group,p.stem,sample,error)
for sample in SPECS:
 folder=sample[:2]+'_window_comparison'
 for p in sorted((R/folder).glob('*.dat')):add(p.relative_to(R),'Non-EFA window and baseline',sample[:2]+' '+p.stem,sample,'Shared baseline propagated within each curve; paired covariance exported separately')
# Distinct non-EFA windows previously fitted by BIFT.
for p in sorted((R/'poster_results').glob('*_bift_window*.npz')):
 a=np.load(p);path='results_comparison/curves/'+p.stem+'_Iq.dat'
 np.savetxt(R/path,np.c_[a['q'],a['I'],a['error']],header='q I sigma_saved_BIFT_input')
 add(path,'Non-EFA BIFT windows',p.stem+' input I(q)',p.name[:9],'Saved selected-window I(q) errors; q covariance unavailable')

# Report-hypothesis rotation: derive conditional profile errors with fixed fitted C.
for p in sorted((R/'results').glob('*report_rotation_*.npz')):
 a=np.load(p);L=np.linalg.pinv(a['C'].T);se=np.sqrt(a['error']**2 @ (L**2))
 for k in range(a['S'].shape[1]):
  out=f'results_comparison/curves/{p.stem}_component{k+1}.dat';np.savetxt(R/out,np.c_[a['q'],a['S'][:,k],se[:,k]],header='q I conditional_sigma_fixed_C_no_cross_frame_covariance')
  add(out,'Report-hypothesis EFA',p.stem+f' component {k+1}',p.name[:9], 'Fixed-C conditional propagation; report-frame covariance unavailable; exploratory rotation')
# Select leading edge WITHOUT optimizing Rg: five consecutive frames before main peak,
# >=10% main peak coefficient, minimum max low-q/wide-q absolute component contribution ratio.
edge_results={}
for sample,spec in SPECS.items():
 q,Y,E,frames,A,I,err=transform(sample)
 f=np.load(R/'component_results'/f'{sample}_fiducial.npz');C=f['C'];S=f['S'];peak=int(np.argmax(C[:,1]))
 ratios=[]
 for lo,hi in [(.012,.030),(.012,.15)]:
  m=(q>=lo-1e-9)&(q<=hi+1e-9);strength=np.trapezoid(np.abs(S[m]),q[m],axis=0)
  z=np.abs(C)*strength[None,:];ratios.append(z[:,0]/np.maximum(z.sum(1),1e-100))
 ratio=np.maximum(*ratios); eligible=(np.arange(len(frames))<peak)&(C[:,1]>=.1*C[:,1].max())
 runs=[]
 for j in range(peak-3):
  ix=np.arange(j,j+5)
  if np.all(eligible[ix]):runs.append((float(np.max(ratio[ix])),int(j)))
 score,j=min(runs); chosen=frames[j:j+5]; windows={'least_overlap_leading':(int(chosen[0]),int(chosen[-1])),'main_apex_5':(int(frames[peak]-2),int(frames[peak]+2))}
 # A strict pre-second-onset window uses a descriptive 5% threshold and may not exist.
 isolated=bool(score<=.05)
 edge_results[sample]=dict(windows=windows,max_other_fraction=score,passes_5percent_rule=isolated,rule='Five consecutive leading frames; main >=10% peak; minimize worst absolute other fraction across q=.012-.030 and .012-.15; 5% descriptive threshold, not biological purity',frames=frames.tolist(),other_fraction=ratio.tolist(),C=C.tolist())
 for baseline in ['linear','pre','post']:
  qt,Y,E,fr,A,I,err=transform(sample,baseline);weights=[]
  for name,(lo,hi) in windows.items():
   w=A[:,(fr>=lo)&(fr<=hi)].mean(1);y=Y@w;e=np.sqrt(E**2@(w*w));weights.append(w)
   p=f'results_comparison/curves/{sample}_{baseline}_{name}_{lo}-{hi}.dat';np.savetxt(R/p,np.c_[q,y,e],header='q I sigma_shared_baseline')
   add(p,'EFA-guided direct edges',f'{sample[:2]} {baseline} {name} {lo}–{hi}',sample,'Exact fixed-baseline linear propagation; edge selection uncertainty excluded',frames=[lo,hi])
  W=np.array(weights).T;cov=np.einsum('qt,ti,tj->qij',E**2,W,W)
  np.savez_compressed(O/f'{sample}_{baseline}_edge_covariance.npz',q=q,covariance=cov)
 # Original four-window covariance for paired comparisons.
 a=np.load(R/(sample[:2]+'_window_comparison')/'profiles.npz')
 for baseline in ['linear','pre','post']:
  W=a[baseline+'_weights'];cov=np.einsum('qt,ti,tj->qij',E**2,W,W)
  np.savez_compressed(O/f'{sample}_{baseline}_window_covariance.npz',q=q,covariance=cov)
(O/'curve_registry.json').write_text(json.dumps(registry,indent=2));(O/'edge_selection.json').write_text(json.dumps(edge_results,indent=2))
print('Registered',len(registry),'curves');print({k:{a:v[a] for a in ['windows','max_other_fraction','passes_5percent_rule']} for k,v in edge_results.items()})
