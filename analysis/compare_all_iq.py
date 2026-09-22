"""Read exported I(q); compare conventional/extended Guinier without rerunning EFA.
All Monte Carlo intervals are conditional on fixed fit ranges and exported errors.
"""
from pathlib import Path
import sys,json,hashlib
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams["figure.dpi"]=85
from scipy.optimize import least_squares
from IPython.display import display,Markdown
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'extended_guinier'))
from extended_guinier import fit,scan,shape,rg_from_nu
OUT=ROOT/'results_comparison'

def conventional(q,y,e,lo=.012,hi=.030):
 m=(q>=lo-1e-9)&(q<=hi+1e-9)&np.isfinite(y)&np.isfinite(e)&(e>0)
 if m.sum()<8 or np.any(y[m]<=0) or np.sum(y[m]>3*e[m])<6:return None
 X=np.c_[np.ones(m.sum()),q[m]**2];sig=e[m]/y[m];cov=np.linalg.inv(X.T@(X/sig[:,None]**2));b=cov@(X.T@(np.log(y[m])/sig**2))
 if b[1]>=0:return None
 rg=np.sqrt(-3*b[1]);return dict(Rg=float(rg),I0=float(np.exp(b[0])),Rg_conditional_se=float(3*np.sqrt(cov[1,1])/(2*rg)),qmin=float(q[m][0]),qmax=float(q[m][-1]),qRgmax=float(q[m][-1]*rg),chi2_reduced=float(np.sum(((np.log(y[m])-X@b)/sig)**2)/(m.sum()-2)),n=int(m.sum()))

def standard_scan(q,y,e,lo=.012):
 rows=[conventional(q,y,e,lo,hi) for hi in q[(q>=lo+.007-1e-9)&(q<=.08)]]
 valid=[r for r in rows if r and r['qRgmax']<=1.1]
 return max(valid,key=lambda r:r['qmax']) if valid else None

def extended_fixed(q,y,e,reference):
 if reference is None:return None
 m=(q>=reference['qmin']-1e-9)&(q<=reference['qmax']+1e-9)
 r=fit(q[m],y[m],e[m]);return r if r and not r['boundary'] and r['qRgmax']<=2 else None

def quick_extended(q,y,e,reference):
 if reference is None:return None
 m=(q>=reference['qmin']-1e-9)&(q<=reference['qmax']+1e-9);q=q[m];y=y[m];e=e[m]
 def residual(x):return (y-np.exp(x[0])*shape(q,x[1]))/e
 r=least_squares(residual,[np.log(reference['I0']),reference['nu']],bounds=([-100,.2],[100,.8]),max_nfev=60)
 return float(rg_from_nu(r.x[1])) if r.success and .201<r.x[1]<.799 else None

def mc(q,y,e,g,x,seed=20260922,n=200,replicates=None):
 rng=np.random.default_rng(seed);draws=y+e*rng.normal(size=(n,len(q))) if replicates is None else replicates
 pairs=[]
 for b in draws:
  gr=conventional(q,b,e,g['qmin'],g['qmax']) if g else None
  xr=quick_extended(q,b,e,x)
  pairs.append([gr['Rg'] if gr else np.nan,xr if xr else np.nan])
 a=np.array(pairs);out={'n':len(a),'valid':np.isfinite(a).sum(0).tolist(),'interval_type':'Conditional 95% noise/replicate interval; ranges fixed; no baseline-model uncertainty'}
 for j,k in enumerate(['conventional95','extended95']):out[k]=np.nanpercentile(a[:,j],[2.5,97.5]).tolist() if np.isfinite(a[:,j]).sum()>=20 else None
 ok=np.all(np.isfinite(a),axis=1);out['extended_minus_conventional95']=np.percentile(a[ok,1]-a[ok,0],[2.5,97.5]).tolist() if ok.sum()>=20 else None
 return out

def analyze(entry,n=200):
 p=ROOT/entry['path'];d=np.loadtxt(p);q,y,e=d[:,:3].T
 assert np.all(np.diff(q)>0) and np.all(np.isfinite(d[:,:3])) and np.all(e>0),p
 g=standard_scan(q,y,e);x,scanrows=scan(q,y,e);same=extended_fixed(q,y,e,g)
 fixed=conventional(q,y,e);fixed_valid=fixed and fixed['qRgmax']<=1.1
 replicates=None
 if entry['group']=='Audited EFA and matched average':
  a=np.load(ROOT/'component_results'/f"{entry['sample']}_audit.npz")
  if 'candidate_main' in entry['path']:replicates=a['bootstrap']
  elif 'frame_average_same' in entry['path']:replicates=a['bootstrap_average']
 sensitivity={str(lo):{'conventional':standard_scan(q,y,e,lo),'extended':scan(q,y,e,lo)[0]} for lo in [.010,.015,.020]}
 out=dict(**entry,sha256=hashlib.sha256(p.read_bytes()).hexdigest(),conventional=g,extended=x,extended_same_range=same,legacy_fixed_012_030=fixed,legacy_fixed_within_limit=bool(fixed_valid),uncertainty=mc(q,y,e,g,x,n=n,replicates=replicates),qmin_sensitivity=sensitivity,negative_fraction=float(np.mean(y<0)),conditional_replicates='Existing paired EFA perturbations' if replicates is not None else 'Independent q Gaussian perturbations')
 out['warnings']=[]
 if 'EFA' in entry['group']:out['warnings'].append('Component identity and physical validity are not established by a converged fit.')
 if g and g['chi2_reduced']>2:out['warnings'].append('Conventional residual chi-square/dof >2.')
 if x and x['chi2_reduced']>2:out['warnings'].append('Extended residual chi-square/dof >2.')
 if x is None:out['warnings'].append('No supported extended fit; do not interpret as a species Rg.')
 return out

def panel(entry,result):
 q,y,e=np.loadtxt(ROOT/entry['path'])[:,:3].T;fig,grid=plt.subplots(2,2,figsize=(11,7));ax=[grid[0,0],grid[1,0],grid[1,1],grid[0,1]]
 m=(q>=.006)&(q<=.25);ax[0].errorbar(q[m],y[m],e[m],fmt='.',ms=2,elinewidth=.4,color='0.35')
 if np.any(y[m]<=0):ax[0].set_yscale('symlog',linthresh=max(float(np.median(e[m])),1e-8));ax[0].axhline(0,color='k',lw=.5)
 else:ax[0].set_yscale('log')
 ax[0].set(xlabel='q (Å⁻¹)',ylabel='I(q), exported scale',title='I(q): all signed points retained')
 m=(q>=.010)&(q<=.085);ax[1].errorbar(q[m]**2,y[m],e[m],fmt='.',ms=2,elinewidth=.4,color='0.55')
 for key,color in [('conventional','C0'),('extended','C1')]:
  r=result[key]
  if not r:continue
  m=(q>=r['qmin']-1e-9)&(q<=r['qmax']+1e-9);pred=r['I0']*(np.exp(-q[m]**2*r['Rg']**2/3) if key=='conventional' else shape(q[m],r['nu']))
  ax[1].plot(q[m]**2,pred,color=color,label=f"{key}: {r['Rg']:.2f} Å\nq={r['qmin']:.3f}–{r['qmax']:.3f}; qRg={r['qRgmax']:.2f}")
  ax[2].plot(q[m],(y[m]-pred)/e[m],'.-',lw=.5,ms=2,color=color,label=key)
 ax[1].set(xlabel='q² (Å⁻²)',ylabel='I(q)',title='Signed intensity: conventional / extended')
 m=(q>=.010)&(q<=.085)&(y>0)
 ax[3].errorbar(q[m]**2,np.log(y[m]),e[m]/y[m],fmt='.',ms=2,elinewidth=.4,color='0.55')
 for key,color in [('conventional','C0'),('extended','C1')]:
  r=result[key]
  if r:
   m=(q>=r['qmin']-1e-9)&(q<=r['qmax']+1e-9);logpred=np.log(r['I0'])+(-q[m]**2*r['Rg']**2/3 if key=='conventional' else np.log(shape(q[m],r['nu'])))
   ax[3].plot(q[m]**2,logpred,color=color,label=key)
 ax[3].set(xlabel='q² (Å⁻²)',ylabel='ln I(q)',title='Guinier plot (positive data shown only)')
 if result['conventional'] or result['extended']:ax[3].legend(fontsize=8)
 if result['conventional'] or result['extended']:ax[1].legend(fontsize=7);ax[2].legend(fontsize=7)
 else:ax[1].text(.1,.8,'No supported fits',transform=ax[1].transAxes)
 ax[2].axhline(0,color='k',lw=.5);ax[2].set(xlabel='q (Å⁻¹)',ylabel='(I − model) / σ',title='Intensity-space residuals')
 fig.suptitle(entry['label'],fontsize=10);fig.tight_layout();return fig

def table(rows):
 def val(r):return 'unsupported' if r is None else f"{r['Rg']:.2f} ± {r['Rg_conditional_se']:.2f}"
 h='| Curve | Conventional Rg ± SE (Å) | Extended, same q range | Extended, qRg≤2 | Conditional extended 95% |\n|---|---:|---:|---:|---|\n'
 for r in rows:
  ci=r['uncertainty']['extended95'];h+=f"|{r['label']}|{val(r['conventional'])}|{val(r['extended_same_range'])}|{val(r['extended'])}|{str(np.round(ci,2).tolist()) if ci else 'unsupported'}|\n"
 display(Markdown(h))

def show_group(groups,n=200):
 registry=json.loads((OUT/'curve_registry.json').read_text());results=[]
 for entry in registry:
  if entry['group'] not in groups:continue
  r=analyze(entry,n=n);results.append(r);fig=panel(entry,r);display(fig);plt.close(fig)
  if r['warnings']:display(Markdown(' '.join(r['warnings'])))
 table(results);return results

def edge_plot():
 data=json.loads((OUT/'edge_selection.json').read_text());fig,axes=plt.subplots(2,2,figsize=(11,7))
 for col,(sample,d) in enumerate(data.items()):
  fr=np.array(d['frames']);C=np.array(d['C']);axes[0,col].plot(fr,C[:,1]/C[:,1].max(),label='Main / main maximum');axes[0,col].plot(fr,d['other_fraction'],label='Other absolute contribution fraction')
  axes[0,col].axhline(.05,color='k',ls=':',label='5% descriptive threshold')
  for label,(lo,hi) in d['windows'].items():axes[0,col].axvspan(lo,hi,alpha=.15,label=label)
  axes[0,col].set(title=sample,xlabel='Frame (zero-based)',ylabel='Relative contribution');axes[0,col].legend(fontsize=7)
  registry=json.loads((OUT/'curve_registry.json').read_text())
  for ent in registry:
   if ent['sample']==sample and ent['group']=='EFA-guided direct edges' and ' linear ' in ent['label']:
    q,y,e=np.loadtxt(ROOT/ent['path'])[:,:3].T;m=(q>=.012)&(q<=.15);norm=np.trapezoid(y[m],q[m]);axes[1,col].errorbar(q[m],y[m]/norm,e[m]/norm,fmt='.',ms=3,label=ent['label'])
  a=np.load(ROOT/'component_results'/f'{sample}_audit.npz');q=a['q'];m=(q>=.012)&(q<=.15);y=a['candidate'];axes[1,col].plot(q[m],y[m]/np.trapezoid(y[m],q[m]),label='EFA main (shape only)')
  axes[1,col].set(xlabel='q (Å⁻¹)',ylabel='Area-normalized I(q)',yscale='log');axes[1,col].legend(fontsize=7)
 fig.tight_layout();display(fig);plt.close(fig)
 return data

def paired_windows(n=300):
 """Draw correlated curve errors using per-q cross-window covariance."""
 rng=np.random.default_rng(20260923);registry=json.loads((OUT/'curve_registry.json').read_text());out=[]
 for sample in ['WT_092025','AA_052026']:
  for baseline in ['linear','pre','post']:
   for kind in ['edge','window']:
    c=np.load(OUT/f'{sample}_{baseline}_{kind}_covariance.npz');q=c['q'];cov=c['covariance']
    if kind=='edge':ents=[r for r in registry if r['sample']==sample and r['group']=='EFA-guided direct edges' and f' {baseline} ' in r['label']]
    else:
     names=['leading','center','trailing','late_tail'];ents=[next(r for r in registry if r['path']==f'{sample[:2]}_window_comparison/{baseline}_{name}.dat') for name in names]
    a=np.stack([np.loadtxt(ROOT/r['path'])[:,1:3] for r in ents],axis=1);Y=a[:,:,0];E=a[:,:,1];refs=[]
    for j in range(Y.shape[1]):
     m=(q>=.012-1e-9)&(q<=.055+1e-9);r=fit(q[m],Y[m,j],E[m,j]);refs.append(r if r and not r['boundary'] and r['qRgmax']<=2 else None)
    ev,U=np.linalg.eigh(cov);L=U*np.sqrt(np.maximum(ev,0))[:,None,:];values=[]
    for rep in range(n):
     noise=np.einsum('qij,qj->qi',L,rng.normal(size=Y.shape));values.append([quick_extended(q,Y[:,j]+noise[:,j],E[:,j],refs[j]) for j in range(Y.shape[1])])
    values=np.array(values,dtype=float)
    for j in range(Y.shape[1]):
     if j==1:continue
     valid=np.all(np.isfinite(values[:,[j,1]]),axis=1);ci=np.percentile(values[valid,j]-values[valid,1],[2.5,97.5]).tolist() if valid.sum()>=20 else None
     out.append(dict(sample=sample,baseline=baseline,kind=kind,comparison=ents[j]['label']+' minus '+ents[1]['label'],extended_delta95=ci,valid_draws=int(valid.sum()),n=n,fit_range=[.012,.055],conditional=True))
 return out

def run_all(n=200):
 registry=json.loads((OUT/'curve_registry.json').read_text());results=[]
 for i,e in enumerate(registry):
  results.append(analyze(e,n));print(f"{i+1}/{len(registry)} {e['label']}",flush=True)
 (OUT/'fit_results.json').write_text(json.dumps(results,indent=2));return results

def gallery(results,groups):
 rows=[r for r in results if r['group'] in groups]
 for r in rows:
  fig=panel(r,r);display(fig);plt.close(fig)
  if r['warnings']:display(Markdown(' '.join(r['warnings'])))
 table(rows)

def summary_plot(results):
 groups=['Original averages','Selected non-EFA','Audited EFA and matched average','EFA-guided direct edges']
 rows=[r for r in results if r['group'] in groups and (r['group']!='EFA-guided direct edges' or ' linear ' in r['label'])]
 fig,ax=plt.subplots(figsize=(11,max(5,len(rows)*.34)))
 for j,r in enumerate(rows):
  for delta,key,color in [(-.12,'conventional','C0'),(.12,'extended','C1')]:
   f=r[key]
   if f:ax.errorbar(f['Rg'],j+delta,xerr=f['Rg_conditional_se'],fmt='o',ms=4,color=color,label=key if j==0 else None)
 ax.set(yticks=range(len(rows)),yticklabels=[r['label'].replace('WT_092025','WT').replace('AA_052026','AA') for r in rows],xlabel='Rg (Å); bars = conditional fit SE',title='Conventional and extended fits on the same exported curves');ax.invert_yaxis();ax.legend();fig.tight_layout();display(fig);fig.savefig(OUT/'Rg_overview.png',dpi=140);plt.close(fig)

def paired_efa(results):
 out=[]
 for sample in ['WT_092025','AA_052026']:
  a=np.load(ROOT/'component_results'/f'{sample}_audit.npz');q=a['q'];vals=[]
  m=(q>=.012-1e-9)&(q<=.055+1e-9)
  ref=[fit(q[m],a[k][m],a[ek][m]) for k,ek in [('average','average_error'),('candidate','candidate_bootstrap_error')]]
  for avg,main in zip(a['bootstrap_average'],a['bootstrap']):
   v=[quick_extended(q,avg,a['average_error'],ref[0]),quick_extended(q,main,a['candidate_bootstrap_error'],ref[1])];vals.append(v)
  v=np.array(vals,dtype=float);ok=np.all(np.isfinite(v),axis=1)
  out.append(dict(sample=sample,comparison='EFA main minus matched average',extended_delta95=np.percentile(v[ok,1]-v[ok,0],[2.5,97.5]).tolist(),valid_draws=int(ok.sum()),qrange=[.012,.055],conditional='Existing paired perturbations refit EFA; baseline/support choice fixed'))
 return out

def pr_check():
 rows=[];fig,ax=plt.subplots(1,2,figsize=(10,3.5))
 for sample,color in [('WT_092025','C0'),('AA_052026','C1')]:
  for folder,kind in [('poster_results','primary'),('component_results','average'),('component_results','candidate')]:
   p=ROOT/folder/(f'{sample}_bift_primary.npz' if folder=='poster_results' else f'{sample}_{kind}_bift.npz')
   a=np.load(p);r=a['r'];pr=a['p'];area=np.trapezoid(pr,r);rg=float(np.sqrt(np.trapezoid(pr*r*r,r)/(2*area)))
   rows.append(dict(sample=sample,curve=folder+'/'+kind,Rg_from_Pr=rg,uncertainty='No new Pr confidence interval: full regularized P(r) covariance unavailable; see original BIFT audit'))
   ax[0 if sample.startswith('WT') else 1].plot(r,pr/area,label=kind if folder=='component_results' else 'selected average')
 for a,s in zip(ax,['WT','AA']):a.set(title=s,xlabel='r (Å)',ylabel='Area-normalized P(r)');a.legend(fontsize=8)
 fig.tight_layout();display(fig);plt.close(fig);return rows
