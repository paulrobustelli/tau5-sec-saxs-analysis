from pathlib import Path
import os,json
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/saxs_mpl')
import numpy as np,h5py
from scipy.optimize import minimize_scalar
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from extended_guinier import fit,scan
OUT=Path(__file__).resolve().parent
SOURCE=Path('/Users/f0044gk/Desktop/SAXS/repository/BNL_SAXS/data/WT_092025.hdf5')
WINDOWS={'leading':(150,159),'center':(160,169),'trailing':(170,179),'late_tail':(180,189)}

def guinier(q,y,e):
 m=(q>=.012-1e-9)&(q<=.030+1e-9)
 if np.any(y[m]<=0):return None
 X=np.c_[np.ones(m.sum()),q[m]**2];s=e[m]/y[m];cov=np.linalg.inv(X.T@(X/s[:,None]**2));b=cov@(X.T@(np.log(y[m])/s**2))
 if b[1]>=0:return None
 rg=np.sqrt(-3*b[1]);return dict(Rg=float(rg),qRgmax=float(q[m][-1]*rg),I0=float(np.exp(b[0])),se=float(3*np.sqrt(cov[1,1])/(2*rg)))

def run():
 with h5py.File(SOURCE) as f:
  q=f['profiles/q'][:];a=np.stack([f[f'profiles/{i:06d}'][:] for i in range(450)],axis=1);Y,E=a[:,:,0],a[:,:,1]
 names=list(WINDOWS);summary={}; saved={}; rng=np.random.default_rng(22102026)
 for baseline in ['linear','pre','post']:
  W=[]
  for lo,hi in WINDOWS.values():
   w=np.zeros(450);w[lo:hi+1]=1/(hi-lo+1)
   frac=np.clip((np.arange(lo,hi+1)-74.5)/(222.0-74.5),0,1).mean()
   if baseline=='pre':frac=0
   if baseline=='post':frac=1
   w[55:95]-=(1-frac)/40;w[210:235]-=frac/25;W.append(w)
  W=np.array(W).T;D=Y@W;var=(E**2)@(W**2);err=np.sqrt(var);saved[baseline]=(D,err,W)
  for j,name in enumerate(names):
   r=dict(frames=WINDOWS[name],standard=guinier(q,D[:,j],err[:,j]))
   m=(q>=.012-1e-9)&(q<=.055+1e-9)
   r['extended_common']=fit(q[m],D[m,j],err[m,j]);r['extended_adaptive']=scan(q,D[:,j],err[:,j])[0]
   summary[baseline+'_'+name]=r
   np.savetxt(OUT/f'{baseline}_{name}.dat',np.c_[q,D[:,j],err[:,j]],header='q_A^-1 I sigma_shared_background')
  # Compare curve shapes allowing intensity scaling, retaining shared baseline covariance.
  for j,name in enumerate(names):
   if j==1:continue
   m=(q>=.012-1e-9)&(q<=.15+1e-9)
   def objective(scale):
    v=(E[m]**2)@((W[:,j]-scale*W[:,1])**2)
    return np.sum((D[m,j]-scale*D[m,1])**2/v)
   z=minimize_scalar(objective,bounds=(.01,5),method='bounded')
   summary[baseline+'_'+name]['shape_vs_center']=dict(scale=float(z.x),chi2_per_dof=float(z.fun/(m.sum()-1)),n=int(m.sum()))
  # Paired noise perturbations of all frames preserve cross-window baseline covariance.
  m=(q>=.012-1e-9)&(q<=.055+1e-9);boots=[]
  for b in range(300):
   pert=D[m]+(rng.normal(size=E[m].shape)*E[m])@W
   vals=[]
   for j in range(4):
    r=fit(q[m],pert[:,j],err[m,j]);vals.append(r['Rg'] if r else np.nan)
   boots.append(vals)
  boots=np.array(boots)
  for j,name in enumerate(names):
   summary[baseline+'_'+name]['extended_common_bootstrap95']=np.nanpercentile(boots[:,j],[2.5,97.5]).tolist()
   summary[baseline+'_'+name]['delta_vs_center_bootstrap95']=np.nanpercentile(boots[:,j]-boots[:,1],[2.5,97.5]).tolist()
 np.savez(OUT/'profiles.npz',q=q,**{k+'_'+n:v[j] for k,v in saved.items() for j,n in enumerate(['I','error','weights'])})
 (OUT/'results.json').write_text(json.dumps(summary,indent=2))
 fig,ax=plt.subplots(2,2,figsize=(11,8));D,err,W=saved['linear'];colors=['C0','C1','C2','C3']
 t=np.arange(120,216);f=(t-74.5)/147.5;traceY=Y[:,t]-Y[:,55:95].mean(1)[:,None]*(1-f)-Y[:,210:235].mean(1)[:,None]*f
 m=(q>=.012)&(q<=.15);ax[0,0].plot(t,np.trapezoid(traceY[m],q[m],axis=0),color='k')
 for j,(name,(lo,hi)) in enumerate(WINDOWS.items()):
  color=colors[j];ax[0,0].axvspan(lo,hi,alpha=.2,color=color,label=f'{name}: {lo}–{hi}')
  r=summary['linear_'+name];scale=r.get('shape_vs_center',{}).get('scale',1);ax[0,1].errorbar(q[m],D[m,j]/scale,err[m,j]/scale,fmt='.',ms=2,color=color,alpha=.65,label=name)
  if j!=1:
   v=(E[m]**2)@((W[:,j]-scale*W[:,1])**2);ax[1,0].plot(q[m],(D[m,j]-scale*D[m,1])/np.sqrt(v),'.',ms=3,color=color,label=name)
 for baseline,marker in [('pre','s'),('linear','o'),('post','^')]:
  vals=[summary[baseline+'_'+n]['extended_common']['Rg'] for n in names];ax[1,1].plot(range(4),vals,marker=marker,label=baseline)
 ax[0,0].set(xlabel='Frame (zero-based)',ylabel='Integrated scattering',title='WT direct averaging windows');ax[0,0].legend(fontsize=8)
 ax[0,1].set(xlabel='q (Å⁻¹)',ylabel='I(q), scaled to center',yscale='log',title='Same-baseline shape comparison');ax[0,1].legend(fontsize=8)
 ax[1,0].axhline(0,color='k',lw=.6);ax[1,0].set(xlabel='q (Å⁻¹)',ylabel='Difference / propagated σ',title='Residuals relative to center');ax[1,0].legend(fontsize=8)
 ax[1,1].set(xticks=range(4),xticklabels=names,ylabel='Extended Rg (Å)',title='Matched fit range: 0.012–0.055 Å⁻¹');ax[1,1].legend()
 fig.tight_layout();fig.savefig(OUT/'WT_window_comparison.png',dpi=170)
 return summary
if __name__=='__main__':
 for k,v in run().items():print(k,round(v['standard']['Rg'],2),round(v['extended_common']['Rg'],2),v.get('shape_vs_center'),v['delta_vs_center_bootstrap95'])
