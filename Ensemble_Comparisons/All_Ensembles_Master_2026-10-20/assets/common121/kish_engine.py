"""Borthakur dual, Kish calibration, one shared SAXS scale per ensemble.
No experimental normalization of individual conformers. Sampling errors remain
fixed at their uniform-prior estimate during each calibration scan.
"""
import ast
from pathlib import Path
import numpy as np
from scipy.optimize import minimize,minimize_scalar
from scipy.special import softmax
R=Path(__file__).resolve().parents[2]
source=R/'work/Borthakur_MaxEnt_IDPs_2024/PaaA2/reweighting_examples/reweight_functions.py'
ns={'np':np}
for node in ast.parse(source.read_text()).body:
 if isinstance(node,ast.FunctionDef) and node.name in ['func_sigma_reg_sigma_md','Weighted_Independent_Blocks','Weighted_Sigma_of_Mean']:
  exec(compile(ast.Module(body=[node],type_ignores=[]),str(source),'exec'),ns)
dual=ns['func_sigma_reg_sigma_md']
def sampling_error(X,kind,time_ns=None):
 if kind=='IID':return np.std(X,axis=1,ddof=1)/np.sqrt(X.shape[1]),{}
 time_ns=np.asarray(time_ns);assert np.all(np.diff(time_ns)>0)
 dt=np.median(np.diff(time_ns));scan={}
 for duration in [10,25,50,100,150,200]:
  dim=max(1,int(round(duration/dt)))
  if X.shape[1]//dim<5:continue
  values=[]
  for row in X:
   a,w=ns['Weighted_Independent_Blocks'](row,np.ones(len(row)),dim)
   values.append(ns['Weighted_Sigma_of_Mean'](a,w))
  scan[str(duration)]=np.array(values)
 assert '100' in scan,'Need at least five100nsblocks for chosenMDerror'
 return scan['100'],scan

def fit(X,y,sigma_reg,sigma_md,saxs_mask=None,raw_saxs=None,initial_scale=1.):
 X=np.array(X,float);y=np.asarray(y,float);sr=np.asarray(sigma_reg);sm=np.asarray(sigma_md)
 mask=np.zeros(len(y),bool) if saxs_mask is None else np.asarray(saxs_mask,bool)
 norm=np.maximum.reduce([np.std(X,axis=1),sr,sm,np.full(len(y),1e-7)])
 warm=np.zeros(len(y));cache={}
 def inner(logratio):
  nonlocal warm
  scale=initial_scale*np.exp(float(logratio));Z=X.copy()
  if mask.any():Z[mask]=raw_saxs*scale
  A=(Z-y[:,None])/norm[:,None]
  opt=minimize(dual,warm,args=(A,np.zeros(len(y)),np.ones(X.shape[1]),sr/norm,sm/norm),jac=True,method='L-BFGS-B',options={'maxiter':10000,'ftol':1e-13,'gtol':1e-8,'maxls':60})
  assert opt.success,opt.message
  warm=opt.x.copy();w=softmax(-opt.x@A)
  return opt,w,scale,A
 if mask.any():
  # Profile the regularized primal objective = minus minimized dual.
  # Optimize one scale shared by every conformer, not per-conformer scales.
  def profile(u):return -inner(u)[0].fun
  optscale=minimize_scalar(profile,bounds=(-np.log(20),np.log(20)),method='bounded',options={'xatol':1e-7,'maxiter':100})
  assert optscale.success
  assert abs(optscale.x)<np.log(20)-.001,'Scale reached search boundary'
  result,w,scale,A=inner(optscale.x)
  scale_gradient=float(np.dot(result.x[mask]/norm[mask],scale*np.einsum('ij,j->i',raw_saxs,w)))
 else:result,w,scale,A=inner(0);scale_gradient=0.
 kish=1/(len(w)*np.dot(w,w));assert np.isfinite(w).all() and np.isclose(w.sum(),1)
 return {'weights':w,'Kish_ratio':float(kish),'scale':float(scale),'dual':float(result.fun),'max_gradient':float(np.max(abs(result.jac))),'logscale_primal_gradient':scale_gradient,'lambda':result.x,'norm':norm,'centered_normalized':A}

def calibrate(X,y,base,sem,label,records,**kwargs):
 def ev(theta):
  r=fit(X,y,base*theta,sem,**kwargs);records.append({'stage':label,'theta':theta,**{k:r[k] for k in ['Kish_ratio','scale','dual','max_gradient','logscale_primal_gradient']}});return r
 lo=.01;r=ev(lo)
 if r['Kish_ratio']>=.1:return lo,r
 hi=1.;r=ev(hi)
 while r['Kish_ratio']<.1:
  hi*=2;r=ev(hi);assert hi<1e8
 for _ in range(20):
  mid=np.sqrt(lo*hi);q=ev(mid)
  if q['Kish_ratio']>=.1:hi=mid;r=q
  else:lo=mid
 return hi,r
