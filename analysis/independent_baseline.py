"""Independent residual-baseline assessment; no Rg or protein fit enters anchor selection."""
from inspect_data import *
from scipy.stats import chi2
from scipy.ndimage import median_filter
import json,csv
OUT=ROOT/'independent_baseline';OUT.mkdir(exist_ok=True)
BANDS=[(.012,.035),(.04,.10),(.10,.25),(.5,1.)]
SEARCH={'WT_092025':[(10,95),(200,235)],'AA_052026':[(10,75),(225,265)]}

def window_stats(q,Y,E,a,b):
 t=np.arange(a,b+1);mid=(a+b)//2;zs=[]
 for lo,hi in BANDS:
  m=(q>=lo)&(q<=hi);y=Y[m,a:b+1].mean(0);v=(E[m,a:b+1]**2).sum(0)/m.sum()**2
  X=np.c_[np.ones(len(t)),t-t.mean()];cov=np.linalg.inv(X.T@((1/v)[:,None]*X));beta=cov@(X.T@(y/v));zs.append(float(beta[1]/np.sqrt(cov[1,1])))
 m=(q>=.012)&(q<=.25);d=Y[m,a:mid+1].mean(1)-Y[m,mid+1:b+1].mean(1)
 v=(E[m,a:mid+1]**2).sum(1)/(mid-a+1)**2+(E[m,mid+1:b+1]**2).sum(1)/(b-mid)**2
 halfms=float(np.mean(d*d/v))
 return dict(start=a,end=b,n=b-a+1,max_band_slope_z=float(max(abs(z) for z in zs)),band_slope_z=zs,half_profile_ms=halfms,passes=bool(max(abs(z) for z in zs)<3 and halfms<1.5))

def find_anchors(name):
 q,Y,E=load(name);rows=[];chosen=[]
 for side,(lo,hi) in zip(['pre','post'],SEARCH[name]):
  rr=[]
  for n in [20,25,30,35,40]:
   for a in range(lo,hi-n+2,5):
    r=window_stats(q,Y,E,a,a+n-1);r.update(sample=name,side=side);rr.append(r)
  valid=[r for r in rr if r['passes']]
  if not valid:raise RuntimeError('No stable anchors: '+name+' '+side)
  # Maximize observations, then use proximity to the target peak. No Rg target.
  best=max(valid,key=lambda r:(r['n'],r['start'] if side=='pre' else -r['end']))
  rows.extend(rr);chosen.append([best['start'],best['end']])
 return chosen,rows

def anchor_means(Y,E,anchors):
 means=[];vars=[];times=[]
 for a,b in anchors:
  means.append(Y[:,a:b+1].mean(1));vars.append((E[:,a:b+1]**2).sum(1)/(b-a+1)**2);times.append((a+b)/2)
 return np.array(means),np.array(vars),np.array(times)

def baseline(Y,E,anchors,kind='linear'):
 means,vs,ts=anchor_means(Y,E,anchors);t=np.arange(Y.shape[1]);f=np.clip((t-ts[0])/(ts[1]-ts[0]),0,1)
 if kind=='pre':f=np.zeros_like(f)
 if kind=='post':f=np.ones_like(f)
 if kind=='integral':
  # q-specific integral correction following the iterative intensity-integral model.
  # Used as a model sensitivity, not proof of capillary fouling. Endpoint change
  # must be nonnegative under this mechanism; negative q offsets retain pre baseline.
  start,end=anchors[0][1],anchors[1][0];delta=np.maximum(means[1]-means[0],0)
  y=Y[:,start:end+1]-means[0,:,None]
  from scipy.signal import savgol_filter
  smooth=savgol_filter(y, min(21, y.shape[1]//2*2-1),3,axis=1)
  b=delta[:,None]*np.linspace(0,1,y.shape[1])[None,:]
  converged=False
  for iteration in range(200):
   # Deposition is nonnegative. Smoothing/positive clipping only defines this
   # alternative baseline; measured scattering profiles themselves are not clipped.
   signal=np.maximum(smooth-b,0);cum=np.cumsum(signal,axis=1);cum-=cum[:,:1]
   frac=np.divide(cum,cum[:,-1:],out=np.zeros_like(cum),where=cum[:,-1:]>0)
   new=delta[:,None]*frac
   change=np.max(np.abs(new-b));b=.5*b+.5*new
   if change<1e-8:converged=True;break
  B=np.broadcast_to(means[0,:,None],Y.shape).copy();B[:,start:end+1]+=b;B[:,end+1:]+=delta[:,None]
  return B,dict(converged=converged,iterations=iteration+1,negative_endpoint_fraction=float(np.mean(means[1]<means[0])))
 B=means[0,:,None]*(1-f)[None,:]+means[1,:,None]*f[None,:]
 return B,dict(converged=True)

def average_independent(name,a,b,anchors,kind='linear'):
 q,Y,E=load(name);B,meta=baseline(Y,E,anchors,kind);I=(Y-B)[:,a:b+1].mean(1);means,vs,ts=anchor_means(Y,E,anchors)
 f=np.clip(((a+b)/2-ts[0])/(ts[1]-ts[0]),0,1)
 if kind=='pre':f=0
 if kind=='post':f=1
 if kind=='integral':
  # Bootstrap variance includes estimated anchor profiles and the nonlinear baseline.
  rng=np.random.default_rng(20260921);sims=[]
  for _ in range(120):
   yn=Y+rng.normal(size=Y.shape)*E;bn,_=baseline(yn,E,anchors,kind);sims.append((yn-bn)[:,a:b+1].mean(1))
  error=np.std(sims,axis=0,ddof=1)
 else:error=np.sqrt((E[:,a:b+1]**2).sum(1)/(b-a+1)**2+(1-f)**2*vs[0]+f*f*vs[1])
 return q,I,error,meta

def fixed_guinier(q,I,E,lo=.012,hi=.030):
 m=(q>=lo-1e-10)&(q<=hi+1e-10)&(I>0)&(E>0);X=np.c_[np.ones(m.sum()),q[m]**2];v=(E[m]/I[m])**2;cov=np.linalg.inv(X.T@((1/v)[:,None]*X));beta=cov@(X.T@(np.log(I[m])/v))
 if beta[1]>=0:return dict(Rg=None)
 rg=float(np.sqrt(-3*beta[1]));return dict(Rg=rg,se=float(1.5*np.sqrt(cov[1,1])/rg),I0=float(np.exp(beta[0])),qmin=float(q[m][0]),qmax=float(q[m][-1]),qRgmax=float(q[m][-1]*rg),chi2=float(np.sum((np.log(I[m])-X@beta)**2/v)/(m.sum()-2)))

if __name__=='__main__':
 selections={};audit=[];results=[]
 for name in SEARCH:
  anchors,rows=find_anchors(name);selections[name]=anchors;audit.extend(rows);print('ANCHORS',name,anchors,flush=True)
 (OUT/'anchor_selection.json').write_text(json.dumps(selections,indent=2));(OUT/'anchor_candidate_audit.json').write_text(json.dumps(audit,indent=2))
 for name,a,b in [('WT_092025',159,175),('AA_052026',169,188)]:
  for kind in ['pre','post','linear','integral']:
   q,I,E,meta=average_independent(name,a,b,selections[name],kind);f=fixed_guinier(q,I,E)
   r=dict(sample=name,frames=[a,b],model=kind,anchors=selections[name],**f,**meta);results.append(r);print(r,flush=True)
   np.savetxt(OUT/f'{name}_{kind}_average.dat',np.c_[q,I,E],header='q_A^-1 corrected_I conditional_sigma; independently selected residual baseline')
 (OUT/'model_comparison.json').write_text(json.dumps(results,indent=2))
