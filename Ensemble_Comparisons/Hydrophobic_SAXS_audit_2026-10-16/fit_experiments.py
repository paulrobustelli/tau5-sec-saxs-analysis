"""Exploratory SAXS-only Borthakur fits; no production weights overwritten.
All WT tests share the frozen original REST 100-ns-block sampling SEM,
unit-per-point base regularization (as in original engine), and a profiled
single intensity scale. Calibration preserves absolute Kish support of
0.1*N_original_pool, including removals/additions. Not held-out validation.
"""
from pathlib import Path
import sys,json,time
import numpy as np,pandas as pd
from analyze import load,O,C
sys.path.insert(0,str(C));from kish_engine import fit,dual
D=O/'fit_experiments';D.mkdir(exist_ok=True)
r,rw,rr,q,y,sy=load('REST_WT');z=np.load(C/'REST_WT/observable_matrix.npz',allow_pickle=True);sem=z['sigma_ensemble'][z['datatype']=='SAXS'];s0=float(z['initial_global_scale'])
compact=np.flatnonzero(r.Rg_CA_A<19);expanded=np.flatnonzero(r.Rg_CA_A>=25);nadd=min(160,len(expanded),len(compact));compact=compact[np.linspace(0,len(compact)-1,nadd,dtype=int)];expanded=expanded[np.linspace(0,len(expanded)-1,nadd,dtype=int)]
np.savez(D/'donor_indices.npz',compact_REST=compact,expanded_REST=expanded)

def run(label,raw,meta,neff,donor_start=None):
 path=D/label;path.mkdir(exist_ok=True)
 if (path/'complete.json').exists():return
 if (path/'started.json').exists():raise RuntimeError('Existing incomplete experiment requires inspection: '+label)
 (path/'started.json').write_text(json.dumps({'time':time.time(),'N':raw.shape[1],'target_absolute_Neff':neff,**meta},indent=2));records=[];N=raw.shape[1];target=neff/N;kwargs=dict(saxs_mask=np.ones(len(q),bool),raw_saxs=raw,initial_scale=s0)
 def ev(t):
  a=fit(raw*s0,y,np.full(len(q),t),sem,**kwargs);records.append({'theta':t,'Kish_ratio':a['Kish_ratio'],'scale':a['scale'],'dual_gradient':a['max_gradient'],'scale_gradient':a['logscale_primal_gradient']});return a
 lo=.01;a=ev(lo)
 if a['Kish_ratio']<target:
  hi=1.;a=ev(hi)
  while a['Kish_ratio']<target:hi*=2;a=ev(hi);assert hi<1e8
  for _ in range(20):
   mid=np.sqrt(lo*hi);b=ev(mid)
   if b['Kish_ratio']>=target:hi=mid;a=b
   else:lo=mid
  theta=hi
 else:theta=lo
 w=a['weights'];pred=a['scale']*(raw@w);res=(pred-y)/sy
 # Independent checks of weights, primal/dual equality and finite differences.
 assert np.isfinite(w).all() and np.isclose(w.sum(),1) and 1/(w@w)>=neff-1e-5
 var=theta**2+sem**2;primal=np.sum(w[w>0]*np.log(w[w>0]*N))+.5*np.sum((pred-y)**2/var)
 assert np.isclose(primal,-a['dual'],atol=1e-5,rtol=1e-5),(primal,a['dual'])
 A=a['centered_normalized'];norm=a['norm'];lam=a['lambda'];grad=dual(lam,A,np.zeros(len(q)),np.ones(N),np.full(len(q),theta)/norm,sem/norm)[1];checks=[]
 for k in [0,len(q)//2,len(q)-1]:
  eps=1e-5;lp=lam.copy();lm=lam.copy();lp[k]+=eps;lm[k]-=eps
  fd=(dual(lp,A,np.zeros(len(q)),np.ones(N),np.full(len(q),theta)/norm,sem/norm)[0]-dual(lm,A,np.zeros(len(q)),np.ones(N),np.full(len(q),theta)/norm,sem/norm)[0])/(2*eps)
  checks.append({'k':k,'analytic':float(grad[k]),'finite_difference':float(fd)});assert abs(fd-grad[k])<1e-4
 assert a['max_gradient']<1e-3 and abs(a['logscale_primal_gradient'])<1e-3
 report={'label':label,'N':N,'Neff':float(1/(w@w)),'theta':theta,'scale':a['scale'],'chi2_all':float(np.mean(res**2)),'chi2_mid':float(np.mean(res[(q>.08)&(q<.25)]**2)),'chi2_low':float(np.mean(res[q<=.08]**2)),'chi2_high':float(np.mean(res[q>=.25]**2)),'max_gradient':a['max_gradient'],'scale_gradient':a['logscale_primal_gradient'],'primal':float(primal),'negative_dual':-a['dual'],'gradient_checks':checks,**meta}
 if donor_start is not None:report['added_REST_weight']=float(w[donor_start:].sum())
 pd.DataFrame({'weight':w}).to_csv(path/'weights.csv',index=False);pd.DataFrame({'q':q,'I_exp':y,'sigma':sy,'prediction':pred}).to_csv(path/'curve.csv',index=False);pd.DataFrame(records).to_csv(path/'calibration.csv',index=False);(path/'complete.json').write_text(json.dumps(report,indent=2));print(label,report['chi2_all'],report['chi2_mid'],report.get('added_REST_weight'),flush=True)

if __name__=='__main__':
 for name in ['REST_WT','IDPForge_WT','bAIes_WT','ASTEROIDS_WT']:
  d,w,raw,qq,yy,ss=load(name);assert np.array_equal(qq,q) and np.array_equal(yy,y)
  neff=.1*len(d);run(name+'_baseline',raw,{'base_pool':name,'experiment':'baseline; SAXS only'},neff)
  if name=='REST_WT':
   ii=np.flatnonzero(d.Rg_CA_A>=19);run(name+'_remove_Rg_lt19',raw[:,ii],{'base_pool':name,'experiment':'Remove Rg<19A','source_indices':ii.tolist()},neff)
   # Equal-population subsets with same 2-A Rg-bin counts, differing in burial.
   low=[];high=[];bins=np.floor(d.Rg_CA_A/2).astype(int)
   for b in sorted(set(bins)):
    ii=np.flatnonzero(bins==b);ii=ii[np.argsort(d.loc[ii,'buried_HSA_nonlocal_A2'].to_numpy(),kind='stable')];half=len(ii)//2
    if half>=3:low.extend(ii[:half]);high.extend(ii[-half:])
   for lab,ii in [('low_burial_Rgmatched',low),('high_burial_Rgmatched',high)]:run(name+'_'+lab,raw[:,ii],{'base_pool':name,'experiment':lab,'source_indices':list(map(int,ii))},neff)
  else:
   for lab,ii in [('add_compact_REST',compact),('add_expanded_REST',expanded)]:run(name+'_'+lab,np.c_[raw,rr[:,ii]],{'base_pool':name,'experiment':lab,'added_N':len(ii),'added_REST_indices':ii.tolist()},neff,donor_start=len(d))
