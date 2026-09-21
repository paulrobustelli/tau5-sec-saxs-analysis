from run_efa import *
from scipy.optimize import minimize_scalar
from scipy.stats import chi2
PLOT=ROOT/'poster_results';PLOT.mkdir(exist_ok=True)

def background_cov(name,Y,E):
 n=Y.shape[1]
 if name.startswith('WT'):
  vv=(E[:,220:237]**2).sum(1)/17**2
  return lambda t1,t2:vv
 idx=np.r_[0:11,439:450];X=np.c_[np.ones(len(idx)),idx]
 covs=np.array([np.linalg.inv(X.T@((1/ee[idx]**2)[:,None]*X)) for ee in E])
 return lambda t1,t2:np.einsum('i,qij,j->q',np.array([1,t1]),covs,np.array([1,t2]))

def average(name,a,b):
 q,Y,E=load(name);y,e=correct(Y,E,name);bcov=background_cov(name,Y,E)
 mean=y[:,a:b+1].mean(1);err=np.sqrt((E[:,a:b+1]**2).sum(1)/(b-a+1)**2+bcov((a+b)/2,(a+b)/2))
 return q,mean,err

def candidates(name,lo,hi):
 q,Y,E=load(name);y,e=correct(Y,E,name);bcov=background_cov(name,Y,E);m=(q>=.012)&(q<=.25)
 trace=np.trapezoid(y[m],q[m],axis=0);peak=lo+np.argmax(trace[lo:hi+1]);rows=[]
 for a in range(lo,peak+1):
  for b in range(peak,hi+1):
   n=b-a+1
   if n<7:continue
   avg=y[:,a:b+1].mean(1);err=np.sqrt((E[:,a:b+1]**2).sum(1)/n**2+bcov((a+b)/2,(a+b)/2));fit=guinier(q,avg,err)
   if fit['Rg'] is None:continue
   cuts=[guinier(q,avg,err,qmin=x) for x in [.010,.012,.015,.018]]
   if any(f['Rg'] is None for f in cuts):continue
   drift=(max(f['Rg'] for f in cuts)-min(f['Rg'] for f in cuts))/fit['Rg']
   mid=(a+b)//2;A=y[:,a:mid+1].mean(1);B=y[:,mid+1:b+1].mean(1);ta=(a+mid)/2;tb=(mid+1+b)/2
   va=(E[:,a:mid+1]**2).sum(1)/(mid-a+1)**2+bcov(ta,ta);vb=(E[:,mid+1:b+1]**2).sum(1)/(b-mid)**2+bcov(tb,tb);cab=bcov(ta,tb)
   def fun(s):return np.mean(((A[m]-s*B[m])**2)/(va[m]+s*s*vb[m]-2*s*cab[m]))
   sc=minimize_scalar(fun,bounds=(.01,10),method='bounded').x;halfms=fun(sc)
   # Rank-one diagnostics with uncertainties used only to evaluate consistency.
   D=y[m,a:b+1]/e[m,a:b+1].mean(1)[:,None];u,s,vt=np.linalg.svd(D,full_matrices=False)
   snr=np.sqrt(np.sum((avg[m]/err[m])**2))
   rows.append(dict(sample=name,start=a,end=b,n=n,peak=int(peak),profile_snr=float(snr),Rg=fit['Rg'],Rg_fit_se=fit['fit_se'],guinier_chi2=fit['chi2'],guinier_qmin=fit['qmin'],guinier_qmax=fit['qmax'],qmin_Rg_relative_span=float(drift),halves_mean_squared_z=float(halfms),s2=float(s[1]),mode2_q_ac=float(np.sum(u[:-1,1]*u[1:,1])),mode2_t_ac=float(np.sum(vt[1,:-1]*vt[1,1:])),eligible=bool(drift<=.05 and fit['chi2']<=1.5 and halfms<=1.5)))
 valid=[r for r in rows if r['eligible']]
 if not valid:raise RuntimeError('No windows passed prespecified quality criteria')
 best=max(valid,key=lambda r:r['profile_snr']);return rows,best
if __name__=='__main__':
 allrows=[];chosen=[]
 for name,lo,hi in [('WT_092025',150,195),('AA_052026',155,215)]:
  rows,best=candidates(name,lo,hi);allrows.extend(rows);chosen.append(best);print('SELECTED',best,flush=True)
  print('TOP5',sorted([r for r in rows if r['eligible']],key=lambda r:-r['profile_snr'])[:5],flush=True)
  q,I,E=average(name,best['start'],best['end']);np.savetxt(PLOT/f'{name}_selected_average.dat',np.c_[q,I,E],header=f"q_A^-1 I_arbitrary sigma_with_shared_residual_background; frames {best['start']}–{best['end']} inclusive")
 with open(PLOT/'window_candidates.csv','w') as f:
  w=csv.DictWriter(f,fieldnames=list(allrows[0]));w.writeheader();w.writerows(allrows)
 (PLOT/'selected_windows.json').write_text(json.dumps(chosen,indent=2))
