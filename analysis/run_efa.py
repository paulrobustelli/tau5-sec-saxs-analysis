"""Reproducible SEC-SAXS EFA using pure numerical functions from BioXTAS RAW.
No notebook is executed and no GUI/runtime stubs are installed.
"""
import os
os.environ.setdefault('MPLCONFIGDIR',str(__import__('pathlib').Path(__file__).resolve().parents[1]/'.mplconfig'))
import ast,json,csv,hashlib,subprocess
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import packaging.version as packv
from inspect_data import load,correct,ROOT
OUT=ROOT/'results'; OUT.mkdir(exist_ok=True)
import raw_efa_core
ns=vars(raw_efa_core)

def guinier(q,I,E,qmin=.012,qmax=.04):
 upper=qmax
 for _ in range(15):
  m=(q>=qmin)&(q<=upper)&(I>0)&(E>0)&np.isfinite(I)&np.isfinite(E)
  if m.sum()<8:return {'Rg':None}
  x=q[m]**2;y=np.log(I[m]);e=E[m]/I[m];X=np.c_[np.ones(len(x)),x];cov=np.linalg.inv(X.T@((1/e**2)[:,None]*X));beta=cov@(X.T@(y/e**2))
  if beta[1]>=0:return {'Rg':None}
  rg=np.sqrt(-3*beta[1]);new=min(qmax,1/rg)
  if upper<=new+1e-10:break
  upper=new
 residual=(y-X@beta)/e
 return dict(Rg=float(rg),fit_se=float(3*np.sqrt(cov[1,1])/(2*rg)),I0=float(np.exp(beta[0])),qmin=float(q[m][0]),qmax=float(q[m][-1]),qRgmax=float(q[m][-1]*rg),n=int(m.sum()),chi2=float(np.sum(residual**2)/(len(x)-2)))

def analyze(name,a,b,qmin=.012,qmax=.25,variant='default',simulate=0,plots=False):
 q,Y,E=load(name);y,e=correct(Y,E,name,variant);m=(q>=qmin)&(q<=qmax);I=y[m,a:b+1];err=e[m,a:b+1];w=err.mean(1);D=I/w[:,None];u,s,vt=np.linalg.svd(D,full_matrices=False)
 au=np.sum(u[:-1,:5]*u[1:,:5],0);av=np.sum(vt[:5,:-1]*vt[:5,1:],1)
 fit1=np.outer(u[:,0]*s[0],vt[0]);res=(D-fit1)*w[:,None]/err
 result=dict(sample=name,frames=[a,b],q=[qmin,qmax],background=variant,sv=s[:5].tolist(),q_ac=au.tolist(),t_ac=av.tolist(),rank1_mean_squared_standardized_residual=float(np.mean(res**2)))
 # Null: rank-one mean + independent input errors, then same background estimator.
 # This preserves covariance due to the NEW residual-background estimate, but cannot
 # recover covariance already present in facility-subtracted input.
 if simulate:
  rng=np.random.default_rng(20260921);ss=[]
  for _ in range(simulate):
   noise=rng.normal(size=Y.shape)*E
   noise,_=correct(noise,E,name,variant)
   sim=fit1+noise[m,a:b+1]/w[:,None]
   ss.append(np.linalg.svd(sim,compute_uv=False)[1])
  result.update(null_s2_99=float(np.quantile(ss,.99)),null_s2_median=float(np.median(ss)),p_rank1=float((1+np.sum(np.array(ss)>=s[1]))/(simulate+1)),null_n=simulate)
 # Rank-one SAXS profile recovered by linear LS in ORIGINAL intensity units.
 c=vt[0].copy()
 if c.sum()<0:c=-c
 c/=c.sum();mult=c/(c@c)
 spectrum=y[:,a:b+1]@mult;se=np.sqrt((e[:,a:b+1]**2)@(mult**2))
 result['rank1_guinier']=guinier(q,spectrum,se)
 label=f'{name}_{a}-{b}_{variant}_{qmin:g}-{qmax:g}'
 if plots:
  frames=np.arange(a,b+1);fw=ns['runEFA'](D)[:4];bw=ns['runEFA'](D,False)[:4,::-1]
  np.savetxt(OUT/(label+'_efa.csv'),np.c_[frames,fw.T,bw.T],delimiter=',',header='frame,forward_s1,forward_s2,forward_s3,forward_s4,backward_s1,backward_s2,backward_s3,backward_s4',comments='')
  np.savez_compressed(OUT/(label+'_matrices.npz'),q=q[m],frames=frames,I=I,error=err,row_error=w,U=u,s=s,Vt=vt,forward=fw,backward=bw,rank1_residual=res)
  np.savetxt(OUT/(label+'_rank1_profile.dat'),np.c_[q,spectrum,se],header='q_A^-1 I_arbitrary conditional_sigma; scale: sum of elution coefficients = 1')
  fig,axs=plt.subplots(2,3,figsize=(14,8));trace=np.trapezoid(y[m],q[m],axis=0)
  axs[0,0].plot(np.arange(450),trace);axs[0,0].axvspan(a,b,alpha=.2,color='green');axs[0,0].set(xlim=(100,230),xlabel='Frame (zero indexed)',ylabel='Integrated I, arbitrary units',title=name)
  axs[0,1].semilogy(np.arange(1,min(15,len(s))+1),s[:15],'o-');axs[0,1].set(xlabel='Singular value index',ylabel='Singular value',title=f'Frames {a}–{b}')
  if simulate:axs[0,1].axhline(result['null_s2_99'],color='grey',ls='--',label='99% rank-one null s2');axs[0,1].legend(fontsize=8)
  axs[0,2].plot(np.arange(1,6),au,'o-',label='q vectors');axs[0,2].plot(np.arange(1,6),av,'s-',label='frame vectors');axs[0,2].axhline(.6,ls=':',c='grey');axs[0,2].set(xlabel='Mode',ylabel='Lag-one dot product',ylim=(-.5,1.05),title='Smoothness, not a species count');axs[0,2].legend()
  for j in range(4):
   axs[1,0].semilogy(frames,fw[j],label=f'Mode {j+1}');axs[1,1].semilogy(frames,bw[j],label=f'Mode {j+1}')
  for ax,t in [(axs[1,0],'Forward EFA'),(axs[1,1],'Backward EFA')]:ax.set(xlabel='Frame',ylabel='Singular value',title=t);ax.legend(fontsize=8)
  im=axs[1,2].imshow(res,origin='lower',aspect='auto',extent=(a-.5,b+.5,0,len(q[m])),vmin=-3,vmax=3,cmap='RdBu_r');axs[1,2].set(xlabel='Frame',ylabel='q-bin index (increasing q)',title='One-component residual / σ');fig.colorbar(im,ax=axs[1,2])
  fig.suptitle(f'{name}: q = {qmin:g}–{qmax:g} Å⁻¹; residual-background method: {variant}',fontsize=12);fig.tight_layout();fig.savefig(OUT/(label+'_diagnostics.png'),dpi=170);plt.close(fig)
 return result,(q,y,e,m,D,u,s,vt)

def rotation(name,a,b,ranges,qmin=.012,qmax=.25,variant='default',save=False):
 result,data=analyze(name,a,b,qmin,qmax,variant);q,y,e,m,D,u,s,vt=data
 rr=np.array(ranges)-a
 with np.errstate(all='ignore'):
  converged,conv,rot=ns['runRotation'](D,y[m,a:b+1],e[m,a:b+1],rr,[True]*len(rr),vt.T,niter=3000,tol=1e-9)
 r=dict(sample=name,frames=[a,b],ranges=ranges,q=[qmin,qmax],background=variant,converged=bool(converged),iterations=conv['iterations'])
 if not converged:return r
 C=rot['M']*rot['C'];mult=np.linalg.pinv(C.T);S=y[:,a:b+1]@mult;SE=np.sqrt(e[:,a:b+1]**2@mult**2)
 order=np.argsort(np.argmax(C,axis=0));C=C[:,order];S=S[:,order];SE=SE[:,order]
 r['components']=[dict(peak_frame=int(a+np.argmax(C[:,k])),guinier=guinier(q,S[:,k],SE[:,k]),negative_q_fraction=float(np.mean(S[m,k]<0))) for k in range(len(rr))]
 R=(y[m,a:b+1]-S[m]@C.T)/e[m,a:b+1];r['residual_ms']=float(np.mean(R**2));r['concentration_condition']=float(np.linalg.cond(C))
 if save:
  label=f'{name}_{a}-{b}_rotation'
  np.savetxt(OUT/(label+'_elution.csv'),np.c_[np.arange(a,b+1),C],delimiter=',',header='frame,'+','.join(f'component{k+1}' for k in range(len(rr))),comments='')
  for k in range(len(rr)):np.savetxt(OUT/(label+f'_component{k+1}.dat'),np.c_[q,S[:,k],SE[:,k]],header='q_A^-1 I_arbitrary conditional_sigma (not total decomposition uncertainty)')
  fig,axs=plt.subplots(2,2,figsize=(11,8))
  for k in range(len(rr)):
   rg=r['components'][k]['guinier']['Rg'];labelk=f'Component {k+1}: apparent Rg {rg:.1f} Å' if rg else f'Component {k+1}: Guinier fit failed'
   axs[0,0].plot(np.arange(a,b+1),C[:,k],label=f'Component {k+1}')
   axs[0,1].plot(q[m],S[m,k]/max(S[m,k]),label=labelk)
   sel=(q>=.012)&(q<=.04)&(S[:,k]>0);axs[1,0].plot(q[sel]**2,np.log(S[sel,k]),'.-',label=f'Component {k+1}')
  axs[0,0].set(xlabel='Frame',ylabel='Elution coefficient (sum = 1)',title='Relative elution, not mole fractions');axs[0,0].legend()
  axs[0,1].set(xlabel='q (Å⁻¹)',ylabel='I / maximum I',yscale='log',title='Recovered curves');axs[0,1].legend(fontsize=8)
  axs[1,0].set(xlabel='q² (Å⁻²)',ylabel='ln I (arbitrary scale)',title='Low-q inspection');axs[1,0].legend()
  im=axs[1,1].imshow(R,aspect='auto',origin='lower',extent=(a-.5,b+.5,0,len(q[m])),vmin=-3,vmax=3,cmap='RdBu_r');axs[1,1].set(xlabel='Frame',ylabel='q-bin index',title='Two-component residual / σ');fig.colorbar(im,ax=axs[1,1]);fig.suptitle('Exploratory WT shoulder separation; component windows are assumptions');fig.tight_layout();fig.savefig(OUT/(label+'.png'),dpi=170);plt.close(fig)
 return r

if __name__=='__main__':
 rows=[]
 for name,rois in [('WT_092025',[(155,180),(125,195),(120,205)]),('AA_052026',[(163,204),(155,215)])]:
  for a,b in rois:
   r,_=analyze(name,a,b,simulate=300,plots=True);rows.append(r);print(json.dumps(r),flush=True)
 (OUT/'primary_results.json').write_text(json.dumps(rows,indent=2))
 sensitivity=[]
 for name,rois in [('WT_092025',[(155,180),(125,195)]),('AA_052026',[(163,204),(155,215)])]:
  for a,b in rois:
   for variant in (['default','early'] if name.startswith('WT') else ['default','buffer']):
    for qmin,qmax in [(.008,.25),(.012,.15),(.012,.25),(.015,.25),(.020,.25),(.012,.35)]:
     r,_=analyze(name,a,b,qmin,qmax,variant,simulate=100);sensitivity.append(r)
 (OUT/'sensitivity_results.json').write_text(json.dumps(sensitivity,indent=2))
 # Candidate windows are subsequently inspected against the EFA curves.
 rotations=[]
 for endearly in [153,158,163,168]:
  for startlate in [145,150,155]:
   r=rotation('WT_092025',120,205,[[120,endearly],[startlate,205]],save=(endearly==163 and startlate==150));rotations.append(r);print('ROT',json.dumps(r),flush=True)
 (OUT/'rotation_sensitivity.json').write_text(json.dumps(rotations,indent=2))
