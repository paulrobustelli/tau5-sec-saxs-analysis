import os
os.environ.setdefault('MPLCONFIGDIR',str(__import__('pathlib').Path(__file__).resolve().parents[1]/'.mplconfig'))
import h5py,numpy as np,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def load(name):
 with h5py.File(ROOT/'BNL_SAXS/data'/f'{name}.hdf5') as f:
  q=f['profiles/q'][:];a=np.stack([f[f'profiles/{i:06d}'][:] for i in range(450)],axis=1)
 return q,a[:,:,0],a[:,:,1]
def correct(y,e,name,variant='default'):
 n=y.shape[1]; t=np.arange(n)
 if name.startswith('WT') or variant=='buffer':
  a,b=(220,236) if name.startswith('WT') else (24,70)
  if variant=='early': a,b=50,75
  bg=y[:,a:b+1].mean(1);be=np.sqrt((e[:,a:b+1]**2).sum(1))/(b-a+1)
  return y-bg[:,None],np.sqrt(e**2+be[:,None]**2)
 idx=np.r_[0:11,439:450];X=np.c_[np.ones(len(idx)),idx];XX=np.c_[np.ones(n),t]
 yc=y.copy();ec=e.copy()
 for j in range(len(y)):
  w=1/e[j,idx]**2;cov=np.linalg.inv(X.T@(w[:,None]*X));beta=cov@(X.T@(w*y[j,idx]));yc[j]-=XX@beta;ec[j]=np.sqrt(e[j]**2+np.einsum('ij,jk,ik->i',XX,cov,XX))
 return yc,ec
if __name__=='__main__':
 fig,ax=plt.subplots(2,2,figsize=(12,8))
 for k,name in enumerate(['WT_092025','AA_052026']):
  q,y,e=load(name);Y,E=correct(y,e,name)
  sel=(q>=.012)&(q<=.25);trace=np.trapezoid(Y[sel],q[sel],axis=0)
  ax[k,0].plot(trace);ax[k,0].set_title(name);ax[k,0].set_xlim(0,449)
  ax[k,1].plot(trace);ax[k,1].set_xlim(100,240);ax[k,1].set_ylim(-.01,max(trace[100:240])*1.1)
  ranges=[(155,180),(145,195),(125,195),(162,173),(175,185)] if k==0 else [(163,204),(155,215),(172,182)]
  for a,b in ranges:
   D=Y[sel,a:b+1]/E[sel,a:b+1].mean(1)[:,None];u,s,v=np.linalg.svd(D,full_matrices=False)
   acu=np.sum(u[:-1,:4]*u[1:,:4],axis=0);acv=np.sum(v[:4,:-1]*v[:4,1:],axis=1)
   print(name,(a,b),'sv',np.round(s[:5],2),'q_ac',np.round(acu,2),'t_ac',np.round(acv,2),'noiseedge',round(np.sqrt(D.shape[0])+np.sqrt(D.shape[1]),2))
  for a,b in [(0,10),(24,70),(50,75),(220,236),(439,449)]:print('bg',name,a,b,np.round(np.mean(y[(q>.04)&(q<.15),a:b+1]),5))
 fig.tight_layout();fig.savefig(ROOT/'results/initial_traces.png',dpi=160)
