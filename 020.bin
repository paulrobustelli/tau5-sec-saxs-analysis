"""WT two-signal EFA; shared-buffer parametric bootstrap; own-peak scaling."""
from pathlib import Path
import sys,json,os
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/saxs_mpl')
REPO=Path(os.environ.get('SAXS_REPO','/Users/f0044gk/Desktop/SAXS/repository'))
sys.path.insert(0,str(REPO/'analysis'))
from component_recovery import load,solve
import numpy as np
OUT=Path(__file__).resolve().parent
q,Y,E=load('WT_092025');fr=np.arange(120,206);A=np.eye(Y.shape[1])[:,fr];A[220:237]-=1/17
I=Y@A;err=np.sqrt(E**2@(A*A));info,o=solve(q,I,err,fr,175,140)
assert o is not None
S=o['S']*o['C'].max(0);C=o['C']/o['C'].max(0)
rng=np.random.default_rng(20260924);bs=[];bc=[]
for k in range(200):
 _,b=solve(q,(Y+rng.normal(size=Y.shape)*E)@A,err,fr,175,140)
 if b is not None:
  bs.append(b['S']*b['C'].max(0));bc.append(b['C']/b['C'].max(0))
 if k%25==0:print(k,flush=True)
bs=np.array(bs);bc=np.array(bc);se=bs.std(0,ddof=1)
np.savez_compressed(OUT/'WT_components.npz',q=q,frames=fr,S=S,C=C,error=se,bootstrap_S=bs,bootstrap_C=bc,I=I,frame_error=err)
for j,name in enumerate(['oligomer_candidate','monomer']):np.savetxt(OUT/f'WT_{name}_Iq.dat',np.c_[q,S[:,j],se[:,j]],header='q_A^-1 I sigma_bootstrap; own-peak component amplitude; buffer frames 220-236')
info.update(bootstrap_success=len(bs),bootstrap_attempts=200,buffer=[220,236],component_order=['oligomer_candidate','monomer'],noise_handling='Rank-two signal reconstruction; remaining modes stay in the residual, not displayed as species.',normalization='Each component concentration profile has maximum 1; profiles are not molecular concentrations.')
(OUT/'decomposition.json').write_text(json.dumps(info,indent=2));print('DONE',flush=True)
