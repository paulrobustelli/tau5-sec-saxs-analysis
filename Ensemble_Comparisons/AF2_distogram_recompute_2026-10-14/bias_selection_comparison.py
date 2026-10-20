from pathlib import Path
import numpy as np,pandas as pd,mdtraj as md
from scipy.special import erf
ROOT=Path(__file__).resolve().parents[2];OUT=Path(__file__).resolve().parent
rows=[]
for v in ['WT','AA']:
 base=ROOT/'outputs/Tau5_bAIes_OpenMM_restart_2026-10-09'/v;top=md.load(str(base/'protein.pdb'))
 for c in ['matrix','8A']:
  run=OUT/f'{v}_{c}_20ns';params=np.loadtxt(run/'baies_params.dat',comments='#');pairs=params[:,1:3].astype(int)-1;mu=params[:,3];sigma=params[:,4]
  # Original serials verified contiguous and equal to topology atom indices+1.
  def energy(xyz):
   r=np.linalg.norm(xyz[:,pairs[:,0]]-xyz[:,pairs[:,1]],axis=-1);delta=r-mu
   ratio=np.empty_like(delta);nz=abs(delta)>1e-10;ratio[nz]=erf(delta[nz]/(np.sqrt(2)*np.broadcast_to(sigma,delta.shape)[nz]))/delta[nz];ratio[~nz]=np.broadcast_to(np.sqrt(2/np.pi)/sigma,delta.shape)[~nz]
   return -2.4789570296*np.log(ratio).sum(axis=1)
  try:t=md.load(str(run/'trajectory.xtc'),top=str(base/'protein.pdb'))
  except (OSError,ValueError):continue
  e=energy(t.xyz);af=float(energy(top.xyz)[0]);hs=md.compute_dssp(t,simplified=False)=='H'
  rows.append({'variant':v,'condition':c,'frames':len(t),'last_time_ns':len(t)*.1,'AF2_bias_energy_kJmol':af,'last_half_mean_bias_energy_kJmol':float(e[len(e)//2:].mean()),'AF2_minus_last_half_bias_kJmol':af-float(e[len(e)//2:].mean()),'last_half_R2_AR394_402_H':float(hs[len(e)//2:,66:75].mean()),'preliminary':not (run/'final_state.xml').exists()})
pd.DataFrame(rows).to_csv(OUT/'bias_selection_comparison.csv',index=False);print(pd.DataFrame(rows).to_string(index=False))
