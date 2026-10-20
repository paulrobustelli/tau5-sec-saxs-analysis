from pathlib import Path
import numpy as np,pandas as pd,mdtraj as md
ROOT=Path(__file__).resolve().parents[2];OUT=Path(__file__).resolve().parent
regions={'R2_main_AR390_404':(390,404),'R2_minor_AR404_410':(404,410)}
rows=[];resrows=[]
def calc(v,label,t,mask=None):
 h=t if isinstance(t,np.ndarray) else md.compute_dssp(t,simplified=False)=='H'
 if mask is not None:h=h[mask]
 for region,(a,b) in regions.items():
  idx=np.arange(a-328,b-328+1);frac=h[:,idx].mean(axis=1)
  rows.append({'variant':v,'ensemble':label,'region':region,'frames_or_conformers':len(h),'alpha_fraction':float(frac.mean()),'alpha_percent':float(frac.mean()*100)})
  for i in idx:resrows.append({'variant':v,'ensemble':label,'AR_residue':i+328,'region':region,'alpha_fraction':float(h[:,i].mean())})
for v in ['WT','AA']:
 base=ROOT/'outputs/Tau5_bAIes_OpenMM_restart_2026-10-09'/v
 calc(v,'AF2_rank1',md.load(str(base/'prediction.pdb')))
 old_h=[]
 cpu=md.load(str(base/'production_2us.xtc'),top=str(base/'protein.pdb'))
 old_h.append(md.compute_dssp(cpu[:2],simplified=False)=='H')
 for chunk in md.iterload(str(base/'production_2us_gpu.xtc'),top=str(base/'protein.pdb'),chunk=1000):old_h.append(md.compute_dssp(chunk,simplified=False)=='H')
 original=np.concatenate(old_h);assert len(original)==20000
 calc(v,'original_2us_production',original)
 for c in ['matrix','8A','8A_sep2']:
  d=OUT/f'{v}_{c}_20ns';t=md.load(str(d/'trajectory.xtc'),top=str(d/'protein.pdb'));calc(v,f'diagnostic_{c}_10_20ns',t,np.arange(len(t))>=100)
 d=ROOT/'outputs/bAIes_8A_seqsep2_2us_2026-10-14'/v;t=md.load(str(d/'production.xtc'),top=str(d/'protein.pdb'));calc(v,f'new_production_partial_first_{len(t)*.1:.1f}ns',t)
 # Independently generated conformers currently available; snapshot, not final.
 d=ROOT/'outputs/IDPForge_D2D_200step_2026-10-14'/v;paths=sorted(d.glob('*_validated.pdb'))
 ts=[]
 for p in paths:
  try:t=md.load(str(p));assert len(list(t.top.residues))==120;ts.append(md.compute_dssp(t,simplified=False)=='H')
  except (OSError,ValueError):pass
 if ts:calc(v,'IDPForge_D2D_200step_partial',np.concatenate(ts))
pd.DataFrame(rows).to_csv(OUT/'helicity_corrected_R2_regions.csv',index=False);pd.DataFrame(resrows).to_csv(OUT/'helicity_corrected_R2_residue_profiles.csv',index=False)
print(pd.DataFrame(rows).to_string(index=False))
