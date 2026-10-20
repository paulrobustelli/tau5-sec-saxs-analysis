from pathlib import Path
import mdtraj as md,numpy as np,pandas as pd,json
OUT=Path(__file__).resolve().parent
rows=[]
for v in ['WT','AA']:
 for c in ['matrix','8A','8A_sep2']:
  d=OUT/f'{v}_{c}_20ns'
  if not (d/'final_state.xml').exists():continue
  t=md.load(str(d/'trajectory.xtc'),top=str(d/'protein.pdb'));hs=md.compute_dssp(t,simplified=False)=='H';times=np.arange(1,len(t)+1)*.1
  data=pd.DataFrame({'time_ns':times,'full_H':hs[:,2:].mean(axis=1),'R2_AR394_402_H':hs[:,66:75].mean(axis=1),'R3_AR430_440_H':hs[:,102:113].mean(axis=1),'Rg_CA_A':md.compute_rg(t.atom_slice(t.top.select('name CA')))*10});data.to_csv(d/'observables.csv',index=False)
  avg=data[data.time_ns>=10].mean(numeric_only=True).to_dict();avg.update(variant=v,condition=c,frames=len(t),window='10-20 ns; diagnostic, not converged');rows.append(avg)
pd.DataFrame(rows).to_csv(OUT/'sensitivity_summary.csv',index=False);(OUT/'sensitivity_complete.json').write_text(json.dumps(rows,indent=2));print(rows)
