from packing_metrics import calc
from pathlib import Path
import numpy as np,pandas as pd,mdtraj as md,json
O=Path(__file__).resolve().parent;R=O.parents[1];A=O/'comparison_analysis';A.mkdir(exist_ok=True)
for method,base in [('harmonic',O),('Jeffreys',R/'outputs/bAIes_8A_seqsep2_2us_2026-10-14')]:
 for v in ['WT','AA']:
  out=A/f'{method}_{v}_packing.csv'
  if out.exists():continue
  rows=[];bur=[];coord=[];count=0
  for t in md.iterload(str(base/v/'production.xtc'),top=str(base/v/'protein.pdb'),chunk=400):
   t=t[19::20];t=t.atom_slice([a.index for a in t.top.atoms if a.element.symbol!='H']);d,b,c,ar=calc(t);d.insert(0,'production_ns',t.time/1000);d.insert(0,'frame_index',np.arange(count+19,count+400,20));rows.append(d);bur.append(b);coord.append(c);count+=400
   if count%4000==0:print(method,v,count,flush=True)
  d=pd.concat(rows,ignore_index=True);assert len(d)==1000 and np.allclose(d.production_ns,np.arange(1,1001)*2)
  d.to_csv(out,index=False);np.savez_compressed(A/f'{method}_{v}_packing_residues.npz',buried_HSA_A2=np.concatenate(bur),coordination=np.concatenate(coord),AR_residues=ar)
  print(method,v,'packing complete',flush=True)
(A/'packing_protocol.json').write_text(json.dumps({'N_per_trajectory':1000,'sampling':'frames19,39,...19999; production2,4,...2000ns','construct':'native120res GP+AR330447 for BOTH likelihoods; no Cys extension','definitions_source':'outputs/Hydrophobic_SAXS_audit_2026-10-16/METHODS.txt','main_contact':'sidechain C/S <5A; |i-j|>6','burial':'SASA local +/-6-residue fragment minus whole structure,240spherepoints,1.4Aprobe','weights':'uniform; no SAXS or CS reweighting'},indent=2))
