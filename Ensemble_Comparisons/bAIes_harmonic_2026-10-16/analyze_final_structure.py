from pathlib import Path
import ast,json,hashlib
import mdtraj as md,numpy as np,pandas as pd
O=Path(__file__).resolve().parent;R=O.parents[1];A=O/'final_region_analysis'
assert (A/'complete.json').exists() and not (A/'structure_complete.json').exists()
source=R/'outputs/Tau5_bAIes_OpenMM_restart_2026-10-09/Analysis/analyze_production.py';ns={'np':np}
for f in ast.parse(source.read_text()).body:
 if isinstance(f,ast.FunctionDef) and f.name in ['rmsd','sa','rg']:exec(compile(ast.Module(body=[f],type_ignores=[]),str(source),'exec'),ns)
ref=md.load(str(source.parent/'R2R3_all_helix.pdb'));ns['refs']=np.array([ref.xyz[0,ref.top.select(f'residue {i} to {i+6} and backbone')] for i in range(1,51)],float)
summary=[]
for v in ['WT','AA']:
 top=md.load(str(O/v/'protein.pdb')).top;idx=np.array([[next(a.index for a in res.atoms if a.name==n) for n in ['N','CA','C','O']] for res in top.residues]);pairs=np.array([(i,j) for i in range(120) for j in range(i+1,120)])
 sums=np.zeros((10,len(pairs)));counts=np.zeros(10,int);rows=[];count=0
 for t in md.iterload(str(O/v/'production.xtc'),top=top,chunk=100):
  sa,crop=ns['sa'](t.xyz[:,idx].astype(float));rg=ns['rg'](t.xyz[:,idx[63:119,1]])
  d,_=md.compute_contacts(t,contacts=pairs,scheme='closest-heavy',periodic=False)
  for k in range(len(t)):
   b=(count+k)//2000;sums[b]+=(d[k]<1.2);counts[b]+=1
   rows.append({'frame_index':count+k,'production_ns':float(t.time[k]/1000),'Sa_full':sa[k],'Sa_AR391_446':crop[k],'Rg_R2R3_CA_A':rg[k]})
  count+=len(t)
  if count%2000==0:print(v,count,flush=True)
 assert count==20000 and np.all(counts==2000)
 df=pd.DataFrame(rows);df.to_csv(A/f'{v}_Sa_timeseries.csv',index=False);df.groupby(df.frame_index//2000).mean().to_csv(A/f'{v}_Sa_200ns_blocks.csv')
 prob=sums.sum(0)/count;matrix=np.zeros((120,120));matrix[pairs[:,0],pairs[:,1]]=prob;matrix+=matrix.T
 np.savetxt(A/f'{v}_contact_map.csv',matrix,delimiter=',');np.savez_compressed(A/f'{v}_block_contacts.npz',probabilities=sums/counts[:,None],counts=counts,pairs=pairs)
 summary.append({'variant':v,'frames':count,'means':df[['Sa_full','Sa_AR391_446','Rg_R2R3_CA_A']].mean().to_dict(),'trajectory_sha256':hashlib.sha256((O/v/'production.xtc').read_bytes()).hexdigest()})
(A/'structure_complete.json').write_text(json.dumps(summary,indent=2))
