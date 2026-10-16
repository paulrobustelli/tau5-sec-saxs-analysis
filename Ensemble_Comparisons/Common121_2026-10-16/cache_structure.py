"""Common121res structural matrices in manifest order, reusable for any weights."""
from pathlib import Path
import json,ast,sys
import numpy as np,pandas as pd,mdtraj as md
R=Path(__file__).resolve().parents[2];O=Path(__file__).resolve().parent
source=R/'outputs/Tau5_bAIes_OpenMM_restart_2026-10-09/Analysis/analyze_production.py';ns={'np':np}
for n in ast.parse(source.read_text()).body:
 if isinstance(n,ast.FunctionDef) and n.name=='rmsd':exec(compile(ast.Module(body=[n],type_ignores=[]),str(source),'exec'),ns)
ref=md.load(str(source.parent/'R2R3_all_helix.pdb'));refs=np.array([ref.xyz[0,ref.top.select(f'residue {i} to {i+6} and backbone')] for i in range(1,51)],float)
def analyze(t):
 assert t.n_residues==121
 idx=np.array([[next(a.index for a in res.atoms if a.name==n) for n in ['N','CA','C','O']] for res in t.top.residues]);bb=t.xyz[:,idx].astype(float);rr=np.repeat(refs[:1],115,axis=0);rr[63:113]=refs;wins=np.stack([bb[:,i:i+7].reshape(len(t),28,3) for i in range(115)],1);q=(ns['rmsd'](wins,rr)/.1)**4;sw=(1+q)/(1+q+q*q)
 def rg(x):return np.sqrt(np.mean(np.sum((x-x.mean(1,keepdims=True))**2,2),1))*10
 vals=np.c_[rg(bb[:,:,1]),rg(bb[:,63:119,1]),sw.sum(1),sw[:,63:113].sum(1)];h=md.compute_dssp(t,simplified=False)=='H';pairs=np.array([(i,j) for i in range(121) for j in range(i+1,121)]);d,p=md.compute_contacts(t,contacts=pairs,scheme='closest-heavy',periodic=False)
 return vals,h,d<1.2,pairs
extension=json.loads((R/'outputs/Cys121_extensions_2026-10-15/manifest.json').read_text());datasets={}
for e in ['IDPForge','bAIes']:
 for v in ['WT','AA']:datasets[e+'_'+v]=[m for m in extension if m['ensemble']==e and m['variant']==v]
datasets['REST_WT']=json.loads((R/'outputs/WT_MD_300K_2026-10-15/manifest.json').read_text())
for name,manifest in datasets.items():
 dest=O/name;dest.mkdir(exist_ok=True)
 if (dest/'structural_cache.npz').exists():continue
 vv=[];hh=[];cc=[]
 for start in range(0,len(manifest),20):
  frames=[]
  for m in manifest[start:start+20]:
   t=md.load(m['file']);t=t.atom_slice([a.index for a in t.top.atoms if a.element.symbol!='H']);frames.append(t)
  t=md.join(frames,check_topology=True);v,h,c,p=analyze(t);vv.append(v);hh.append(h);cc.append(c)
  if (start+20)%200==0:print(name,start+20,flush=True)
 vals=np.concatenate(vv);np.savez_compressed(dest/'structural_cache.npz',values=vals,columns=['Rg_CA_A','Rg_AR391446_CA_A','Sa_full','Sa_AR391446'],helix=np.concatenate(hh),contacts=np.concatenate(cc),contact_pairs=p,files=np.array([m['file'] for m in manifest]))
 pd.DataFrame(vals,columns=['Rg_CA_A','Rg_AR391446_CA_A','Sa_full','Sa_AR391446']).to_csv(dest/'structural_values.csv',index=False);(dest/'manifest.json').write_text(json.dumps(manifest,indent=2));print(name,'complete',flush=True)
