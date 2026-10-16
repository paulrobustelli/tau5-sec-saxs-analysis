from pathlib import Path
import hashlib,json
import numpy as np,pandas as pd,mdtraj as md
from cache_structure import analyze
R=Path(__file__).resolve().parents[2];O=Path(__file__).resolve().parent;old=Path('/Users/paulrobustelli/Documents/Codex/2026-09-28/is-th');F=old/'outputs/Ensemble-Selection/final_ensembles';E=R/'outputs/IDPForge_D2D_MaxEnt_2026-10-15/experimental';reports=[]
for v,case in [('WT','WT'),('AA','W397A_W433A')]:
 D=O/('ASTEROIDS_'+v);D.mkdir(exist_ok=True);src=F/(v+'_shifts_plus_SAXS');selected=pd.read_csv(src/'selected.csv');saxsdir=old/'work'/f'refined-{case}-saxs';raw=pd.read_csv(saxsdir/'predicted_dro5.csv').set_index('conformer_id');assert len(selected)==500
 q,iy,sy=np.loadtxt(E/f'{v}_SAXS.dat',unpack=True);nativeq=np.array([float(x.split('=')[1]) for x in raw.columns]);assert q.min()>=nativeq.min()-1e-6 and q.max()<=nativeq.max()+1e-6
 values=raw.loc[selected.conformer_id].to_numpy();w=selected.weight.to_numpy();w/=w.sum();curves=np.array([np.interp(q,nativeq,x) for x in values]);prior=np.mean(raw.to_numpy(),axis=0);prior=np.interp(q,nativeq,prior);chosen=np.einsum('i,ij->j',w,curves)
 def scale(z):return np.sum(z*iy/sy**2)/np.sum((z/sy)**2)
 ps=scale(prior);ss=scale(chosen);pd.DataFrame({'q_Ainv':q,'I_exp':iy,'sigma_exp':sy,'I_prior_dro5':prior*ps,'I_selected_dro5':chosen*ss}).to_csv(D/'matched_saxs_fit.csv',index=False)
 if not (D/'structural_cache.npz').exists():
  t=md.load(str(src/'ensemble.pdb.gz'));assert len(t)==500 and t.n_residues==121;vv=[];hh=[];cc=[]
  for i in range(0,len(t),20):
   z=t[i:i+20];z=z.atom_slice([a.index for a in z.top.atoms if a.element.symbol!='H']);val,h,c,p=analyze(z);vv.append(val);hh.append(h);cc.append(c)
  np.savez_compressed(D/'structural_cache.npz',values=np.concatenate(vv),columns=['Rg_CA_A','Rg_AR391446_CA_A','Sa_full','Sa_AR391446'],helix=np.concatenate(hh),contacts=np.concatenate(cc),contact_pairs=p,weights=w,conformer_ids=selected.conformer_id.to_numpy())
 for stage,file in [('prior','prior_predictions.csv'),('selected','predictions.csv')]:
  df=pd.read_csv(src/file);df.to_csv(D/(stage+'_archived_predictions.csv'),index=False)
 report={'variant':v,'N_selected':len(selected),'N_prior':len(raw),'SAXS_prior_chi2_mean_common_dro5':float(np.mean(((prior*ps-iy)/sy)**2)),'SAXS_selected_chi2_mean_common_dro5':float(np.mean(((chosen*ss-iy)/sy)**2)),'scale_prior':float(ps),'scale_selected':float(ss),'selection_provenance':'Existing ASTEROIDS final selection preserved. Re-evaluated with canonicaldro5 and shared experimentalgrid, independently optimized global ensemble scale/nooffset; originalselectionoptimizedhydration maydiffer. Not a newly rerun ASTEROIDS selection or KishMaxEnt.','source_SAXS_sha256':hashlib.sha256((saxsdir/'predicted_dro5.csv').read_bytes()).hexdigest(),'source_ensemble_sha256':hashlib.sha256((src/'ensemble.pdb.gz').read_bytes()).hexdigest()};(D/'reference_provenance.json').write_text(json.dumps(report,indent=2));reports.append(report);print(v,'reference ready',flush=True)
(O/'ASTEROIDS_reference_complete.json').write_text(json.dumps(reports,indent=2))
