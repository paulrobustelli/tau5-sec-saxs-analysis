"""Audit normalized support, concentration, chronological blocks, and PPII populations."""
from pathlib import Path
import json
import numpy as np,pandas as pd
O=Path(__file__).resolve().parent;R=O.parents[1]
for name in ['REST_WT','bAIes_WT','bAIes_AA','IDPForge_WT','IDPForge_AA']:
 d=O/name
 if not (d/'fit_complete.json').exists():continue
 frame=pd.read_csv(d/'weights.csv');w=frame.combined.to_numpy();assert np.isfinite(w).all() and np.all(w>=0) and abs(w.sum()-1)<1e-8
 n=len(w);kish=1/(n*np.einsum('i,i->',w,w));assert kish>=.1
 order=np.argsort(w)[::-1];cum=np.cumsum(w[order]);report={'N':n,'Kish_ratio':kish,'Neff':kish*n,'max_weight':float(w.max()),'top10_weight':float(w[order[:10]].sum()),'N_for_50percent':int(np.searchsorted(cum,.5)+1),'N_for_90percent':int(np.searchsorted(cum,.9)+1),'Shannon_effective_N':float(np.exp(-np.sum(w*np.log(np.maximum(w,1e-300))))) }
 manifest=json.loads((d/'manifest.json').read_text())
 if not name.startswith('IDPForge'):
  t=np.array([m.get('production_ns',m.get('time_ns')) for m in manifest]);assert np.all(np.diff(t)>0)
  block=np.floor((t-t.min())/100).astype(int);b=pd.DataFrame({'time_ns':t,'weight':w,'block100ns':block}).groupby('block100ns').agg(N=('weight','size'),weight=('weight','sum'),first_ns=('time_ns','min'),last_ns=('time_ns','max'));b.to_csv(d/'weight_chronological_blocks.csv');report['largest_100ns_block_weight']=float(b.weight.max());report['block_note']='Last block may be shorter. REST saved ordering is not physical kinetics; block weights assess concentration only.'
 src='MD_WT' if name=='REST_WT' else name
 a=np.load(R/'outputs/PPII_dihedrals_2026-10-15'/f'{src}_angles.npz');angles=a['angles_degrees']
 if len(angles)!=n:
  # Full bAIes trajectories map to the frozen chronological predictor snapshot.
  indices=np.array([m['frame_index'] for m in manifest]);angles=angles[indices]
 assert len(angles)==n
 omega=angles[:,:,2];pd.DataFrame({'AR_residue':a['AR_residues'],'prior_cis_percent':100*(abs(omega)<30).mean(0),'MaxEnt_cis_percent':100*np.einsum('i,ij->j',w,abs(omega)<30),'prior_other_percent':100*((abs(omega)>=30)&(abs(omega)<=150)).mean(0),'MaxEnt_other_percent':100*np.einsum('i,ij->j',w,(abs(omega)>=30)&(abs(omega)<=150))}).to_csv(d/'PPII_cis_before_after.csv',index=False)
 (d/'weight_support_audit.json').write_text(json.dumps(report,indent=2));print(name,report,flush=True)
