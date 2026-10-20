from pathlib import Path
import sys,json,time
import numpy as np,pandas as pd
R=Path(__file__).resolve().parents[2];O=Path(__file__).resolve().parent;H=R/'outputs/Hydrophobic_SAXS_audit_2026-10-16';sys.path.insert(0,str(H))
from analyze import load,C
sys.path.insert(0,str(C));from kish_engine import fit,dual
D=O/'motif_ablations';D.mkdir(exist_ok=True)
_,_,_,q,y,sy=load('REST_WT');z=np.load(C/'REST_WT/observable_matrix.npz',allow_pickle=True);sem=z['sigma_ensemble'][z['datatype']=='SAXS'];s0=float(z['initial_global_scale'])
source=(H/'fit_experiments.py').read_text();exec(source[source.index('def run('):source.index("if __name__=='__main__':")])
# These are geometric masks defined without SAXS residuals. Controls remove identical counts per 2-A Rg bin.
name='ASTEROIDS_WT';d=pd.read_csv(O/(name+'_metrics.csv'));_,_,raw,_,_,_=load(name)
masks={'broad_reversals_low_burial':(d.broad_antiparallel_pairs>=2)&(d.buried_HSA_nonlocal_A2<400),'high_burial':d.buried_HSA_nonlocal_A2>=800}
records=[]
for label,mask in masks.items():
 removed=np.flatnonzero(mask);keep=np.flatnonzero(~mask);records.append({'class':label,'removed':len(removed),'remaining':len(keep)})
 if len(removed)==0 or len(keep)<200:continue
 run('AST_WT_remove_'+label,raw[:,keep],{'base_pool':name,'experiment':'geometric motif removal','motif':label,'source_indices':keep.tolist()},200)
 bins=np.floor(d.Rg_CA_A/2).astype(int).to_numpy()
 for seed in [17,29]:
  rng=np.random.default_rng(seed);ctrl=[]
  for b in np.unique(bins[removed]):
   pool=np.flatnonzero(bins==b);ctrl.extend(rng.choice(pool,np.sum(bins[removed]==b),replace=False).tolist())
  keep=np.setdiff1d(np.arange(len(d)),ctrl)
  run('AST_WT_control_'+label+'_'+str(seed),raw[:,keep],{'base_pool':name,'experiment':'Rg-bin-matched random removal','motif':label,'seed':seed,'source_indices':keep.tolist()},200)
pd.DataFrame(records).to_csv(O/'ablation_masks.csv',index=False)
