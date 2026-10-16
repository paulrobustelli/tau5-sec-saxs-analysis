"""Cache original 2000-member 121-residue ASTEROIDS priors without selection changes."""
from pathlib import Path
import tarfile,tempfile,json,hashlib
import numpy as np,pandas as pd,mdtraj as md
from cache_structure import analyze
O=Path(__file__).resolve().parent
B=Path('/Users/paulrobustelli/Documents/Codex/2026-09-28/is-th')
for variant,case in [('WT','WT'),('AA','W397A_W433A')]:
 d=O/('ASTEROIDS_'+variant);out=d/'prior_structural_cache.npz'
 if out.exists():continue
 archive=B/'outputs/Ensemble-Selection/refined_reproducibility'/case/'prepared_pool.tar.gz'
 ids=pd.read_csv(B/'work'/f'refined-{case}-saxs/predicted_dro5.csv').conformer_id.astype(str).tolist()
 vv=[];hh=[];cc=[]
 with tarfile.open(archive) as tf,tempfile.TemporaryDirectory(prefix='tau5_ast_prior_') as tmp:
  members={Path(m.name).stem:m for m in tf.getmembers() if m.isfile() and m.name.endswith('.pdb')}
  assert len(ids)==2000 and all(x in members for x in ids)
  for start in range(0,len(ids),20):
   frames=[]
   for cid in ids[start:start+20]:
    p=Path(tmp)/(cid+'.pdb');p.write_bytes(tf.extractfile(members[cid]).read());t=md.load(str(p));t=t.atom_slice([a.index for a in t.top.atoms if a.element.symbol!='H']);frames.append(t);p.unlink()
   t=md.join(frames,check_topology=True);assert t.n_residues==121
   v,h,c,pairs=analyze(t);vv.append(v);hh.append(h);cc.append(c)
   if (start+20)%200==0:print(variant,start+20,flush=True)
 np.savez_compressed(out,values=np.concatenate(vv),columns=['Rg_CA_A','Rg_AR391446_CA_A','Sa_full','Sa_AR391446'],helix=np.concatenate(hh),contacts=np.concatenate(cc),contact_pairs=pairs,conformer_ids=np.array(ids),weights=np.full(len(ids),1/len(ids)))
 (d/'prior_structure_provenance.json').write_text(json.dumps({'N':len(ids),'residues':121,'archive':str(archive),'archive_sha256':hashlib.sha256(archive.read_bytes()).hexdigest(),'order':'predicted_dro5.csv conformer_id','selection':'uniform original prior'},indent=2))
 print(variant,'complete',flush=True)
