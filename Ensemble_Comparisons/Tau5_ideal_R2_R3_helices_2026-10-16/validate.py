from pathlib import Path
import json,hashlib
import mdtraj as md,numpy as np
O=Path(__file__).resolve().parent;seqs=json.loads((O/'sequences.json').read_text());rows=[]
for p in sorted(O.glob('*.pdb')):
 t=md.load(str(p));v='AA' if 'AA' in p.stem else 'WT';full=p.stem.startswith('Tau5');a,b=(328,448) if full else ((390,410) if '_R2_' in p.stem else (430,440))
 assert t.n_residues==b-a+1 and t.top.to_fasta()[0]==seqs[v][a-328:b-327]
 assert np.isfinite(t.xyz).all()
 errors={}
 for label,fn,target,center in [('phi',md.compute_phi,-57,1),('psi',md.compute_psi,-47,1)]:
  ix,angles=fn(t);ars=np.array([t.top.atom(i[center]).residue.resSeq for i in ix]);mask=((ars>=390)&(ars<=410))|((ars>=430)&(ars<=440));mask &= (ars!=390) if label=="phi" else True;err=np.abs((np.degrees(angles[0,mask])-target+180)%360-180);errors[label]=float(err.max());assert err.max()<1, (p,label,err.max())
 bonds=np.array([(x.index,y.index) for x,y in t.top.bonds]);maxbond=float(md.compute_distances(t,bonds).max());assert maxbond<.20,(p,maxbond)
 rows.append(dict(file=p.name,residues=t.n_residues,sequence=t.top.to_fasta()[0],max_torsion_error_degrees=errors,max_bond_nm=maxbond,sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
(O/'validation.json').write_text(json.dumps(rows,indent=2));print(json.dumps(rows,indent=2))
