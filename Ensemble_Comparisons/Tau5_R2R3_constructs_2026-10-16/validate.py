from pathlib import Path
import json,hashlib
import mdtraj as md,numpy as np
O=Path(__file__).resolve().parent;seqs=json.loads((O/'sequences.json').read_text());records=[]
for p in sorted(O.glob('*.pdb')):
 t=md.load(str(p));v='AA' if '_AA_' in p.name else 'WT';mega='mega' in p.name;two='two_helices' in p.name
 a,b=(391,446) if mega or two else ((391,413) if p.stem.endswith('_R2') else (426,446))
 assert t.top.to_fasta()[0]==seqs[v][a-328:b-327];assert t.n_residues==b-a+1;assert np.isfinite(t.xyz).all()
 errors={};prolines={}
 for label,fn,target in [('phi',md.compute_phi,-57),('psi',md.compute_psi,-47)]:
  ix,ang=fn(t);values=np.degrees(ang[0]);rs=np.array([t.top.atom(x[1]).residue.resSeq for x in ix]);mask=np.ones(len(rs),bool)
  if two:mask=((rs>=391)&(rs<=413))|((rs>=426)&(rs<=446))
  if label=='phi':
   for r in [420,425]:
    if (rs==r).any():prolines[str(r)]=float(values[rs==r][0])
   mask &= ~np.isin(rs,[420,425])
  err=np.abs((values[mask]-target+180)%360-180);assert err.max()<1
  errors[label]=float(err.max())
 _,omega=md.compute_omega(t);dev=float(np.max(abs(abs(np.degrees(omega))-180)));assert dev<1
 bonds=np.array([(x.index,y.index) for x,y in t.top.bonds]);maxbond=float(md.compute_distances(t,bonds).max());assert maxbond<.2
 records.append(dict(file=p.name,AR_range=[a,b],residues=t.n_residues,sequence=t.top.to_fasta()[0],max_target_torsion_error_degrees=errors,proline_phi_degrees=prolines,max_trans_omega_deviation_degrees=dev,max_bond_nm=maxbond,sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
(O/'validation.json').write_text(json.dumps(records,indent=2));print('PASS: 8 sequences, numbering, helical torsions, trans peptide bonds, finite coordinates, bond lengths')
