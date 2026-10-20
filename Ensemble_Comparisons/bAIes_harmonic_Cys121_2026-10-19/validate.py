from pathlib import Path
import json,hashlib,numpy as np,pandas as pd,mdtraj as md
O=Path(__file__).resolve().parent;R=O.parents[1];m=json.loads((O/'manifest.json').read_text());seq=json.loads((R/'outputs/IDPForge_2026-10-13/sequences.json').read_text());report=[]
for r in m:
 t=md.load(r['file']);assert t.n_residues==121 and t.top.to_fasta()==[seq[r['variant']]+'C'];assert np.isfinite(t.xyz).all();assert hashlib.sha256(Path(r['file']).read_bytes()).hexdigest()==r['sha256']
 j=Path(r['job']);assert hashlib.sha256((j/'input.pdb').read_bytes()).hexdigest()==r['predictor_input_sha256'];assert not (j/'sparta_ready').exists()
 if (j/'done').exists():
  assert (j/'done').read_text().strip()=='0';x=np.loadtxt(j/'prediction.dat');assert x.shape==(601,8) and np.isfinite(x).all() and np.all(x[:,1]>0) and np.allclose(x[:,0],np.linspace(0,.6,601));report.append(r['id'])
v=pd.read_csv(O/'extension_validation.csv');assert len(v)==2000 and v.min_nonbonded_extension_core_nm.min()>=.16 and v.Cys_phi_deg.between(-180,-20).all()
(O/'validation.json').write_text(json.dumps({'validated_structures':len(m),'validated_native_SAXS':len(report),'complete_SAXS_ids':report,'minimum_nonbonded_distance_nm':float(v.min_nonbonded_extension_core_nm.min()),'terminal_Pro_carbonyl_rotations_needed':int((abs(v.terminal_Pro_psi_rotation_deg)>0).sum()),'CS':'New native120 CS in source_job; no appended-Cys CS'},indent=2));print(len(m),'structures verified;',len(report),'native SAXS verified')
