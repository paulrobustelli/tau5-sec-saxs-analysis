from pathlib import Path
import subprocess,numpy as np,json
R=Path(__file__).resolve().parents[2];O=Path(__file__).resolve().parent;rows=[]
def read(p):
 l=p.read_text().splitlines();i=next(i for i,x in enumerate(l) if x.startswith('ITEM: ATOMS'));return np.loadtxt(l[i+1:])[:,4:]*41.84
for case in ['initial','perturbed','expanded']:
 c=O/'cases'/case;s=(c/'reference.in').read_text()
 s+='\nthermo_style custom step pe\nunfix drycmap\nrun 0\nwrite_dump all custom no_cmap.dat id x y z fx fy fz modify sort id format float %.15g\n'
 s+='dihedral_style multi/harmonic\ndihedral_coeff * 0 0 0 0 0\nrun 0\nwrite_dump all custom no_torsion.dat id x y z fx fy fz modify sort id format float %.15g\n'
 s+='pair_coeff * * 0.0 1.0 1.122462048309373\nrun 0\nwrite_dump all custom no_nonbonded.dat id x y z fx fy fz modify sort id format float %.15g\n'
 (c/'isolate.in').write_text(s);subprocess.run([str(R/'work/lammps-build/lmp'),'-in','isolate.in','-log','isolate.log'],cwd=c,check=True,stdout=(c/'isolate_stdout.log').open('w'))
 z=np.load(c/'force_comparison.npz');a=read(c/'reference_forces.dat');b=read(c/'no_cmap.dat');d=read(c/'no_torsion.dat');e=read(c/'no_nonbonded.dat')
 for term,ref in [('CMAPTorsionForce',a-b),('PeriodicTorsionForce',b-d),('CustomNonbondedForce',d-e),('bond_plus_angle',e)]:
  om=z[term] if term!='bond_plus_angle' else z['HarmonicBondForce']+z['HarmonicAngleForce'];dif=om-ref
  rows.append(dict(case=case,term=term,RMS_difference=float(np.sqrt(np.mean(dif*dif))),max_difference=float(np.max(abs(dif))),reference_RMS=float(np.sqrt(np.mean(ref*ref)))))
(O/'isolated_force_comparison.json').write_text(json.dumps(rows,indent=2));print(json.dumps(rows,indent=2))
