from pathlib import Path
import json
import pymol
pymol.finish_launching(['pymol','-cq'])
from pymol import cmd
O=Path(__file__).resolve().parent
seqs=json.loads((O/'sequences.json').read_text())
def atom(obj,r,n):return f'{obj} and resi {r} and name {n}'
for v,s in seqs.items():
 name=f'Tau5_{v}_R2_R3_ideal_helices'
 cmd.fab(s,name,resi=328,chain='A',ss=4,hydro=0)
 for a,b in [(390,410),(430,440)]:
  for r in range(a,b+1):
   cmd.set_dihedral(atom(name,r-1,'C'),atom(name,r,'N'),atom(name,r,'CA'),atom(name,r,'C'),-57)
   cmd.set_dihedral(atom(name,r,'N'),atom(name,r,'CA'),atom(name,r,'C'),atom(name,r+1,'N'),-47)
 cmd.save(str(O/f'{name}.pdb'),name)
 cmd.delete('all')
 for region,a,b in [('R2',390,410),('R3',430,440)]:
  name=f'{v}_{region}_{a}_{b}_ideal_helix';cmd.fab(s[a-328:b-327],name,resi=a,chain='A',ss=1,hydro=0)
  cmd.save(str(O/f'{name}.pdb'),name);cmd.delete('all')
print('Six PDBs built',flush=True)
cmd.quit()
