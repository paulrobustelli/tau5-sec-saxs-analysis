from pathlib import Path
import json
import pymol
pymol.finish_launching(['pymol','-cq'])
from pymol import cmd
O=Path(__file__).resolve().parent
seqs=json.loads((O/'sequences.json').read_text())
def atom(obj,r,n):return f'{obj} and resi {r} and name {n}'
for v,s in seqs.items():
 for label,a,b,ss in [('R2R3_two_helices',391,446,4),('R2',391,413,1),('R3',426,446,1),('R2R3_mega_helix',391,446,1)]:
  name=f'Tau5_{v}_{label}';cmd.fab(s[a-328:b-327],name,resi=a,chain='A',ss=ss,hydro=0)
  if label=='R2R3_two_helices':
   for start,end in [(391,413),(426,446)]:
    for r in range(start,end+1):
     if r>a:cmd.set_dihedral(atom(name,r-1,'C'),atom(name,r,'N'),atom(name,r,'CA'),atom(name,r,'C'),-57)
     if r<b:cmd.set_dihedral(atom(name,r,'N'),atom(name,r,'CA'),atom(name,r,'C'),atom(name,r+1,'N'),-47)
  cmd.save(str(O/f'{name}.pdb'),name);cmd.delete('all')
print('Built 8 PDBs',flush=True);cmd.quit()
