from pathlib import Path
import json,hashlib,numpy as np,mdtraj as md
R=Path(__file__).resolve().parents[2];O=Path(__file__).resolve().parent;S=R/'outputs/bAIes_harmonic_2026-10-16'
assert not (O/'manifest.json').exists()
seqs=json.loads((R/'outputs/IDPForge_2026-10-13/sequences.json').read_text());records=[]
for v in ['WT','AA']:
 t=md.load(str(S/v/'production.xtc'),top=str(S/v/'protein.pdb'));assert len(t)==20000 and t.top.to_fasta()==[seqs[v]];assert np.allclose(t.time,np.arange(1,20001)*100)
 ii=np.arange(19,20000,20);heavy=t.top.select('element != H');sel=t[ii].atom_slice(heavy);dest=O/v;dest.mkdir(exist_ok=True)
 for k,i in enumerate(ii):
  ident=f'HARMONIC120_{v}_{k+1:04d}';p=dest/(ident+'.pdb');sel[k].save_pdb(str(p));lines=[]
  for l in p.read_text().splitlines():
   if l.startswith(('ATOM','HETATM')):
    if l[17:20] in ['HID','HIE','HIP']:l=l[:17]+'HIS'+l[20:]
    lines.append(l)
  data=('\n'.join(lines)+'\nEND\n').encode();p.write_bytes(data);job=R/'work/maxent_vm_share/jobs'/('CS_'+ident);job.mkdir(exist_ok=False);(job/'input.pdb').write_bytes(data);(job/'sparta_ready').touch()
  records.append(dict(variant=v,id=ident,file=str(p),job=str(job),frame_index=int(i),production_ns=float(t.time[i]/1000),absolute_ns=float(t.time[i]/1000+20),sha256=hashlib.sha256(data).hexdigest(),predictor_input_sha256=hashlib.sha256(data).hexdigest()))
 print(v,'1000 extracted',flush=True)
(O/'manifest.json').write_text(json.dumps(records,indent=2));(O/'sampling_protocol.json').write_text(json.dumps({'N_per_variant':1000,'frame_indices':'19,39,...19999','production_ns':'2,4,...2000','absolute_ns':'22,24,...2020','CS_construct':'native120 GP+AR330-447; no shifts recomputed for appended Cys','SAXS_construct':'121 with terminal Cys appended using previously validated empirical geometry procedure','SAXS_protocol':'Pepsi3.0 -ms0.6 -ns601 --dro5 -n80 --noSmearing; default radius1.61A; raw intensities; no fitting at this stage','reweighting':'not launched'},indent=2))
