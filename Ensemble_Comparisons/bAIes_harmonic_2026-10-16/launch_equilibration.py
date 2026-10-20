from pathlib import Path
import json,shutil,subprocess,hashlib,datetime
import numpy as np
R=Path(__file__).resolve().parents[2];O=Path(__file__).resolve().parent
assert json.loads((O/'harmonic_validation.json').read_text())['passed']
for v in ['WT','AA']:
 if (O/v/'launch.json').exists():raise RuntimeError('Launch already recorded: '+v)
for v in ['WT','AA']:
 src=R/'outputs/AF2_distogram_recompute_2026-10-14'/f'{v}_8A_sep2_20ns';dst=O/v;dst.mkdir(exist_ok=True)
 for f in ['protein.pdb','baies_params.dat']:shutil.copy2(src/f,dst/f)
 cmd=json.loads((src/'command.json').read_text());cmd[2]=str(O/'code/bAIes_openMM.py')
 changes={'-pdb':str(dst/'protein.pdb'),'-baies':str(dst/'baies_params.dat'),'-xtc':str(dst/'equilibration.xtc'),'-checkpoint':str(dst/'equilibration.chk'),'-state_out':str(dst/'equilibration_final.xml'),'-seed':str(2026101600+(v=='AA'))}
 for flag,value in changes.items():cmd[cmd.index(flag)+1]=value
 assert '-restart_state' not in cmd and '-restart' not in cmd
 params=np.atleast_2d(np.loadtxt(dst/'baies_params.dat'));assert len(params)==(45 if v=='WT' else 37)
 with (dst/'run.log').open('xb') as log:
  p=subprocess.Popen(cmd,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
 record=dict(pid=p.pid,command=cmd,started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),duration_ns=20,restraints=len(params),selection='8 Angstrom cutoff; seqsep2 includes i,i+3',energy='0.5*kBT*((r-mu)/sigma)^2',sigma='Fixed AF2 fitted standard deviation in nm',starting_structure='Original diagnostic input PDB; fresh minimization and velocities; no Jeffreys restart',hashes={f:hashlib.sha256((dst/f).read_bytes()).hexdigest() for f in ['protein.pdb','baies_params.dat']})
 (dst/'launch.json').write_text(json.dumps(record,indent=2))
 c=subprocess.Popen(['/usr/bin/caffeinate','-i','-w',str(p.pid)],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL,start_new_session=True)
 print(v,p.pid,'caffeinate',c.pid,flush=True)
