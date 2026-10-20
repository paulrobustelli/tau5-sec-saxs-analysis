from pathlib import Path
import json,hashlib,subprocess,datetime
import mdtraj as md,numpy as np
O=Path(__file__).resolve().parent;checks={}
for v in ['WT','AA']:
 d=O/v
 assert not (d/'production_launch.json').exists(),'Production launch already recorded'
 log=(d/'run.log').read_text();assert 'Completed 10000000 steps' in log
 t=md.load(str(d/'equilibration.xtc'),top=str(d/'protein.pdb'))
 assert t.n_frames==200 and np.isfinite(t.xyz).all()
 assert np.allclose(t.time,np.arange(1,201)*100,atol=.1)
 rows=np.array([[float(x) for x in line.split('\t')[:3]] for line in log.splitlines() if line and line[0].isdigit() and '\t' in line])
 assert rows.shape==(200,3) and np.isfinite(rows).all()
 assert np.all((rows[:,2]>200)&(rows[:,2]<400))
 h=md.compute_dssp(t,simplified=False);ar=np.arange(t.n_residues)+328
 checks[v]={'frames':200,'production_time_ns':20,'temperature_mean_K':float(rows[:,2].mean()),'temperature_range_K':[float(rows[:,2].min()),float(rows[:,2].max())],'potential_energy_range_kJmol':[float(rows[:,1].min()),float(rows[:,1].max())],'last10ns_R2_main_H_percent':float((h[100:,(ar>=390)&(ar<=404)]=='H').mean()*100),'last10ns_R3_H_percent':float((h[100:,(ar>=430)&(ar<=440)]=='H').mean()*100),'final_state_sha256':hashlib.sha256((d/'equilibration_final.xml').read_bytes()).hexdigest(),'passed':True}
(O/'equilibration_validation.json').write_text(json.dumps(checks,indent=2));print(json.dumps(checks),flush=True)
for v in ['WT','AA']:
 d=O/v;cmd=json.loads((d/'launch.json').read_text())['command']
 for flag,val in {'-xtc':str(d/'production.xtc'),'-Nsteps':'1000000000','-checkpoint':str(d/'production.chk'),'-state_out':str(d/'production_final.xml')}.items():cmd[cmd.index(flag)+1]=val
 cmd+=['-restart_state',str(d/'equilibration_final.xml')]
 with (d/'production.log').open('xb') as log:p=subprocess.Popen(cmd,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
 (d/'production_launch.json').write_text(json.dumps({'pid':p.pid,'command':cmd,'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'additional_ns':2000,'initial_step':10000000,'expected_final_step':1010000000,'method':'harmonic fixed-sigma bAIes-N likelihood; 8A seqsep2 includes i,i+3','equilibration_validation':checks[v]},indent=2))
 subprocess.Popen(['/usr/bin/caffeinate','-i','-w',str(p.pid)],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL,start_new_session=True)
 print(v,'production PID',p.pid,flush=True)
