from pathlib import Path
import subprocess,sys,time,json,fcntl,datetime,hashlib
import numpy as np
R=Path(__file__).resolve().parents[2];O=Path(__file__).resolve().parent;E=R/'outputs/bAIes_harmonic_Cys121_2026-10-19'
lock=(O/'vm.lock').open('w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
assert not (O/'backcalculation_complete.json').exists()
logfile=O/'vm.log';offset=logfile.stat().st_size if logfile.exists() else 0
log=logfile.open('ab');p=subprocess.Popen([sys.executable,str(R/'work/maxent_vm.py')],stdin=subprocess.PIPE,stdout=log,stderr=subprocess.STDOUT)
(O/'vm_pid.json').write_text(json.dumps({'vm_pid':p.pid,'driver_pid':__import__('os').getpid(),'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'workers':'6 SAXS +2 CS, separate unique-prefix queues'}))
awake=subprocess.Popen(['/usr/bin/caffeinate','-i','-w',str(p.pid)])
for _ in range(180):
 if p.poll() is not None:raise RuntimeError('VM exited during boot')
 with logfile.open('rb') as f:f.seek(offset);ready=b'VM_READY' in f.read()
 if ready:break
 time.sleep(1)
else:raise RuntimeError('VM readiness timeout')
p.stdin.write(b'for i in 1 2 3 4 5 6; do sh /share/worker_harmonic_saxs.sh "$i" > /share/harmonic-saxs-$i.log 2>&1 & done\nfor i in 1 2; do sh /share/worker_harmonic_sparta.sh "$i" > /share/harmonic-sparta-$i.log 2>&1 & done\nwait\n');p.stdin.flush();print('WORKERS_STARTED',flush=True)
cs=json.loads((O/'manifest.json').read_text());seq=json.loads((R/'outputs/IDPForge_2026-10-13/sequences.json').read_text());verified={};keys={}
while True:
 if p.poll() is not None:raise RuntimeError('VM exited before completion')
 sa=json.loads((E/'manifest.json').read_text()) if (E/'complete.json').exists() else []
 counts={}
 for typ,records,flag in [('CS',cs,'sparta_done'),('SAXS',sa,'done')]:
  for m in records:
   ident=typ+m['id'];j=Path(m['job'])
   if ident in verified or not (j/flag).exists():continue
   assert (j/flag).read_text().strip()=='0',j
   assert hashlib.sha256((j/'input.pdb').read_bytes()).hexdigest()==m['predictor_input_sha256']
   if typ=='SAXS':
    f=j/'prediction.dat';x=np.loadtxt(f);assert x.shape==(601,8) and np.isfinite(x).all() and (x[:,1]>0).all() and np.allclose(x[:,0],np.linspace(0,.6,601))
   else:
    f=j/'pred.tab';lines=f.read_text().splitlines();s=''.join(''.join(l.split()[2:]) for l in lines if l.startswith('DATA SEQUENCE'));assert s==seq[m['variant']]
    rows=[l.split() for l in lines if len(l.split())==9 and l.split()[0].isdigit()];kk=[(int(a[0]),a[1],a[2]) for a in rows];assert len(kk)>500 and len(set(kk))==len(kk);assert all(1<=i<=len(s) and s[i-1]==a for i,a,b in kk);assert np.isfinite(np.array([a[3:] for a in rows],float)).all()
    if m['variant'] in keys:assert kk==keys[m['variant']]
    else:keys[m['variant']]=kk
   verified[ident]={'type':typ,'variant':m['variant'],'id':m['id'],'output_sha256':hashlib.sha256(f.read_bytes()).hexdigest()}
  counts[typ]={v:sum(a['type']==typ and a['variant']==v for a in verified.values()) for v in ['WT','AA']}
 state={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'verified':counts,'reweighting':'not launched'};(O/'status.json').write_text(json.dumps(state,indent=2));(O/'native_validation.json').write_text(json.dumps(verified,indent=2))
 if len(verified)==4000:
  (O/'backcalculation_complete.json').write_text(json.dumps(state,indent=2));p.terminate();p.wait(timeout=30);print('ALL_4000_VERIFIED_VM_STOPPED',flush=True);break
 time.sleep(30)
