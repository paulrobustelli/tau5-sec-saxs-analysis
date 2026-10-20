"""Launch each completed native dataset once; no repeated failed-fit launches."""
from pathlib import Path
import subprocess,sys,json,time,datetime,os
O=Path(__file__).resolve().parent;R=O.parents[1]
lock=O/'watch_ready_fits.lock'
try:lock.mkdir()
except FileExistsError:raise SystemExit('Watcher already registered; check its PID before considering restart')
(lock/'pid').write_text(str(os.getpid()))
names=['REST_WT','bAIes_WT','bAIes_AA','IDPForge_WT','IDPForge_AA']
while True:
 pending=False
 for name in names:
  D=O/name
  if (D/'fit_complete.json').exists() or (D/'fit_launch.json').exists():continue
  pending=True;manifest=json.loads((D/'manifest.json').read_text())
  if not all((Path(m['job'])/'done').exists() for m in manifest):continue
  # run_fit.py independently verifies exit statuses, inputs, q grids and CS sequences.
  f=(D/'fit.log').open('ab');p=subprocess.Popen([sys.executable,'-u',str(O/'run_fit.py'),name],stdout=f,stderr=subprocess.STDOUT,start_new_session=True)
  (D/'fit_launch.json').write_text(json.dumps({'pid':p.pid,'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'mode':'readiness-triggered, never duplicate/retry without inspection'},indent=2));print(name,'launched',p.pid,flush=True)
 if not pending:break
 time.sleep(30)
