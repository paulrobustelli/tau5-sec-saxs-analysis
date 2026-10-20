from pathlib import Path
import json,subprocess,shutil
ROOT=Path(__file__).resolve().parents[2];audit=Path(__file__).resolve().parent;out=ROOT/'outputs/bAIes_8A_seqsep2_2us_2026-10-14';out.mkdir(exist_ok=True);pids={}
for v in ['WT','AA']:
 src=audit/f'{v}_8A_sep2_20ns';dst=out/v;dst.mkdir(exist_ok=True)
 if (dst/'launch.json').exists():raise RuntimeError('Already launched; use recorded job status before attempting resume')
 for name in ['protein.pdb','baies_params.dat','final_state.xml']:shutil.copyfile(src/name,dst/('start_state.xml' if name=='final_state.xml' else name))
 cmd=json.load(open(src/'command.json'))
 replacements={'-pdb':str(dst/'protein.pdb'),'-baies':str(dst/'baies_params.dat'),'-xtc':str(dst/'production.xtc'),'-Nsteps':'1000000000','-checkpoint':str(dst/'production.chk'),'-state_out':str(dst/'production_final.xml')}
 for flag,val in replacements.items():cmd[cmd.index(flag)+1]=val
 cmd.extend(['-restart_state',str(dst/'start_state.xml')])
 log=(dst/'run.log').open('ab');p=subprocess.Popen(cmd,stdout=log,stderr=subprocess.STDOUT,start_new_session=True);pids[v]=p.pid
 (dst/'launch.json').write_text(json.dumps({'pid':p.pid,'command':cmd,'initial_step':10000000,'initial_time_ns':20,'additional_time_ns':2000,'method_variant':'constant8A cutoff, seqsep2 includes i,i+3, original Jeffreys likelihood/kBT weight and OpenMM coil port','restraints':45 if v=='WT' else 37},indent=2))
 subprocess.Popen(['/usr/bin/caffeinate','-i','-w',str(p.pid)],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL,start_new_session=True)
(out/'pids.json').write_text(json.dumps(pids,indent=2));print(pids)
(audit/'production_decision.json').write_text(json.dumps({'decision':'launch exploratory selection-variant 2us WT/AA','reason':'matched20ns tests increase R2 H to21.12%WT/21.89%AA with8A+seqsep2; exact equilibrium/experimental correctness not established','not_exact_paper_reproduction':True,'extra_restraint_multiplier':False,'new_production_launched':True,'output':str(out),'pids':pids},indent=2))
