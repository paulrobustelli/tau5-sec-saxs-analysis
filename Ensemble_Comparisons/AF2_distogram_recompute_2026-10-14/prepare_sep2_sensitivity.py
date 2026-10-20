from pathlib import Path
import importlib.util,pickle,shutil,json,subprocess
from scipy.special import softmax
ROOT=Path(__file__).resolve().parents[2];OUT=Path(__file__).resolve().parent
s=importlib.util.spec_from_file_location('pre',ROOT/'outputs/Tau5_bAIes_OpenMM_restart_2026-10-09/code/preprocess_bAIes.py');p=importlib.util.module_from_spec(s);s.loader.exec_module(p)
run=ROOT/'outputs/Tau5_bAIes_OpenMM_restart_2026-10-09/code/bAIes_openMM.py';code=run.parent
print([x.name for x in code.glob('*')])
for v in ['WT','AA']:
 base=ROOT/'outputs/Tau5_bAIes_OpenMM_restart_2026-10-09'/v;f=next(f for f in (ROOT/'work/af2_original_dropbox').rglob(f'*{v}*all_rank_001*.pickle') if not f.name.startswith('._'))
 d=pickle.load(f.open('rb'))['distogram'];edges=d['bin_edges'].astype(float);bins=((edges-(edges[1]-edges[0])/2)*.1).tolist();prob=softmax(d['logits'].astype(float),axis=-1)
 res=p.read_input(str(base/'prediction.pdb'));mdres=p.read_input(str(base/'protein.pdb'))
 fit=p.select_and_fit_distograms(bins,prob,res,mdres,8.,seqsep=2)
 path=OUT/f'{v}_8A_sep2_params.dat';p.write_output(fit,output_file=str(path),ndx_file=str(OUT/f'{v}_8A_sep2_atoms.ndx'));print(v,'8A',len(fit))
 # Original matrix condition is preserved exactly, rather than relabeled as an implementation fix.
 for condition in ['8A_sep2']:
  dest=OUT/f'{v}_{condition}_20ns';dest.mkdir(exist_ok=True);shutil.copyfile(base/'protein.pdb',dest/'protein.pdb')
  params=base/'baies_params.dat' if condition=='matrix' else path;shutil.copyfile(params,dest/'baies_params.dat')
  ff=next(code.glob('coil*.xml'));cm=next(code.glob('cmap*.cmap'))
  cmd=[str(ROOT/'work/baies_env_compatible/bin/python'),'-u',str(run),'-ff',str(ff),'-pdb',str(dest/'protein.pdb'),'-cmap',str(cm),'-baies',str(dest/'baies_params.dat'),'-xtc',str(dest/'trajectory.xtc'),'-Nsteps','10000000','-save_xtc','50000','-save_stdout','50000','-checkpoint',str(dest/'checkpoint.chk'),'-state_out',str(dest/'final_state.xml'),'-platform','OpenCL','-seed',str(2026101400+(v=='AA'))]
  (dest/'command.json').write_text(json.dumps(cmd,indent=2))
  log=(dest/'run.log').open('ab');job=subprocess.Popen(cmd,stdout=log,stderr=subprocess.STDOUT,start_new_session=True);(dest/'pid').write_text(str(job.pid));subprocess.Popen(['/usr/bin/caffeinate','-i','-w',str(job.pid)],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL,start_new_session=True)
  print(v,condition,job.pid,flush=True)
