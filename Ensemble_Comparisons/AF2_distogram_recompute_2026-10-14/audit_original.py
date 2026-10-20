from pathlib import Path
import pickle,json,hashlib,importlib.util,re
import numpy as np,pandas as pd,mdtraj as md
from scipy.special import softmax
from scipy.optimize import curve_fit
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[2];OUT=Path(__file__).resolve().parent
s=importlib.util.spec_from_file_location('preprocess',ROOT/'outputs/Tau5_bAIes_OpenMM_restart_2026-10-09/code/preprocess_bAIes.py');p=importlib.util.module_from_spec(s);s.loader.exec_module(p)
summary=[];allrows=[]
for v in ['WT','AA']:
 files=sorted(f for f in (ROOT/'work/af2_original_dropbox').rglob(f'*{v}*all_rank_*.pickle') if not f.name.startswith('._'))
 base=ROOT/'outputs/Tau5_bAIes_OpenMM_restart_2026-10-09'/v
 old=np.loadtxt(base/'baies_params.dat',comments='#');mddata=p.read_input(str(base/'protein.pdb'));byatom={a[0]:i for i,a in mddata.items()};oldpairs={(byatom[int(r[1])],byatom[int(r[2])]):r[3:] for r in old}
 for f in files:
  rank=int(re.search('rank_(\d+)',f.name).group(1));d=pickle.load(f.open('rb'))['distogram'];prob=softmax(d['logits'].astype(float),axis=-1);edges=d['bin_edges'].astype(float);centers=(edges-(edges[1]-edges[0])/2)/10
  pdb=f.with_name(f.name.replace('_all_','_unrelaxed_').replace('.pickle','.pdb'));t=md.load(str(pdb));ss=md.compute_dssp(t,simplified=False)[0];assert t.top.to_fasta()[0]==md.load(str(base/'protein.pdb')).top.to_fasta()[0]
  indices=[next(a.index for a in r.atoms if a.name==('CA' if r.name=='GLY' else 'CB')) for r in t.top.residues]
  xyz=t.xyz[0,indices];res=p.read_input(str(pdb));selected=[];local=[]
  for i in range(120):
   for j in range(i+1,120):
    q=prob[i,j];m=int(q.argmax());mode=centers[m] if m<63 else float('inf');types=(res[i+1][1],res[j+1][1]);key=types if types in p.cutoff_matrix else types[::-1];cut=float(p.cutoff_matrix[key][0])/10
    keep=(j-i>3 and mode<cut)
    if keep:selected.append((i+1,j+1))
    if j-i in [3,4,5,6] and 63<=i<80:
     row={'variant':v,'rank':rank,'model_i':i+1,'model_j':j+1,'AR_i':i+328,'AR_j':j+328,'seqsep':j-i,'DSSP_i':ss[i],'DSSP_j':ss[j],'coordinate_CB_nm':float(np.linalg.norm(xyz[i]-xyz[j])),'mode_nm':mode,'matrix_cutoff_nm':cut,'prob_below_cutoff':float(q[:63][centers<cut].sum()),'overflow_probability':float(q[-1]),'selected':keep,'in_original_restraints':(i+1,j+1) in oldpairs,'reason':'seqsep <= 3' if j-i<=3 else 'mode >= cutoff' if mode>=cut else 'retained'}
     allrows.append(row);local.append((row,q))
  item={'variant':v,'rank':rank,'pickle_sha256':hashlib.sha256(f.read_bytes()).hexdigest(),'selected_count':len(selected),'R2_AR394_402_DSSP':''.join(ss[66:75]),'selected_R2_pairs':[x for x in selected if 67<=x[0]<=75 and x[1]<=75]}
  if rank==1:
   assert np.array_equal(t.xyz,md.load(str(base/'prediction.pdb')).xyz),'rank1 coordinates differ'
   item['original_pair_set_exact_match']=set(selected)==set(oldpairs)
   bins,pr=p.read_pkl(str(f));fit=p.select_and_fit_distograms(bins,pr,res,mddata,'matrix',seqsep=3)
   p.write_output(fit,output_file=str(OUT/f'{v}_recomputed_original_params.dat'),ndx_file=str(OUT/f'{v}_atoms.ndx'))
   new=np.loadtxt(OUT/f'{v}_recomputed_original_params.dat',comments='#');item['max_mu_sigma_delta_nm']=float(np.max(abs(new[:,3:]-old[:,3:])))
   fig,axs=plt.subplots(3,3,figsize=(12,9),layout='constrained')
   for ax,(row,q) in zip(axs.flat,[x for x in local if x[0]['seqsep']==4][:9]):
    ax.plot(centers,q[:63]);ax.axvline(row['coordinate_CB_nm'],color='green',label='AF2 structure');ax.axvline(row['matrix_cutoff_nm'],color='red',ls=':',label='cutoff');ax.set_title(f"AR{row['AR_i']}–{row['AR_j']}: {row['reason']}");ax.set_xlim(.3,1.4)
   axs.flat[0].legend(fontsize=7);fig.savefig(OUT/f'{v}_R2_i_i4_distograms.png',dpi=150);plt.close(fig)
  summary.append(item)
pd.DataFrame(allrows).to_csv(OUT/'R2_contact_selection_all_models.csv',index=False);(OUT/'recompute_summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))
