from pathlib import Path
import numpy as np,pandas as pd,mdtraj as md,json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
OUT=Path(__file__).resolve().parent
fig,axs=plt.subplots(2,1,figsize=(11,6),layout='constrained');rows=[];energy=[]
for ax,v in zip(axs,['WT','AA']):
 for c in ['matrix','8A','8A_sep2']:
  d=OUT/f'{v}_{c}_20ns';obs=pd.read_csv(d/'observables.csv');ax.plot(obs.time_ns,obs.R2_AR394_402_H,label=c,alpha=.8)
  pars=np.loadtxt(d/'baies_params.dat',comments='#');t=md.load(str(d/'trajectory.xtc'),top=str(d/'protein.pdb'));pairs=pars[:,1:3].astype(int)-1;dist=md.compute_distances(t,pairs);last=dist[100:]
  df=pd.DataFrame({'atom_i':pars[:,1].astype(int),'atom_j':pars[:,2].astype(int),'AR_i':[t.top.atom(i).residue.index+328 for i in pairs[:,0]],'AR_j':[t.top.atom(j).residue.index+328 for j in pairs[:,1]],'mu_nm':pars[:,3],'sigma_nm':pars[:,4],'mean_distance_nm':last.mean(axis=0),'std_distance_nm':last.std(axis=0),'within_mu_sigma':np.mean(abs(last-pars[:,3])<=pars[:,4],axis=0)})
  df.to_csv(d/'restraint_distance_audit.csv',index=False)
  h=obs[obs.time_ns>10].R2_AR394_402_H.to_numpy();blocks=h.reshape(5,20).mean(axis=1)
  rows.append({'variant':v,'condition':c,'mean_R2_H_10_20ns':float(h.mean()),'2ns_block_means':blocks.tolist(),'block_SD_descriptive_not_SEM':float(blocks.std(ddof=1)),'retained_pairs_inside_AR394_402':int((df.AR_i.between(394,402)&df.AR_j.between(394,402)).sum())})
  records=[]
  for line in (d/'run.log').read_text().splitlines():
   fields=line.split('\t')
   if len(fields)==5 and fields[0].isdigit():records.append([float(x) for x in fields[:3]])
  a=np.array(records);assert len(a)==200 and a[-1,0]==10000000 and np.isfinite(a).all()
  energy.append({'variant':v,'condition':c,'completed_steps':int(a[-1,0]),'last_half_temperature_K':float(a[100:,2].mean()),'last_half_energy_range_kJmol':[float(a[100:,1].min()),float(a[100:,1].max())]})
 ax.set_title(v+' R2 AR394–402');ax.set_ylabel('DSSP alpha fraction');ax.legend();ax.set_ylim(-.02,1.02)
axs[-1].set_xlabel('Time (ns)');fig.savefig(OUT/'R2_helicity_selection_sensitivity.png',dpi=160)
(OUT/'block_variation.json').write_text(json.dumps(rows,indent=2));(OUT/'energy_stability_checks.json').write_text(json.dumps(energy,indent=2));print(json.dumps(rows,indent=2))
