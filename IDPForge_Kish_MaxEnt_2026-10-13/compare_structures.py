from pathlib import Path
import ast,json,numpy as np,pandas as pd,mdtraj as md
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(__file__).resolve().parents[2];O=Path(__file__).resolve().parent
src=R/'outputs/Tau5_bAIes_OpenMM_restart_2026-10-09/Analysis/analyze_production.py';tr=ast.parse(src.read_text());ns={'np':np}
for n in tr.body:
 if isinstance(n,ast.FunctionDef) and n.name in ['rmsd','sa','rg']:exec(compile(ast.Module(body=[n],type_ignores=[]),str(src),'exec'),ns)
ref=md.load(str(R/'outputs/Tau5_bAIes_OpenMM_restart_2026-10-09/Analysis/R2R3_all_helix.pdb'));ns['refs']=np.array([ref.xyz[0,ref.top.select(f'residue {i} to {i+6} and backbone')] for i in range(1,51)],float)
rows=[];profiles={}
for v in ['WT','AA']:
 ws=pd.read_csv(O/f'{v}_weights.csv');traj=md.load(ws.file.tolist());h=md.compute_dssp(traj,simplified=False)=='H';top=traj.top
 idx=np.array([[next(a.index for a in res.atoms if a.name==n) for n in ['N','CA','C','O']] for res in top.residues]);bb=traj.xyz[:,idx].astype(float);sa,cropsa=ns['sa'](bb);fullrg=ns['rg'](traj.xyz[:,idx[:,1]]);croprg=ns['rg'](traj.xyz[:,idx[63:119,1]])
 pairs=np.array([(i,j) for i in range(120) for j in range(i+1,120)]);d,_=md.compute_contacts(traj,contacts=pairs,scheme='closest-heavy',periodic=False);c=d<1.2
 obs=pd.DataFrame({'file':ws.file,'Sa_full':sa,'Sa_AR391_446':cropsa,'Rg_CA_A':fullrg,'Rg_R2R3_CA_A':croprg});obs.to_csv(O/f'{v}_structural_observables.csv',index=False)
 for label,w in [('uniform',ws.uniform.to_numpy()),('Kish_MaxEnt',ws.combined.to_numpy())]:
  profiles[v+' '+label]=w@h[:,2:]
  prob=np.einsum('i,ij->j',w,c);m=np.zeros((120,120));m[pairs[:,0],pairs[:,1]]=prob;m+=m.T;np.savetxt(O/f'{v}_{label}_contact_map.csv',m,delimiter=',')
  for scope,rg,salpha,hh in [('full',fullrg,sa,h[:,2:]),('R2R3',croprg,cropsa,h[:,63:119])]:
   rows.append({'ensemble':v+' IDPForge '+label,'view':scope,'n':200,'effective_size':float(1/np.dot(w,w)),'Rg_mean_A':float(w@rg),'Sa_mean':float(w@salpha),'DSSP_alpha_fraction':float(w@hh.mean(1))})
# Reuse frozen structural summaries for historical ensembles; label unavailable values explicitly.
old=pd.read_csv(R/'work/asteroids_reference/ensemble_summary.csv')
for row in old.itertuples():
 rows.append({'ensemble':row.ensemble,'view':row.view,'n':row.n,'effective_size':None,'Rg_mean_A':row.Rg_mean_nm*10,'Sa_mean':row.Sa_mean,'DSSP_alpha_fraction':None,'source':'archived ASTEROIDS master summary'})
for v in ['WT','AA']:
 b=pd.read_csv(R/f'outputs/Tau5_bAIes_OpenMM_restart_2026-10-09/Analysis/{v}_timeseries.csv')
 for scope,rg,salpha in [('full','Rg_CA_A','Sa_full'),('R2R3','Rg_R2R3_CA_A','Sa_AR391_446')]:
  rows.append({'ensemble':v+' bAIes 2us MD','view':scope,'n':len(b),'effective_size':None,'Rg_mean_A':b[rg].mean(),'Sa_mean':b[salpha].mean(),'DSSP_alpha_fraction':b.DSSP_alpha_fraction.mean() if scope=='full' else None,'source':'previously computed bAIes trajectory observables'})
pd.DataFrame(rows).to_csv(O/'structural_ensemble_comparison.csv',index=False)
fig,axes=plt.subplots(2,1,figsize=(11,7),layout='constrained')
for ax,v in zip(axes,['WT','AA']):
 for label in ['uniform','Kish_MaxEnt']:ax.plot(np.arange(330,448),profiles[v+' '+label],label='IDPForge '+label)
 for pool in ['base','expanded','aggressive']:
  p=R/f'work/asteroids_reference/data/{v}_{pool}_ens1_helicity.csv'
  if p.exists():
   f=pd.read_csv(p);f=f[f.definition=='DSSP_H'];ax.plot(f.AR_residue,f.selected,label='ASTEROIDS '+pool,alpha=.7)
 b=pd.read_csv(R/f'outputs/Tau5_bAIes_OpenMM_restart_2026-10-09/Analysis/{v}_helicity.csv');ax.plot(b.AR_number,b.alpha_propensity,label='bAIes MD',ls='--')
 ax.set(title=v,xlabel='AR residue',ylabel='DSSP alpha propensity');ax.legend(fontsize=8)
fig.savefig(O/'helicity_ensemble_comparison.png',dpi=180);plt.close(fig)
fig,axes=plt.subplots(1,2,figsize=(11,4),layout='constrained')
for ax,v in zip(axes,['WT','AA']):
 ws=pd.read_csv(O/f'{v}_weights.csv');obs=pd.read_csv(O/f'{v}_structural_observables.csv')
 for label,w in [('uniform',ws.uniform),('MaxEnt',ws.combined)]:ax.hist(obs.Rg_CA_A,bins=np.linspace(10,70,31),weights=w,histtype='step',label='IDPForge '+label)
 b=pd.read_csv(R/f'outputs/Tau5_bAIes_OpenMM_restart_2026-10-09/Analysis/{v}_timeseries.csv');ax.hist(b.Rg_CA_A,bins=np.linspace(10,70,31),weights=np.ones(len(b))/len(b),histtype='step',label='bAIes MD');ax.set(title=v,xlabel='CA Rg (Å)',ylabel='Probability per bin');ax.legend(fontsize=8)
fig.savefig(O/'Rg_ensemble_comparison.png',dpi=180);plt.close(fig)
print('Structural comparisons complete')
