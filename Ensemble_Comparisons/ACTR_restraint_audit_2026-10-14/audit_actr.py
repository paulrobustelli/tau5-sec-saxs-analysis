from pathlib import Path
import numpy as np,pandas as pd,json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[2];OUT=Path(__file__).resolve().parent
base=ROOT/'work/baies-upstream/benchmark/bAIes/ACTR';atoms={}
for line in (base/'ACTR.pdb').read_text().splitlines():
 if line.startswith('ATOM'):atoms[int(line[6:11])]=(int(line[22:26]),line[17:20],line[12:16].strip())
a=np.loadtxt(base/'baies_gauss_matrix.dat',comments='#');rows=[]
for row in a:
 i,ri,ai=atoms[int(row[1])];j,rj,aj=atoms[int(row[2])]
 assert ai==('CA' if ri=='GLY' else 'CB') and aj==('CA' if rj=='GLY' else 'CB')
 rows.append({'Id':int(row[0]),'residue_i':i,'residue_j':j,'type_i':ri,'type_j':rj,'seqsep':abs(j-i),'mu_A':row[3]*10,'sigma_A':row[4]*10})
df=pd.DataFrame(rows);df.to_csv(OUT/'ACTR_restraints_residue_mapped.csv',index=False)
sumry={'total':len(df),'minimum_separation':int(df.seqsep.min()),'separation_counts':{str(k):int(v) for k,v in df.seqsep.value_counts().sort_index().items()},'local_i_i4':df[df.seqsep==4][['mu_A','sigma_A']].describe().to_dict(),'other_pairs':df[df.seqsep>4][['mu_A','sigma_A']].describe().to_dict(),'sequence':' '.join(f'{i}:{r}' for i,r in sorted({(x[0],x[1]) for x in atoms.values()}))}
(OUT/'summary.json').write_text(json.dumps(sumry,indent=2));print(json.dumps(sumry,indent=2));print(df[df.seqsep==4].to_string(index=False))
fig,axs=plt.subplots(1,2,figsize=(11,4),layout='constrained');local=df[df.seqsep==4]
axs[0].errorbar(local.residue_i,local.mu_A,yerr=local.sigma_A,fmt='o',capsize=2);axs[0].set(xlabel='ACTR construct residue i',ylabel='Fitted mean ± sigma (Å)',title='Published ACTR i→i+4 Gaussian restraints')
for ax,subset,label,color in [(axs[1],local,'i→i+4','tab:blue'),(axs[1],df[df.seqsep>4],'separation >4','tab:orange')]:ax.scatter(subset.mu_A,subset.sigma_A,label=label,color=color,alpha=.8)
axs[1].set(xlabel='Fitted mean (Å)',ylabel='Fitted sigma (Å)',title='Gaussian widths');axs[1].legend();fig.savefig(OUT/'ACTR_restraint_parameters.png',dpi=160)
# Gaussian curves reconstructed from published fitted parameters, not raw logits.
examples=local.sort_values('sigma_A').iloc[[0,len(local)//2,-1]];fig,axs=plt.subplots(1,3,figsize=(11,3),layout='constrained')
for ax,(_,r) in zip(axs,examples.iterrows()):
 x=np.linspace(max(0,r.mu_A-4*r.sigma_A),r.mu_A+4*r.sigma_A,400);y=np.exp(-.5*((x-r.mu_A)/r.sigma_A)**2)/(r.sigma_A*np.sqrt(2*np.pi));ax.plot(x,y);ax.set_title(f'{int(r.residue_i)}–{int(r.residue_j)}: μ={r.mu_A:.2f}, σ={r.sigma_A:.2f} Å');ax.set_xlabel('Distance (Å)')
axs[0].set_ylabel('Reconstructed Gaussian density');fig.savefig(OUT/'ACTR_example_fitted_gaussians.png',dpi=160)
