from pathlib import Path
import numpy as np,pandas as pd,json,re
R=Path(__file__).resolve().parents[2];O=Path(__file__).resolve().parent;S=R/'outputs/Ensemble_Comparison_121_2026-10-15'
d=pd.read_csv(S/'matched_saxs_curves.csv');rows=[]
for n,g in d.groupby('ensemble'):
 for stage,col in [('prior','I_uniform'),('after','I_MaxEnt')]:
  q=g.q_Ainv.to_numpy();y=g.I_exp.to_numpy();e=g.sigma_exp.to_numpy();p=g[col].to_numpy()
  for band,mask in [('low',q<=.08),('middle',(q>.08)&(q<.25)),('high',q>=.25),('tail',q>=.4)]:
   a=np.sum(p[mask]*y[mask]/e[mask]**2)/np.sum((p[mask]/e[mask])**2)
   rows.append(dict(ensemble=n,stage=stage,band=band,N=int(mask.sum()),mean_ratio=float(np.mean(p[mask]/y[mask])),mean_z=float(np.mean((p[mask]-y[mask])/e[mask])),chi2mean=float(np.mean(((p[mask]-y[mask])/e[mask])**2)),band_only_scale_multiplier=float(a),band_only_rescaled_chi2mean=float(np.mean(((a*p[mask]-y[mask])/e[mask])**2))))
pd.DataFrame(rows).to_csv(O/'residual_bands.csv',index=False)
print(pd.DataFrame(rows).query('band=="high"').to_string(index=False))
components=[];counts=[]
for n in ['IDPForge_WT','IDPForge_AA','bAIes_WT','bAIes_AA','REST_WT']:
 manifest=json.loads((S/n/'manifest.json').read_text());mats=[];formulas=[]
 for m in manifest:
  job=Path(m['job']);p=np.loadtxt(job/'prediction.dat');assert p.shape==(601,8);mats.append(p[[0,100,250,400,500],1:]);log=(job/'run.log').read_text()
  formulas.append(tuple(line.split(':',1)[-1].strip() for line in log.splitlines() if line.startswith(('Number of atoms read','Chemical formula','r0:'))))
 vals=np.array(mats);weights=pd.read_csv(S/n/'weights.csv').combined.to_numpy()
 for stage,w in [('prior',np.ones(len(vals))/len(vals)),('after',weights)]:
  mean=np.einsum('i,ijk->jk',w,vals)
  for k,q in enumerate([0,.1,.25,.4,.5]):components.append(dict(ensemble=n,stage=stage,q=q,**dict(zip(['I','Iat','Iev','Ihs','AatAev','AatAhs','AevAhs'],mean[k]))))
 counts.append(dict(ensemble=n,unique_atom_formula_records=list(set(formulas))))
pd.DataFrame(components).to_csv(O/'native_component_means.csv',index=False);(O/'native_formula_audit.json').write_text(json.dumps(counts,indent=2))
