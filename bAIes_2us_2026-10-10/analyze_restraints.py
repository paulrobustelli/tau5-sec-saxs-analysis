import argparse,json
from pathlib import Path
import numpy as np,pandas as pd,mdtraj as md
p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--repo',type=Path,required=True);a=p.parse_args();out=a.root/'Analysis';rows=[]
for v,start in [('WT',4695000),('AA',4675000)]:
 params=np.atleast_2d(np.loadtxt(a.root/v/'baies_params.dat'));pairs=params[:,1:3].astype(int)-1;dist=[]
 for filename,offset in [('production_2us.xtc',4650000),('production_2us_gpu.xtc',4750000)]:
  for t in md.iterload(str(a.root/v/filename),top=str(a.root/v/'protein.pdb'),chunk=500):
   absps=t.time+offset*.002;mask=(absps>start*.002)&(absps<= (4750000 if filename=='production_2us.xtc' else start+1000000000)*.002+1e-3)
   if mask.any():dist.append(md.compute_distances(t[mask],pairs,periodic=False))
 d=np.concatenate(dist);np.savez_compressed(out/f'{v}_restraint_distances.npz',distances_nm=d,parameters=params)
 for i,row in enumerate(params):
  x=d[:,i];rows.append(dict(variant=v,restraint_id=int(row[0]),atom_i=int(row[1]),atom_j=int(row[2]),mu_nm=row[3],sigma_nm=row[4],mean_nm=x.mean(),SD_nm=x.std(ddof=1),mean_minus_mu_nm=x.mean()-row[3],fraction_within_one_sigma=np.mean(abs(x-row[3])<row[4])))
pd.DataFrame(rows).to_csv(out/'restraint_summary.csv',index=False)
source=json.load(open(a.repo/'Ensemble_Comparisons/Sa_Rg/summary.json'));comparison=[]
for label,values in source.items():
 for scope,d in values.items():comparison.append({'ensemble':label,'scope':scope,**d})
for v in ['WT','AA']:
 df=pd.read_csv(out/f'{v}_timeseries.csv')
 for scope,sa,rg in [('full','Sa_full','Rg_CA_A'),('R2R3','Sa_AR391_446','Rg_R2R3_CA_A')]:
  comparison.append(dict(ensemble=v+' bAIes 2us',scope=scope,n=len(df),mean_Salpha=df[sa].mean(),mean_Rg_nm=df[rg].mean()/10,Sa_Rg_correlation=df[sa].corr(df[rg]),fraction_Rg_lt_1p3_Sa_gt_6=float(((df[rg]<13)&(df[sa]>6)).mean()) if scope=='R2R3' else None))
pd.DataFrame(comparison).to_csv(out/'ensemble_comparison.csv',index=False)
# Correct default MDTraj Rg label: default assumes equal atom weights.
for file in out.glob('*.csv'):
 text=file.read_text();file.write_text(text.replace('Rg_allatom_mass_A','Rg_allatom_equalweight_A'))
