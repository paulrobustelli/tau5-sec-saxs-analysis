from pathlib import Path
import numpy as np,pandas as pd,mdtraj as md
O=Path(__file__).resolve().parent;A=O/'final_region_analysis'
for v in ['WT','AA']:
 p=np.loadtxt(O/v/'baies_params.dat');top=md.load(str(O/v/'protein.pdb')).top;pairs=p[:,1:3].astype(int)-1
 distances=[]
 for t in md.iterload(str(O/v/'production.xtc'),top=top,chunk=200):distances.append(md.compute_distances(t,pairs,periodic=False))
 d=np.concatenate(distances);assert len(d)==20000
 data={'id':p[:,0].astype(int),'atom_i':p[:,1].astype(int),'atom_j':p[:,2].astype(int),'AR_i':[top.atom(i).residue.index+328 for i in pairs[:,0]],'AR_j':[top.atom(i).residue.index+328 for i in pairs[:,1]],'mu_nm':p[:,3],'sigma_nm':p[:,4],'mean_nm':d.mean(0),'SD_nm':d.std(0),'fraction_within_mu_sigma':(abs(d-p[:,3])<=p[:,4]).mean(0)}
 pd.DataFrame(data).to_csv(A/f'{v}_restraint_agreement.csv',index=False)
 np.savez_compressed(A/f'{v}_restraint_distance_blocks.npz',block_mean_nm=d.reshape(10,2000,-1).mean(1),block_within_mu_sigma=(abs(d-p[:,3])<=p[:,4]).reshape(10,2000,-1).mean(1))
 print(v,len(p),'restraints checked',flush=True)
