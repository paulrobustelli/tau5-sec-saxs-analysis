from pathlib import Path
import numpy as np,pandas as pd,mdtraj as md
O=Path(__file__).resolve().parent;R=O.parents[1]
regions={'R1':(350,360),'R2_main':(390,404),'R2_minor':(404,410),'R3':(430,440)};rows=[]
for v in ['WT','AA']:
 top=md.load(str(R/f'outputs/bAIes_8A_seqsep2_2us_2026-10-14/{v}/protein.pdb')).top
 for condition in ['matrix','8A','8A_sep2']:
  params=np.loadtxt(O/f'{v}_{condition}_20ns/baies_params.dat',comments='#')
  for p in params:
   i,j=map(int,p[1:3]);ai,aj=top.atom(i-1),top.atom(j-1);ri,rj=ai.residue.index+328,aj.residue.index+328
   assert ai.serial==i and aj.serial==j
   row=dict(variant=v,condition=condition,atom_i=i,atom_j=j,AR_i=ri,AR_j=rj,residue_i=ai.residue.name,residue_j=aj.residue.name,sequence_separation=abs(rj-ri),mu_nm=p[3],sigma_nm=p[4])
   for a,(lo,hi) in regions.items():
    for b,(lo2,hi2) in regions.items():
     if a>b:continue
     row[a+'__'+b]=bool((lo<=ri<=hi and lo2<=rj<=hi2) or (lo<=rj<=hi and lo2<=ri<=hi2))
   rows.append(row)
f=pd.DataFrame(rows);f.to_csv(O/'region_pair_selection_detail.csv',index=False)
cols=[c for c in f if '__' in c];f.groupby(['variant','condition'])[cols].sum().to_csv(O/'region_pair_selection_counts.csv')
print(f.groupby(['variant','condition'])[cols].sum().to_string())
