from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
O=Path(__file__).resolve().parent/'final_region_analysis'
fig,axes=plt.subplots(1,4,figsize=(18,13),width_ratios=[1.5,1,1.5,1],layout='constrained')
for k,v in enumerate(['WT','AA']):
 d=pd.read_csv(O/f'{v}_restraint_agreement.csv'); b=np.load(O/f'{v}_restraint_distance_blocks.npz'); y=np.arange(len(d));a,c=axes[2*k:2*k+2]
 assert b['block_mean_nm'].shape==(10,len(d))
 a.errorbar(d.mu_nm*10,y,xerr=d.sigma_nm*10,fmt='o',color='gray',alpha=.7,label='Target μ ± σ',capsize=2,ms=3)
 a.scatter(d.mean_nm*10,y,color='#0072b2',s=18,label='Production mean',zorder=3)
 a.set_yticks(y,[f'{i}–{j}' for i,j in zip(d.AR_i,d.AR_j)]);a.invert_yaxis();a.set(xlabel='Distance (Å)',ylabel='AR residue pair',title=f'{v}: {len(d)} restraints');a.legend(fontsize=9)
 p=d.fraction_within_mu_sigma.to_numpy()*100; blocks=b['block_within_mu_sigma']*100
 c.barh(y,p,color='#0072b2',alpha=.65)
 c.errorbar(p,y,xerr=np.vstack([p-blocks.min(0),blocks.max(0)-p]),fmt='none',ecolor='black',lw=.7,capsize=2)
 c.set(yticks=y,yticklabels=[],xlim=(0,100),xlabel='Frames within μ ± σ (%)',title=f'{v}: mean occupancy {p.mean():.1f}%');c.invert_yaxis();c.grid(axis='x',alpha=.2)
fig.suptitle('Harmonic bAIes production restraint agreement — full 2 μs per variant\n20,000 frames; occupancy whiskers: range across ten 200 ns blocks\nSoft restraints: μ ± σ occupancy is descriptive, not a hard pass/fail test',fontsize=14)
for ext in ['png','pdf']:fig.savefig(O/f'production_restraint_satisfaction.{ext}',dpi=180)
print(O/'production_restraint_satisfaction.png')
