from pathlib import Path
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
O=Path(__file__).resolve().parent
d=pd.read_csv(O/'REST_WT_cis_counts_vs_time.csv');x=d.all_peptide_bonds_cis.to_numpy();k=np.arange(x.max()+1);n=np.bincount(x,minlength=len(k));pct=n/len(x)*100
pd.DataFrame({'cis_bonds_per_frame':k,'frames':n,'percent':pct}).to_csv(O/'REST_WT_cis_count_histogram.csv',index=False)
fig,ax=plt.subplots(figsize=(7.5,4.8),layout='constrained');ax.bar(k,pct,width=.78,color='#0072b2')
for i,p,c in zip(k,pct,n):ax.text(i,p+.7,f'{p:.1f}%\n({c})',ha='center',fontsize=10)
ax.axvline(x.mean(),color='#d55e00',ls='--',lw=1.5,label=f'Mean = {x.mean():.2f}')
ax.set(xticks=k,xlabel='Number of cis peptide bonds per frame',ylabel='Frames (%)',title='WT REST 300 K — cis-bond count distribution',ylim=(0,pct.max()+8));ax.legend(frameon=False);ax.spines[['top','right']].set_visible(False)
fig.supxlabel(f'All 120 peptide bonds; cis: |ω| < 30°. {len(x):,} unweighted frames, 0–746 ns.',fontsize=10)
for ext in ['png','pdf']:fig.savefig(O/f'REST_WT_cis_count_histogram.{ext}',dpi=180)
print(list(zip(k.tolist(),n.tolist(),pct.tolist())))
