from pathlib import Path
import os,json,numpy as np
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/saxs_mpl')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
p=Path(__file__).resolve().parent
r=json.loads((p/'summary_histograms.json').read_text())
fig,axes=plt.subplots(1,2,figsize=(11,4.2),sharex=True,sharey=True)
for ax,n in zip(axes,[118,120]):
    for label,color in [('WT','#2471a3'),('AA','#d97706')]:
        d=r[f'{label}_{n}'];counts=np.array(d['histogram_counts']);edges=np.array(d['edges_A'])
        ax.stairs(counts/(d['n']*np.diff(edges)),edges,color=color,linewidth=2,label=f"{label}: mean {d['mean_Rg_A']:.1f} Å")
        ax.axvline(d['mean_Rg_A'],color=color,ls='--',alpha=.65)
    ax.set(xlim=(14,60),xlabel='Conformational Rg (Å)',title=f'{n} residues');ax.legend(frameon=False)
axes[0].set_ylabel('Probability density (Å⁻¹)')
fig.suptitle('STARLING 2.0.2 • 400 conformations per construct • 188 mM',fontsize=12)
fig.tight_layout();fig.savefig(p/'Rg_histograms.png',dpi=180);fig.savefig(p/'Rg_histograms.pdf')
