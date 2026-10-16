from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
O=Path(__file__).resolve().parent;df=pd.read_csv(O/'helicity_profiles_before_after.csv');names=list(df.ensemble.unique());fig,axes=plt.subplots(4,2,figsize=(13,12),sharex=True,sharey=True);axes=axes.ravel()
for ax,name in zip(axes,names):
 for stage,g in df[(df.ensemble==name)&df.AR_residue.notna()].groupby('stage',sort=False):
  ax.plot(g.AR_residue,g.DSSP_H_percent,label=stage,color='#777777' if stage=='prior' else '#d95f02',lw=1.4)
 for lo,hi,color in [(390,404,'#56b4e9'),(404,410,'#009e73'),(430,440,'#cc79a7')]:ax.axvspan(lo-.5,hi+.5,color=color,alpha=.13)
 ax.set(title=name.replace('_',' '),xlim=(330,448),ylim=(0,100),xlabel='AR residue',ylabel='DSSP H (%)');ax.legend(frameon=False,fontsize=8)
for ax in axes[len(names):]:ax.set_visible(False)
fig.suptitle('R2 main390–404 (blue), minor404–410 (green), R3430–440 (purple)',fontsize=12);fig.tight_layout();fig.savefig(O/'master_figures/helicity_before_after.png',dpi=180);fig.savefig(O/'master_figures/helicity_before_after.pdf');plt.close(fig)
