"""Display anchor locations and baseline mismatch relative to protein signal."""
from pathlib import Path
import numpy as np,json,matplotlib.pyplot as plt
from inspect_data import load
from component_recovery import SPECS
ROOT=Path(__file__).resolve().parents[1]
def plot_baseline_anchors():
 fig,axes=plt.subplots(2,2,figsize=(11,7))
 for col,(sample,spec) in enumerate(SPECS.items()):
  q,Y,E=load(sample);anchors=spec['anchors'];bs=[Y[:,lo:hi+1].mean(1) for lo,hi in anchors];mid=[(lo+hi)/2 for lo,hi in anchors];m=(q>=.012-1e-9)&(q<=.15+1e-9)
  ax=axes[0,col];ax.plot(np.arange(Y.shape[1]),np.trapezoid(Y[m],q[m],axis=0),color='0.25',lw=1,label='Provided input trace')
  for (lo,hi),label,color in zip(anchors,['pre anchor','post anchor'],['C0','C2']):ax.axvspan(lo,hi,color=color,alpha=.3,label=f'{label}: {lo}–{hi}')
  main=(160,169) if sample.startswith('WT') else (170,179);tail=(180,189) if sample.startswith('WT') else (190,199)
  ax.axvspan(*main,color='C1',alpha=.3,label=f'center: {main[0]}–{main[1]}');ax.set(xlim=(0,280),xlabel='Frame (zero-based)',ylabel='Integrated input intensity',title=sample[:2]+' — residual-background anchor frames');ax.legend(fontsize=7)
  m=(q>=.012-1e-9)&(q<=.055+1e-9)
  for (lo,hi),label,color in [(main,'center','C1'),(tail,'late tail','C3')]:
   frac=np.clip(((lo+hi)/2-mid[0])/(mid[1]-mid[0]),0,1);signal=Y[:,lo:hi+1].mean(1)-((1-frac)*bs[0]+frac*bs[1]);ratio=100*(bs[1]-bs[0])/signal
   axes[1,col].plot(q[m],ratio[m],'.-',color=color,label=label)
  axes[1,col].axhline(0,color='0.5',lw=.5);axes[1,col].set(xlabel='q (Å⁻¹)',ylabel='100 × (Bpost − Bpre) / corrected I(q)',title='Same background mismatch, larger relative effect in tail');axes[1,col].legend(fontsize=8)
 fig.tight_layout();out=ROOT/'results_comparison/baseline_anchor_explanation.png';fig.savefig(out,dpi=150);return fig
