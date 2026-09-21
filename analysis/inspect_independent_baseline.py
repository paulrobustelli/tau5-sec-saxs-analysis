from inspect_data import *
from scipy.ndimage import median_filter
out=ROOT/'independent_baseline'
fig,axes=plt.subplots(2,4,figsize=(17,8))
for row,name in enumerate(['WT_092025','AA_052026']):
 q,Y,E=load(name)
 for j,(lo,hi) in enumerate([(.012,.035),(.04,.1),(.1,.25),(.5,1.)]):
  m=(q>=lo)&(q<=hi);t=Y[m].mean(0);ax=axes[row,j]
  ax.plot(t,alpha=.4,lw=.6);ax.plot(median_filter(t,size=7),lw=1)
  ax.set(title=f'{name} q={lo}–{hi}',xlabel='Frame',ylabel='Mean facility-subtracted I')
  relevant=np.r_[30:115,210:250];v=t[relevant];pad=max(np.std(v)*2,.02);ax.set_ylim(np.median(v)-pad,np.median(v)+pad);ax.set_xlim(0,449)
  print(name,(lo,hi),[(a,b,round(t[a:b+1].mean(),5),round(t[a:b+1].std(),5)) for a,b in [(20,49),(50,79),(80,109),(100,119),(210,229),(230,249),(320,349),(400,429)]])
fig.tight_layout();fig.savefig(out/'raw_baseline_inventory.png',dpi=170)
