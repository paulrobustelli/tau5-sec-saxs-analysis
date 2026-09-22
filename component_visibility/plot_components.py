from pathlib import Path
import os
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/saxs_mpl')
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path(__file__).resolve().parent
fig,axes=plt.subplots(2,2,figsize=(11,8))
for row,s in enumerate(['WT_092025','AA_052026']):
 for col,suffix in enumerate(['EFA_candidate_main','EFA_other_component']):
  q,y,e=np.loadtxt(P/f'{s}_{suffix}.dat',unpack=True);m=(q>=.012)&(q<=.25);a=axes[row,col]
  a.errorbar(q[m],y[m],e[m],fmt='.',ms=3,elinewidth=.6)
  if row==1 and col==1:a.set_yscale('symlog',linthresh=.005);a.axhline(0,color='k',lw=.7)
  else:a.set_yscale('log')
  a.set(xlabel='q (Å⁻¹)',ylabel='I(q), component scale',title=s[:2]+(' main component' if col==0 else ' earlier/secondary component'))
fig.suptitle('Actual EFA component I(q): no smoothing or positivity clipping');fig.tight_layout();fig.savefig(P/'all_component_Iq.png',dpi=170)
q,y,e=np.loadtxt(P/'WT_092025_EFA_other_component.dat',unpack=True)
fig,ax=plt.subplots(1,2,figsize=(11,4));m=(q>=.012)&(q<=.25)
ax[0].errorbar(q[m],y[m],e[m],fmt='.',ms=3,lw=.6);ax[0].set(yscale='log',xlabel='q (Å⁻¹)',ylabel='WT shoulder I(q)',title='WT earlier component: positive, structured signal')
m=(q>=.01)&(q<=.06);ax[1].errorbar(q[m]**2,y[m],e[m],fmt='.',ms=4,lw=.7)
nu=.6148307989724572;rg=40.86640256343632;i0=.128224660766131
x=np.linspace(.012,.048,150);ax[1].plot(x*x,i0*np.exp(-(x*rg)**2/3+.0479*(nu-.212)*(x*rg)**4),label='Extended fit: Rg 40.9 Å',color='C1');ax[1].set(xlabel='q² (Å⁻²)',ylabel='I(q)',title='Exploratory single-chain fit; low-q excess remains');ax[1].legend()
fig.tight_layout();fig.savefig(P/'WT_shoulder_Iq_extended.png',dpi=170)
