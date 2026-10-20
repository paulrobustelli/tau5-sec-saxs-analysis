from pathlib import Path
import pandas as pd,numpy as np
from scipy.signal import fftconvolve
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
O=Path(__file__).resolve().parent/'final_region_analysis';rows=[]
fig,axes=plt.subplots(3,2,figsize=(11,9),layout='constrained')
for vi,v in enumerate(['WT','AA']):
 d=pd.read_csv(O/f'{v}_timeseries.csv').merge(pd.read_csv(O/f'{v}_Sa_timeseries.csv'),on=['frame_index','production_ns'],validate='one_to_one')
 for c in ['Rg_CA_A','R2_main_DSSP_H','R2_minor_DSSP_H','R3_DSSP_H','Sa_full','Sa_AR391_446']:
  x=d[c].to_numpy();blocks=x.reshape(10,2000).mean(1);se=blocks.std(ddof=1)/np.sqrt(10);xc=x-x.mean();ac=fftconvolve(xc,xc[::-1],mode='full')[len(x)-1:]/np.arange(len(x),0,-1)
  if ac[0]>0:
   ac/=ac[0];cut=next((i for i in range(1,len(ac)) if ac[i]<=0),len(ac));g=max(1,1+2*ac[1:cut].sum())
  else:g=1
  rows.append({'variant':v,'observable':c,'mean':x.mean(),'SD':x.std(ddof=1),'first_1us_mean':x[:10000].mean(),'last_1us_mean':x[10000:].mean(),'block_200ns_SE':se,'approx_autocorrelation_ESS':len(x)/g,'interpretation':'within-trajectory diagnostics; no independent replicate or equilibrium convergence claim'})
 for i,c in enumerate(['Rg_CA_A','R2_main_DSSP_H','Sa_AR391_446']):
  axes[i,vi].plot(d.production_ns,d[c].rolling(100,min_periods=1).mean(),lw=.6);axes[i,vi].plot(np.arange(100,2000,200),d[c].to_numpy().reshape(10,2000).mean(1),'o-',color='black',ms=3)
  axes[i,vi].set(xlabel='Production time (ns)',ylabel=c,title=v)
pd.DataFrame(rows).to_csv(O/'convergence_summary.csv',index=False);fig.savefig(O/'convergence.png',dpi=180);fig.savefig(O/'convergence.pdf');plt.close(fig)
print(pd.DataFrame(rows)[['variant','observable','first_1us_mean','last_1us_mean','approx_autocorrelation_ESS']].to_string(index=False))
