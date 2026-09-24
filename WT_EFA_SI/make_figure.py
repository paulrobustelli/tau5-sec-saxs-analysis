from pathlib import Path
import os,sys,json
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/saxs_mpl')
import numpy as np
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
REPO=Path(os.environ.get('SAXS_REPO','/Users/f0044gk/Desktop/SAXS/repository'));sys.path.insert(0,str(REPO/'analysis'))
from compare_all_iq import standard_scan,scan,conventional,fit,shape
OUT=Path(__file__).resolve().parent
a=np.load(OUT/'WT_components.npz');q=a['q'];results={};colors=['#c45100','#000080'];names=['Oligomer candidate','Monomer / main peak']
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none'})
fig=plt.figure(figsize=(12,12));gs=fig.add_gridspec(3,2,height_ratios=[.95,1.15,1.15],hspace=.40,wspace=.32)
ax=fig.add_subplot(gs[0,0]);iq=fig.add_subplot(gs[0,1])
for j in [1,0]:
 c=colors[j];ax.plot(a['frames'],a['C'][:,j],color=c,lw=2,label=names[j]);lo,hi=np.percentile(a['bootstrap_C'][:,:,j],[2.5,97.5],axis=0);ax.fill_between(a['frames'],lo,hi,color=c,alpha=.16)
 y=a['S'][:,j];e=a['error'][:,j];m=(q>=.012)&(q<=.25);iq.errorbar(q[m],y[m],e[m],fmt='o',ms=2.5,lw=.55,color=c,alpha=.8,label=names[j])
ax.set(title='A  WT EFA elution profiles',xlabel='SEC frame (zero-based)',ylabel='Component amplitude / own maximum',xlim=(120,205));ax.legend(frameon=False,fontsize=9)
iq.set(title='B  Extracted component scattering',xlabel='q (Å⁻¹)',ylabel='I(q), own-peak scale (a.u.)',yscale='symlog',ylim=(-.025,3),xlim=(.005,.255));iq.set_yscale('symlog',linthresh=.01);iq.axhline(0,color='.7',lw=.5);iq.legend(frameon=False,fontsize=9)
for row,j in enumerate([1,0],1):
 y=a['S'][:,j];e=a['error'][:,j];lo=.012 if j==1 else .015
 g=standard_scan(q,y,e,lo);x=scan(q,y,e,lo)[0]
 results[names[j]]={'conventional':g,'extended_single_chain_diagnostic' if j==0 else 'extended':x,'qmin_sensitivity':{str(l):{'guinier':standard_scan(q,y,e,l),'extended':scan(q,y,e,l)[0]} for l in [.010,.012,.015,.020]}}
 for col,(typ,f) in enumerate([('Guinier',g),('Extended Guinier',x)]):
  sub=gs[row,col].subgridspec(2,1,height_ratios=[3,1],hspace=.06);top=fig.add_subplot(sub[0]);res=fig.add_subplot(sub[1],sharex=top)
  inside=(q>=f['qmin']-1e-9)&(q<=f['qmax']+1e-9);ctx=(q>=.006)&(q<=f['qmax']*1.17)&(y>0);outside=ctx&~inside;c=colors[j]
  top.errorbar(q[outside]**2,np.log(y[outside]),e[outside]/y[outside],fmt='o',ms=3,color=c,alpha=.18,lw=.6)
  top.errorbar(q[inside]**2,np.log(y[inside]),e[inside]/y[inside],fmt='o',ms=3,color=c,lw=.65,capsize=1.5)
  pred=lambda z:f['I0']*(np.exp(-z*z*f['Rg']**2/3) if col==0 else shape(z,f['nu'],120))
  z=np.linspace(q[inside][0],q[inside][-1],160);top.plot(z*z,np.log(pred(z)),color=c,lw=2)
  vals=[]
  for b in a['bootstrap_S'][:,:,j]:
   ff=conventional(q,b,e,f['qmin'],f['qmax']) if col==0 else fit(q[inside],b[inside],e[inside],120)
   if ff and (col==0 or not ff['boundary']):vals.append(ff['Rg'])
  ci=np.percentile(vals,[2.5,97.5]);f['bootstrap95']=ci.tolist();f['bootstrap_valid']=len(vals)
  top.text(.98,.97,f"Rg = {f['Rg']:.1f} Å\n95% interval: {ci[0]:.1f}–{ci[1]:.1f} Å\nqRg max = {f['qRgmax']:.2f}",ha='right',va='top',transform=top.transAxes,fontsize=9,bbox=dict(facecolor='white',edgecolor='none',alpha=.8))
  label='CDEF'[(row-1)*2+col];title=f'{label}  {names[j]} — {typ}'
  if j==0 and col==1:title+='*'
  top.set(title=title,ylabel='ln I(q)');top.tick_params(labelbottom=False)
  residual=(y[inside]-pred(q[inside]))/e[inside];res.errorbar(q[inside]**2,residual,yerr=1,fmt='o',ms=2.5,lw=.5,color=c,alpha=.7);res.axhline(0,color='.25',lw=.8);res.set(xlabel='q² (Å⁻²)',ylabel='ΔI / σ');res.xaxis.set_major_locator(MaxNLocator(4))
fig.suptitle('WT SEC–SAXS: EFA component separation',fontsize=16,y=.985)
fig.text(.08,.025,'Shading: conditional 95% bootstrap intervals; error bars: ±1σ. Faint points are excluded from fitting.\n*Single-chain extended model applied diagnostically to the oligomer candidate; stoichiometry is unassigned.',fontsize=9)
fig.subplots_adjust(top=.935,bottom=.10);fig.savefig(OUT/'WT_EFA_SI.png',dpi=240);fig.savefig(OUT/'WT_EFA_SI.pdf');fig.savefig(OUT/'WT_EFA_SI.svg');plt.close(fig)
(OUT/'fit_results.json').write_text(json.dumps(results,indent=2))
print(json.dumps({k:{kk:vv for kk,vv in v.items() if kk!='qmin_sensitivity'} for k,v in results.items()},indent=2))
