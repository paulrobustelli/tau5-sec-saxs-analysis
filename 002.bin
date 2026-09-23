"""Natalie-style display of corrected data, explicit errors and transparent excluded q points."""
from paper_figures import *
from IPython.display import display,Markdown
STYLE=OUT/'natalie_style';STYLE.mkdir(exist_ok=True)
COL={'WT_092025':'navy','AA_052026':'orangered'}
def output(fig,name):
 fig.tight_layout();fig.savefig(STYLE/(name+'.png'),dpi=180);fig.savefig(STYLE/(name+'.pdf'));plt.close(fig)
def calculate():
 rows=json.loads((OUT/'curve_results.json').read_text());summary=[];high=[]
 for sample in SPECS:
  r=next(r for r in rows if r['sample']==sample and r['kind']=='chosen');q,y,e=np.loadtxt(OUT/(r['key']+'.dat')).T
  for lo,hi in [(.15,.25),(.25,.35),(.35,.5),(.5,.6)]:
   m=(q>=lo)&(q<=hi);high.append(dict(sample=sample,q=[lo,hi],median_I_over_sigma=float(np.median(y[m]/e[m])),n=int(m.sum()),nonpositive=int(np.sum(y[m]<=0))))
  result=dict(sample=sample,conventional=r['conventional'],extended=r['extended'],extended95=r['extended95'],bift=[])
  for hi in [.25,.35,.5,.6]:
   m=(q>=.012-1e-9)&(q<=hi+1e-9);np.random.seed(20260923)
   o=BIFT.doBift(q[m],y[m],e[m],sample,40,150,1e10,12,50,220,15,80,single_proc=True,nprocs=1)
   if o is not None:
    p=o.getAllParameters();result['bift'].append(dict(qmax=hi,**{k:float(p[k]) for k in ['rg','rger','dmax','dmaxer','i0','i0er','chisq']}))
    np.savez_compressed(STYLE/f'{sample}_Pr_{hi}.npz',r=o.r,p=o.p,error=o.err,q=o.q_orig,I=o.i_orig,sigma=o.err_orig,fit=o.i_fit)
  summary.append(result)
 (STYLE/'uncertainty_and_qmax.json').write_text(json.dumps(summary,indent=2));(STYLE/'high_q_signal.json').write_text(json.dumps(high,indent=2))
 return summary

def make():
 rows=json.loads((OUT/'curve_results.json').read_text());summary=json.loads((STYLE/'uncertainty_and_qmax.json').read_text())
 plt.rcParams.update({'axes.spines.top':True,'axes.spines.right':True,'font.size':11,'pdf.fonttype':42})
 fig,ax=plt.subplots(figsize=(6,4));figs,axs=plt.subplots(figsize=(6,4));fk,ak=plt.subplots(figsize=(6,6));fp,ap=plt.subplots(figsize=(6,4));fh,ah=plt.subplots(2,1,figsize=(7,7))
 for sample in SPECS:
  c=COL[sample];name=sample[:2];r=next(r for r in rows if r['sample']==sample and r['kind']=='chosen');u=next(z for z in summary if z['sample']==sample);q,y,e=np.loadtxt(OUT/(r['key']+'.dat')).T;m=(q>=.005)&(q<=.6);good=m&(y>0);f=r['extended'];g=r['conventional'];bands=np.load(OUT/f'{sample}_normalized_bands.npz')
  ax.plot(q[good],y[good],color=c,lw=1,label=name);ax.errorbar(q[good],y[good],e[good],fmt='none',ecolor=c,elinewidth=.7,capsize=2,alpha=.7);ax.fill_between(q[m],np.where(y[m]-e[m]>0,y[m]-e[m],np.nan),y[m]+e[m],color=c,alpha=.25)
  axs.errorbar(q[good],y[good]/f['I0'],e[good]/f['I0'],fmt='.',ms=3,capsize=2,color=c,label=name);axs.fill_between(q[m],*bands['normalized95'][:,m],color=c,alpha=.2)
  x=q[m]*f['Rg'];ky=x*x*y[m]/f['I0'];ak.plot(x,ky,color=c,label=name);ak.errorbar(x,ky,x*x*e[m]/f['I0'],fmt='none',ecolor=c,capsize=2,elinewidth=.6,alpha=.65)
  # Extend dimensionless bands to full displayed x=8 from paired original perturbations.
  boot=np.load(OUT/f'{sample}_bootstrap.npz')['direct'][:,:,0];xx=np.linspace(.4,8,160);draws=[]
  fm=(q>=f['qmin']-1e-9)&(q<=f['qmax']+1e-9)
  for b in boot:
   ff=fit(q[fm],b[fm],e[fm])
   if ff and not ff['boundary']:draws.append(np.interp(xx,q*ff['Rg'],(q*ff['Rg'])**2*b/ff['I0']))
  ak.fill_between(xx,*np.percentile(draws,[2.5,97.5],axis=0),color=c,alpha=.2)
  b=u['bift'][0];p=np.load(STYLE/f'{sample}_Pr_0.25.npz');ap.plot(p['r'],p['p']/b['i0'],color=c,label=f"{name}: Rg={b['rg']:.2f} ± {b['rger']:.2f} Å");ap.errorbar(p['r'],p['p']/b['i0'],p['error']/b['i0'],fmt='none',capsize=2,elinewidth=.7,color=c);ap.fill_between(p['r'],(p['p']-p['error'])/b['i0'],(p['p']+p['error'])/b['i0'],color=c,alpha=.25)
  ah[0].errorbar(q[m],y[m],e[m],fmt='.',ms=3,capsize=2,color=c,label=name);ah[1].errorbar(q[m],q[m]**2*y[m]/f['I0'],q[m]**2*e[m]/f['I0'],fmt='.',ms=3,capsize=2,color=c)
  # Separate conventional and extended panels: all context shown, only fit range opaque.
  gf,ga=plt.subplots(2,2,figsize=(12,6),gridspec_kw={'height_ratios':[3,1]},sharex='col')
  for j,(key,ff) in enumerate([('Conventional Guinier',g),('Extended Guinier',f)]):
   ctx=(q<=.085)&(y>0);inside=ctx&(q>=ff['qmin']-1e-9)&(q<=ff['qmax']+1e-9);outside=ctx&~inside
   ga[0,j].errorbar(q[outside]**2,np.log(y[outside]),e[outside]/y[outside],fmt='o',ms=4,color=c,alpha=.15,capsize=2,label='Outside fit range')
   ga[0,j].errorbar(q[inside]**2,np.log(y[inside]),e[inside]/y[inside],fmt='o',ms=4,color=c,capsize=2,label='Fitted points ±1σ')
   pred=ff['I0']*(np.exp(-q[inside]**2*ff['Rg']**2/3) if j==0 else shape(q[inside],ff['nu']))
   ga[0,j].plot(q[inside]**2,np.log(pred),color=c,lw=2,label=f"Rg={ff['Rg']:.2f} ± {ff['Rg_conditional_se']:.2f} Å (fit SE)")
   ga[1,j].errorbar(q[inside]**2,(np.log(y[inside])-np.log(pred))/(e[inside]/y[inside]),yerr=np.ones(inside.sum()),fmt='o',ms=3,capsize=2,color=c)
   ga[1,j].axhline(0,color='red',lw=1);ga[0,j].set(title=key,ylabel='ln(I(q))');ga[1,j].set(xlabel='q² (Å⁻²)',ylabel='Residual / σ');ga[0,j].legend(fontsize=8)
  gf.suptitle(f'{name} average: frames {r["frames"][0]}–{r["frames"][1]}');output(gf,name+'_Guinier')
 for a,title,yl in [(ax,'Average scattering — measured range','Average intensity (a.u.)'),(axs,'WT–AA scattering shape','I(q) / I(0)')]:a.set(xscale='log',yscale='log',xlim=(.005,.6),xlabel='q (Å⁻¹)',ylabel=yl,title=title);a.legend()
 ak.axvline(np.sqrt(3),color='grey',ls='--');ak.axhline(3/np.e,color='grey',ls='--');ak.axhline(2,color='teal');ak.grid(True,alpha=.4);ak.set(xlim=(0,8),ylim=(0,3.2),xlabel='q Rg',ylabel='(q Rg)² I(q) / I(0)',title='Dimensionless Kratky');ak.legend()
 ap.set(xlabel='r (Å)',ylabel='P(r) / I(0)',title='Pair-distance distribution — BIFT');ap.legend(fontsize=9)
 ah[0].set(yscale='symlog',xlabel='q (Å⁻¹)',ylabel='I(q)',title='Signed I(q), ±1σ');ah[0].legend();ah[1].set(xlabel='q (Å⁻¹)',ylabel='q² I(q) / I(0)',title='High-q Kratky — ±1σ at fixed normalization')
 output(fig,'Iq');output(figs,'Iq_normalized');output(fk,'Kratky');output(fp,'Pr');output(fh,'high_q')
 fc,ac=plt.subplots(1,2,figsize=(11,4))
 for j,s in enumerate(summary):
  for b in s['bift']:ac[j].errorbar(b['qmax'],b['rg'],yerr=b['rger'],fmt='o',capsize=4,color=COL[s['sample']])
  ac[j].set(xlabel='P(r) fit qmax (Å⁻¹)',ylabel='P(r) Rg ± local BIFT σ (Å)',title=s['sample'][:2]+' q-range sensitivity')
 output(fc,'Pr_qmax')
 return summary

def table():
 s=json.loads((STYLE/'uncertainty_and_qmax.json').read_text());h='|Sample|Guinier Rg ± fit SE|Extended Rg ± fit SE|Extended bootstrap95%|P(r) Rg ± local BIFT σ|Dmax ± local BIFT σ|\n|---|---:|---:|---|---:|---:|\n'
 for r in s:
  g=r['conventional'];f=r['extended'];b=r['bift'][0];h+=f"|{r['sample'][:2]}|{g['Rg']:.2f} ± {g['Rg_conditional_se']:.2f}|{f['Rg']:.2f} ± {f['Rg_conditional_se']:.2f}|{np.round(r['extended95'],2).tolist()}|{b['rg']:.2f} ± {b['rger']:.2f}|{b['dmax']:.1f} ± {b['dmaxer']:.1f}|\n"
 display(Markdown(h+'\nAll distances in Å. Local/conditional errors exclude baseline-model uncertainty.'))
if __name__=='__main__':calculate();make()

def efa_bands():
 choices=json.loads((OUT/'choices.json').read_text());sample='WT_092025';q,Y,E=load(sample);spec=SPECS[sample];frames=np.arange(spec['roi'][0],spec['roi'][1]+1);A=matrix(Y.shape[1],frames,choices[sample]['primary_model']);err=np.sqrt(E**2@(A*A));rng=np.random.default_rng(20260923);draws=[]
 for _ in range(100):
  Yn=Y+rng.normal(size=Y.shape)*E;_,o=solve(q,Yn@A,err,frames,*spec['fid'])
  if o is not None:draws.append(o['C']/o['C'].max(0))
 np.savez_compressed(STYLE/'WT_elution95.npz',frames=frames,interval=np.percentile(draws,[2.5,97.5],axis=0),success=len(draws))

def efa_plot():
 if not (STYLE/'WT_elution95.npz').exists():efa_bands()
 b=np.load(STYLE/'WT_elution95.npz');d=np.load(OUT/'WT_092025_post17_efa.npz');fig,ax=plt.subplots(1,2,figsize=(12,4))
 for j,(kind,label,color) in enumerate([('EFA_shoulder_or_other','Early shoulder','grey'),('EFA_main','Main component','navy')]):
  q,y,e=np.loadtxt(OUT/f'WT_092025_{kind}.dat').T;m=(q>=.006)&(q<=.25);ax[0].errorbar(q[m],y[m],e[m],fmt='.',ms=3,capsize=2,color=color,label=label);ax[0].fill_between(q[m],y[m]-e[m],y[m]+e[m],color=color,alpha=.2)
  ax[1].plot(d['frames'],d['C'][:,j]/d['C'][:,j].max(),color=color,label=label);ax[1].fill_between(b['frames'],b['interval'][0,:,j],b['interval'][1,:,j],color=color,alpha=.25)
 ax[0].set_yscale('symlog',linthresh=.003);ax[0].set(xlabel='q (Å⁻¹)',ylabel='Component contribution I(q)',title='WT EFA curves ±1σ');ax[1].set(xlabel='Frame',ylabel='Elution / own maximum',title='WT EFA — conditional95% bands')
 for a in ax:a.legend()
 output(fig,'WT_EFA')

def reference_guinier():
 """Match supplied Natalie screenshots; faint context is not used in regression."""
 rows=json.loads((OUT/'curve_results.json').read_text())
 for sample in SPECS:
  r=next(z for z in rows if z['sample']==sample and z['kind']=='chosen');q,y,e=np.loadtxt(OUT/(r['key']+'.dat')).T;c='navy' if sample.startswith('WT') else 'orangered';name='WT' if sample.startswith('WT') else 'W397A/W433A'
  for key,title in [('conventional','Guinier'),('extended','Extended Guinier')]:
   f=r[key];fig,ax=plt.subplots(2,1,figsize=(6,6),gridspec_kw={'height_ratios':[3,1]});fig.subplots_adjust(hspace=.25)
   ctx=(q<=.06)&(y>0);inside=ctx&(q>=f['qmin']-1e-9)&(q<=f['qmax']+1e-9);outside=ctx&~inside;inds=np.flatnonzero(inside)
   ax[0].errorbar(q[outside]**2,np.log(y[outside]),e[outside]/y[outside],fmt='o',ms=6,color='lightgray',alpha=.65,elinewidth=.5,label='Data outside fit')
   ax[0].errorbar(q[inside]**2,np.log(y[inside]),e[inside]/y[inside],fmt='o',ms=5,color=c,alpha=.85,elinewidth=.7,capsize=2,label=f"Guinier points\nidxmin={inds[0]}, idxmax={inds[-1]}\nI(0)={f['I0']:.2f}")
   pred=lambda z: f['I0']*(np.exp(-z*z*f['Rg']**2/3) if key=='conventional' else shape(z,f['nu']))
   ax[0].plot(q[inside]**2,np.log(pred(q[inside])),color=c,lw=2,label=f"Rg={f['Rg']:.2f} ± {f['Rg_conditional_se']:.2f} Å")
   # Same visual continuation as reference, dashed to identify extrapolation.
   for mm in [ctx&(q<f['qmin']-1e-9),ctx&(q>f['qmax']+1e-9)]:
    ax[0].plot(q[mm]**2,np.log(pred(q[mm])),'--',color=c,lw=1,alpha=.45)
   resid=(np.log(y[inside])-np.log(pred(q[inside])))/(e[inside]/y[inside]);ax[1].errorbar(q[inside]**2,resid,yerr=np.ones(len(resid)),fmt='o',ms=5,color=c,alpha=.55,elinewidth=.6,capsize=2);ax[1].axhline(0,color='red',lw=1)
   ax[0].set(title=f'{title} Fit for {name} Average Profile',ylabel='ln(I(q))',xlabel='q² (Å⁻²)');ax[0].legend(fontsize=9);ax[1].set(xlabel='q² (Å⁻²)',ylabel='Δln(I) / σln(I)');output(fig,sample[:2]+'_'+key+'_reference')
