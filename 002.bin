"""Four-panel WT/AA summary from existing fits; no new fitting or window selection."""
from natalie_style_figures import *
def composite():
 rows=json.loads((OUT/'curve_results.json').read_text());summary=json.loads((STYLE/'uncertainty_and_qmax.json').read_text())
 plt.rcParams.update({'font.size':10,'axes.spines.top':True,'axes.spines.right':True,'pdf.fonttype':42})
 fig=plt.figure(figsize=(11,9),layout='constrained');gs=fig.add_gridspec(2,2,wspace=.08,hspace=.10)
 a=fig.add_subplot(gs[0,0]);sub=gs[0,1].subgridspec(2,1,height_ratios=[3,1],hspace=.04);gax=fig.add_subplot(sub[0]);rax=fig.add_subplot(sub[1],sharex=gax);k=fig.add_subplot(gs[1,0]);pax=fig.add_subplot(gs[1,1])
 for sample in SPECS:
  color=COL[sample];name=sample[:2];r=next(z for z in rows if z['sample']==sample and z['kind']=='chosen');u=next(z for z in summary if z['sample']==sample);q,y,e=np.loadtxt(OUT/(r['key']+'.dat')).T;f=r['extended'];g=r['conventional'];m=(q<=.6)&(y>0);b=np.load(OUT/f'{sample}_normalized_bands.npz')
  a.errorbar(q[m],y[m]/f['I0'],e[m]/f['I0'],fmt='.',ms=3,elinewidth=.5,capsize=1,color=color,label=name);a.fill_between(q[m],*b['normalized95'][:,m],color=color,alpha=.17)
  ctx=(q<=.06)&(y>0);fitmask=ctx&(q>=f['qmin']-1e-9)&(q<=f['qmax']+1e-9);outside=ctx&~fitmask
  gax.errorbar(q[outside]**2,np.log(y[outside]/f['I0']),e[outside]/y[outside],fmt='o',ms=3,color='lightgray',alpha=.5,elinewidth=.5)
  gax.errorbar(q[fitmask]**2,np.log(y[fitmask]/f['I0']),e[fitmask]/y[fitmask],fmt='o',ms=3,color=color,elinewidth=.5,capsize=1)
  for ff,ls,marker,label in [(g,'--','x','G'),(f,'-','o','EG')]:
   mm=(q>=ff['qmin']-1e-9)&(q<=ff['qmax']+1e-9);pred=ff['I0']*(np.exp(-q[mm]**2*ff['Rg']**2/3) if label=='G' else shape(q[mm],ff['nu']))
   gax.plot(q[mm]**2,np.log(pred/f['I0']),ls,color=color,lw=1.5,label=f"{name} {label}: {ff['Rg']:.2f} ± {ff['Rg_conditional_se']:.2f} Å")
   rax.plot(q[mm]**2,(np.log(y[mm])-np.log(pred))/(e[mm]/y[mm]),marker,ms=3,color=color,alpha=.7)
  x=q*f['Rg'];mm=(x<=8)&(q>=.012);k.errorbar(x[mm],x[mm]**2*y[mm]/f['I0'],x[mm]**2*e[mm]/f['I0'],fmt='.-',ms=2,lw=.8,elinewidth=.5,capsize=1,color=color,label=name);k.fill_between(b['x'],*b['dimensionless95'],color=color,alpha=.2)
  pp=np.load(STYLE/f'{sample}_Pr_0.25.npz');bf=u['bift'][0];pax.plot(pp['r'],pp['p']/bf['i0'],color=color,lw=1.5,label=f"{name}: {bf['rg']:.2f} ± {bf['rger']:.2f} Å");pax.errorbar(pp['r'],pp['p']/bf['i0'],pp['error']/bf['i0'],fmt='none',elinewidth=.6,capsize=2,color=color);pax.fill_between(pp['r'],(pp['p']-pp['error'])/bf['i0'],(pp['p']+pp['error'])/bf['i0'],color=color,alpha=.2)
 a.set(xscale='log',yscale='log',xlim=(.005,.6),xlabel='q (Å⁻¹)',ylabel='I(q) / I(0)',title='A   Scattering');a.legend(frameon=False)
 gax.set(ylabel='ln[I(q) / I(0)]',title='B   Guinier and extended Guinier');gax.tick_params(labelbottom=False);gax.legend(fontsize=7,loc='lower left',framealpha=.9)
 rax.axhline(0,color='red',lw=.7);rax.axhspan(-1,1,color='gray',alpha=.12);rax.set(xlabel='q² (Å⁻²)',ylabel='Residual / σ');rax.ticklabel_format(axis='x',style='sci',scilimits=(-3,-3))
 k.axvline(np.sqrt(3),color='gray',ls='--',lw=.7);k.axhline(3/np.e,color='gray',ls='--',lw=.7);k.axhline(2,color='teal',lw=.7);k.set(xlim=(0,8),ylim=(0,3.2),xlabel='q Rg',ylabel='(q Rg)² I(q) / I(0)',title='C   Dimensionless Kratky');k.legend(frameon=False)
 pax.set(xlabel='r (Å)',ylabel='P(r) / I(0)',title='D   Pair-distance distribution');pax.legend(fontsize=8,frameon=False)
 for ext in ['png','pdf']:fig.savefig(STYLE/('WT_AA_four_panel.'+ext),dpi=240)
 plt.close(fig)
if __name__=='__main__':composite()
