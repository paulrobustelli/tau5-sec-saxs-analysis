from select_windows import *
from matplotlib.ticker import ScalarFormatter
from matplotlib.lines import Line2D
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.titlesize':12,'axes.labelsize':11,'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.8,'xtick.direction':'out','ytick.direction':'out','legend.frameon':False,'svg.fonttype':'none','savefig.facecolor':'white'})
COLORS={'WT_092025':'#243B82','AA_052026':'#D95F02'}
LABELS={'WT_092025':'WT','AA_052026':'W397A/W433A'}
chosen=json.loads((PLOT/'selected_windows.json').read_text());items=[]
for s in chosen:
 name=s['sample'];q,I,E=average(name,s['start'],s['end']);g=guinier(q,I,E);bp=json.loads((PLOT/f'{name}_bift_primary.json').read_text());bd=np.load(PLOT/f'{name}_bift_primary.npz');items.append((s,q,I,E,g,bp,bd))

def save(fig,name):
 fig.savefig(PLOT/(name+'.png'),dpi=300,bbox_inches='tight');fig.savefig(PLOT/(name+'.svg'),bbox_inches='tight');plt.close(fig)

def iq(ax,normalized=True):
 for s,q,I,E,g,bp,bd in items:
  c=COLORS[s['sample']];m=(q>=.008)&(q<=.30);norm=g['I0'] if normalized else 1
  ax.errorbar(q[m],I[m]/norm,yerr=E[m]/norm,fmt='.',ms=3,elinewidth=.5,alpha=.65,color=c)
  ax.plot(bd['q'],bd['fit']/norm,color=c,lw=1.5,label=LABELS[s['sample']])
 ax.set(xscale='log',yscale='log',xlabel='q (Å⁻¹)',ylabel='I(q) / I(0)' if normalized else 'I(q) (arbitrary units)',title='Scattering profiles' if normalized else 'Scattering profiles, measured scale')
 ax.legend();ax.text(.03,.05,'Points: averaged data   Lines: BIFT fit',transform=ax.transAxes,fontsize=8,color='#555555')

def pr(ax):
 for s,q,I,E,g,bp,bd in items:
  r=bd['r'];p=bd['p']/bp['i0'];pe=bd['p_error']/bp['i0'];c=COLORS[s['sample']]
  ax.fill_between(r,np.maximum(0,p-pe),p+pe,color=c,alpha=.13,lw=0)
  ax.plot(r,p,color=c,lw=1.8,label=f"{LABELS[s['sample']]}: Rg = {bp['rg']:.1f} Å, Dmax ≈ {bp['dmax']:.0f} Å")
 ax.set(xlabel='Pair distance r (Å)',ylabel='P(r) / I(0) (Å⁻¹)',title='Pair-distance distributions',ylim=(0,None));ax.legend(fontsize=8.5)

def guinier_panel(fig,slot,item):
 s,q,I,E,g,bp,bd=item;c=COLORS[s['sample']];sub=slot.subgridspec(2,1,height_ratios=[3.2,1],hspace=.06);ax=fig.add_subplot(sub[0]);res=fig.add_subplot(sub[1],sharex=ax)
 valid=(q>=.006)&(q<=.060)&(I>0);fit=(q>=g['qmin'])&(q<=g['qmax'])&(I>0)
 ax.errorbar(q[valid]**2,np.log(I[valid]),yerr=E[valid]/I[valid],fmt='.',ms=3.5,color='#C1C1C1',elinewidth=.6,alpha=.75)
 ax.errorbar(q[fit]**2,np.log(I[fit]),yerr=E[fit]/I[fit],fmt='o',ms=3,color=c,elinewidth=.7)
 xx=q[fit]**2;pred=np.log(g['I0'])-xx*g['Rg']**2/3;ax.plot(xx,pred,color='black',lw=1.2)
 ax.set(ylabel='ln I(q)',title=f"{LABELS[s['sample']]} Guinier fit · frames {s['start']}–{s['end']}")
 ax.text(.98,.97,f"Rg = {g['Rg']:.1f} ± {g['fit_se']:.1f} Å\nqRg ≤ {g['qRgmax']:.2f}",ha='right',va='top',transform=ax.transAxes,fontsize=10)
 z=(np.log(I[fit])-pred)/(E[fit]/I[fit]);res.axhline(0,color='#777777',lw=.7);res.axhline(2,color='#C1C1C1',lw=.5,ls=':');res.axhline(-2,color='#C1C1C1',lw=.5,ls=':');res.plot(xx,z,'o',ms=3,color=c)
 res.set(ylabel='Δln I / σ',xlabel='q² (Å⁻²)',ylim=(-3,3));res.set_yticks([-2,0,2]);ax.tick_params(labelbottom=False)
 fmt=ScalarFormatter(useMathText=True);fmt.set_powerlimits((-3,-3));res.xaxis.set_major_formatter(fmt);ax.set_xlim(0,.0037)
 return ax,res

# Main panel deliberately contains the requested three analysis types.
fig=plt.figure(figsize=(12,9.8));gs=fig.add_gridspec(2,2,height_ratios=[1,1.15],hspace=.34,wspace=.28)
a=fig.add_subplot(gs[0,0]);iq(a);b=fig.add_subplot(gs[0,1]);pr(b);c,_=guinier_panel(fig,gs[1,0],items[0]);d,_=guinier_panel(fig,gs[1,1],items[1])
for ax,l in zip([a,b,c,d],'ABCD'):ax.text(-.14,1.08,l,transform=ax.transAxes,weight='bold',fontsize=14)
fig.suptitle('Tau-5* SEC-SAXS · selected peak windows',fontsize=16,y=.98)
save(fig,'poster_saxs_panel')
for name,fn in [('Iq_normalized',lambda ax:iq(ax,True)),('Iq_measured_scale',lambda ax:iq(ax,False)),('Pr_comparison',pr)]:
 fig,ax=plt.subplots(figsize=(6.2,4.8));fn(ax);fig.tight_layout();save(fig,name)
for item in items:
 fig=plt.figure(figsize=(6.2,5.6));gs=fig.add_gridspec(1,1);guinier_panel(fig,gs[0],item);fig.subplots_adjust(left=.16,right=.97,bottom=.13,top=.90);save(fig,item[0]['sample']+'_Guinier')
# SEC window documentation with recomputed five-frame apparent Rg.
fig,axs=plt.subplots(1,2,figsize=(12,4.4))
for ax,item in zip(axs,items):
 s,q,I,E,g,bp,bd=item;name=s['sample'];q,Y,EE=load(name);y,e=correct(Y,EE,name);m=(q>=.012)&(q<=.25);trace=np.trapezoid(y[m],q[m],axis=0);c=COLORS[name]
 ax.plot(np.arange(450),trace,color=c,lw=1.5);ax.axvspan(s['start']-.5,s['end']+.5,color=c,alpha=.13);ax.set(xlim=(125,215) if name.startswith('WT') else (145,220),xlabel='Frame (zero indexed)',ylabel='Integrated I, q = 0.012–0.25 Å⁻¹',title=f"{LABELS[name]} · selected {s['start']}–{s['end']}")
 twin=ax.twinx();frames=[];rgs=[];ers=[]
 for frame in range(145,211):
  if trace[frame]<.12*max(trace[145:211]):continue
  qq,ii,ee=average(name,frame-2,frame+2);gg=guinier(qq,ii,ee)
  if gg['Rg'] is not None:frames.append(frame);rgs.append(gg['Rg']);ers.append(gg['fit_se'])
 twin.errorbar(frames,rgs,yerr=ers,fmt='o',ms=2.7,color='#8B2525',elinewidth=.6,alpha=.75);twin.set(ylabel='Apparent Rg (Å), 5-frame averages',ylim=(15,55));twin.spines['right'].set_visible(True);twin.tick_params(axis='y',labelcolor='#8B2525')
fig.tight_layout();save(fig,'SEC_selected_windows')
# Diagnostics supporting final figure: BIFT residuals and nearby window Rg.
fig,axs=plt.subplots(2,2,figsize=(12,7.5))
for j,item in enumerate(items):
 s,q,I,E,g,bp,bd=item;c=COLORS[s['sample']];z=(bd['I']-bd['fit'])/bd['error'];axs[0,j].plot(bd['q'],z,'.-',ms=3,lw=.6,color=c);axs[0,j].axhline(0,c='grey',lw=.6);axs[0,j].set(xlabel='q (Å⁻¹)',ylabel='(I − I fit) / σ',title=LABELS[s['sample']]+' BIFT residuals')
 candidates_list=list(csv.DictReader(open(PLOT/'window_candidates.csv')));near=[r for r in candidates_list if r['sample']==s['sample'] and r['eligible']=='True' and abs(int(r['start'])-s['start'])<=2 and abs(int(r['end'])-s['end'])<=2]
 near.sort(key=lambda r:(int(r['start']),int(r['end'])));xs=np.arange(len(near));axs[1,j].errorbar(xs,[float(r['Rg']) for r in near],yerr=[float(r['Rg_fit_se']) for r in near],fmt='o',ms=3,color=c);axs[1,j].axhline(g['Rg'],c='grey',ls='--',lw=.8);axs[1,j].set(xlabel='Nearby eligible window',ylabel='Guinier Rg (Å)',title='Bounds shifted by up to 2 frames');axs[1,j].set_xticks(xs[::max(1,len(xs)//5)],[f"{near[k]['start']}–{near[k]['end']}" for k in xs[::max(1,len(xs)//5)]],rotation=35,ha='right')
fig.tight_layout();save(fig,'fit_and_window_diagnostics')
# Optional companion from the poster; not used to select windows.
fig,ax=plt.subplots(figsize=(6.2,4.8))
for s,q,I,E,g,bp,bd in items:
 m=(q>=.012)&(q<=.22);x=q[m]*g['Rg'];y=x*x*I[m]/g['I0'];ax.plot(x,y,'.-',color=COLORS[s['sample']],ms=3,lw=.7,label=LABELS[s['sample']])
ax.set(xlabel='q Rg',ylabel='(q Rg)² I(q) / I(0)',title='Dimensionless Kratky comparison',ylim=(0,4));ax.legend();fig.tight_layout();save(fig,'Kratky_companion')
(PLOT/'plot_summary.json').write_text(json.dumps([dict(sample=s['sample'],frames=[s['start'],s['end']],guinier=g,bift=bp) for s,q,I,E,g,bp,bd in items],indent=2))
print('Saved poster panel, individual PNG/SVG figures, window and fit diagnostics.',flush=True)
