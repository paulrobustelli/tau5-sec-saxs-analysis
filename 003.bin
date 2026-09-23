"""Reproducible paper-figure candidates; choices use stationarity/shape, never target Rg."""
import os
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/saxs_mpl');os.environ.setdefault('NUMBA_CACHE_DIR','/private/tmp/saxs_numba')
from pathlib import Path
import sys,json,itertools
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar
from inspect_data import load,ROOT
from component_recovery import SPECS,solve,transform
from independent_baseline import window_stats
from compare_all_iq import conventional,standard_scan,fit,shape,quick_extended
sys.path.insert(0,str(ROOT/'bioxtasraw-source'))
from bioxtasraw import BIFT
OUT=ROOT/'paper_figures';OUT.mkdir(exist_ok=True)
MODELS={'early25':(25,35),'early60':(60,70),'early47':(24,70),'post11':(210,220),'post17':(220,236)}
COLORS={'WT_092025':'#2563a6','AA_052026':'#d55e00'}

def save(obj,name):
 (OUT/name).write_text(json.dumps(obj,indent=2,default=lambda x:x.item() if isinstance(x,np.generic) else x.tolist()))
def matrix(n,frames,model):
 A=np.eye(n)[:,frames]
 if model=='linear':
  lo,hi=24,70;f=np.clip((frames-47)/(228-47),0,1)
  A[lo:hi+1]-=(1-f)[None,:]/47;A[220:237]-=f[None,:]/17
 else:
  lo,hi=MODELS[model];A[lo:hi+1]-=1/(hi-lo+1)
 return A

def getfits(q,y,e):
 g=standard_scan(q,y,e);m=(q>=.012-1e-9)&(q<=.055+1e-9);x=fit(q[m],y[m],e[m]);x=x if x and not x['boundary'] and x['qRgmax']<=2 else None
 return {'conventional':g,'extended':x}

def shape_difference(q,Y,E,w1,w2):
 m=(q>=.012)&(q<=.15);a=Y[m]@w1;b=Y[m]@w2;v1=E[m]**2@(w1*w1);v2=E[m]**2@(w2*w2);c=E[m]**2@(w1*w2)
 def fun(s):return np.mean((a-s*b)**2/np.maximum(v1+s*s*v2-2*s*c,1e-30))
 r=minimize_scalar(fun,bounds=(.05,20),method='bounded');return float(r.fun)

def bift(q,y,e,key,dmax=(50,220)):
 m=(q>=.012-1e-9)&(q<=.25+1e-9);np.random.seed(20260923)
 o=BIFT.doBift(q[m],y[m],e[m],key,40,150,1e10,12,*dmax,15,80,single_proc=True,nprocs=1)
 if o is None:return None
 p=o.getAllParameters();np.savez_compressed(OUT/(key+'_Pr.npz'),r=o.r,p=o.p,error=o.err,q=o.q_orig,I=o.i_orig,sigma=o.err_orig,fit=o.i_fit)
 return {k:float(p[k]) for k in ['rg','dmax','chisq']}

def run(nboot=100):
 """Full EFA baseline audit, window screen, paired uncertainty and BIFT."""
 diagnostics=[];windowrows=[];choices={};curves=[];supportgrid=[]
 for si,(sample,spec) in enumerate(SPECS.items()):
  q,Y,E=load(sample);frames=np.arange(*[spec['roi'][0],spec['roi'][1]+1]);m=(q>=.012)&(q<=.25);modeldata={}
  for model in list(MODELS)+['linear']:
   A=matrix(Y.shape[1],frames,model);I=Y@A;err=np.sqrt(E**2@(A*A));info,o=solve(q,I,err,frames,*spec['fid']);checks=[window_stats(q,Y,E,*MODELS[k]) for k in (['early47','post17'] if model=='linear' else [model])]
   D=I[m]/err[m].mean(1)[:,None];U,s,Vt=np.linalg.svd(D,full_matrices=False);rank1=(U[:,:1]*s[:1])@Vt[:1];rank2=(U[:,:2]*s[:2])@Vt[:2]
   info.update(sample=sample,model=model,anchor_pass=all(c['passes'] for c in checks),anchor_checks=checks,s2_over_s1=float(s[1]/s[0]),mode2_q_ac=float(U[:-1,1]@U[1:,1]),mode2_t_ac=float(Vt[1,:-1]@Vt[1,1:]),unconstrained_rank2_ms=float(np.mean(((D-rank2)*err[m].mean(1)[:,None]/err[m])**2)))
   if o is not None:
    scale=o['C'].mean(0);S=o['S']*scale;L=A@o['L']*scale;se=np.sqrt(E**2@(L*L));band=(q>=.012)&(q<=.15)
    info['strong_negative_fraction']=np.mean(S[m]<-3*se[m],axis=0).tolist();info['component_band_snr']=np.sqrt(np.mean((S[band]/se[band])**2,axis=0)).tolist()
    info['relative_other_absolute_contribution']=float(np.sum(np.abs(o['S'][band,0,None]*o['C'][None,:,0]))/np.sum(np.abs(o['S'][band]@o['C'].T)))
   modeldata[model]=(A,I,err,info,o);diagnostics.append(info)
   np.savez_compressed(OUT/f'{sample}_{model}_efa.npz',q=q,frames=frames,I=I,error=err,C=o['C'] if o else np.empty((0,2)),S=o['S'] if o else np.empty((len(q),0)),singular_values=s,U=U,Vt=Vt)
   print('EFA',sample,model,info['converged'],flush=True)
  # Primary buffers follow the documented source-notebook reference choices.
  # Both pass stationarity. This is an explicit reference convention, not proof
  # of optimal subtraction; matched early/post/linear alternatives remain visible.
  primary='post17' if sample.startswith('WT') else 'early47';A,I,err,info,o=modeldata[primary];trace=np.trapezoid(I[m],q[m],axis=0);peak=int(frames[np.argmax(trace)])
  candidates=[]
  for a in range(peak-12,peak+1):
   for b in range(peak,peak+16):
    if not 7<=b-a+1<=21:continue
    half=(a+b)//2;metrics=[];snrs=[]
    for model in ['early25','early60','early47','post17','linear']:
     aa,ii,ee,_,_=modeldata[model];w=aa[:,(frames>=a)&(frames<=b)].mean(1);w1=aa[:,(frames>=a)&(frames<=half)].mean(1);w2=aa[:,(frames>half)&(frames<=b)].mean(1)
     metrics.append(shape_difference(q,Y,E,w1,w2));snrs.append(float(np.sqrt(np.sum((Y[m]@w)**2/(E[m]**2@(w*w))))))
    r=dict(sample=sample,start=a,end=b,n=b-a+1,peak=peak,max_halves_ms=max(metrics),halves_ms=metrics,primary_snr=snrs[3] if sample.startswith('WT') else snrs[2],eligible=max(metrics)<=1.5)
    if sample.startswith('WT'):
     fractions=[]
     for mod in ['early47','post17','linear']:
      oo=modeldata[mod][4];mm=(q>=.012)&(q<=.15);cc=oo['C'][(frames>=a)&(frames<=b)];ss=oo['S'][mm];fractions.append(float(np.sum(abs(ss[:,0,None]*cc[None,:,0]))/np.sum(abs(ss@cc.T))))
     r['max_model_other_fraction']=max(fractions);r['eligible'] &= max(fractions)<=.05
    candidates.append(r)
  valid=[r for r in candidates if r['eligible']]
  if not valid:raise RuntimeError('No shape-consistent window; do not choose by Rg')
  best=max(valid,key=lambda r:r['primary_snr']);choices[sample]=dict(primary_model=primary,buffer=list(MODELS[primary]),window=[best['start'],best['end']],selection=best,alternatives=sorted(valid,key=lambda r:-r['primary_snr'])[:5]);windowrows.extend(candidates)
  print('WINDOW',sample,best,flush=True)
  windows={'chosen':choices[sample]['window'],'leading':([161,165] if sample.startswith('WT') else [165,169]),'apex':[peak-2,peak+2],'trailing':[peak+6,peak+10],'previous':spec['avg']}
  if sample.startswith('WT'):windows['shoulder']=[135,145]
  rng=np.random.default_rng(20260923+si)
  # Perturb every source frame once per draw; all windows/baselines see same draw.
  boot={k:[] for k in ['direct','efa']};W=np.stack([matrix(Y.shape[1],np.arange(a,b+1),primary).mean(1) for a,b in windows.values()],axis=1)
  for rep in range(nboot):
   Yn=Y+rng.normal(size=Y.shape)*E;boot['direct'].append(Yn@W);_,ob=solve(q,Yn@A,err,frames,*spec['fid'])
   boot['efa'].append(ob['S']*ob['C'][(frames>=best['start'])&(frames<=best['end'])].mean(0) if ob else np.full((len(q),2),np.nan))
  bd=np.array(boot['direct']);be=np.array(boot['efa']);good=np.all(np.isfinite(be),axis=(1,2));choices[sample]['efa_bootstrap_success']=int(good.sum())
  np.savez_compressed(OUT/f'{sample}_bootstrap.npz',direct=bd,efa=be,q=q,window_names=np.array(list(windows)),W=W)
  for j,(kind,win) in enumerate(windows.items()):
   y=Y@W[:,j];e=np.sqrt(E**2@(W[:,j]**2));key=f'{sample}_{kind}';np.savetxt(OUT/(key+'.dat'),np.c_[q,y,e],header=f'q_A^-1 I sigma; buffer{MODELS[primary]}; frames{win}')
   f=getfits(q,y,e);r=dict(sample=sample,kind=kind,frames=win,buffer=list(MODELS[primary]),key=key,**f)
   v=[quick_extended(q,d,e,f['extended']) for d in bd[:,:,j]];v=np.array(v,dtype=float);r['extended95']=np.nanpercentile(v,[2.5,97.5]).tolist() if np.isfinite(v).sum()>20 else None
   if kind in ['chosen','leading','apex','trailing','previous']:r['bift']=bift(q,y,e,key)
   curves.append(r)
  if o is not None and good.sum()>=20:
   scale=o['C'][(frames>=best['start'])&(frames<=best['end'])].mean(0);S=o['S']*scale
   for j,kind in enumerate(['EFA_shoulder_or_other','EFA_main']):
    y=S[:,j];e=be[good,:,j].std(0,ddof=1);key=f'{sample}_{kind}';np.savetxt(OUT/(key+'.dat'),np.c_[q,y,e],header='q_A^-1 I sigma; component contribution to chosen window; refitted-bootstrap sigma')
    f=getfits(q,y,e);r=dict(sample=sample,kind=kind,key=key,**f,bootstrap_success=int(good.sum()))
    vals=np.array([quick_extended(q,d,e,f['extended']) for d in be[good,:,j]],dtype=float);r['extended95']=np.nanpercentile(vals,[2.5,97.5]).tolist() if np.isfinite(vals).sum()>20 else None
    if kind=='EFA_main':r['bift']=bift(q,y,e,key)
    # A single-chain extended model cannot establish a multimer identity.
    if kind=='EFA_shoulder_or_other':r['interpretation']='Diagnostic fit only; no oligomer identity or quantitative Rg accepted without a stable Guinier regime.'
    curves.append(r)
  # Baseline spread on exactly the chosen window, independent of frame selection.
  for model,(aa,ii,ee,inf,oo) in modeldata.items():
   w=aa[:,(frames>=best['start'])&(frames<=best['end'])].mean(1);y=Y@w;e=np.sqrt(E**2@(w*w));key=f'{sample}_chosen_{model}';np.savetxt(OUT/(key+'.dat'),np.c_[q,y,e]);curves.append(dict(sample=sample,kind='baseline_sensitivity',model=model,key=key,**getfits(q,y,e)))
  # Support/q-range robustness, fixed grid; all failures recorded.
  for model,end,start,qmax in itertools.product(['early47','post17','linear'],[spec['fid'][0]-5,spec['fid'][0],spec['fid'][0]+5],[spec['fid'][1]-5,spec['fid'][1],spec['fid'][1]+5],[.15,.25]):
   aa,ii,ee,_,_=modeldata[model];inf,oo=solve(q,ii,ee,frames,end,start,qmax=qmax);inf.update(sample=sample,model=model)
   if oo:
    lw=aa@oo['L'];se=np.sqrt(E**2@(lw*lw));inf['other_guinier']=standard_scan(q,oo['S'][:,0],se[:,0]);inf['main_fits']=getfits(q,oo['S'][:,1],se[:,1])
   supportgrid.append(inf)
  save(diagnostics,'efa_diagnostics.json');save(windowrows,'window_screen.json');save(choices,'choices.json');save(curves,'curve_results.json');save(supportgrid,'support_grid.json')
 return choices

def emit(fig,name):
 fig.tight_layout();fig.savefig(OUT/(name+'.png'),dpi=220);fig.savefig(OUT/(name+'.pdf'));plt.close(fig)

def figures():
 rows=json.loads((OUT/'curve_results.json').read_text());choices=json.loads((OUT/'choices.json').read_text());diag=json.loads((OUT/'efa_diagnostics.json').read_text())
 plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
 fig,ax=plt.subplots(2,3,figsize=(14,8))
 for sample in SPECS:
  color=COLORS[sample];r=next(r for r in rows if r['sample']==sample and r['kind']=='chosen');q,y,e=np.loadtxt(OUT/(r['key']+'.dat')).T;f=r['extended'];g=r['conventional'];m=(q>=.012)&(q<=.25)
  ax[0,0].errorbar(q[m],y[m]/f['I0'],e[m]/f['I0'],fmt='.',ms=3,elinewidth=.5,color=color,label=sample[:2]);ax[0,0].set(yscale='log',xlabel='q (Å⁻¹)',ylabel='I(q) / I(0)',title='A  Measured scattering')
  qq=q[(q>=.012)&(q<=.055)];yy=y[(q>=.012)&(q<=.055)];ee=e[(q>=.012)&(q<=.055)]
  ax[0,1].errorbar(qq**2,np.log(yy/f['I0']),ee/yy,fmt='.',ms=3,elinewidth=.5,color=color)
  ax[0,1].plot(qq**2,np.log(shape(qq,f['nu'])),color=color,label=f"{sample[:2]} extended {f['Rg']:.2f} Å")
  mm=(q>=g['qmin'])&(q<=g['qmax']);ax[0,1].plot(q[mm]**2,np.log(g['I0']/f['I0'])-q[mm]**2*g['Rg']**2/3,'--',color=color,label=f"{sample[:2]} Guinier {g['Rg']:.2f} Å")
  ax[0,1].set(xlabel='q² (Å⁻²)',ylabel='ln[I(q) / I(0)]',title='B  Conventional and extended Guinier')
  ax[1,0].errorbar(q[m],q[m]**2*y[m]/f['I0'],q[m]**2*e[m]/f['I0'],fmt='.',ms=3,elinewidth=.5,color=color);ax[1,0].set(xlabel='q (Å⁻¹)',ylabel='q² I(q) / I(0) (Å⁻²)',title='D  Normalized Kratky')
  ax[1,1].errorbar(q[m]*f['Rg'],(q[m]*f['Rg'])**2*y[m]/f['I0'],(q[m]*f['Rg'])**2*e[m]/f['I0'],fmt='.',ms=3,elinewidth=.5,color=color);ax[1,1].set(xlabel='q Rg',ylabel='(q Rg)² I(q) / I(0)',title='E  Dimensionless Kratky')
  p=np.load(OUT/(r['key']+'_Pr.npz'));area=np.trapezoid(p['p'],p['r']);ax[0,2].plot(p['r'],p['p']/area,color=color,label=sample[:2]);ax[0,2].fill_between(p['r'],(p['p']-p['error'])/area,(p['p']+p['error'])/area,color=color,alpha=.15);ax[0,2].set(xlabel='r (Å)',ylabel='P(r) / area (Å⁻¹)',title='C  Pair-distance distribution')
  bandpath=OUT/f'{sample}_normalized_bands.npz'
  if bandpath.exists():
   band=np.load(bandpath);ax[0,0].fill_between(q[m],*band['normalized95'][:,m],color=color,alpha=.13);ax[1,0].fill_between(q[m],*band['kratky95'][:,m],color=color,alpha=.13);ax[1,1].fill_between(band['x'],*band['dimensionless95'],color=color,alpha=.13)
  pred=f['I0']*shape(qq,f['nu']);ax[1,2].plot(qq,(yy-pred)/ee,'.',color=color,label=sample[:2]);ax[1,2].set(xlabel='q (Å⁻¹)',ylabel='(I − extended fit) / σ',title='F  Fit residuals')
 for a in [ax[0,0],ax[0,1],ax[0,2]]:a.legend(fontsize=8)
 ax[1,2].axhline(0,color='0.5',lw=.7);emit(fig,'01_WT_AA_comparison')
 # Window comparison: include curves BEFORE fit summaries.
 fig,ax=plt.subplots(2,3,figsize=(14,8))
 for k,sample in enumerate(SPECS):
  q,Y,E=load(sample);frames=np.arange(120,221);A=matrix(Y.shape[1],frames,choices[sample]['primary_model']);I=Y@A;m=(q>=.012)&(q<=.25);tr=np.trapezoid(I[m],q[m],axis=0);ax[k,0].plot(frames,tr,color='0.3')
  rr=[r for r in rows if r['sample']==sample and r['kind'] in ['chosen','leading','apex','trailing','previous']]
  for j,r in enumerate(rr):
   q,y,e=np.loadtxt(OUT/(r['key']+'.dat')).T;f=r['extended'];color=f'C{j}';a,b=r['frames'];ax[k,0].axvspan(a,b,alpha=.13,color=color,label=f"{r['kind']} {a}–{b}");mm=(q>=.012)&(q<=.15);ax[k,1].errorbar(q[mm],y[mm]/f['I0'],e[mm]/f['I0'],fmt='.',ms=2,elinewidth=.4,color=color);ci=r['extended95'];ax[k,2].errorbar(f['Rg'],j,xerr=np.array([[f['Rg']-ci[0]],[ci[1]-f['Rg']]]),fmt='o',color=color)
  ax[k,0].set(title=sample[:2]+' SEC and tested windows',xlabel='Frame (zero based)',ylabel='Integrated scattering');ax[k,0].legend(fontsize=7);ax[k,1].set(yscale='log',xlabel='q (Å⁻¹)',ylabel='I / I0',title='Normalized measured curves');ax[k,2].set(yticks=range(len(rr)),yticklabels=[r['kind'] for r in rr],xlabel='Extended Rg (Å), conditional 95%',title='Window sensitivity')
 emit(fig,'02_window_selection')
 fig,ax=plt.subplots(2,3,figsize=(14,8))
 for k,sample in enumerate(SPECS):
  for model in list(MODELS)+['linear']:
   d=np.load(OUT/f'{sample}_{model}_efa.npz');q=d['q'];m=(q>=.012)&(q<=.15);s=d['singular_values'];ax[k,0].plot(range(1,5),s[:4],'o-',label=model)
   if d['C'].size:
    S=d['S'];C=d['C'];scale=C.mean(0);y=S[:,0]*scale[0];ax[k,1].plot(q[m],y[m],lw=1,label=model)
   r=next(r for r in rows if r['sample']==sample and r['kind']=='baseline_sensitivity' and r['model']==model);ax[k,2].plot(model,r['extended']['Rg'],'o')
  ax[k,0].set(yscale='log',xlabel='Singular mode',ylabel='Weighted singular value',title=sample[:2]+' rank evidence');ax[k,0].legend(fontsize=7);ax[k,1].set(xlabel='q (Å⁻¹)',ylabel='Other component contribution',title='Signed other I(q), mean ROI scale');ax[k,1].set_yscale('symlog',linthresh=.005);ax[k,1].axhline(0,color='k',lw=.5);ax[k,1].legend(fontsize=7);ax[k,2].set(ylabel='Chosen-window extended Rg (Å)',title='Buffer-model sensitivity');ax[k,2].tick_params(axis='x',rotation=45)
 emit(fig,'03_buffer_effect_on_EFA')
 fig,ax=plt.subplots(2,2,figsize=(11,8));sample='WT_092025';d=np.load(OUT/f'{sample}_post17_efa.npz');q=d['q'];C=d['C'];S=d['S'];fr=d['frames'];m=(q>=.012)&(q<=.25)
 for j,label in enumerate(['Early shoulder — identity unassigned','Main component']):
  key=f'{sample}_'+(['EFA_shoulder_or_other','EFA_main'][j]);qq,y,e=np.loadtxt(OUT/(key+'.dat')).T;ax[0,0].errorbar(qq,y,e,fmt='.',ms=2,elinewidth=.4,label=label)
  ax[0,1].plot(fr,C[:,j]/C[:,j].max(),label=label)
 ax[0,0].set(xlim=(.006,.25),xlabel='q (Å⁻¹)',ylabel='I(q): contribution in selected window',title='A  Both recovered component curves');ax[0,0].set_yscale('symlog',linthresh=.003);ax[0,0].legend(fontsize=7);ax[0,1].set(xlabel='Frame',ylabel='Elution / own maximum',title='B  EFA elution profiles');ax[0,1].legend(fontsize=7)
 for kind in ['shoulder','chosen','EFA_shoulder_or_other','EFA_main']:
  r=next(r for r in rows if r['sample']==sample and r['kind']==kind);q,y,e=np.loadtxt(OUT/(r['key']+'.dat')).T;mm=(q>=.012)&(q<=.15);norm=np.trapezoid(y[mm],q[mm]);ax[1,0].errorbar(q[mm],y[mm]/norm,e[mm]/abs(norm),fmt='.',ms=2,elinewidth=.4,label=kind)
 ax[1,0].set(xlabel='q (Å⁻¹)',ylabel='Area-normalized measured I(q)',title='C  Direct shoulder and main versus EFA');ax[1,0].set_yscale('symlog',linthresh=.1);ax[1,0].legend(fontsize=7)
 resid=(d['I'][m]-S[m]@C.T)/d['error'][m];im=ax[1,1].imshow(resid,aspect='auto',origin='lower',extent=[fr[0],fr[-1],q[m][0],q[m][-1]],vmin=-3,vmax=3,cmap='RdBu_r');fig.colorbar(im,ax=ax[1,1],label='Residual / σ');ax[1,1].set(xlabel='Frame',ylabel='q (Å⁻¹)',title='D  Reconstruction residuals')
 emit(fig,'04_WT_EFA')

if __name__=='__main__':run();figures()
