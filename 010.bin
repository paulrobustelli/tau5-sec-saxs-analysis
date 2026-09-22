"""Controlled buffer-window sensitivity from input HDF5 profiles, subtracting ONCE."""
from pathlib import Path
import sys,os,json,hashlib
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/saxs_mpl');os.environ.setdefault('NUMBA_CACHE_DIR','/private/tmp/saxs_numba')
import numpy as np
import matplotlib.pyplot as plt
from inspect_data import load
from component_recovery import SPECS,solve
from independent_baseline import window_stats
from compare_all_iq import conventional,fit,scan,shape,quick_extended
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'bioxtasraw-source'))
from bioxtasraw import BIFT
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'buffer_window_test';OUT.mkdir(exist_ok=True)
REQUESTED={'early_25_35':(25,35),'early_60_70':(60,70),'post_210_220':(210,220)}

def fits(q,y,e):
 g=conventional(q,y,e,.012,.030);m=(q>=.012-1e-9)&(q<=.055+1e-9);x=fit(q[m],y[m],e[m])
 if x and (x['boundary'] or x['qRgmax']>2):x=None
 return dict(conventional=g,conventional_within_IDP_limit=bool(g and g['qRgmax']<=1.1),extended=x,extended_adaptive=scan(q,y,e)[0])

def run(nboot=100,do_bift=True):
 results=[];checks=[];paired=[]
 for sample,spec in SPECS.items():
  q,Y,E=load(sample);frames=np.arange(spec['roi'][0],spec['roi'][1]+1);anchors=dict(REQUESTED);anchors['Natalie_buffer_reference']=(220,236) if sample.startswith('WT') else (24,70)
  windows={'selected_average':spec['avg'],'leading_5':([161,165] if sample.startswith('WT') else [165,169]),'apex_5':([165,169] if sample.startswith('WT') else [175,179])}
  data={}
  for label,(lo,hi) in anchors.items():
   A=np.eye(Y.shape[1])[:,frames];A[lo:hi+1]-=1/(hi-lo+1);I=Y@A;err=np.sqrt(E**2@(A*A));av=(frames>=spec['avg'][0])&(frames<=spec['avg'][1]);W=np.array([A[:,(frames>=a)&(frames<=b)].mean(1) for a,b in windows.values()]).T
   ys=Y@W;es=np.sqrt(E**2@(W*W));info,o=solve(q,I,err,frames,*spec['fid']);entry=dict(A=A,I=I,error=err,W=W,ys=ys,es=es,info=info,rotation=o,boot_main=[],boot_other=[],boot_direct=[]);data[label]=entry
   checks.append(dict(sample=sample,buffer=label,**window_stats(q,Y,E,lo,hi),note='Stationarity only, not proof of protein-free buffer. 210–220 is explicitly a tail-contamination stress test.'))
   if o is not None:
    scale=o['C'][av].mean(0);entry['components']=o['S']*scale[None,:];LW=A@o['L']*scale[None,:];entry['component_error']=np.sqrt(E**2@(LW*LW))
   print(sample,label,'rotation',info['converged'],flush=True)
  # Same perturbation is used for every buffer choice and sample window.
  rng=np.random.default_rng(20260925 if sample.startswith('WT') else 20260926)
  for rep in range(nboot):
   Yn=Y+rng.normal(size=Y.shape)*E
   for label,d in data.items():
    d['boot_direct'].append(Yn@d['W']);ii=Yn@d['A'];info,o=solve(q,ii,d['error'],frames,*spec['fid'])
    if o is None:d['boot_main'].append(np.full_like(q,np.nan));d['boot_other'].append(np.full_like(q,np.nan))
    else:
     comps=o['S']*o['C'][av].mean(0)[None,:];d['boot_main'].append(comps[:,1]);d['boot_other'].append(comps[:,0])
  paired_rg={}
  for label,d in data.items():
   allcurves={name:(d['ys'][:,j],d['es'][:,j],np.array(d['boot_direct'])[:,:,j]) for j,name in enumerate(windows)}
   if d['rotation'] is not None:
    bm=np.array(d['boot_main']);bo=np.array(d['boot_other']);valid=np.all(np.isfinite(bm),axis=1);count=int(valid.sum());d['info']['bootstrap_converged']=count;d['info']['bootstrap_requested']=nboot
    if count>=20:
     allcurves['EFA_main']=(d['components'][:,1],np.std(bm[valid],axis=0,ddof=1),bm)
     allcurves['EFA_other']=(d['components'][:,0],np.std(bo[valid],axis=0,ddof=1),bo)
   for kind,(y,e,draws) in allcurves.items():
    key=f'{sample}_{label}_{kind}';np.savetxt(OUT/(key+'.dat'),np.c_[q,y,e],header='q_A^-1 I sigma; input profiles minus ONE averaged buffer; no additional baseline')
    f=fits(q,y,e);r=dict(sample=sample,buffer=label,buffer_frames=anchors[label],curve=kind,sample_frames=windows.get(kind,spec['avg']),path='buffer_window_test/'+key+'.dat',**f,rotation=d['info'] if kind.startswith('EFA') else None)
    br=[];bg=[]
    for b in draws:
     if not np.all(np.isfinite(b)):br.append(np.nan);bg.append(np.nan);continue
     v=quick_extended(q,b,e,f['extended']);g=conventional(q,b,e,.012,.030);br.append(v if v else np.nan);bg.append(g['Rg'] if g else np.nan)
    br=np.array(br);bg=np.array(bg);paired_rg[label,kind]=br
    r['extended_bootstrap95']=np.nanpercentile(br,[2.5,97.5]).tolist() if np.isfinite(br).sum()>=20 else None
    r['conventional_bootstrap95']=np.nanpercentile(bg,[2.5,97.5]).tolist() if np.isfinite(bg).sum()>=20 else None
    r['uncertainty']='Paired input-frame perturbations include shared buffer noise and refitted EFA; fixed supports/windows/q ranges. No source q covariance or buffer-purity uncertainty.'
    if do_bift and kind in ['selected_average','EFA_main']:
     m=(q>=.012-1e-9)&(q<=.25+1e-9);np.random.seed(20260925)
     obj=BIFT.doBift(q[m],y[m],e[m],key,40,150,1e10,12,50,220,15,80,single_proc=True,nprocs=1)
     if obj is not None:
      pars=obj.getAllParameters();r['bift']={k:float(pars[k]) for k in ['rg','dmax','chisq']};np.savez_compressed(OUT/(key+'_Pr.npz'),r=obj.r,p=obj.p,p_error=obj.err,q=obj.q_orig,I=obj.i_orig,error=obj.err_orig,fit=obj.i_fit)
      r['bift']['uncertainty']='Conditional RAW BIFT settings; Dmax regularization/model sensitivity remains separate.'
     else:r['bift']=None
    results.append(r);print(key,'extended',None if f['extended'] is None else round(f['extended']['Rg'],3),flush=True)
   np.savez_compressed(OUT/f'{sample}_{label}_propagation.npz',q=q,frames=frames,buffer=np.array(anchors[label]),W=d['W'],boot_direct=np.array(d['boot_direct']),boot_main=np.array(d['boot_main']),boot_other=np.array(d['boot_other']),C=d['rotation']['C'] if d['rotation'] is not None else np.empty((0,2)))
  for (label,kind),v in paired_rg.items():
   if label=='Natalie_buffer_reference':continue
   ref=paired_rg.get(('Natalie_buffer_reference',kind))
   if ref is None:continue
   delta=v-ref;good=np.isfinite(delta);paired.append(dict(sample=sample,buffer=label,curve=kind,reference='Natalie_buffer_reference',extended_delta95=np.percentile(delta[good],[2.5,97.5]).tolist() if good.sum()>=20 else None,successful_pairs=int(good.sum())))
  (OUT/'results.json').write_text(json.dumps(results,indent=2));(OUT/'anchor_checks.json').write_text(json.dumps(checks,indent=2));(OUT/'paired_differences.json').write_text(json.dumps(paired,indent=2))
 return results

def show_results():
 from IPython.display import display,Markdown
 rows=json.loads((OUT/'results.json').read_text())
 h='|Sample|Buffer frames|Curve|Conventional Rg|Extended Rg|Extended95%|BIFT Rg|\n|---|---|---|---:|---:|---|---:|\n'
 for r in rows:
  if r['curve']=='EFA_other':continue
  val=lambda x:'unsupported' if x is None else f"{x['Rg']:.2f}"
  h+=f"|{r['sample'][:2]}|{r['buffer_frames']}|{r['curve']}|{val(r['conventional'])}|{val(r['extended'])}|{r['extended_bootstrap95']}|{None if not r.get('bift') else round(r['bift']['rg'],2)}|\n"
 display(Markdown(h))
 for sample in SPECS:
  for kind in ['selected_average','leading_5','apex_5','EFA_main','EFA_other']:
   rr=[r for r in rows if r['sample']==sample and r['curve']==kind];fig,ax=plt.subplots(2,2,figsize=(11,7))
   for color_index,r in enumerate(rr):
    color=plt.cm.tab10.colors[color_index]
    q,y,e=np.loadtxt(ROOT/r['path']).T;m=(q>=.012)&(q<=.25);ax[0,0].errorbar(q[m],y[m],e[m],fmt='.',ms=2,elinewidth=.3,color=color,label=r['buffer']);f=r['extended']
    if r['conventional']:
     g=r['conventional'];mg=(q>=.012-1e-9)&(q<=.030+1e-9);ax[0,1].plot(q[mg]**2,g['I0']*np.exp(-q[mg]**2*g['Rg']**2/3),'--',lw=.8,alpha=.7,color=color)
    if f:
     m=(q>=.012-1e-9)&(q<=.055+1e-9);line=ax[0,1].plot(q[m]**2,y[m],'.',ms=3,color=color,label=r['buffer'])[0];ax[0,1].plot(q[m]**2,f['I0']*shape(q[m],f['nu']),color=line.get_color());m=(q>=.012)&(q<=.25);ax[1,0].plot(q[m]*f['Rg'],(q[m]*f['Rg'])**2*y[m]/f['I0'],label=r['buffer'])
    if r.get('bift'):
     a=np.load(OUT/(Path(r['path']).stem+'_Pr.npz'));area=np.trapezoid(a['p'],a['r']);ax[1,1].plot(a['r'],a['p']/area,label=r['buffer'])
   ax[0,0].set_yscale('symlog',linthresh=.01);ax[0,0].set(xlabel='q (Å⁻¹)',ylabel='I(q)',title='Input − one buffer average');ax[0,0].legend(fontsize=7)
   ax[0,1].set(xlabel='q² (Å⁻²)',ylabel='I(q)',title='Guinier dashed; extended solid (fixed q ranges)')
   ax[1,0].set(xlabel='q Rg',ylabel='(q Rg)² I / I0',title='Dimensionless Kratky: measured data')
   ax[1,1].set(xlabel='r (Å)',ylabel='Normalized P(r)',title='Same BIFT settings (average and EFA main only)')
   fig.suptitle(sample[:2]+' — '+kind+' — buffer-window sensitivity');fig.tight_layout();fig.savefig(OUT/f'{sample}_{kind}_comparison.png',dpi=130);display(fig);plt.close(fig)
 return rows
if __name__=='__main__':run()

def compare_samples():
 rows=json.loads((OUT/'results.json').read_text());diffs=[]
 for label in list(REQUESTED)+['Natalie_buffer_reference']:
  for kind in ['selected_average','EFA_main']:
   vals=[];point=[]
   for sample in SPECS:
    rr=[r for r in rows if r['sample']==sample and r['buffer']==label and r['curve']==kind]
    if not rr or rr[0]['extended'] is None:break
    r=rr[0];q,y,e=np.loadtxt(ROOT/r['path']).T;a=np.load(OUT/f'{sample}_{label}_propagation.npz');draws=a['boot_direct'][:,:,0] if kind=='selected_average' else a['boot_main'];v=[]
    for d in draws:
     f=quick_extended(q,d,e,r['extended']) if np.all(np.isfinite(d)) else None;v.append(f if f is not None else np.nan)
    vals.append(np.array(v));point.append(r['extended']['Rg'])
   if len(vals)==2:
    d=vals[1]-vals[0];good=np.isfinite(d);diffs.append(dict(buffer=label,curve=kind,AA_minus_WT=point[1]-point[0],conditional95=np.percentile(d[good],[2.5,97.5]).tolist(),valid_pairs=int(good.sum()),note='Independent sample perturbations; not biological replicate variability. Natalie-reference buffers differ between samples.'))
 (OUT/'AA_minus_WT.json').write_text(json.dumps(diffs,indent=2));return diffs

def show_anchors():
 from IPython.display import display
 fig,ax=plt.subplots(2,2,figsize=(11,7))
 for col,sample in enumerate(SPECS):
  q,Y,E=load(sample);m=(q>=.012)&(q<=.035);v=Y[m].mean(0);e=np.sqrt((E[m]**2).sum(0))/m.sum();t=np.arange(len(v))
  ax[0,col].plot(t,v,color='0.3',lw=1)
  for (label,(lo,hi)),color in zip(REQUESTED.items(),['C0','C1','C2']):ax[0,col].axvspan(lo,hi,alpha=.3,color=color,label=label)
  ax[0,col].set(xlim=(0,250),xlabel='Frame',ylabel='Mean input I(q), q=.012–.035',title=sample[:2]+' buffer candidates');ax[0,col].legend(fontsize=7)
  m=(t>=195)&(t<=235);ax[1,col].errorbar(t[m],v[m],e[m],fmt='.-',ms=3,lw=.5,color='0.3');ax[1,col].axvspan(210,220,alpha=.3,color='C2');ax[1,col].set(xlabel='Frame',ylabel='Mean input I(q)',title='Post-peak zoom: 210–220 shaded')
 fig.tight_layout();fig.savefig(OUT/'buffer_anchor_traces.png',dpi=140);display(fig);plt.close(fig)
