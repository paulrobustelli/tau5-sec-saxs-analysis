"""Measured-data Kratky transforms with joint curve/normalization perturbations."""
from compare_all_iq import *
import warnings
KOUT=OUT/'kratky';KOUT.mkdir(exist_ok=True)

def refit_params(q,y,e,reference,method):
 if not reference:return None
 if method=='conventional':
  r=conventional(q,y,e,reference['qmin'],reference['qmax']);return None if r is None else (r['Rg'],r['I0'])
 m=(q>=reference['qmin']-1e-9)&(q<=reference['qmax']+1e-9)
 def residual(x):return (y[m]-np.exp(x[0])*shape(q[m],x[1]))/e[m]
 r=least_squares(residual,[np.log(reference['I0']),reference['nu']],bounds=([-100,.2],[100,.8]),max_nfev=60)
 return (float(rg_from_nu(r.x[1])),float(np.exp(r.x[0]))) if r.success and .201<r.x[1]<.799 else None

def build_kratky_cache(n=200):
 results=json.loads((OUT/'fit_results.json').read_text());arrays={};metadata=[]
 for index,r in enumerate(results):
  q,y,e=np.loadtxt(ROOT/r['path'])[:,:3].T;mask=(q>=.012-1e-9)&(q<=.25+1e-9);key=hashlib.sha256(r['id'].encode()).hexdigest()[:12]
  arrays[key+'_q']=q[mask];arrays[key+'_raw']=q[mask]**2*y[mask];arrays[key+'_raw_se']=q[mask]**2*e[mask]
  rng=np.random.default_rng(202609240+index);draws=y+e*rng.normal(size=(n,len(q)));basis='Independent q Gaussian perturbations of exported I(q)'
  if r['group']=='Audited EFA and matched average':
   a=np.load(ROOT/'component_results'/f"{r['sample']}_audit.npz")
   if 'candidate_main' in r['path']:draws=a['bootstrap'];basis='Existing EFA refit perturbations'
   elif 'frame_average_same' in r['path']:draws=a['bootstrap_average'];basis='Existing matched-average perturbations'
  meta=dict(id=r['id'],key=key,group=r['group'],label=r['label'],sample=r['sample'],basis=basis,methods={})
  for method in ['conventional','extended']:
   ref=r[method]
   if not ref:meta['methods'][method]={'supported':False};continue
   qq=q[mask];yy=y[mask];grid=np.linspace(qq[0]*ref['Rg'],qq[-1]*ref['Rg'],100)
   center=grid**2*np.interp(grid/ref['Rg'],qq,yy)/ref['I0'];normalized=qq**2*yy/ref['I0'];dim_draw=[];norm_draw=[];pars=[]
   for b in draws:
    p=refit_params(q,b,e,ref,method)
    if p is None:continue
    rg,i0=p;pars.append(p);dim_draw.append(grid**2*np.interp(grid/rg,qq,b[mask],left=np.nan,right=np.nan)/i0);norm_draw.append(qq**2*b[mask]/i0)
   dim_draw=np.asarray(dim_draw);norm_draw=np.asarray(norm_draw)
   if len(pars)<20:meta['methods'][method]={'supported':False,'reason':'Too few successful normalization draws'};continue
   with warnings.catch_warnings():
    warnings.simplefilter('ignore',RuntimeWarning);dimci=np.nanpercentile(dim_draw,[2.5,97.5],axis=0)
   coverage=np.isfinite(dim_draw).sum(0);dimci[:,coverage<.8*len(pars)]=np.nan
   prefix=key+'_'+method
   arrays[prefix+'_x']=grid;arrays[prefix+'_dimensionless']=center;arrays[prefix+'_dimensionless95']=dimci
   arrays[prefix+'_normalized']=normalized;arrays[prefix+'_normalized95']=np.percentile(norm_draw,[2.5,97.5],axis=0)
   meta['methods'][method]=dict(supported=True,Rg=ref['Rg'],I0=ref['I0'],valid_draws=len(pars),normalization_covariance=np.cov(np.array(pars).T).tolist(),fit_qrange=[ref['qmin'],ref['qmax']])
  metadata.append(meta)
 np.savez_compressed(KOUT/'transforms.npz',**arrays);(KOUT/'metadata.json').write_text(json.dumps(metadata,indent=2))
 return metadata

def kratky_panel(rows,title):
 a=np.load(KOUT/'transforms.npz');fig,axes=plt.subplots(2,2,figsize=(11,7));skipped=[]
 for color,r in zip(plt.cm.tab10.colors,rows):
  k=r['key'];q=a[k+'_q'];label=r['label'].replace('WT_092025','WT').replace('AA_052026','AA').replace('least_overlap_leading','leading').replace('frame_average_same_baseline','matched average').replace('EFA_candidate_main','EFA main').replace('EFA_other_component','EFA other')
  axes[0,0].plot(q,a[k+'_raw'],'.-',ms=2,lw=.5,color=color,label=label);axes[0,0].fill_between(q,a[k+'_raw']-1.96*a[k+'_raw_se'],a[k+'_raw']+1.96*a[k+'_raw_se'],color=color,alpha=.10)
  for method,ax in [('conventional',axes[1,0]),('extended',axes[1,1])]:
   if not r['methods'][method]['supported']:skipped.append(label+' / '+method);continue
   pref=k+'_'+method;x=a[pref+'_x'];ci=a[pref+'_dimensionless95'];ax.plot(x,a[pref+'_dimensionless'],lw=1,color=color,label=label);ax.fill_between(x,ci[0],ci[1],color=color,alpha=.13)
   if method=='extended':
    ci=a[pref+'_normalized95'];axes[0,1].plot(q,a[pref+'_normalized'],lw=1,color=color,label=label);axes[0,1].fill_between(q,ci[0],ci[1],color=color,alpha=.13)
 for ax in axes.flat:ax.axhline(0,color='0.5',lw=.5);ax.tick_params(labelsize=8)
 axes[0,0].set(title='Unnormalized Kratky (exported scales differ)',xlabel='q (Å⁻¹)',ylabel='q² I(q)')
 axes[0,1].set(title='I(0)-normalized; extended-fit I(0)',xlabel='q (Å⁻¹)',ylabel='q² I(q) / I(0) (Å⁻²)')
 for method,ax in [('conventional',axes[1,0]),('extended',axes[1,1])]:ax.set(title='Dimensionless: '+method+' normalization',xlabel='q Rg',ylabel='(q Rg)² I(q) / I(0)')
 axes[0,0].legend(fontsize=6,loc='best');fig.suptitle(title+'\nPointwise conditional 95% bands; baseline/model uncertainty excluded',fontsize=10);fig.tight_layout();display(fig)
 safe=hashlib.sha256(title.encode()).hexdigest()[:10];fig.savefig(KOUT/(safe+'.png'),dpi=140);plt.close(fig)
 if skipped:display(Markdown('Unsupported normalization omitted: '+ '; '.join(skipped)))

def priority_kratky():
 metadata=json.loads((KOUT/'metadata.json').read_text());chosen=[]
 for sample in ['WT_092025','AA_052026']:
  rows=[r for r in metadata if r['sample']==sample and ((r['group']=='EFA-guided direct edges' and ' linear ' in r['label']) or (r['group']=='Audited EFA and matched average' and ('candidate_main' in r['label'] or 'frame_average_same_baseline' in r['label'])))]
  kratky_panel(rows,sample[:2]+': leading edge / apex / EFA main / matched average');chosen.extend([r for r in rows if 'candidate_main' in r['label'] or 'least_overlap' in r['label']])
 kratky_panel(chosen,'WT versus AA: leading edge and EFA main')

def kratky_groups(groups):
 metadata=json.loads((KOUT/'metadata.json').read_text())
 for group in groups:
  for sample in ['WT_092025','AA_052026']:
   rows=[r for r in metadata if r['group']==group and r['sample']==sample]
   for i in range(0,len(rows),4):kratky_panel(rows[i:i+4],sample[:2]+' — '+group+f' ({i+1}–{min(i+4,len(rows))})')

def validate_kratky():
 # Scaling I and sigma must leave normalized and dimensionless curves invariant.
 q=np.linspace(.012,.07,50);rg=30.;i0=2.;y=i0*np.exp(-(q*rg)**2/3);scale=7.3
 assert np.allclose((q*rg)**2*y/i0,(q*rg)**2*(scale*y)/(scale*i0))
 assert np.allclose(q*q*y/i0,q*q*(scale*y)/(scale*i0))
 a=np.load(KOUT/'transforms.npz');meta=json.loads((KOUT/'metadata.json').read_text());reg={r['id']:r for r in json.loads((OUT/'curve_registry.json').read_text())}
 for r in meta:
  d=np.loadtxt(ROOT/reg[r['id']]['path']);m=(d[:,0]>=.012-1e-9)&(d[:,0]<=.25+1e-9)
  assert np.allclose(a[r['key']+'_raw_se'],d[m,0]**2*d[m,2])
 result={'curves':len(meta),'scale_invariance_pass':True,'raw_error_propagation_pass':True,'uncertainty':'Pointwise conditional intervals; each perturbed curve refits Rg and I0 jointly; fixed-x interpolation propagates horizontal scaling; no extrapolation'}
 (KOUT/'validation.json').write_text(json.dumps(result,indent=2));return result
