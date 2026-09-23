"""Additional controls for Paper Figures; no frame or baseline tuning to an Rg target."""
from paper_figures import *

def run_controls():
 rows=json.loads((OUT/'curve_results.json').read_text());choices=json.loads((OUT/'choices.json').read_text());stats=[];pr=[]
 for si,sample in enumerate(SPECS):
  q,Y,E=load(sample);rng=np.random.default_rng(403+si)
  for model in list(MODELS)+['linear']:
   d=np.load(OUT/f'{sample}_{model}_efa.npz');frames=d['frames'];A=matrix(Y.shape[1],frames,model);m=(q>=.012)&(q<=.25);I=d['I'][m];err=d['error'][m];scale=err.mean(1);D=I/scale[:,None];U,s,Vt=np.linalg.svd(D,full_matrices=False);null=(U[:,:1]*s[:1])@Vt[:1]*scale[:,None];sv=[]
   for rep in range(60):
    noise=(rng.normal(size=E.shape)*E)@A;sv.append(np.linalg.svd((null+noise[m])/scale[:,None],compute_uv=False)[1])
   stats.append(dict(sample=sample,model=model,observed_s2=float(s[1]),rank1_noise_s2_95=float(np.percentile(sv,95)),null_exceedances=int(np.sum(np.array(sv)>=s[1])),n=60,note='Conditional rank-one noise diagnostic; shared-buffer covariance retained; no q covariance or baseline-model uncertainty.'))
  rr=next(r for r in rows if r['sample']==sample and r['kind']=='chosen');q,y,e=np.loadtxt(OUT/(rr['key']+'.dat')).T;a=np.load(OUT/f'{sample}_bootstrap.npz');draws=a['direct'][:,:,0];mm=(q>=.012-1e-9)&(q<=.055+1e-9);vals=[];ni=[];kr=[];dk=[];xgrid=np.linspace(.4,6,120)
  for b in draws:
   f=fit(q[mm],b[mm],e[mm])
   if f is None or f['boundary'] or f['qRgmax']>2:continue
   vals.append([f['Rg'],f['I0']]);ni.append(b/f['I0']);kr.append(q*q*b/f['I0']);dk.append(np.interp(xgrid,q*f['Rg'],(q*f['Rg'])**2*b/f['I0'],left=np.nan,right=np.nan))
  np.savez_compressed(OUT/f'{sample}_normalized_bands.npz',q=q,normalized95=np.percentile(ni,[2.5,97.5],axis=0),kratky95=np.percentile(kr,[2.5,97.5],axis=0),x=xgrid,dimensionless95=np.nanpercentile(dk,[2.5,97.5],axis=0),fit_draws=np.array(vals))
  # BIFT under alternative admissible backgrounds, plus deliberately varied initial Dmax search upper limits.
  for model in ['early47','post17','linear']:
   key=f'{sample}_chosen_{model}';q,y,e=np.loadtxt(OUT/(key+'.dat')).T;f=bift(q,y,e,key);pr.append(dict(sample=sample,model=model,key=key,cap=220,bift=f))
  for cap in [160,280]:
   key=rr['key']+f'_Dmax{cap}';q,y,e=np.loadtxt(OUT/(rr['key']+'.dat')).T;f=bift(q,y,e,key,(50,cap));pr.append(dict(sample=sample,model=choices[sample]['primary_model'],key=key,cap=cap,bift=f))
 save(stats,'rank1_noise_controls.json');save(pr,'Pr_sensitivity.json')
 return stats

def control_figures():
 rows=json.loads((OUT/'curve_results.json').read_text());pr=json.loads((OUT/'Pr_sensitivity.json').read_text());fig,ax=plt.subplots(2,2,figsize=(11,8))
 for k,sample in enumerate(SPECS):
  for r in pr:
   if r['sample']!=sample or not r['bift']:continue
   p=np.load(OUT/(r['key']+'_Pr.npz'));area=np.trapezoid(p['p'],p['r']);ax[k,0].plot(p['r'],p['p']/area,label=f"{r['model']}, grid≤{r['cap']}: Rg{r['bift']['rg']:.1f}")
  ax[k,0].set(xlabel='r (Å)',ylabel='Normalized P(r)',title=sample[:2]+' baseline / Dmax sensitivity');ax[k,0].legend(fontsize=7)
  for r in rows:
   if r['sample']!=sample or r['kind'] not in ['chosen','leading','apex','trailing','previous']:continue
   p=np.load(OUT/(r['key']+'_Pr.npz'));area=np.trapezoid(p['p'],p['r']);ax[k,1].plot(p['r'],p['p']/area,label=r['kind'])
  ax[k,1].set(xlabel='r (Å)',ylabel='Normalized P(r)',title='Fixed primary buffer: frame-window comparison');ax[k,1].legend(fontsize=7)
 emit(fig,'05_Pr_controls')
 # Full-frame baseline overview with the anchors used in the primary figure.
 fig,ax=plt.subplots(2,1,figsize=(11,7));choices=json.loads((OUT/'choices.json').read_text())
 for k,sample in enumerate(SPECS):
  q,Y,E=load(sample);m=(q>=.012)&(q<=.035);tr=Y[m].mean(0);ax[k].plot(np.arange(250),tr[:250],color='0.3',lw=1);a,b=choices[sample]['buffer'];ax[k].axvspan(a,b,color='C2',alpha=.25,label=f'Primary buffer {a}–{b}');a,b=choices[sample]['window'];ax[k].axvspan(a,b,color=COLORS[sample],alpha=.25,label=f'Selected average {a}–{b}');ax[k].set(xlabel='Frame (zero based)',ylabel='Mean input I, q .012–.035',title=sample[:2]+' source profiles and selected regions');ax[k].legend()
 emit(fig,'06_buffer_and_SEC')
if __name__=='__main__':run_controls();control_figures()

def efa_progression():
 fig,ax=plt.subplots(1,2,figsize=(11,4));choices=json.loads((OUT/'choices.json').read_text())
 for k,sample in enumerate(SPECS):
  d=np.load(OUT/f"{sample}_{choices[sample]['primary_model']}_efa.npz");m=(d['q']>=.012)&(d['q']<=.25);D=d['I'][m]/d['error'][m].mean(1)[:,None];t=d['frames'];f=[];b=[]
  for j in range(3,len(t)-2):
   f.append(np.linalg.svd(D[:,:j+1],compute_uv=False)[:3]);b.append(np.linalg.svd(D[:,j:],compute_uv=False)[:3])
  for j in range(3):
   ax[k].plot(t[3:-2],np.array(f)[:,j],color=f'C{j}',label=f'Mode{j+1} forward');ax[k].plot(t[3:-2],np.array(b)[:,j],'--',color=f'C{j}',label=f'Mode{j+1} backward')
  ax[k].set(yscale='log',xlabel='Frame',ylabel='Weighted singular value',title=sample[:2]+' forward / backward EFA');ax[k].legend(fontsize=7)
 emit(fig,'07_forward_backward_EFA')

def window_kratky():
 rows=json.loads((OUT/'curve_results.json').read_text());fig,ax=plt.subplots(2,2,figsize=(11,8))
 for k,sample in enumerate(SPECS):
  for r in rows:
   if r['sample']!=sample or r['kind'] not in ['chosen','leading','apex','trailing','previous']:continue
   q,y,e=np.loadtxt(OUT/(r['key']+'.dat')).T;m=(q>=.012)&(q<=.25);f=r['extended'];ax[k,0].plot(q[m],q[m]**2*y[m]/f['I0'],'.-',ms=2,lw=.7,label=r['kind']);ax[k,1].plot(q[m]*f['Rg'],(q[m]*f['Rg'])**2*y[m]/f['I0'],'.-',ms=2,lw=.7,label=r['kind'])
  ax[k,0].set(xlabel='q (Å⁻¹)',ylabel='q² I / I0 (Å⁻²)',title=sample[:2]+' window Kratky comparison');ax[k,1].set(xlabel='q Rg',ylabel='(q Rg)² I / I0',title='Dimensionless Kratky');ax[k,0].legend(fontsize=8)
 emit(fig,'08_window_Kratky')
