from run_efa import *
import shutil,importlib.metadata
# Validate numerical primitives on known synthetic data, independent of measurements.
rng=np.random.default_rng(17);q=np.linspace(.008,.15,100);t=np.arange(60)
C=np.c_[np.maximum(0,1-abs(t-18)/15),np.maximum(0,1-abs(t-40)/15)];C/=C.sum(0)
S=np.c_[np.exp(-q*q*30**2/3),np.exp(-q*q*20**2/3)];D=S@C.T;u,s,v=np.linalg.svd(D,full_matrices=False)
fw=ns['runEFA'](D);bw=ns['runEFA'](D,False)
assert np.allclose(fw[:len(s),-1],s)
assert np.allclose(bw[:len(s),-1],s)
ok,_,rot=ns['runRotation'](D,D,np.full_like(D,.001),np.array([[3,33],[25,55]]),[True,True],v.T,niter=3000,tol=1e-12)
assert ok
recon=rot['int']@(rot['M']*rot['C']).T
rel=np.linalg.norm(D-recon)/np.linalg.norm(D);assert rel<1e-5
fit=guinier(q,np.exp(-q*q*25**2/3),np.full(len(q),.001));assert abs(fit['Rg']-25)<1e-7
checks=dict(synthetic_reconstruction_relative_error=float(rel),synthetic_guinier_Rg=fit['Rg'],forward_backward_singular_values_match=True)
(OUT/'validation.json').write_text(json.dumps(checks,indent=2));print('checks',checks)
# Bootstrap clean-profile apparent Rg with re-estimation of both background and SVD.
boot=[];scans=[]
for name,a,b in [('WT_092025',155,180),('AA_052026',163,204)]:
 q,Y,E=load(name);y,e=correct(Y,E,name);m=(q>=.012)&(q<=.25);w=e[m,a:b+1].mean(1)
 for variant in (['default','early'] if name.startswith('WT') else ['default','buffer']):
  rr,data=analyze(name,a,b,variant=variant);qq,yy,ee,mm,DD,uu,ss,vv=data;c=vv[0];c=c/c.sum();mult=c/(c@c);spec=yy[:,a:b+1]@mult;se=np.sqrt(ee[:,a:b+1]**2@mult**2)
  for qm in [.008,.010,.012,.015,.018,.020]:scans.append(dict(sample=name,background=variant,qmin=qm,fit=guinier(q,spec,se,qmin=qm)))
 # Parametric perturbation around the observed data, not a species-count null.
 values=[]
 for j in range(200):
  noise=rng.normal(size=Y.shape)*E;noise,_=correct(noise,E,name)
  yy=y+noise;_,_,vt=np.linalg.svd(yy[m,a:b+1]/w[:,None],full_matrices=False);c=vt[0];c/=c.sum();mult=c/(c@c)
  fit=guinier(q,yy[:,a:b+1]@mult,np.sqrt(e[:,a:b+1]**2@mult**2))
  if fit['Rg'] is not None:values.append(fit['Rg'])
 # Difference-based check of supplied uncertainty scale in pre-peak residual background.
 idx=np.arange(30,76);dz=np.diff(Y[m][:,idx],axis=1)/np.sqrt(E[m][:,idx[:-1]]**2+E[m][:,idx[1:]]**2)
 scale=float(np.median(abs(dz-np.median(dz)))/.67448975)
 br=dict(sample=name,bootstrap_n=200,valid=len(values),perturbation_percentiles=np.quantile(values,[.025,.5,.975]).tolist(),background_difference_noise_scale=scale)
 boot.append(br);print('bootstrap',br)
(OUT/'guinier_sensitivity.json').write_text(json.dumps(scans,indent=2));(OUT/'rg_bootstrap.json').write_text(json.dumps(boot,indent=2))
# Export compact summary figure.
fig,ax=plt.subplots(2,2,figsize=(12,8))
for j,(name,a,b) in enumerate([('WT_092025',155,180),('AA_052026',163,204)]):
 q,Y,E=load(name);y,e=correct(Y,E,name);m=(q>=.012)&(q<=.25);trace=np.trapezoid(y[m],q[m],axis=0)
 ax[0,j].plot(trace,color=['navy','orangered'][j]);ax[0,j].axvspan(a,b,color='green',alpha=.16,label=f'Clean region: {a}–{b}');ax[0,j].set(xlim=(110,225),xlabel='Frame (zero indexed)',ylabel='Integrated I (arbitrary)',title=name);ax[0,j].legend()
 r,data=analyze(name,a,b);q,y,e,m,D,u,s,vt=data
 c=vt[0]/vt[0].sum();mult=c/(c@c);spec=y[:,a:b+1]@mult;se=np.sqrt(e[:,a:b+1]**2@mult**2);fit=guinier(q,spec,se)
 qm=(q>=.008)&(q<=.05)&(spec>0);ax[1,j].errorbar(q[qm]**2,np.log(spec[qm]),yerr=se[qm]/spec[qm],fmt='.',color=['navy','orangered'][j],alpha=.7)
 fitm=(q>=fit['qmin'])&(q<=fit['qmax']);ax[1,j].plot(q[fitm]**2,np.log(fit['I0'])-q[fitm]**2*fit['Rg']**2/3,c='black',label=f"Apparent Rg = {fit['Rg']:.1f} Å")
 ax[1,j].set(xlabel='q² (Å⁻²)',ylabel='ln I (arbitrary scale)',title='Dominant clean-region profile; qRg ≤ 1');ax[1,j].legend()
fig.suptitle('Clean peaks: one dominant scattering pattern; additional oligomers not resolved');fig.tight_layout();fig.savefig(OUT/'clean_regions_summary.png',dpi=170);plt.close(fig)
# Provenance for copied source and original inputs.
git='/Users/f0044gk/.cache/codex-runtimes/codex-primary-runtime/dependencies/bin/fallback/git'
prov=dict(input_repository='https://github.com/Natalie-Loui/BNL_SAXS',input_commit=subprocess.check_output([git,'-C',str(ROOT/'BNL_SAXS'),'rev-parse','HEAD'],text=True).strip(),raw_repository='https://github.com/jbhopkins/bioxtasraw',raw_commit=subprocess.check_output([git,'-C',str(ROOT/'bioxtasraw-source'),'rev-parse','HEAD'],text=True).strip(),python=__import__('sys').version,packages={n:importlib.metadata.version(n) for n in ['numpy','scipy','h5py','matplotlib','nbformat']},input_sha256={name:hashlib.sha256((ROOT/'BNL_SAXS/data'/name).read_bytes()).hexdigest() for name in ['WT_092025.hdf5','AA_052026.hdf5']})
(OUT/'provenance.json').write_text(json.dumps(prov,indent=2))
(ROOT/'requirements.txt').write_text('\n'.join(f'{k}=={v}' for k,v in prov['packages'].items())+'\npackaging=='+importlib.metadata.version('packaging')+'\n')
