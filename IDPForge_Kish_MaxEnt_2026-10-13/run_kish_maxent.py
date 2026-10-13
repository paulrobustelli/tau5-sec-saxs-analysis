"""Borthakur repository dual objective + per-datatype/global Kish calibration.
Backcalculation collector waits for the actual native predictor jobs.
"""
from pathlib import Path
import ast,json,time,sys,hashlib
import numpy as np,pandas as pd,scipy as sp,mdtraj as md
from scipy.optimize import minimize
from scipy.special import softmax
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(__file__).resolve().parents[2];O=Path(__file__).resolve().parent
source=R/'work/Borthakur_MaxEnt_IDPs_2024/PaaA2/reweighting_examples/reweight_functions.py'
tree=ast.parse(source.read_text());fn=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='func_sigma_reg_sigma_md')
code=ast.get_source_segment(source.read_text(),fn);(O/'borthakur_dual_objective.py').write_text(code+'\n');exec(compile(ast.Module(body=[fn],type_ignores=[]),str(source),'exec'))
manifest=json.load(open(O/'manifest.json'));experiment=pd.read_csv(O/'experimental/chemical_shifts.csv')
while True:
 done=[m for m in manifest if (Path(m['job'])/'done').exists()]
 failed=[m for m in done if int((Path(m['job'])/'done').read_text())!=0]
 if failed:raise RuntimeError('Native predictor failures: '+','.join(m['id'] for m in failed))
 if len(done)==len(manifest):break
 print('Backcalc',len(done),'/',len(manifest),flush=True);time.sleep(30)

def parse_sparta(path):
 rows={}
 for line in path.read_text().splitlines():
  f=line.split()
  if len(f)==9 and f[0].isdigit():
   atom='H' if f[2] in ['HN','NH'] else f[2];rows[(int(f[0]),atom)]=(float(f[4]),float(f[5]))
 return rows

def fit_at(X,y,base,sem,theta):
 # Centering subtracts a constant per measurement; it leaves the original dual objective and weights unchanged.
 A=(X-y[:,None]);norm=np.maximum(np.maximum(np.std(A,axis=1),base),1e-6);An=A/norm[:,None];sr=theta*base/norm;sm=sem/norm
 opt=minimize(func_sigma_reg_sigma_md,np.zeros(len(y)),args=(An,np.zeros(len(y)),np.ones(X.shape[1]),sr,sm),jac=True,method='L-BFGS-B',options={'maxiter':10000,'ftol':1e-13,'gtol':1e-8,'maxls':50})
 if not opt.success:raise RuntimeError(str(opt.message))
 w=softmax(-opt.x@An);K=1/(len(w)*np.dot(w,w))
 return w,float(K),float(opt.fun),float(np.max(abs(opt.jac)))

def calibrate(X,y,base,sem,label,records):
 def ev(t):
  w,k,f,g=fit_at(X,y,base,sem,t);records.append(dict(stage=label,theta=t,Kish_ratio=k,Kish_Neff=k*len(w),objective=f,max_gradient=g));return w,k
 w0,k0=ev(.01)
 if k0>=.1:return .01,w0,k0
 lo=.01;hi=1.;wh,kh=ev(hi)
 while kh<.1:
  hi*=2;wh,kh=ev(hi)
  if hi>1e8:raise RuntimeError('No Kish bracket')
 for _ in range(20):
  mid=np.sqrt(lo*hi);wm,km=ev(mid)
  if km>=.1:hi=mid;wh,kh=wm,km
  else:lo=mid
 return hi,wh,kh

summaries=[];comparisons=[]
for v in ['WT','AA']:
 members=[m for m in manifest if m['variant']==v];n=len(members);assert n==200
 tables=[parse_sparta(Path(m['job'])/'pred.tab') for m in members]
 obs=experiment[experiment.variant==v].copy();retained=[];excluded=[]
 for row in obs.to_dict('records'):
  key=(row['model_residue'],row['atom'])
  (retained if all(key in t for t in tables) else excluded).append(row)
 obs=pd.DataFrame(retained);pd.DataFrame(excluded).to_csv(O/f'{v}_excluded_shifts.csv',index=False)
 Xcs=np.array([[t[(r.model_residue,r.atom)][0] for t in tables] for r in obs.itertuples()]);rc=np.array([tables[0][(r.model_residue,r.atom)][1] for r in obs.itertuples()]);ycs=obs.value.to_numpy();groups=obs.atom.to_numpy()
 q,iy,sy=np.loadtxt(O/f'experimental/{v}_SAXS.dat',unpack=True)
 curves=[]
 for m in members:
  d=np.loadtxt(Path(m['job'])/'prediction.dat');assert d.shape==(501,8);assert np.isfinite(d).all();curves.append(np.interp(q,d[:,0],d[:,1]))
 raw=np.array(curves);scales=(raw/sy)@(iy/sy)/np.sum((raw/sy)**2,axis=1)
 # Match calc_saxs.py fixed-parameter native profile fitting (verified by explicit native experimental-profile pilot).
 saxs=raw*scales[:,None];Xs=saxs.T
 X=np.concatenate([Xcs,Xs]);y=np.r_[ycs,iy];types=np.r_[groups,np.repeat('SAXS',len(q))];sem=np.std(X,axis=1,ddof=1)/np.sqrt(n)
 np.savez_compressed(O/f'{v}_observable_matrix.npz',predicted=X,observed=y,datatype=types,sigma_ensemble=sem,raw_saxs=raw,pepsi_scales=scales,q=q,intensity_observed=iy,intensity_sigma=sy,chemical_random_coil=rc)
 # Gradient check of the verbatim repository objective.
 A=X-y[:,None];rng=np.random.default_rng(16);z=rng.normal(0,1e-4,len(y));direction=rng.normal(size=len(y));direction/=np.linalg.norm(direction);eps=1e-6
 f,g=func_sigma_reg_sigma_md(z,A,np.zeros(len(y)),np.ones(n),np.ones(len(y)),sem)
 fd=(func_sigma_reg_sigma_md(z+eps*direction,A,np.zeros(len(y)),np.ones(n),np.ones(len(y)),sem)[0]-func_sigma_reg_sigma_md(z-eps*direction,A,np.zeros(len(y)),np.ones(n),np.ones(len(y)),sem)[0])/(2*eps)
 assert np.isclose(fd,g@direction,rtol=1e-5,atol=1e-6)
 scan=[];theta={};weights={'uniform':np.full(n,1/n)};base=np.ones(len(y))
 for group in np.unique(types):
  mask=types==group;tt,ww,kk=calibrate(X[mask],y[mask],np.ones(mask.sum()),sem[mask],group,scan);theta[group]=tt;base[mask]=tt;weights[group+'_only']=ww
 global_t,w,k=calibrate(X,y,base,sem,'combined',scan);weights['combined']=w
 assert abs(w.sum()-1)<1e-12 and w.min()>=0 and k>=.1
 pd.DataFrame(scan).to_csv(O/f'{v}_regularization_scan.csv',index=False)
 pd.DataFrame({'file':[m['file'] for m in members],**weights}).to_csv(O/f'{v}_weights.csv',index=False)
 helices=[]
 for member in members:
  t=md.load(member['file']);helices.append((md.compute_dssp(t,simplified=False)=='H')[0])
 helices=np.array(helices);pd.DataFrame({'AR_residue':np.arange(330,448),'uniform':helices[:,2:].mean(0),'Kish_MaxEnt':w@helices[:,2:]}).to_csv(O/f'{v}_helicity.csv',index=False)
 shifts=obs.copy();shifts['random_coil_ppm']=rc;shifts['unweighted_ppm']=Xcs.mean(1);shifts['reweighted_ppm']=Xcs@w;shifts['experimental_secondary_ppm']=ycs-rc;shifts['unweighted_secondary_ppm']=Xcs.mean(1)-rc;shifts['reweighted_secondary_ppm']=Xcs@w-rc;shifts.to_csv(O/f'{v}_shift_fit.csv',index=False)
 fit=pd.DataFrame({'q_Ainv':q,'I_exp':iy,'sigma_exp':sy,'I_uniform':saxs.mean(0),'I_MaxEnt':w@saxs});fit.to_csv(O/f'{v}_saxs_fit.csv',index=False)
 geom=pd.read_csv(R/f'outputs/IDPForge_2026-10-13/{v}_observables.csv');lookup=geom.set_index('file');geom=lookup.loc[[Path(m['file']).name for m in members]]
 metrics={}
 for label,ww in [('uniform',weights['uniform']),('Kish_MaxEnt',w)]:
  metrics[label]={col:float(ww@geom[col]) for col in geom.columns}
  metrics[label]['SAXS_chi2_mean']=float(np.mean(((np.einsum("i,ij->j",ww,saxs)-iy)/sy)**2));metrics[label]['Kish_ratio']=float(1/(n*np.dot(ww,ww)));metrics[label]['Kish_Neff']=float(1/np.dot(ww,ww))
  for atom in np.unique(groups):metrics[label]['RMSE_'+atom+'_ppm']=float(np.sqrt(np.mean((np.einsum("ij,j->i",Xcs[groups==atom],ww)-ycs[groups==atom])**2)))
  comparisons.append(dict(variant=v,method='IDPForge '+label,**metrics[label]))
 report={'variant':v,'method':'Borthakur per-datatype and global Kish-calibrated maximum entropy','source_objective_sha256':hashlib.sha256(code.encode()).hexdigest(),'Kish_min':.1,'N':n,'single_datatype_sigma_reg':theta,'global_factor':global_t,'final_Kish_ratio':k,'final_Neff':float(1/np.dot(w,w)),'largest_weight':float(w.max()),'ensemble_error':'IID standard error over independently generated conformers, replacing MD time-block error','metrics':metrics,'independent_validation':False}
 (O/f'{v}_fit_summary.json').write_text(json.dumps(report,indent=2));summaries.append(report)
 fig,ax=plt.subplots(2,2,figsize=(11,8),layout='constrained')
 ax[0,0].errorbar(q,iy,yerr=sy,fmt='.',alpha=.35,label='Experimental');ax[0,0].plot(q,fit.I_uniform,label='IDPForge uniform');ax[0,0].plot(q,fit.I_MaxEnt,label='Kish MaxEnt');ax[0,0].set(yscale='log',xlabel='q (Å⁻¹)',ylabel='Intensity',title=v);ax[0,0].legend()
 ax[1,0].plot(q,(fit.I_uniform-iy)/sy,label='Uniform');ax[1,0].plot(q,(fit.I_MaxEnt-iy)/sy,label='MaxEnt');ax[1,0].axhline(0,color='gray');ax[1,0].set(xlabel='q (Å⁻¹)',ylabel='SAXS residual / experimental σ')
 ax[0,1].hist(geom.Rg_CA_A,bins=25,weights=weights['uniform'],histtype='step',label='Uniform');ax[0,1].hist(geom.Rg_CA_A,bins=25,weights=w,histtype='step',label='MaxEnt');ax[0,1].set(xlabel='CA Rg (Å)',ylabel='Weighted probability');ax[0,1].legend()
 ax[1,1].plot(np.sort(w)[::-1],'.');ax[1,1].set(xlabel='Weight rank',ylabel='Weight',title=f'K={k:.3f}, Neff={1/np.dot(w,w):.1f}')
 fig.savefig(O/f'{v}_fit.png',dpi=180);plt.close(fig)
 print(v,json.dumps(report),flush=True)
# Compare the archived ASTEROIDS outputs on exactly overlapping measured nuclei/residues.
for v in ['WT','AA']:
 new=pd.read_csv(O/f'{v}_shift_fit.csv')
 for p in (R/'work/asteroids_reference/data').glob(f'{v}_*_shifts.csv'):
  old=pd.read_csv(p);old.atom=old.atom.replace({'NH':'H','HN':'H'});joined=new.merge(old,on=['AR_residue','atom'])
  if not len(joined):continue
  metric={'variant':v,'method':'ASTEROIDS '+p.stem[len(v)+1:].removesuffix('_shifts'),'common_shift_points':len(joined)}
  for atom,g in joined.groupby('atom'):
   metric['RMSE_'+atom+'_ppm']=float(np.sqrt(np.mean((g.selected_secondary_ppm-g.experimental_secondary_ppm_y)**2)))
  saxsfile=p.with_name(p.name.replace('_shifts.csv','_saxs.csv'));summaryfile=p.with_name(p.name.replace('_shifts.csv','_summary.json'))
  if saxsfile.exists():
   f=pd.read_csv(saxsfile);metric['SAXS_chi2_mean']=float(np.mean(((f.I_selected-f.I_exp)/f.sigma_exp)**2))
  hfile=p.with_name(p.name.replace('_shifts.csv','_helicity.csv'))
  if hfile.exists():
   h=pd.read_csv(hfile);h=h[h.definition=='DSSP_H'];metric['DSSP_alpha_fraction']=float(h.selected.mean())
  structural=pd.read_csv(R/'work/asteroids_reference/ensemble_summary.csv');key=p.stem.removesuffix('_shifts');match=structural[(structural.ensemble==key)&(structural.view=='full')]
  if len(match):metric['Rg_CA_A']=float(match.iloc[0].Rg_mean_nm*10)
  if summaryfile.exists():
   s=json.load(open(summaryfile));metric['selected_count']=s.get('n_selected')
  comparisons.append(metric)
pd.DataFrame(comparisons).to_csv(O/'ensemble_comparison.csv',index=False)
(O/'fit_complete.json').write_text(json.dumps({'status':'complete','results':summaries},indent=2))

report=['# IDPForge Kish-restrained MaxEnt and ASTEROIDS comparison','', 'Uses the verbatim Borthakur repository dual objective, single-datatype calibration followed by global calibration, minimum Kish ratio 0.10. No ASTEROIDS selection was run.','', 'The prior is 200 validated IDPForge conformers per variant. Sampling-error regularization is estimated with IID standard errors, because these are independently generated conformers rather than a chronological MD trajectory. This is an explicit adaptation of the paper’s MD block-error calculation. No fixed SPARTA model-error tolerances from the earlier exploratory fit are imposed.','', 'Native SPARTA+ and Pepsi-SAXS predictions were used. Pepsi parameters follow the Borthakur calc_saxs.py values (dro=3.34, r0 min/max factor=1.68). Per-conformer scale fits reproduce the native fixed-parameter experimental-profile calculation; its intensity profile is used in the linear MaxEnt matrix. No additive background was fitted.','', 'The historical ASTEROIDS ensembles used different SAXS forward-model parameters and optimization scales. Comparisons are descriptive; no independent validation or matched rerun of the historical ensembles is claimed. Chemical-shift comparisons use overlapping residue/nucleus pairs.','', 'The regularization scan, retained/excluded observations, weights, predictions, effective sizes and fit summaries are saved as CSV/JSON/NPZ. With K=0.10 and N=200, the final effective size may be only 20; this does not establish adequate ensemble support.','', '![WT fit](WT_fit.png)','', '![AA fit](AA_fit.png)','', 'See ensemble_comparison.csv for the raw/MaxEnt and archived ASTEROIDS metrics.']
(O/'README.md').write_text('\n'.join(report)+'\n')
