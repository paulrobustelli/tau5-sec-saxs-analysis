"""Collect validated native predictions; per-datatype then globalKish0.10 fit."""
from pathlib import Path
import json,hashlib,sys
import numpy as np,pandas as pd
from kish_engine import sampling_error,calibrate,dual
R=Path(__file__).resolve().parents[2];O=Path(__file__).resolve().parent;name=sys.argv[1];D=O/name
assert (O/'engine_validation.json').exists();assert not (D/'fit_complete.json').exists()
manifest=json.loads((D/'manifest.json').read_text());variant=name.rsplit('_',1)[1];n=len(manifest);tables=[];native=[]
seq=json.loads((R/'outputs/IDPForge_2026-10-13/sequences.json').read_text())[variant]
for m in manifest:
 j=Path(m['job']);assert (j/'done').exists() and (j/'done').read_text().strip()=='0',f'SAXS pending: {m["id"]}'
 assert hashlib.sha256((j/'input.pdb').read_bytes()).hexdigest()==m['predictor_input_sha256']
 sj=Path(m.get('source_job',m['job']));flag='sparta_done' if name.startswith(('REST','bAIes')) else 'done';assert (sj/flag).read_text().strip()=='0'
 lines=(sj/'pred.tab').read_text().splitlines();s=''.join(''.join(l.split()[2:]) for l in lines if l.startswith('DATA SEQUENCE'));assert s==seq+('C' if name.startswith('REST') else '')
 table={}
 for l in lines:
  f=l.split()
  if len(f)==9 and f[0].isdigit():table[(int(f[0]),'H' if f[2] in ['HN','NH'] else f[2])]=float(f[4])
 assert np.isfinite(list(table.values())).all();tables.append(table)
 x=np.loadtxt(j/'prediction.dat');assert x.shape==(601,8) and np.isfinite(x).all() and np.all(x[:,1]>0) and np.allclose(x[:,0],np.linspace(0,.6,601));native.append(x[:,1])
expdir=R/'outputs/IDPForge_D2D_MaxEnt_2026-10-15/experimental';obs=pd.read_csv(expdir/'chemical_shifts.csv');obs=obs[obs.variant==variant].copy();keep=obs.apply(lambda r:all((r.model_residue,r.atom) in t for t in tables),axis=1);obs[~keep].to_csv(D/'excluded_observations.csv',index=False);obs=obs[keep]
cs=np.array([[t[(r.model_residue,r.atom)] for t in tables] for r in obs.itertuples()]);yc=obs.value.to_numpy();groups=obs.atom.to_numpy();q,iy,sy=np.loadtxt(expdir/f'{variant}_SAXS.dat',unpack=True);raw=np.array([np.interp(q,np.linspace(0,.6,601),z) for z in native]).T
uniform=raw.mean(1);s0=float(np.sum(uniform*iy/sy**2)/np.sum((uniform/sy)**2));X=np.r_[cs,raw*s0];y=np.r_[yc,iy];types=np.r_[groups,np.repeat('SAXS',len(q))];saxsmask=types=='SAXS'
kind='IID' if name.startswith('IDPForge') else 'MDblocks';times=None if kind=='IID' else np.array([m.get('production_ns',m.get('time_ns')) for m in manifest]);sem,blocks=sampling_error(X,kind,times)
np.savez_compressed(D/'observable_matrix.npz',predicted_reference_scale=X,observed=y,datatype=types,sigma_ensemble=sem,raw_saxs=raw,q=q,intensity_sigma=sy,initial_global_scale=s0,chemical_shifts=cs,time_ns=np.array([]) if times is None else times)
if blocks:np.savez_compressed(D/'sampling_error_block_scan.npz',**blocks)
scan=[];base=np.ones(len(y));theta={};weights={'uniform':np.ones(n)/n}
for group in np.unique(types):
 mask=types==group;kwargs={} if group!='SAXS' else {'saxs_mask':np.ones(mask.sum(),bool),'raw_saxs':raw,'initial_scale':s0}
 tt,result=calibrate(X[mask],y[mask],np.ones(mask.sum()),sem[mask],group,scan,**kwargs);base[mask]=tt;theta[group]=tt;weights[group+'_only']=result['weights'];pd.DataFrame(scan).to_csv(D/'regularization_scan_in_progress.csv',index=False);print(name,group,tt,result['Kish_ratio'],flush=True)
global_t,r=calibrate(X,y,base,sem,'combined',scan,saxs_mask=saxsmask,raw_saxs=raw,initial_scale=s0);w=r['weights'];assert r['Kish_ratio']>=.1;weights['combined']=w
# Independently validate the real-data optimum, not only a synthetic case.
actualX=np.r_[cs,raw*r['scale']];actualmean=np.einsum('ij,j->i',actualX,w);var=(base*global_t)**2+sem**2
primal=float(np.einsum('i,i->',w,np.log(np.maximum(w,1e-300)*n))+.5*np.sum((actualmean-y)**2/var))
assert np.isclose(primal,-r['dual'],rtol=1e-4,atol=1e-5),(primal,r['dual'])
assert r['max_gradient']<1e-3 and abs(r['logscale_primal_gradient'])<1e-3,r
rng=np.random.default_rng(1729);gradchecks=[];args=(r['centered_normalized'],np.zeros(len(y)),np.ones(n),base*global_t/r['norm'],sem/r['norm'])
z=r['lambda']+rng.normal(0,.0001,len(y));ff,gg=dual(z,*args)
for ep in [1e-5,1e-6,1e-7]:
 direction=rng.normal(size=len(y));direction/=np.linalg.norm(direction);fd=(dual(z+ep*direction,*args)[0]-dual(z-ep*direction,*args)[0])/(2*ep);analytic=float(np.einsum('i,i->',gg,direction));assert np.isclose(fd,analytic,rtol=1e-4,atol=1e-5);gradchecks.append({'epsilon':ep,'finite_difference':float(fd),'analytic':analytic})
(D/'fit_validation.json').write_text(json.dumps({'primal':primal,'negative_dual':-r['dual'],'gradient_checks':gradchecks,'Kish_ratio':r['Kish_ratio'],'logscale_primal_gradient':r['logscale_primal_gradient'],'finite_weights':bool(np.isfinite(w).all())},indent=2))
pd.DataFrame({'id':[m['id'] for m in manifest],'file':[m['file'] for m in manifest],**weights}).to_csv(D/'weights.csv',index=False);pd.DataFrame(scan).to_csv(D/'regularization_scan.csv',index=False)
mean=np.einsum('ij,j->i',raw,w)*r['scale'];pd.DataFrame({'q_Ainv':q,'I_exp':iy,'sigma_exp':sy,'I_uniform':uniform*s0,'I_MaxEnt':mean}).to_csv(D/'saxs_fit.csv',index=False);obs['unweighted_ppm']=cs.mean(1);obs['reweighted_ppm']=np.einsum('ij,j->i',cs,w);obs.to_csv(D/'chemical_shift_fit.csv',index=False)
metrics={}
z=np.load(D/'structural_cache.npz');assert list(z['files'])==[m['file'] for m in manifest]
for label,ww,scale in [('uniform',weights['uniform'],s0),('MaxEnt',w,r['scale'])]:
 met={col:float(np.einsum('i,i->',ww,z['values'][:,i])) for i,col in enumerate(z['columns'])};met['Kish_ratio']=float(1/(n*np.einsum('i,i->',ww,ww)));met['Neff']=float(1/np.einsum('i,i->',ww,ww));met['SAXS_chi2_mean']=float(np.mean(((np.einsum('ij,j->i',raw,ww)*scale-iy)/sy)**2))
 for atom in np.unique(groups):
  m=groups==atom;met['RMSE_'+atom+'_ppm']=float(np.sqrt(np.mean((np.einsum('ij,j->i',cs[m],ww)-yc[m])**2)))
 for reg,(lo,hi) in {'R2_main':(390,404),'R2_minor':(404,410),'R3':(430,440)}.items():met[reg+'_DSSP_H']=float(np.einsum('i,i->',ww,z['helix'][:,lo-328:hi-327].mean(1)))
 metrics[label]=met
report={'ensemble':name,'N':n,'datatype_theta':theta,'global_theta':global_t,'scale_uniform':s0,'scale_MaxEnt':r['scale'],'Kish_ratio':r['Kish_ratio'],'max_dual_gradient':r['max_gradient'],'logscale_primal_gradient':r['logscale_primal_gradient'],'metrics':metrics,'sampling_error':kind+';100ns blocks for correlated ensembles, fixed uniform-prior sampling error during calibration','CS_reuse':'Original120res predictions retained for terminallyextendedIDPForge/bAIes by explicituserinstruction; comparisonobservationscommonAR330446','SAXS':'IdenticalcanonicalPepsi3.0dro5defaultFraserradius,oneglobalensemblescale/nooffset,noindividualscales','independent_validation':False}
(D/'fit_complete.json').write_text(json.dumps(report,indent=2));print(json.dumps(report),flush=True)
