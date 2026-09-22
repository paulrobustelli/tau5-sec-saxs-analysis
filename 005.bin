from component_recovery import *
import sys
sys.path.insert(0,str(ROOT/'bioxtasraw-source'))
from bioxtasraw import BIFT

def normalized(q,I):
 m=(q>=.06)&(q<=.10);s=I[m].mean()
 return I/s if s>0 else np.full_like(I,np.nan)

def audit(name,bootstrap=100,null_reps=200):
 spec=SPECS[name];q,Y,E,frames,A,I,err=transform(name);info,out=evaluate(name,'linear',*spec['fid'])
 assert out is not None
 av=(frames>=spec['avg'][0])&(frames<=spec['avg'][1]);avw=A[:,av].mean(1);avg=Y@avw;avge=np.sqrt((E**2)@(avw**2))
 cscale=out['C'][av].mean(0);components=out['S']*cscale[None,:];ce=out['error']*cscale[None,:]
 rng=np.random.default_rng(20260922);boot=[];rgs=[];bootavg=[];diffs=[];fails=0
 for j in range(bootstrap):
  Yn=Y+rng.normal(size=Y.shape)*E;In=Yn@A;ri,oi=solve(q,In,err,frames,*spec['fid'])
  if oi is None:fails+=1;continue
  curve=oi['S'][:,1]*oi['C'][av,1].mean();boot.append(curve);bootavg.append(In[:,av].mean(1));f=fixed_guinier(q,curve,ce[:,1]);fa=fixed_guinier(q,In[:,av].mean(1),avge)
  if f['Rg'] is not None and fa['Rg'] is not None:rgs.append(f['Rg']);diffs.append(f['Rg']-fa['Rg'])
 boot=np.array(boot);bootavg=np.array(bootavg);bs=np.std(boot,axis=0,ddof=1)
 rgs=[];diffs=[]
 for curve,avgdraw in zip(boot,bootavg):
  f=fixed_guinier(q,curve,bs);fa=fixed_guinier(q,avgdraw,avge)
  if f['Rg'] is not None and fa['Rg'] is not None:rgs.append(f['Rg']);diffs.append(f['Rg']-fa['Rg'])
 m=out['mask'];D=I[m]/err[m].mean(1)[:,None];u,s,vt=np.linalg.svd(D,full_matrices=False);rankone=np.outer(u[:,0]*s[0],vt[0]);null=[]
 for _ in range(null_reps):
  noise=(rng.normal(size=Y.shape)*E)@A;ds=rankone+noise[m]/err[m].mean(1)[:,None];null.append(np.linalg.svd(ds,compute_uv=False)[1])
 info.update(bootstrap_requested=bootstrap,bootstrap_converged=len(boot),bootstrap_failed=fails,bootstrap_Rg_95=np.quantile(rgs,[.025,.975]).tolist(),bootstrap_delta_Rg_vs_average_95=np.quantile(diffs,[.025,.975]).tolist(),average_guinier=fixed_guinier(q,avg,avge),candidate_guinier=fixed_guinier(q,components[:,1],bs),candidate_guinier_fixedC=fixed_guinier(q,components[:,1],ce[:,1]),null_s2_99=float(np.quantile(null,.99)),s2_p_rank1=float((1+sum(v>=s[1] for v in null))/(null_reps+1)),bootstrap_note='Parametric perturbations of supplied errors; residual-baseline covariance and refitted C included. Conditional on anchors/supports/model. Source facility covariance unknown.')
 # A rank-one estimate is exported separately; it is NOT a purified component.
 cc=vt[0].copy();cc*=1 if cc.sum()>0 else -1;cc/=cc.sum();ll=cc/(cc@cc);rank1=I@ll*cc[av].mean();rw=A@ll*cc[av].mean();rank1e=np.sqrt((E**2)@(rw**2))
 np.savetxt(OUT/f'{name}_EFA_candidate_main.dat',np.c_[q,components[:,1],bs],header='q_A^-1 I_candidate_main bootstrap_sigma; CANDIDATE ONLY, not validated monomer; scaled to contribution in selected average')
 np.savetxt(OUT/f'{name}_frame_average_same_baseline.dat',np.c_[q,avg,avge],header='q_A^-1 I_mean sigma_with_shared_baseline; same input and independent linear baseline as candidate')
 np.savetxt(OUT/f'{name}_EFA_other_component.dat',np.c_[q,components[:,0],ce[:,0]],header='q_A^-1 I_other_component conditional_fixed_C_sigma; exploratory, no species assignment')
 np.savetxt(OUT/f'{name}_rank1_dominant_not_purified.dat',np.c_[q,rank1,rank1e],header='q_A^-1 I_rank1 conditional_sigma; rank-one dominant signal NOT a cleaned monomer')
 np.savetxt(OUT/f'{name}_EFA_elution.csv',np.c_[frames,out['C']],delimiter=',',header='frame,other_coefficient,main_candidate_coefficient',comments='')
 np.savez_compressed(OUT/f'{name}_audit.npz',q=q,average=avg,average_error=avge,candidate=components[:,1],candidate_fixedC_error=ce[:,1],candidate_bootstrap_error=bs,other=components[:,0],other_error=ce[:,0],bootstrap=boot,bootstrap_average=bootavg,rank1=rank1,rank1_error=rank1e,frames=frames,C=out['C'],residual=out['residual'],efa_forward=efa.runEFA(D)[:4],efa_backward=efa.runEFA(D,False)[:4,::-1])
 (OUT/f'{name}_audit.json').write_text(json.dumps(info,indent=2));print('AUDIT',name,json.dumps(info),flush=True)
 # Conditional BIFT visualization of each mean, with the SAME settings.
 for tag,curve,error in [('candidate',components[:,1],bs),('average',avg,avge)]:
  sel=(q>=.012)&(q<=.25);np.random.seed(20260922)
  obj=BIFT.doBift(q[sel],curve[sel],error[sel],tag,40,150,1e10,12,50,220,15,100,single_proc=True,nprocs=1)
  if obj is None:continue
  pars=obj.getAllParameters();pars={k:float(v) if isinstance(v,np.floating) else v for k,v in pars.items()}
  pars['warning']='Conditional BIFT display; candidate is not a validated monomer. Bootstrap diagonal errors omit q covariance; Dmax is model dependent.'
  (OUT/f'{name}_{tag}_bift.json').write_text(json.dumps(pars,indent=2))
  np.savez_compressed(OUT/f'{name}_{tag}_bift.npz',r=obj.r,p=obj.p,p_error=obj.err,q=obj.q_orig,fit=obj.i_fit)
  np.savetxt(OUT/f'{name}_{tag}_Pr.dat',np.c_[obj.r,obj.p,obj.err],header='r_A P(r) conditional_BIFT_error; exploratory component/average comparison')
 return info

def validate():
 # Exact two-species recovery in a known identifiable synthetic case.
 q=np.linspace(.012,.25,100);t=np.arange(86)+120
 c0=np.exp(-.5*((t-139)/8)**2);c0[t>160]=0
 c1=np.exp(-.5*((t-176)/10)**2);c1[t<150]=0
 S=np.c_[np.exp(-q*q*45**2/3),np.exp(-q*q*25**2/3)];C=np.c_[c0,c1];I=S@C.T;E=np.ones_like(I)*.001
 r,o=solve(q,I,E,t,160,150)
 assert r['converged'];recon=o['S']@o['C'].T
 assert np.max(abs(recon-I))<1e-7
 assert np.max(abs(normalized(q,o['S'][:,1])-normalized(q,S[:,1])))<1e-6
 # Baseline transform recovers a known linear baseline and analytic variance.
 q,Y,E,t,A,_,_=transform('WT_092025');x=np.arange(450);line=2+.01*x
 assert np.max(abs(line@A))<1e-12
 rng=np.random.default_rng(73);v=A[:,20];draws=rng.normal(size=(20000,450))@v
 assert abs(np.var(draws)/(v@v)-1)<.04
 return {'known_two_component_curve_recovered':True,'linear_background_removed_exactly':True,'shared_background_variance_MC_check':True}
if __name__=='__main__':
 (OUT/'validation.json').write_text(json.dumps(validate(),indent=2))
 for name in SPECS:audit(name)
