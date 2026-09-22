"""Recover and audit actual EFA component scattering profiles, without assigning species."""
from inspect_data import *
from run_efa import guinier
from independent_baseline import fixed_guinier
import raw_efa_core as efa
import json, itertools
OUT=ROOT/'component_results';OUT.mkdir(exist_ok=True)
SPECS={'WT_092025':dict(roi=[120,205],avg=[159,175],anchors=[[55,94],[210,234]],ends=[165,175,185,195],starts=[130,140,150],fid=[175,140]),
       'AA_052026':dict(roi=[145,215],avg=[169,188],anchors=[[35,74],[225,254]],ends=[175,185,195,205],starts=[150,160,170],fid=[185,160])}

def transform(name,kind='linear',anchors=None):
 q,Y,E=load(name);spec=SPECS[name];a,b=spec['roi'];frames=np.arange(a,b+1);anchors=anchors or spec['anchors'];t=np.array([(i+j)/2 for i,j in anchors])
 f=np.clip((frames-t[0])/(t[1]-t[0]),0,1)
 if kind=='pre':f[:]=0
 if kind=='post':f[:]=1
 A=np.eye(Y.shape[1])[:,frames]
 for h,(lo,hi) in enumerate(anchors):A[lo:hi+1,:]-=((1-f) if h==0 else f)[None,:]/(hi-lo+1)
 I=Y@A;err=np.sqrt((E**2)@(A**2))
 return q,Y,E,frames,A,I,err

def solve(q,I,err,frames,ends,starts,qmin=.012,qmax=.25):
 m=(q>=qmin)&(q<=qmax);scale=err[m].mean(1);D=I[m]/scale[:,None];U,s,Vt=np.linalg.svd(D,full_matrices=False)
 supports=np.array([[frames[0],ends],[starts,frames[-1]]]);rel=supports-frames[0]
 with np.errstate(all='ignore'):
  ok,conv,r=efa.runRotation(D,I[m],err[m],rel,[True,True],Vt.T,niter=3000,tol=1e-9)
 info=dict(converged=bool(ok),iterations=conv['iterations'],supports=supports.tolist(),q=[qmin,qmax],singular_values=s[:4].tolist())
 if not ok:return info,None
 C=r['M']*r['C'];L=np.linalg.pinv(C.T);S=I@L;res=(I[m]-S[m]@C.T)/err[m]
 info.update(residual_ms=float(np.mean(res**2)),rank1_residual_ms=float(np.mean(((D-np.outer(U[:,0]*s[0],Vt[0]))*scale[:,None]/err[m])**2)),condition=float(np.linalg.cond(C)),peaks=(frames[np.argmax(C,axis=0)]).tolist(),negative_fraction=np.mean(S[m]<0,axis=0).tolist())
 return info,dict(C=C,S=S,L=L,residual=res,mask=m,U=U,s=s,Vt=Vt)

def evaluate(name,kind,end,start,qmax=.25,anchors=None):
 q,Y,E,frames,A,I,err=transform(name,kind,anchors);info,out=solve(q,I,err,frames,end,start,qmax=qmax)
 info.update(sample=name,baseline=kind)
 if out is None:return info,None
 # Exact propagation through fixed linear background and FIXED estimated elution C.
 # Bootstrap below separately accounts for uncertainty in C.
 W=A@out['L'];SE=np.sqrt((E**2)@(W**2));out.update(q=q,error=SE,A=A,I=I,frame_error=err,frames=frames)
 info['main_guinier']=fixed_guinier(q,out['S'][:,1],SE[:,1])
 info['main_guinier_adaptive']=guinier(q,out['S'][:,1],SE[:,1])
 m=out['mask'];info['strong_negative_fraction'] = np.mean(out['S'][m]<-3*SE[m],axis=0).tolist()
 # Not biological validity: only a numerical/physical screening rule for the main curve.
 info['main_screen']=bool(info['condition']<30 and info['main_guinier'].get('Rg') is not None and info['strong_negative_fraction'][1]==0 and info['residual_ms']<=1.1*info['rank1_residual_ms'])
 return info,out

if __name__=='__main__':
 rows=[]
 for name,spec in SPECS.items():
  for kind,end,start,qmax in itertools.product(['linear','pre','post'],spec['ends'],spec['starts'],[.15,.25]):
   r,_=evaluate(name,kind,end,start,qmax);rows.append(r)
  (OUT/'rotation_grid.json').write_text(json.dumps(rows,indent=2))
  primary,out=evaluate(name,'linear',*spec['fid']);print('FIDUCIAL',json.dumps(primary),flush=True)
  good=[r for r in rows if r['sample']==name and r.get('main_screen')];print('SCREENED',name,len(good),'of',sum(r['sample']==name for r in rows),'Rg range', [min(r['main_guinier']['Rg'] for r in good),max(r['main_guinier']['Rg'] for r in good)] if good else [],flush=True)
  if out:
   np.savez_compressed(OUT/f'{name}_fiducial.npz',**out)
  (OUT/f'{name}_fiducial.json').write_text(json.dumps(primary,indent=2))
