from pathlib import Path
import json,numpy as np
from scipy.special import softmax
from kish_engine import fit,dual
rng=np.random.default_rng(287);n=80;nc=5;ns=10
cs=rng.normal(size=(nc,n));raw=np.exp(-np.linspace(.05,2,ns)[:,None]*(1+.15*rng.normal(size=n)[None,:]));true=softmax(.8*cs[0]);y=np.r_[cs@true,2.3*np.einsum('ij,j->i',raw,true)];s0=2.;X=np.r_[cs,raw*s0];mask=np.arange(nc+ns)>=nc;sr=np.r_[np.ones(nc)*.2,np.ones(ns)*.05];sm=np.ones(nc+ns)*.01
r=fit(X,y,sr,sm,saxs_mask=mask,raw_saxs=raw,initial_scale=s0);A=r['centered_normalized'];norm=r['norm'];l=r['lambda'];w=r['weights'];z=l+rng.normal(size=len(l))*.03;direction=rng.normal(size=len(l));direction/=np.linalg.norm(direction);eps=1e-6
args=(A,np.zeros(len(y)),np.ones(n),sr/norm,sm/norm);f,g=dual(z,*args);fd=(dual(z+eps*direction,*args)[0]-dual(z-eps*direction,*args)[0])/(2*eps);assert np.isclose(fd,np.einsum('i,i->',g,direction),rtol=1e-5,atol=1e-7)
pred=np.r_[cs,np.asarray(raw)*r['scale']];mean=np.einsum('ij,j->i',pred,w);primal=np.einsum('i,i->',w,np.log(w*n))+.5*np.sum((mean-y)**2/(sr*sr+sm*sm));assert np.isclose(primal,-r['dual'],rtol=1e-5,atol=1e-6),(primal,r['dual'])
# At fixed lambda the profile envelope derivative equals lambda dot scaledSAXS / norm.
def logscale_primal_envelope(u):
 z=pred.copy();z[mask]=raw*r['scale']*np.exp(u);aa=(z-y[:,None])/norm[:,None];return -dual(l,aa,np.zeros(len(y)),np.ones(n),sr/norm,sm/norm)[0]
fdscale=(logscale_primal_envelope(eps)-logscale_primal_envelope(-eps))/(2*eps);assert np.isclose(fdscale,r['logscale_primal_gradient'],atol=1e-6)
assert abs(r['logscale_primal_gradient'])<1e-4
out={k:r[k] for k in ['scale','Kish_ratio','dual','max_gradient','logscale_primal_gradient']};out.update(primal=primal,dual_gradient_fd=fd,scale_gradient_fd=fdscale,checks='dual finite differences; sharedscale envelope finite differences; independent primal/dual equality; scale stationary')
Path(__file__).with_name('engine_validation.json').write_text(json.dumps(out,indent=2));print(out)
