from select_windows import *
import sys,threading
os.environ.setdefault('NUMBA_CACHE_DIR',str(ROOT/'.numba_cache'))
sys.path.insert(0,str(ROOT/'bioxtasraw-source'))
from bioxtasraw import BIFT

def fit_bift(name,a,b,qmin=.012,qmax=.25,suffix='primary',mc=150):
 q,I,E=average(name,a,b);m=(q>=qmin)&(q<=qmax)
 np.random.seed(20260921)
 obj=BIFT.doBift(q[m],I[m],E[m],name,40,150,1e10,12,50,220,15,mc,single_proc=True,nprocs=1)
 if obj is None:raise RuntimeError('BIFT failed: '+name)
 par=obj.getAllParameters();par={k:(float(v) if isinstance(v,np.floating) else v) for k,v in par.items()};par.update(sample=name,frames=[a,b],uncertainty='RAW BIFT local alpha/Dmax sampling; excludes background/window/q-range systematic uncertainty')
 np.savez_compressed(PLOT/f'{name}_bift_{suffix}.npz',r=obj.r,p=obj.p,p_error=obj.err,q=obj.q_orig,I=obj.i_orig,error=obj.err_orig,fit=obj.i_fit)
 (PLOT/f'{name}_bift_{suffix}.json').write_text(json.dumps(par,indent=2))
 np.savetxt(PLOT/f'{name}_pr_{suffix}.dat',np.c_[obj.r,obj.p,obj.err],header='r_A P(r) conditional_RAW_BIFT_error; I(q)=4pi integral P(r) sinc(qr) dr')
 print('BIFT',name,suffix,par,flush=True);return par
if __name__=='__main__':
 sels=json.loads((PLOT/'selected_windows.json').read_text());rows=[]
 for s in sels:rows.append(fit_bift(s['sample'],s['start'],s['end']))
 for s in sels:
  for qm,qx,suff in [(.010,.25,'lowq'),(.015,.25,'trimlowq'),(.012,.20,'qmax020'),(.012,.30,'qmax030')]:rows.append(fit_bift(s['sample'],s['start'],s['end'],qm,qx,suff,mc=80))
 (PLOT/'bift_sensitivity.json').write_text(json.dumps(rows,indent=2))
