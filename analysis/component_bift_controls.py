from audit_component_recovery import *
def controls():
 rows=[]
 for name in SPECS:
  d=np.load(OUT/f'{name}_audit.npz');q=d['q'];m=(q>=.012)&(q<=.25)
  for tag,curve,error in [('average_matched_candidate_errors',d['average'],d['candidate_bootstrap_error']),('candidate_matched_average_errors',d['candidate'],d['average_error'])]:
   np.random.seed(20260922);o=BIFT.doBift(q[m],curve[m],error[m],tag,40,150,1e10,12,50,220,15,100,single_proc=True,nprocs=1)
   p=o.getAllParameters();r=dict(sample=name,control=tag,Rg=float(p['rg']),Dmax=float(p['dmax']),chi2=float(p['chisq']));rows.append(r)
   np.savez_compressed(OUT/f'{name}_{tag}_bift.npz',r=o.r,p=o.p,p_error=o.err)
 (OUT/'bift_weight_controls.json').write_text(json.dumps(rows,indent=2));print(rows)
if __name__=='__main__':controls()
