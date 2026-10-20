from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd,mdtraj as md
O=Path(__file__).resolve().parent;R=O.parents[1];A=O/'comparison_analysis';A.mkdir(exist_ok=True)
records=[];cis_summary=[]
for method,base in [('harmonic',O),('Jeffreys',R/'outputs/bAIes_8A_seqsep2_2us_2026-10-14')]:
 for v in ['WT','AA']:
  d=base/v;p=np.loadtxt(d/'baies_params.dat');assert np.isfinite(p).all() and (p[:,4]>0).all()
  if method=='harmonic':assert np.array_equal(p,np.loadtxt(R/'outputs/bAIes_8A_seqsep2_2us_2026-10-14'/v/'baies_params.dat'))
  times=[];states=[];count=0;maxbond=0
  for t in md.iterload(str(d/'production.xtc'),top=str(d/'protein.pdb'),chunk=500):
   assert np.allclose(t.time,np.arange(count+1,count+len(t)+1)*100)
   ix,om=md.compute_omega(t,periodic=False);states.append(np.where(abs(om)<np.pi/6,1,np.where(abs(om)>5*np.pi/6,0,-1)))
   if method=='harmonic':
    bonds=np.array([(a.index,b.index) for a,b in t.top.bonds if a.element.symbol!='H' and b.element.symbol!='H']);maxbond=max(maxbond,float(md.compute_distances(t,bonds,periodic=False).max()))
   times.extend(t.time/1000);count+=len(t)
  assert count==20000
  x=np.concatenate(states);ars=np.array([t.top.atom(int(i[2])).residue.index+328 for i in ix]);pro=np.array([t.top.atom(int(i[2])).residue.name=='PRO' for i in ix]);poly=(ars>=372)&(ars<=379)
  out=A/f'{method}_{v}';np.savez_compressed(str(out)+'_cis_states.npz',states=x,AR_following=ars,is_proline=pro,time_ns=times)
  pd.DataFrame({'production_ns':times,'cis_all':(x==1).sum(1),'cis_XPro':(x[:,pro]==1).sum(1),'cis_Pro372379':(x[:,poly]==1).sum(1)}).to_csv(str(out)+'_cis_counts.csv',index=False)
  for k,ar in enumerate(ars):
   z=x[:,k];cis_summary.append({'method':method,'variant':v,'following_residue':'tag_Pro' if ar==329 else str(ar),'is_proline':bool(pro[k]),'cis_percent':float((z==1).mean()*100),'other_percent':float((z==-1).mean()*100),'observed_cis_entry_exit_count':int(np.sum((z[1:]==1)!=(z[:-1]==1)))})
  if method=='harmonic':
   lines=(d/'production.log').read_text().splitlines();assert lines[-1]=='Completed 1000000000 steps';rows=[list(map(float,l.split('\t')[:4])) for l in lines if l and l[0].isdigit() and '\t' in l];z=pd.DataFrame(rows,columns=['step','potential_kJmol','temperature_K','speed_ns_day']);assert len(z)==20000 and np.isfinite(z.to_numpy()).all();assert np.array_equal(z.step,10000000+np.arange(1,20001)*50000)
   z['production_ns']=(z.step-10000000)*.000002;z.to_csv(A/f'{v}_energy_temperature.csv',index=False);z.groupby(np.arange(len(z))//2000).mean().to_csv(A/f'{v}_energy_temperature_blocks.csv')
   records.append({'variant':v,'frames':count,'initial_step':10000000,'final_step':int(z.step.iloc[-1]),'production_ns':2000,'absolute_end_ns':2020,'maximum_heavy_bond_nm':maxbond,'temperature_mean_K':z.temperature_K.mean(),'temperature_SD_K':z.temperature_K.std(),'potential_mean_kJmol':z.potential_kJmol.mean(),'final_state_exists':(d/'production_final.xml').exists(),'same_restraint_parameters_as_Jeffreys':True,'hashes':{f:hashlib.sha256((d/f).read_bytes()).hexdigest() for f in ['production.xtc','protein.pdb','baies_params.dat','production_final.xml']}})
  print(method,v,'validated and cis complete',flush=True)
pd.DataFrame(cis_summary).to_csv(A/'cis_residue_summary.csv',index=False);(A/'production_validation.json').write_text(json.dumps(records,indent=2))
