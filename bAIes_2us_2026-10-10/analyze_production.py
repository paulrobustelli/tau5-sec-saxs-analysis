"""Chunked production analysis. Run with --root <restart directory>."""
import argparse,json,hashlib,os
from pathlib import Path
os.environ.setdefault('OMP_NUM_THREADS','4')
import numpy as np,pandas as pd,mdtraj as md
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.signal import fftconvolve
p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);a=p.parse_args();root=a.root;out=root/'Analysis';out.mkdir(exist_ok=True)
ref=md.load(str(out/'R2R3_all_helix.pdb'))
refs=np.array([ref.xyz[0,ref.top.select(f'residue {i} to {i+6} and backbone')] for i in range(1,51)],float)
assert refs.shape==(50,28,3)
def rmsd(x,y):
 x=x-x.mean(-2,keepdims=True);y=y-y.mean(-2,keepdims=True)
 u,s,v=np.linalg.svd(np.einsum('...ni,...nj->...ij',x,y));sg=np.where(np.linalg.det(u@v)<0,-1.,1.)
 return np.sqrt(np.maximum(0,(np.sum(x*x,(-2,-1))+np.sum(y*y,(-2,-1))-2*(s[...,0]+s[...,1]+sg*s[...,2]))/28))
def sa(bb):
 w=np.stack([bb[:,i:i+7].reshape(len(bb),28,3) for i in range(114)],1)
 rr=np.repeat(refs[:1],114,axis=0);rr[63:113]=refs
 q=(rmsd(w,rr)/.1)**4;sw=(1+q)/(1+q+q*q)
 return sw.sum(1),sw[:,63:113].sum(1)
def rg(x):return np.sqrt(np.mean(np.sum((x-x.mean(1,keepdims=True))**2,2),1))*10
summ=[];manifest=[];data={};maps={}
for variant,start in [('WT',4695000),('AA',4675000)]:
 top=md.load(str(root/variant/'protein.pdb')).topology
 assert top.n_residues==120
 bbidx=np.array([[next(x.index for x in r.atoms if x.name==n) for n in ['N','CA','C','O']] for r in top.residues]);ca=bbidx[:,1]
 pairs=np.array([(i,j) for i in range(120) for j in range(i+1,120)])
 contactsum=np.zeros(len(pairs));blockcontacts=[];helices=[];records=[];segments=[]
 # Preserve actual times; retain each timestamp once, preferring earlier segment.
 last=-np.inf
 for fname,kind in [('trajectory_partial.xtc','equilibration'),('trajectory_continuation.xtc','equilibration'),('production_2us.xtc','production_cpu'),('production_2us_gpu.xtc','production_gpu')]:
  path=root/variant/fname;count=0;first=None;end=None
  offset={'trajectory_partial.xtc':0,'trajectory_continuation.xtc':4670000 if variant=='WT' else 4650000,'production_2us.xtc':4650000,'production_2us_gpu.xtc':4750000}[fname]
  cutoff={'trajectory_partial.xtc':4670000 if variant=='WT' else 4650000,'trajectory_continuation.xtc':start,'production_2us.xtc':4750000,'production_2us_gpu.xtc':start+1000000000}[fname]
  for t in md.iterload(str(path),top=top,chunk=200):
   t.time=t.time.astype(float)+offset*.002
   t=t[t.time<=cutoff*.002+1e-3]
   if not np.isfinite(t.xyz).all():raise ValueError('nonfinite coordinates')
   keep=t.time>last+1e-3;t=t[keep]
   if not len(t):continue
   last=float(t.time[-1]);count+=len(t);first=float(t.time[0]) if first is None else first;end=last
   production=t.time>start*.002+1e-3
   if not production.any():continue
   t=t[production];bb=t.xyz[:,bbidx,:].astype(float)
   h=md.compute_dssp(t,simplified=False);alpha=h=='H';anyh=np.isin(h,['H','G','I']);helices.append(alpha)
   fullsa,cropsa=sa(bb)
   d,_=md.compute_contacts(t,contacts=pairs,scheme='closest-heavy',periodic=False);cs=d<1.2;contactsum+=cs.sum(0)
   # 20 consecutive ~100 ns production blocks; aggregate by actual time.
   bid=np.minimum(19,((t.time-start*.002)/100000).astype(int))
   for k in range(20):
    if len(blockcontacts)<=k:blockcontacts.append([np.zeros(len(pairs)),0])
    m=bid==k;blockcontacts[k][0]+=cs[m].sum(0);blockcontacts[k][1]+=int(m.sum())
   frame=pd.DataFrame({'time_ns':t.time/1000,'production_ns':(t.time-start*.002)/1000,'segment':kind,'Rg_CA_A':rg(t.xyz[:,ca]),'Rg_allatom_equalweight_A':md.compute_rg(t)*10,'Rg_R2R3_CA_A':rg(t.xyz[:,ca[63:119]]),'DSSP_alpha_fraction':alpha[:,2:].mean(1),'DSSP_any_helix_fraction':anyh[:,2:].mean(1),'Sa_full':fullsa,'Sa_AR391_446':cropsa,'R2_alpha':alpha[:,67:78].mean(1),'R3_alpha':alpha[:,102:113].mean(1)})
   records.append(frame)
  segments.append({'file':fname,'kind':kind,'unique_frames':count,'first_time_ps':first,'last_time_ps':end,'sha256':hashlib.file_digest(path.open('rb'),'sha256').hexdigest()})
  print(variant,fname,count,flush=True)
 df=pd.concat(records,ignore_index=True);assert df.time_ns.is_monotonic_increasing;assert not df.time_ns.duplicated().any()
 df.to_csv(out/f'{variant}_timeseries.csv',index=False);hh=np.concatenate(helices)
 pd.DataFrame({'construct_index':np.arange(1,121),'AR_number':[None,None]+list(range(330,448)),'alpha_propensity':hh.mean(0)}).to_csv(out/f'{variant}_helicity.csv',index=False)
 cmap=np.zeros((120,120));cmap[pairs[:,0],pairs[:,1]]=contactsum/len(df);cmap+=cmap.T;maps[variant]=cmap
 np.savetxt(out/f'{variant}_contact_map.csv',cmap,delimiter=',')
 bc=np.array([s/n if n else np.full_like(s,np.nan) for s,n in blockcontacts]);np.savez_compressed(out/f'{variant}_block_contacts.npz',probabilities=bc,counts=[n for s,n in blockcontacts],pairs=pairs)
 blocks=df.assign(block=np.minimum(19,(df.production_ns/100).astype(int))).groupby('block').mean(numeric_only=True);blocks.to_csv(out/f'{variant}_100ns_blocks.csv')
 for col in ['Rg_CA_A','Rg_allatom_equalweight_A','Rg_R2R3_CA_A','DSSP_alpha_fraction','DSSP_any_helix_fraction','Sa_full','Sa_AR391_446','R2_alpha','R3_alpha']:
  x=df[col].to_numpy();b=blocks[col].to_numpy();se=b.std(ddof=1)/np.sqrt(len(b));xc=x-x.mean();ac=fftconvolve(xc,xc[::-1],mode='full')[len(x)-1:];ac/=np.arange(len(x),0,-1);ac=ac/ac[0] if ac[0]>0 else ac
  cut=next((i for i in range(1,len(ac)) if ac[i]<=0),len(ac));g=max(1,1+2*ac[1:cut].sum())
  summ.append({'variant':variant,'observable':col,'frames':len(x),'mean':x.mean(),'SD':x.std(ddof=1),'block_SE':se,'block95_low':x.mean()-2.093*se,'block95_high':x.mean()+2.093*se,'first_half_mean':x[:len(x)//2].mean(),'second_half_mean':x[len(x)//2:].mean(),'approx_ESS':len(x)/g})
 data[variant]=df;manifest.append({'variant':variant,'production_start_step':start,'segments':segments,'analyzed_frames':len(df),'final_time_ns':float(df.time_ns.iloc[-1])})
 print(variant,'done',flush=True)
s=pd.DataFrame(summ);s.to_csv(out/'summary.csv',index=False)
json.dump(manifest,(out/'trajectory_provenance.json').open('w'),indent=2)
fig,axes=plt.subplots(3,2,figsize=(12,10),layout='constrained')
for v,c in [('WT','#2868ac'),('AA','#db782c')]:
 df=data[v]
 for row,col in enumerate(['Rg_CA_A','DSSP_alpha_fraction','Sa_AR391_446']):
  axes[row,0].plot(df.production_ns,df[col].rolling(100,min_periods=1).mean(),label=v,color=c);axes[row,0].set(xlabel='Production time (ns)',ylabel=col)
  axes[row,1].hist(df[col],bins=60,density=True,histtype='step',label=v,color=c);axes[row,1].set(xlabel=col,ylabel='Density')
for ax in axes.flat:ax.legend()
fig.savefig(out/'observables.png',dpi=180);fig.savefig(out/'observables.pdf');plt.close(fig)
fig,ax=plt.subplots(figsize=(10,4),layout='constrained')
for v in data:
 h=pd.read_csv(out/f'{v}_helicity.csv');ax.plot(h.AR_number,h.alpha_propensity,label=v)
ax.set(xlabel='AR residue',ylabel='DSSP alpha propensity',ylim=(0,1));ax.legend();fig.savefig(out/'helicity.png',dpi=180);plt.close(fig)
fig,ax=plt.subplots(1,3,figsize=(15,5),layout='constrained')
for i,v in enumerate(['WT','AA','AA − WT']):
 m=maps[v] if i<2 else maps['AA']-maps['WT'];im=ax[i].imshow(m,origin='lower',vmin=0 if i<2 else -1,vmax=1,cmap='viridis' if i<2 else 'coolwarm');ax[i].set(title=v,xlabel='Construct residue index (zero based)',ylabel='Construct residue index (zero based)');fig.colorbar(im,ax=ax[i])
fig.savefig(out/'contacts.png',dpi=180);fig.savefig(out/'contacts.pdf');plt.close(fig)
# Comparison uncertainties describe within-trajectory time blocks, not independent replicates.
rows=[]
for col in s.observable.unique():
 w=s[(s.variant=='WT')&(s.observable==col)].iloc[0];m=s[(s.variant=='AA')&(s.observable==col)].iloc[0];se=np.hypot(w.block_SE,m.block_SE)
 rows.append({'observable':col,'AA_minus_WT':m['mean']-w['mean'],'approx95_low':m['mean']-w['mean']-2.1*se,'approx95_high':m['mean']-w['mean']+2.1*se})
pd.DataFrame(rows).to_csv(out/'WT_AA_comparison.csv',index=False)
print(s.to_string(index=False),flush=True)
