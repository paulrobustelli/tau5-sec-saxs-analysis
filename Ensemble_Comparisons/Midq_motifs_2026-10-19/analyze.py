from pathlib import Path
import sys,json,tarfile,io,hashlib
import numpy as np,pandas as pd,mdtraj as md
from scipy.stats import spearmanr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(__file__).resolve().parents[2];O=Path(__file__).resolve().parent
sys.path.insert(0,str(R/'outputs/Hydrophobic_SAXS_audit_2026-10-16'))
from analyze import load,C,NAMES
pairs=np.array([(i,j) for i in range(121) for j in range(i+1,121)]);sep=np.diff(pairs,axis=1)[:,0];edges=np.arange(0,251,1.);cent=(edges[:-1]+edges[1:])/2
classes=[('local_1_6',sep<=6),('intermediate_7_20',(sep>6)&(sep<=20)),('long_range_gt20',sep>20)]
def getxyz(name):
 z=np.load(C/name/('prior_structural_cache.npz' if name.startswith('AST') else 'structural_cache.npz'),allow_pickle=True)
 if name.startswith('AST'):
  ids=z['conformer_ids'].astype(str);p=json.loads((C/name/'prior_structure_provenance.json').read_text());tf=tarfile.open(p['archive']);mem={Path(m.name).stem:m for m in tf.getmembers() if m.name.endswith('.pdb')};out=[]
  for cid in ids:
   lines=tf.extractfile(mem[cid]).read().decode().splitlines();a=[[float(s[30:38]),float(s[38:46]),float(s[46:54])] for s in lines if s.startswith(('ATOM  ','HETATM')) and s[12:16].strip()=='CA'];assert len(a)==121;out.append(a)
  tf.close();return np.array(out)
 return np.array([md.load(f).xyz[0,md.load(f).top.select('name CA')]*10 for f in z['files'].astype(str)])
allrows=[];corr=[];groupcurves=[]
for name in NAMES:
 d,w,raw,q,y,sy=load(name);p=O/(name+'.npz')
 if not p.exists():
  xyz=getxyz(name);dist=np.linalg.norm(xyz[:,pairs[:,0]]-xyz[:,pairs[:,1]],axis=2);assert dist.max()<edges[-1]
  hist=np.array([[np.histogram(v[m],edges)[0] for _,m in classes] for v in dist],dtype=np.int32)
  centered=xyz-xyz.mean(1,keepdims=True);ev=np.linalg.eigvalsh(np.einsum('nri,nrj->nij',centered,centered)/121);rg=np.sqrt(ev.sum(1));assert np.allclose(rg,d.Rg_CA_A,atol=2e-4)
  anis=1-3*(ev[:,0]*ev[:,1]+ev[:,0]*ev[:,2]+ev[:,1]*ev[:,2])/ev.sum(1)**2
  windows=np.array([np.sqrt(((xyz[:,i:i+20]-xyz[:,i:i+20].mean(1,keepdims=True))**2).sum(2).mean(1)) for i in range(2,102,5)]).T
  # Anti-parallel 10-residue segment chords with close segment centers: broad bends, not atom contacts.
  starts=np.arange(2,112,5);chords=xyz[:,starts+9]-xyz[:,starts];chords/=np.linalg.norm(chords,axis=2,keepdims=True);midpoints=(xyz[:,starts+9]+xyz[:,starts])/2
  bends=np.zeros(len(xyz),int)
  for i in range(len(starts)):
   for j in range(i+1,len(starts)):
    if starts[j]-starts[i]>=20:bends+=((chords[:,i]*chords[:,j]).sum(1)<-.7)&(np.linalg.norm(midpoints[:,i]-midpoints[:,j],axis=1)<20)
  np.savez_compressed(p,xyz=xyz,hist=hist,anisotropy=anis,window_rg_min=windows.min(1),window_rg_cv=windows.std(1)/windows.mean(1),broad_antiparallel_pairs=bends,rg=rg)
 z=np.load(p);df=d.copy()
 for f in ['anisotropy','window_rg_min','window_rg_cv','broad_antiparallel_pairs']:df[f]=z[f]
 df['ensemble']=name;df['combined_weight']=w.combined
 base=raw.mean(1);scale=np.sum(base*y/sy**2)/np.sum((base/sy)**2);res=(scale*base-y)/sy;mid=(q>.08)&(q<.25)
 # Infinitesimal enrichment benefit with parent scale fixed, per unit mixture fraction.
 delta=(scale*(raw-base[:,None]))/sy[:,None];df['midq_enrichment_benefit']=-2*np.mean(res[mid,None]*delta[mid],axis=0)
 # Equal-weight CA Debye decomposition is a geometric diagnostic, not the native solvent-aware forward model.
 kernel=np.sinc(q[:,None]*cent[None,:]/np.pi)
 debye=2*np.einsum('ncr,qr->ncq',z['hist'],kernel)/121**2
 np.savez_compressed(O/(name+'_distance_scattering.npz'),q=q,r_A=cent,classes=np.array([a for a,b in classes]),geometric_CA_pair_terms=debye,self_term=1/121,native_Pepsi=raw,conformer_ids=d.conformer_id.to_numpy(dtype=str))
 for ci,(label,m) in enumerate(classes):
  df[label+'_pairs_10_30A']=z['hist'][:,ci,10:30].sum(1)
 for f in ['buried_HSA_nonlocal_A2','hydrophobic_pairs_sep6_50','anisotropy','window_rg_cv','broad_antiparallel_pairs']+[a+'_pairs_10_30A' for a,b in classes]:
  # within 2A Rg-bin centered ranks remove gross size dependence without pooling sequences.
  a=df.groupby(np.floor(df.Rg_CA_A/2))[f].rank(pct=True);b=df.groupby(np.floor(df.Rg_CA_A/2)).midq_enrichment_benefit.rank(pct=True)
  corr.append(dict(ensemble=name,feature=f,within_Rg_bin_rank_correlation=spearmanr(a,b).statistic))
 df.to_csv(O/(name+'_metrics.csv'),index=False);allrows.append(df);print(name,'complete',flush=True)
pd.concat(allrows).to_csv(O/'all_metrics.csv',index=False);pd.DataFrame(corr).to_csv(O/'feature_correlations.csv',index=False)
# Match low-burial AST WT to high-burial REST WT within 1 A, using normalized native curve shape, no experiment used in match selection.
r=allrows[0];a=next(x for x in allrows if x.ensemble.iloc[0]=='ASTEROIDS_WT');_,_,rr,q,y,sy=load('REST_WT');_,_,ar,_,_,_=load('ASTEROIDS_WT');mid=(q>.08)&(q<.25);k=np.argmin(abs(q-.03));matches=[]
for j in np.flatnonzero((a.buried_HSA_nonlocal_A2<400)&(a.Rg_CA_A.between(24,32))):
 ii=np.flatnonzero((r.buried_HSA_nonlocal_A2>800)&(abs(r.Rg_CA_A-a.Rg_CA_A.iloc[j])<1))
 if not len(ii):continue
 err=np.sqrt(np.mean((np.log(rr[mid][:,ii]/rr[k,ii])-np.log(ar[mid,j:j+1]/ar[k,j]))**2,axis=0));i=ii[err.argmin()];matches.append(dict(REST_index=int(i),AST_index=int(j),log_midq_shape_RMS=float(err.min()),REST_Rg=r.Rg_CA_A.iloc[i],AST_Rg=a.Rg_CA_A.iloc[j],REST_HSA=r.buried_HSA_nonlocal_A2.iloc[i],AST_HSA=a.buried_HSA_nonlocal_A2.iloc[j]))
m=pd.DataFrame(matches).sort_values('log_midq_shape_RMS');m.to_csv(O/'matched_REST_ASTEROIDS.csv',index=False)
chosen=[];used=set()
for _,row in m.iterrows():
 if int(row.REST_index) in used:continue
 chosen.append(row);used.add(int(row.REST_index))
 if len(chosen)==3:break
fig=plt.figure(figsize=(16,12));zs=[np.load(O/'REST_WT.npz'),np.load(O/'ASTEROIDS_WT.npz')]
for row,pair in enumerate(chosen):
 ids=[int(pair.REST_index),int(pair.AST_index)]
 for c,(name,z,i) in enumerate(zip(['REST','ASTEROIDS'],zs,ids)):
  ax=fig.add_subplot(3,4,row*4+c+1,projection='3d');x=z['xyz'][i];x=x-x.mean(0);u,s,v=np.linalg.svd(x,full_matrices=False);x=x@v.T;ax.plot(*x.T,lw=1);ax.scatter(*x[::10].T,c=np.arange(0,121,10),cmap='viridis',s=15);ax.set_box_aspect((1,1,1));ax.set(xlim=(-60,60),ylim=(-60,60),zlim=(-60,60),title=f'{name} #{i}; Rg {z["rg"][i]:.1f} Å');ax.set_axis_off()
 ax=fig.add_subplot(3,4,row*4+3)
 for name,z,i in zip(['REST','ASTEROIDS'],zs,ids):
  ax.plot(cent,z['hist'][i,2]/sum(z['hist'][i,2]),label=name)
 ax.set(xlim=(0,100),xlabel='CA distance (Å)',ylabel='Long-range pair probability');ax.legend()
 ax=fig.add_subplot(3,4,row*4+4);ax.plot(q,rr[:,ids[0]]/rr[k,ids[0]],label='REST');ax.plot(q,ar[:,ids[1]]/ar[k,ids[1]],label='ASTEROIDS');ax.axvspan(.08,.25,alpha=.1,color='grey');ax.set(yscale='log',xlabel='q (Å⁻¹)',ylabel='Native I(q) / I(0.03)');ax.legend()
fig.suptitle('Packed REST versus low-burial ASTEROIDS: matched Rg and native mid-q shape\nCA traces; geometric long-range distances; shape normalization only, not experimental fitting');fig.tight_layout()
for ext in ['png','pdf']:fig.savefig(O/f'matched_motifs.{ext}',dpi=160)
(O/'complete.json').write_text(json.dumps({'ensembles':NAMES,'matched_candidates':len(m),'classes':[a for a,b in classes],'note':'CA geometric pair histogram and equal-weight Debye diagnostic; native Pepsi retained, no solvent decomposition claimed'},indent=2))
