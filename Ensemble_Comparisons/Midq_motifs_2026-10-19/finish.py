from pathlib import Path
import sys,json,tarfile,hashlib,zipfile
import numpy as np,pandas as pd,mdtraj as md
from scipy.spatial.distance import pdist
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(__file__).resolve().parents[2];O=Path(__file__).resolve().parent;H=R/'outputs/Hydrophobic_SAXS_audit_2026-10-16';sys.path.insert(0,str(H));from analyze import load,C,NAMES
m=pd.read_csv(O/'matched_REST_ASTEROIDS.csv');chosen=[];used=set()
for _,row in m.iterrows():
 if int(row.REST_index) in used:continue
 chosen.append(row);used.add(int(row.REST_index))
 if len(chosen)==3:break
rows=[];edges=np.arange(0,251,1.);cent=(edges[:-1]+edges[1:])/2;out=O/'representatives';out.mkdir(exist_ok=True)
rc=np.load(C/'REST_WT/structural_cache.npz',allow_pickle=True);ac=np.load(C/'ASTEROIDS_WT/prior_structural_cache.npz',allow_pickle=True);prov=json.loads((C/'ASTEROIDS_WT/prior_structure_provenance.json').read_text());tf=tarfile.open(prov['archive']);members={Path(x.name).stem:x for x in tf.getmembers() if x.name.endswith('.pdb')}
fig,axes=plt.subplots(3,3,figsize=(14,12),layout='constrained');pdbprov=[]
for k,row in enumerate(chosen):
 for name,i in [('REST_WT',int(row.REST_index)),('ASTEROIDS_WT',int(row.AST_index))]:
  f=out/f'pair{k+1}_{name}_{i}.pdb'
  if name=='REST_WT':src=Path(rc['files'][i]);blob=src.read_bytes();source=str(src)
  else:cid=str(ac['conformer_ids'][i]);blob=tf.extractfile(members[cid]).read();source=prov['archive']+'::'+members[cid].name
  f.write_bytes(blob);pdbprov.append(dict(file=f.name,source=source,sha256=hashlib.sha256(blob).hexdigest()))
  t=md.load(str(f));t=t.atom_slice([a.index for a in t.top.atoms if a.element.symbol!='H']);atoms=list(t.top.atoms);ri=np.array([a.residue.index for a in atoms]);ii,jj=np.triu_indices(len(atoms),1);ds=pdist(t.xyz[0]*10);sep=abs(ri[ii]-ri[jj]);Z=np.array([a.element.atomic_number for a in atoms]);ww=Z[ii]*Z[jj];assert ds.max()<250
  h=np.array([a.residue.name in {'ALA','VAL','ILE','LEU','MET','PHE','TRP','TYR','PRO','CYS'} and a.element.symbol in {'C','S'} and a.name not in {'CA','C'} for a in atoms]);hyd=h[ii]&h[jj]
  for label,mask in [('same_residue',sep==0),('local_1_6',(sep>0)&(sep<=6)),('intermediate_7_20',(sep>6)&(sep<=20)),('long_range_gt20',sep>20),('hydrophobic_nonlocal',hyd&(sep>6))]:
   hist=np.histogram(ds[mask],edges,weights=ww[mask])[0];hist=hist/ww.sum() # common all-pair denominator, retains pair mass
   rows.extend(dict(pair=k+1,ensemble=name,source_index=i,category=label,r_A=r,value=v) for r,v in zip(cent,hist))
   if label in ['local_1_6','intermediate_7_20','long_range_gt20']:axes[k,0].plot(cent,hist,label=name.split('_')[0]+' '+label,ls='-' if name=='REST_WT' else '--')
  ca=t.xyz[0,t.top.select('name CA')]*10;dist=np.linalg.norm(ca[:,None]-ca[None,:],axis=2);col=1 if name=='REST_WT' else 2;im=axes[k,col].imshow(dist,vmin=0,vmax=50,cmap='viridis_r',origin='lower');axes[k,col].set(title=f'{name} CA distances; HSA {row.REST_HSA if name=="REST_WT" else row.AST_HSA:.0f} Å²',xlabel='Model residue index (GP + AR330–448)',ylabel='Model residue index');fig.colorbar(im,ax=axes[k,col],label='Å',shrink=.7)
 axes[k,0].set(xlim=(0,100),xlabel='Heavy-atom distance (Å)',ylabel='ZᵢZⱼ pair mass / Å');axes[k,0].legend(fontsize=6)
fig.suptitle('Representative distance fingerprints: electron-count weighted heavy atoms\nStructural diagnostic: q=0 atomic electron weights; excludes hydrogen, solvent and q-dependent form factors')
for ext in ['png','pdf']:fig.savefig(O/f'representative_distance_fingerprints.{ext}',dpi=160)
pd.DataFrame(rows).to_csv(O/'representative_heavy_atom_Pr.csv',index=False);(out/'provenance.json').write_text(json.dumps(pdbprov,indent=2));tf.close()
# Matched-Rg pair-distance / native curve comparisons: equal populations per 2A bin.
fig,axes=plt.subplots(7,3,figsize=(15,25),layout='constrained');summary=[]
for row,name in enumerate(NAMES):
 d=pd.read_csv(O/(name+'_metrics.csv'));z=np.load(O/(name+'.npz'));_,w,raw,q,y,sy=load(name);groups=[[],[]];bins=np.floor(d.Rg_CA_A/2).astype(int)
 for b in np.unique(bins):
  ii=np.flatnonzero(bins==b);ii=ii[np.argsort(d.buried_HSA_nonlocal_A2.iloc[ii],kind='stable')];half=len(ii)//2
  if half>=3:groups[0].extend(ii[:half]);groups[1].extend(ii[-half:])
 parent=raw.mean(1);scale=np.sum(parent*y/sy**2)/np.sum((parent/sy)**2)
 for label,ii in zip(['low burial','high burial'],groups):
  for c in [1,2]:
   hist=z['hist'][ii,c].mean(0);hist=hist/hist.sum();axes[row,c-1].plot(cent,hist,label=label)
  axes[row,2].plot(q,scale*raw[:,ii].mean(1)/y,label=label)
  summary.append(dict(ensemble=name,group=label,N=len(ii),Rg=d.Rg_CA_A.iloc[ii].mean(),HSA=d.buried_HSA_nonlocal_A2.iloc[ii].mean()))
 for c in [0,1]:axes[row,c].set(xlim=(0,100),title=name+(' • seqsep 7–20' if c==0 else ' • seqsep >20'),xlabel='CA distance (Å)',ylabel='Pair probability / Å');axes[row,c].legend(fontsize=7)
 axes[row,2].axhline(1,color='k',lw=.6);axes[row,2].axvspan(.08,.25,alpha=.1,color='grey');axes[row,2].set(ylim=(0,2.5),title=name+' • fixed parent scale',xlabel='q (Å⁻¹)',ylabel='Native I(q) / experiment');axes[row,2].legend(fontsize=7)
for ext in ['png','pdf']:fig.savefig(O/f'burial_distance_midq.{ext}',dpi=150)
pd.DataFrame(summary).to_csv(O/'matched_burial_groups.csv',index=False)
# Rank correlation heatmap.
c=pd.read_csv(O/'feature_correlations.csv').pivot(index='ensemble',columns='feature',values='within_Rg_bin_rank_correlation');fig,ax=plt.subplots(figsize=(14,5),layout='constrained');im=ax.imshow(c,vmin=-1,vmax=1,cmap='coolwarm');ax.set_xticks(range(len(c.columns)),c.columns,rotation=25,ha='right',fontsize=8);ax.set_yticks(range(len(c)),c.index);fig.colorbar(im,ax=ax,label='Within-Rg-bin rank correlation with enrichment benefit');fig.savefig(O/'feature_correlations.png',dpi=160)
print(m.head(8).to_string(index=False));print(c.to_string())
