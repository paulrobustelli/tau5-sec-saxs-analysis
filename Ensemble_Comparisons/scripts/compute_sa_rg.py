"""Repository-matched backbone S-alpha and CA Rg for Tau-5 pools and 2022 apo MD."""
import csv,json,time
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import mdtraj as md
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'outputs/Tau5_joint_refinement';OUT=ROOT/'outputs/Tau5_ASTEROIDS/Sa_Rg';OUT.mkdir(exist_ok=True)
ARCH=BASE/'reference_2022/archive/Zhu_et_al_NatComms_9.27.22/AndrogenReceptor_LigandBinding_Trajectory_Files_9.27.22'
ref=md.load(str(ARCH/'R2R3_all_helix.pdb'))
refs=np.array([ref.xyz[0,ref.top.select(f'residue {i} to {i+6} and backbone')] for i in range(1,51)],float)
assert refs.shape==(50,28,3)
def rmsd(x,y):
 x=np.asarray(x,float);y=np.asarray(y,float);x=x-x.mean(axis=-2,keepdims=True);y=y-y.mean(axis=-2,keepdims=True)
 cov=np.einsum('...ni,...nj->...ij',x,y);u,s,v=np.linalg.svd(cov);sgn=np.where(np.linalg.det(u@v)<0,-1.,1.)
 ss=s[...,0]+s[...,1]+sgn*s[...,2]
 return np.sqrt(np.maximum(0,(np.sum(x*x,axis=(-2,-1))+np.sum(y*y,axis=(-2,-1))-2*ss)/x.shape[-2]))
def switch(r):
 # Algebraically cancels removable 0/0 at r/r0=1; limit=2/3.
 q=(r/.10)**4
 return (1+q)/(1+q+q*q)
def calculate(bb):
 nres=bb.shape[1];refwindows=np.repeat(refs[:1],nres-6,axis=0)
 if nres==120:refwindows[63:113]=refs
 elif nres==56:refwindows=refs
 windows=np.stack([bb[:,i:i+7].reshape(len(bb),28,3) for i in range(nres-6)],axis=1)
 sa=switch(rmsd(windows,refwindows)).sum(axis=1)
 ca=bb[:,:,1];rg=np.sqrt(np.mean(np.sum((ca-ca.mean(axis=1,keepdims=True))**2,axis=2),axis=1))
 return sa,rg
# Check rigid-body invariance, singular switching point, and agreement with MDTraj.
rng=np.random.default_rng(42);rot=np.linalg.qr(rng.normal(size=(3,3)))[0]
x=refs[0]@rot+3;assert rmsd(x,refs[0])<1e-7;assert abs(switch(np.array(.1))-2/3)<1e-12
inds=ref.top.select('residue 1 to 7 and backbone');tr=ref.atom_slice(inds);pert=md.Trajectory((refs[0]+rng.normal(0,.02,refs[0].shape))[None].astype('float32'),tr.top)
assert np.isclose(rmsd(pert.xyz[0],refs[0]),md.rmsd(pert,tr,parallel=False)[0],atol=1e-6)
# PDB coordinates are Angstrom; use nm as in the source notebook.
def backbone(path):
 arr=np.full((120,4,3),np.nan);names={'N':0,'CA':1,'C':2,'O':3}
 for line in path.read_text().splitlines():
  if line.startswith('ATOM') and line[12:16].strip() in names:
   i=int(line[22:26])-1
   arr[i,names[line[12:16].strip()]]=[float(line[30:38])/10,float(line[38:46])/10,float(line[46:54])/10]
 if not np.isfinite(arr).all():raise ValueError(f'Incomplete backbone {path}')
 return arr
snap=json.loads((ROOT/'outputs/Tau5_ASTEROIDS/base_pool_helicity/snapshot.json').read_text());data={};rows=[]
for sample in ['WT','AA']:
 cache=OUT/f'{sample}_observables.npz';paths=snap['source_pdbs'][f'{sample}_1']
 if cache.exists():
  z=np.load(cache);fullsa,fullrg,sa,rg=[z[k] for k in ['full_sa','full_rg','r2r3_sa','r2r3_rg']]
 else:
  bb=np.array([backbone(ROOT/p) for p in paths]);fullsa,fullrg=calculate(bb);sa,rg=calculate(bb[:,63:119]);np.savez_compressed(cache,full_sa=fullsa,full_rg=fullrg,r2r3_sa=sa,r2r3_rg=rg)
 data[sample]={'full':np.c_[fullsa,fullrg],'R2R3':np.c_[sa,rg]}
 for p,a,b,c,d in zip(paths,fullsa,fullrg,sa,rg):rows.append([sample,p,a,b,c,d])
 n=1000 if sample=='WT' else 500
 sel={r['source_pdb'] for r in csv.DictReader((ROOT/f'outputs/Tau5_ASTEROIDS/posthoc_d2d/{sample}_pool1_{n}_members.csv').open())};mask=np.array([p in sel for p in paths]);data[sample+'_calibrated']={k:v[mask] for k,v in data[sample].items()}
 print(sample,len(paths),'complete',flush=True)
with (OUT/'pool_observables.csv').open('w') as f:
 w=csv.writer(f);w.writerow(['protein','source_pdb','full_Salpha','full_CA_Rg_nm','AR391_446_Salpha','AR391_446_CA_Rg_nm']);w.writerows(rows)
# Recalculate the original apo MD using precisely the 50 original reference windows.
mdfile=ARCH/'Tau5R2R3_apo.xtc';mdcache=OUT/'original_apo_observables.npz'
if mdcache.exists():
 z=np.load(mdcache);data['2022 apo MD']={'R2R3':np.c_[z['sa'],z['rg']]}
elif mdfile.exists():
 top=md.load(str(ARCH/'Tau5R2R3_apo.pdb')).top;idx=top.select('residue 1 to 56 and backbone');assert len(idx)==224
 values=[]
 for chunk in md.iterload(str(mdfile),top=top,chunk=1000,atom_indices=idx):
  bb=chunk.xyz.reshape(chunk.n_frames,56,4,3);sa,rg=calculate(bb);values.append(np.c_[sa,rg])
 a=np.concatenate(values);data['2022 apo MD']={'R2R3':a};np.savez_compressed(mdcache,sa=a[:,0],rg=a[:,1]);print('2022 apo MD frames',len(a),flush=True)
summary={}
for label,views in data.items():
 summary[label]={}
 for view,a in views.items():summary[label][view]={'n':len(a),'mean_Salpha':float(a[:,0].mean()),'mean_Rg_nm':float(a[:,1].mean()),'Sa_Rg_correlation':float(np.corrcoef(a.T)[0,1]),'fraction_Rg_lt_1p3_Sa_gt_6':float(np.mean((a[:,1]<1.3)&(a[:,0]>6))) if view=='R2R3' else None}
(OUT/'summary.json').write_text(json.dumps(summary,indent=2))
def plot(name,keys,view,title):
 arrays=[data[k][view] for k in keys];xmax=max(25.,np.ceil(max(a[:,0].max() for a in arrays)/5)*5);ymin=np.floor(min(a[:,1].min() for a in arrays)*10)/10;ymax=np.ceil(max(a[:,1].max() for a in arrays)*10)/10
 xe=np.linspace(0,xmax,31);ye=np.linspace(ymin,ymax,31);hists=[np.histogram2d(a[:,0],a[:,1],bins=[xe,ye])[0]/len(a) for a in arrays];vmax=max(h.max() for h in hists)*100
 fig,axes=plt.subplots(1,len(keys),figsize=(5*len(keys)+1,4.8),sharex=True,sharey=True,layout='constrained')
 cmap=plt.colormaps['magma'].copy();cmap.set_bad('white')
 for ax,key,a,h in zip(np.atleast_1d(axes),keys,arrays,hists):
  im=ax.pcolormesh(xe,ye,np.ma.masked_equal(100*h.T,0),cmap=cmap,vmin=0,vmax=vmax,rasterized=True)
  ax.set(title=f'{key}\nn={len(a):,}; mean Rg={a[:,1].mean():.2f} nm',xlabel='Sα',xlim=(0,xmax),ylim=(ymin,ymax));ax.plot(a[:,0].mean(),a[:,1].mean(),'o',mfc='none',mec='cyan',ms=8,mew=1.5)
  if view=='R2R3':ax.axvline(6,color='#27b8ce',ls='--',lw=.7);ax.axhline(1.3,color='#27b8ce',ls='--',lw=.7)
 np.atleast_1d(axes)[0].set_ylabel('Cα Rg (nm)');fig.colorbar(im,ax=axes,label='Structures per bin (%)',shrink=.8)
 fig.suptitle(title+'\nUniform populations; cyan circle = mean. Sampling distributions, not thermodynamic free energies.',fontsize=11)
 for ext in ['png','pdf']:fig.savefig(OUT/f'{name}.{ext}',dpi=160)
 plt.close(fig)
plot('full_chain_WT_AA',['WT','AA'],'full','Full Tau-5*: GP + AR330–447 (120 residues; 114 seven-residue windows)')
keys=['WT','AA']+(['2022 apo MD'] if '2022 apo MD' in data else [])
plot('R2R3_side_by_side',keys,'R2R3','Matching Tau-5 R2–R3 residues: AR391–446 (56 residues; 50 seven-residue windows)')
plot('R2R3_calibrated_side_by_side',['WT_calibrated','AA_calibrated']+keys[2:],'R2R3','AR391–446: D2D-calibrated subsets versus original apo MD')
(OUT/'README.md').write_text('S-alpha follows executable cell 16 in Zhu_AR_Ligands_Apo.9.27.22.ipynb: 7-residue N/CA/C/O windows, optimally superposed RMSD to R2R3_all_helix.pdb, r0=0.10 nm, switch (1-(r/r0)^8)/(1-(r/r0)^12), summed over windows. The unused calc_SA helper differs by selecting CA only; we use the actual analysis cell. AR391–446 maps to model residues 64–119 and uses the original 50 reference windows exactly. For the full 120-residue chain, windows lying entirely within AR391–446 use their corresponding original references; all other windows use the first seven-residue reference as a generic helical template. Full score has 114 terms, cropped score 50, so totals are not length-normalized. Rg is geometric CA-only Rg (nm), as in the notebook. Original apo trajectory is unweighted, uncapped residues 391–446 selected geometrically; the simulation itself retains its original end caps. Current segment coordinates are extracted from full-chain conformers, not separately simulated cropped constructs. Current pool 1 has 2000 structures per protein, supplements excluded. Calibration subsets use 1000 WT and 500 AA. Histograms show fractions, not free energies for these generated pools. Cyan threshold lines reproduce the study Rg<1.3 nm and S-alpha>6 definition. Structural data retrieved from Zenodo 7120845.\n')
print(OUT)
