from pathlib import Path
import numpy as np,mdtraj as md,csv,json,shutil
ROOT=Path(__file__).resolve().parents[2]; OUT=ROOT/'outputs/ASTEROIDS_master';DATA=OUT/'data';DATA.mkdir(parents=True,exist_ok=True)
ARCH=ROOT/'outputs/Tau5_joint_refinement/reference_2022/archive/Zhu_et_al_NatComms_9.27.22/AndrogenReceptor_LigandBinding_Trajectory_Files_9.27.22'
ref=md.load(str(ARCH/'R2R3_all_helix.pdb'));refs=np.array([ref.xyz[0,ref.top.select(f'residue {i} to {i+6} and backbone')] for i in range(1,51)],float)
def rmsd(x,y):
 x=x-x.mean(axis=-2,keepdims=True);y=y-y.mean(axis=-2,keepdims=True);u,s,v=np.linalg.svd(np.einsum('...ni,...nj->...ij',x,y));sign=np.where(np.linalg.det(u@v)<0,-1.,1.)
 return np.sqrt(np.maximum(0,(np.sum(x*x,axis=(-2,-1))+np.sum(y*y,axis=(-2,-1))-2*(s[...,0]+s[...,1]+sign*s[...,2]))/x.shape[-2]))
def calc(bb):
 nr=bb.shape[1];rr=np.repeat(refs[:1],nr-6,axis=0)
 if nr==120:rr[63:113]=refs
 else:rr=refs
 w=np.stack([bb[:,i:i+7].reshape(len(bb),28,3) for i in range(nr-6)],axis=1).astype(float);q=(rmsd(w,rr)/.1)**4;sa=((1+q)/(1+q+q*q)).sum(1);ca=bb[:,:,1];rg=np.sqrt(((ca-ca.mean(1,keepdims=True))**2).sum(2).mean(1));return sa,rg
meta={'seeds':[20260929,20261029,20261129],'conditions':{},'definitions':{'contacts':'closest-heavy <1.2nm; diagonal0; adjacent pairs included','helix':'DSSP H; maximal consecutive runs; R2R3 runs clipped at segment edges','Sa':'7-residue N/CA/C/O RMSD switch r0=.1nm, exponents8/12; original AR executable cell','Rg':'geometric C-alpha Rg in nm'}}
for protein in ['WT','AA']:
 src=ROOT/'outputs/Tau5_ASTEROIDS'/protein/'replicate_1';records=list(csv.DictReader((src/'expanded/candidates.csv').open()));paths=[r['source_pdb'] for r in records];cache=DATA/f'{protein}_structural.npz'
 if not cache.exists():
  traj=md.load(str(ROOT/f'output/ensembles/Tau5_DCD_ensembles/{protein}_expanded.dcd'),top=str(ROOT/f'output/ensembles/Tau5_DCD_ensembles/{protein}.pdb'))
  assert len(traj)==len(records)
  # Match original PDB reader: convert decimal angstrom coordinates to nm before float32 casting.
  traj.xyz=np.array([[[float(l[30:38])/10,float(l[38:46])/10,float(l[46:54])/10] for l in (ROOT/r['source_pdb']).read_text().splitlines() if l.startswith(('ATOM  ','HETATM'))] for r in records],dtype=np.float32)
  inds=np.array([[next(a.index for a in res.atoms if a.name==n) for n in ['N','CA','C','O']] for res in traj.top.residues]);assert inds.shape==(120,4)
  pairs=np.array([(i,j) for i in range(120) for j in range(i+1,120)]);contacts=[];values=[]
  for start in range(0,len(traj),100):
   t=traj[start:start+100];bb=t.xyz[:,inds];a,b=calc(bb);c,d=calc(bb[:,63:119]);values.append(np.c_[a,b,c,d]);dist,pr=md.compute_contacts(t,contacts=pairs,scheme='closest-heavy',periodic=False);assert np.array_equal(pairs,pr);contacts.append(dist<1.2)
   if start%500==0:print(protein,start,'/',len(traj),flush=True)
  h=np.array([json.loads((ROOT/r['protonated_pdb']).with_name('result.json').read_text())['dssp'] for r in records])=='H';assert h.shape==(len(records),120)
  v=np.concatenate(values);np.savez_compressed(cache,values=v,helix=h,contacts=np.concatenate(contacts),pairs=pairs)
 # record pool/selected mappings; copy experimental/fit numerical evidence
 for mode in ['base','expanded']:
  key=f'{protein}_{mode}';poolrows=list(csv.DictReader((src/mode/'candidates.csv').open()));lookup={p:i for i,p in enumerate(paths)};cfg={'protein':protein,'mode':mode,'pool':[lookup[r['source_pdb']] for r in poolrows],'selected':[]}
  for j,seed in enumerate(meta['seeds'],1):
   run=src/mode/f'seed_{seed}';sel=list(csv.DictReader((run/'selected.csv').open()));idx=[lookup[r['source_pdb']] for r in sel];assert len(set(idx))==500;assert set(idx)<=set(cfg['pool']);cfg['selected'].append(idx)
   for f in ['shifts.csv','saxs.csv','summary.json','helicity.csv']:
    shutil.copyfile(run/f,DATA/f'{key}_ens{j}_{f}')
  meta['conditions'][key]=cfg
 shutil.copyfile(ROOT/f'outputs/Tau5_joint_refinement/{protein}/d2d_populations.csv',DATA/f'{protein}_d2d.csv')
 (DATA/f'{protein}_candidate_paths.json').write_text(json.dumps(paths))
mdcache=DATA/'MD_structural.npz'
if not mdcache.exists():
 old=np.load(ROOT/'outputs/Tau5_ASTEROIDS/Sa_Rg/original_apo_observables.npz'); hs=[]
 for i,t in enumerate(md.iterload(str(ARCH/'Tau5R2R3_apo.xtc'),top=str(ARCH/'Tau5R2R3_apo.pdb'),chunk=1000)):
  core=t.atom_slice(t.top.select('residue 1 to 56')); hs.append(md.compute_dssp(core,simplified=False)=='H')
  if i%10==0:print('MD DSSP',i*1000,flush=True)
 h=np.concatenate(hs);assert len(h)==len(old['sa'])==57144
 cm=np.load(ROOT/'outputs/publish_snapshot/Ensemble_Comparisons/Contact_Maps/MD_R2R3.npz');print('MD contact keys',cm.files,flush=True)
 mat=np.loadtxt(ROOT/'outputs/publish_snapshot/Ensemble_Comparisons/Contact_Maps/MD_R2R3.csv',delimiter=',');assert mat.shape==(56,56)
 np.savez_compressed(mdcache,values=np.c_[old['sa'],old['rg']],helix=h,contacts=mat)
for f in ['fit_metrics.csv','convergence.json','Rg_convergence.json']:
 shutil.copyfile(ROOT/'outputs/publish_snapshot/Ensemble_Comparisons/CS_SAXS_Selection'/f,DATA/f)
(DATA/'metadata.json').write_text(json.dumps(meta,indent=2));print('COMPLETE',flush=True)
