from pathlib import Path
import json,os,tarfile,tempfile,time,hashlib
import numpy as np,pandas as pd,mdtraj as md
R=Path(__file__).resolve().parents[2];O=Path(__file__).resolve().parent;C=R/'outputs/Ensemble_Comparison_121_2026-10-15';OLD=Path('/Users/paulrobustelli/Documents/Codex/2026-09-28/is-th')
HYD={'ALA','VAL','ILE','LEU','MET','PHE','TRP','TYR','PRO','CYS'};ARO={'PHE','TRP','TYR'}
def calc(t,points=240):
 top=t.top;res=list(top.residues);n=len(t)
 assert len(res)==121
 ha={r.index:[a.index for a in r.atoms if a.element.symbol in ['C','S'] and a.name not in ['N','CA','C','O','OXT']] for r in res if r.name in HYD and r.index>=2}
 ha={i:a for i,a in ha.items() if a};hids=sorted(ha);targets=np.array([a for i in hids for a in ha[i]])
 full=md.shrake_rupley(t,probe_radius=.14,n_sphere_points=points,atom_indices=targets)[:,targets]
 # Area loss caused ONLY by atoms >6 residues away: compare full structure
 # with the unchanged local +/-6 residue fragment for each hydrophobic residue.
 burial=[]; exposed=[]
 for i in hids:
  atoms=[a.index for a in top.atoms if abs(a.residue.index-i)<=6]
  loc=t.atom_slice(atoms);lookup={a:j for j,a in enumerate(atoms)};ii=[lookup[a] for a in ha[i]]
  local=md.shrake_rupley(loc,probe_radius=.14,n_sphere_points=points,atom_indices=ii)[:,ii].sum(1)
  pos=[int(np.where(targets==a)[0][0]) for a in ha[i]];f=full[:,pos].sum(1)
  assert np.min(local-f)>-1e-4
  burial.append(np.maximum(local-f,0)*100);exposed.append(f*100)
 burial=np.array(burial).T;exposed=np.array(exposed).T
 pairs=np.array([(i,j) for i in hids for j in hids if j-i>3]);ds=[]
 for i,j in pairs:
  ap=np.array([(a,b) for a in ha[i] for b in ha[j]])
  ds.append(md.compute_distances(t,ap,periodic=False).min(1))
 ds=np.array(ds).T
 coords={}
 for sep in [3,6,10]:
  for cutoff in [.45,.5,.6]:
   use=(pairs[:,1]-pairs[:,0]>sep);ct=(ds<cutoff)&use[None,:];tag=f'sep{sep}_{int(cutoff*100):02d}'
   coords['hydrophobic_pairs_'+tag]=ct.sum(1)
   aro=np.array([res[i].name in ARO or res[j].name in ARO for i,j in pairs]);aa=np.array([res[i].name in ARO and res[j].name in ARO for i,j in pairs])
   coords['aromatic_hydrophobic_pairs_'+tag]=(ct&aro).sum(1);coords['aromatic_pairs_'+tag]=(ct&aa).sum(1)
 main=(ds<.5)&(pairs[:,1]-pairs[:,0]>6)[None,:];coord=np.zeros((n,len(hids)),int)
 for k,(i,j) in enumerate(pairs):coord[:,hids.index(i)]+=main[:,k];coord[:,hids.index(j)]+=main[:,k]
 ca=t.xyz[:,top.select('name CA')];rg=np.sqrt(((ca-ca.mean(1,keepdims=True))**2).sum(2).mean(1))*10
 coords.update(Rg_CA_A=rg,buried_HSA_nonlocal_A2=burial.sum(1),exposed_HSA_A2=exposed.sum(1),buried_aromatic_HSA_A2=burial[:,[res[i].name in ARO for i in hids]].sum(1),hydrophobic_mean_coordination=coord.mean(1),hydrophobic_max_coordination=coord.max(1),hydrophobic_residues_ge3_contacts=(coord>=3).sum(1))
 return pd.DataFrame(coords),burial,coord,np.array(hids)+328

def load_heavy(file):
 t=md.load(str(file));return t.atom_slice([a.index for a in t.top.atoms if a.element.symbol!='H'])

def main():
 for name in ['REST_WT','IDPForge_WT','bAIes_WT','ASTEROIDS_WT','IDPForge_AA','bAIes_AA','ASTEROIDS_AA']:
  d=O/name;d.mkdir(exist_ok=True)
  if (d/'complete.json').exists():continue
  src=C/name; ast=name.startswith('ASTEROIDS');tmp=None;tf=None
  if ast:
   cache=np.load(src/'prior_structural_cache.npz',allow_pickle=True);ids=cache['conformer_ids'].astype(str).tolist();prov=json.loads((src/'prior_structure_provenance.json').read_text());tf=tarfile.open(prov['archive']);members={Path(m.name).stem:m for m in tf.getmembers() if m.isfile() and m.name.endswith('.pdb')};tmp=tempfile.TemporaryDirectory();files=[]
   for cid in ids:
    f=Path(tmp.name)/(cid+'.pdb');f.write_bytes(tf.extractfile(members[cid]).read());files.append(str(f))
   sel=np.load(src/'structural_cache.npz',allow_pickle=True);after=np.zeros(len(ids))
   for cid,w in zip(sel['conformer_ids'].astype(str),sel['weights']):after[ids.index(cid)]+=w
   w=pd.DataFrame({'uniform':np.full(len(ids),1/len(ids)),'combined':after})
  else:
   cache=np.load(src/'structural_cache.npz',allow_pickle=True);files=cache['files'].astype(str).tolist();ids=[Path(f).stem for f in files];w=pd.read_csv(src/'weights.csv')[['uniform','combined']]
  w.to_csv(d/'weights.csv',index=False);rows=[];bur=[];co=[];starttime=time.time()
  for start in range(0,len(files),20):
   cp=d/f'chunk_{start:05d}.npz'
   if cp.exists():z=np.load(cp);frame=pd.DataFrame(z['values'],columns=z['columns']);b=z['burial'];c=z['coordination'];ars=z['AR_residues']
   else:
    t=md.join([load_heavy(f) for f in files[start:start+20]],check_topology=True);frame,b,c,ars=calc(t)
    np.savez_compressed(cp,values=frame.to_numpy(),columns=frame.columns.to_numpy(dtype=str),burial=b,coordination=c,AR_residues=ars)
   rows.append(frame);bur.append(b);co.append(c)
   if start%200==0:print(name,start,'seconds',round(time.time()-starttime),flush=True)
  result=pd.concat(rows,ignore_index=True);assert np.allclose(result.Rg_CA_A,cache['values'][:,0],atol=1e-4)
  result.insert(0,'conformer_id',ids);result.to_csv(d/'metrics.csv',index=False);np.savez_compressed(d/'per_residue.npz',buried_HSA_A2=np.concatenate(bur),coordination=np.concatenate(co),AR_residues=ars)
  (d/'complete.json').write_text(json.dumps({'N':len(files),'seconds':time.time()-starttime,'Rg_matches_existing_cache':True,'source_cache':str(src/('prior_structural_cache.npz' if ast else 'structural_cache.npz')),'source_cache_sha256':hashlib.sha256((src/('prior_structural_cache.npz' if ast else 'structural_cache.npz')).read_bytes()).hexdigest()},indent=2))
  if tf:tf.close();tmp.cleanup()
  print(name,'COMPLETE',flush=True)
if __name__=='__main__':main()
