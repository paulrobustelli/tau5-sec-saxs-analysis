from compute import *
from scipy.stats import spearmanr
checks=[]
for name in ['REST_WT','IDPForge_WT','bAIes_WT','ASTEROIDS_WT','IDPForge_AA','bAIes_AA','ASTEROIDS_AA']:
 dest=O/name
 if (dest/'tertiary_contacts.csv').exists():continue
 ast=name.startswith('ASTEROIDS');src=C/name
 with tempfile.TemporaryDirectory() as tmp:
  if ast:
   ids=pd.read_csv(dest/'metrics.csv').conformer_id.astype(str);prov=json.loads((src/'prior_structure_provenance.json').read_text());files=[]
   with tarfile.open(prov['archive']) as tf:
    members={Path(m.name).stem:m for m in tf.getmembers() if m.isfile() and m.name.endswith('.pdb')}
    for cid in ids:
     f=Path(tmp)/(cid+'.pdb');f.write_bytes(tf.extractfile(members[cid]).read());files.append(str(f))
  else:files=np.load(src/'structural_cache.npz')['files'].astype(str).tolist()
  rows=[];pairs=np.array([(i,j) for i in range(2,121) for j in range(i+7,121)])
  for start in range(0,len(files),40):
   t=md.join([load_heavy(f) for f in files[start:start+40]],check_topology=True);dist,_=md.compute_contacts(t,contacts=pairs,scheme='closest-heavy',periodic=False)
   rows.append(pd.DataFrame({f'all_tertiary_pairs_sep6_{int(c*100)}':(dist<c).sum(1) for c in [.5,.8,1.2]}))
  pd.concat(rows,ignore_index=True).to_csv(dest/'tertiary_contacts.csv',index=False)
  # Fixed, evenly spaced 12-frame discretization sensitivity test, same coords.
  ii=np.linspace(0,len(files)-1,12,dtype=int);t=md.join([load_heavy(files[i]) for i in ii],check_topology=True);high,*_=calc(t,points=960);low=pd.read_csv(dest/'metrics.csv').iloc[ii]
  a=low.buried_HSA_nonlocal_A2.to_numpy();b=high.buried_HSA_nonlocal_A2.to_numpy()
  checks.append({'ensemble':name,'frames':ii.tolist(),'HSA_240':a.tolist(),'HSA_960':b.tolist(),'mean_abs_difference_A2':float(abs(a-b).mean()),'max_abs_difference_A2':float(abs(a-b).max()),'spearman':float(spearmanr(a,b).statistic)})
  print(name,'tertiary + SASA check complete',flush=True)
(O/'SASA_resolution_check.json').write_text(json.dumps(checks,indent=2))
