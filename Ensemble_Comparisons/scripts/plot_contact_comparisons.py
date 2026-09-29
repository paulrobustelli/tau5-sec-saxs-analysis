"""AR_ligand_binding notebook cells 26/27, batched equivalent calculation."""
from pathlib import Path
import json, hashlib
import numpy as np
import mdtraj as md
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'outputs/publish_snapshot/Ensemble_Comparisons';OUT.mkdir(exist_ok=True)
DATA=OUT/'Contact_Maps';DATA.mkdir(exist_ok=True)
nb=json.loads((ROOT/'outputs/Tau5_joint_refinement/reference_2022/Zhu_AR_Ligands_Apo.9.27.22.ipynb').read_text())
(DATA/'original_contact_cells.py').write_text('\n\n'.join(''.join(nb['cells'][i]['source']) for i in (26,27)))
snap=json.loads((ROOT/'outputs/Tau5_ASTEROIDS/base_pool_helicity/snapshot.json').read_text())
seqs=json.loads((ROOT/'outputs/Tau5_bAIes/manifest.json').read_text())['samples']
pairs=np.array([(i,j) for i in range(120) for j in range(i+1,120)])
maps={};meta={'definition':'MDTraj closest-heavy < 1.2 nm; diagonal zero; all off-diagonal sequence separations included','source_notebook':'paulrobustelli/AR_ligand_binding: Zhu_AR_Ligands_Apo.9.27.22.ipynb cells 26,27','samples':{}}
for sample in ('WT','AA'):
 paths=snap['source_pdbs'][sample+'_1']; cache=DATA/f'{sample}_full.npy'
 if cache.exists(): matrix=np.load(cache)
 else:
  counts=np.zeros(len(pairs),dtype=np.int64)
  for start in range(0,len(paths),50):
   trajs=[md.load(str(ROOT/p)) for p in paths[start:start+50]]
   for t in trajs:
    assert t.top.to_fasta()[0]==seqs[sample]['sequence']
   trj=md.join(trajs,check_topology=True)
   d,pr=md.compute_contacts(trj,pairs,scheme='closest-heavy',periodic=False)
   assert np.array_equal(pr,pairs)
   counts+=(d<1.2).sum(axis=0)
   if start%500==0: print(sample,start,flush=True)
  matrix=np.zeros((120,120));matrix[pairs[:,0],pairs[:,1]]=counts/len(paths);matrix+=matrix.T
  np.save(cache,matrix)
 assert np.allclose(matrix,matrix.T) and np.all(np.diag(matrix)==0) and matrix.min()>=0 and matrix.max()<=1
 maps[sample]=matrix
 np.savetxt(DATA/f'{sample}_full.csv',matrix,delimiter=',',fmt='%.6f')
 np.savetxt(DATA/f'{sample}_R2R3.csv',matrix[63:119,63:119],delimiter=',',fmt='%.6f')
 meta['samples'][sample]={'n':len(paths),'members':paths,'sequence':seqs[sample]['sequence'],'uniform_weights':True,'supplements_included':False,'fitted':False}
for view,sl in [('full',slice(None)),('R2R3',slice(63,119))]:
 fig,axes=plt.subplots(1,2,figsize=(15,6.5),layout='constrained')
 for ax,sample in zip(axes,('WT','AA')):
  mat=maps[sample][sl,sl];im=sns.heatmap(mat,cmap='jet',ax=ax,vmin=0,vmax=1,cbar_kws={'label':'Contact probability'},rasterized=True)
  ax.grid(which='both',alpha=.5);ax.invert_yaxis()
  if view=='full':
   ticks=np.array([0,12,32,52,72,92,112,119]);labels=['GP-G','340','360','380','400','420','440','447']
  else:
   ticks=np.arange(0,56,5);labels=[str(391+i) for i in ticks]
  ax.set_xticks(ticks+.5,labels,rotation=45);ax.set_yticks(ticks+.5,labels,rotation=45)
  ax.set(xlabel='AR residue',ylabel='AR residue',title=f'{sample} base pool (n={len(meta["samples"][sample]["members"]):,})')
 fig.suptitle(('Full Tau-5*: GP + AR330–447' if view=='full' else 'R2–R3 zoom: AR391–446, extracted from full Tau-5*')+'\nClosest heavy-atom distance < 1.2 nm; uniform base pools',fontsize=14)
 for ext in ('png','pdf'):fig.savefig(DATA/f'{view}_WT_AA.{ext}',dpi=160)
 plt.close(fig)
(DATA/'manifest.json').write_text(json.dumps(meta,indent=2))
print('Contact maps complete',flush=True)
