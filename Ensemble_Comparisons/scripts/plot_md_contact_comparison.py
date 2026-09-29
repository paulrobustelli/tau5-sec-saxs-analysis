"""Original apo MD contacts using AR_ligand_binding cells 26–27."""
from pathlib import Path
import json,hashlib
import numpy as np
import mdtraj as md
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'outputs/publish_snapshot/Ensemble_Comparisons';D=OUT/'Contact_Maps'
pdb=next((ROOT/'outputs/Tau5_joint_refinement/reference_2022/archive').rglob('Tau5R2R3_apo.pdb'));xtc=pdb.with_suffix('.xtc')
top=md.load(str(pdb)).top
seq=top.to_fasta()[0]
assert seq==json.loads((ROOT/'outputs/Tau5_bAIes/manifest.json').read_text())['samples']['WT']['sequence'][63:119]
assert top.residue(0).name=='ACE' and top.residue(57).name=='NH2'
pairs=np.array([(i,j) for i in range(1,57) for j in range(i+1,57)])
cache=D/'MD_R2R3.npz'
if cache.exists():
 z=np.load(cache);matrix=z['matrix'];n=int(z['n'])
else:
 counts=np.zeros(len(pairs),dtype=np.int64);n=0
 for trj in md.iterload(str(xtc),top=str(pdb),chunk=1000):
  # Same defaults as original md.compute_contacts(trj, [[i,j]]).
  dist,pr=md.compute_contacts(trj,pairs,scheme='closest-heavy',periodic=True)
  assert np.array_equal(pr,pairs)
  counts+=(dist<1.2).sum(axis=0);n+=len(trj)
  if n%10000==0:print(n,flush=True)
 matrix=np.zeros((56,56));matrix[pairs[:,0]-1,pairs[:,1]-1]=counts/n;matrix+=matrix.T
 np.savez_compressed(cache,matrix=matrix,n=n)
assert n==57144 and np.allclose(matrix,matrix.T) and np.all(np.diag(matrix)==0)
np.savetxt(D/'MD_R2R3.csv',matrix,delimiter=',',fmt='%.6f')
provenance={'frames':n,'sequence':seq,'AR_residues':[391,446],'topology_residue_indices':[1,56],'caps_excluded':['ACE','NH2'],'cutoff_nm':1.2,'scheme':'closest-heavy','periodic':True,'weighting':'uniform over all saved frames, no reweighting or burn-in exclusion','source':'https://zenodo.org/records/7120845','local_files':{str(p.relative_to(ROOT)):{'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in [pdb,xtc]}}
(D/'MD_provenance.json').write_text(json.dumps(provenance,indent=2))
def draw(ax,m,title):
 sns.heatmap(m,cmap='jet',ax=ax,vmin=0,vmax=1,cbar_kws={'label':'Contact probability'},rasterized=True)
 ax.grid(which='both',alpha=.5);ax.invert_yaxis();ticks=np.arange(0,56,5)
 ax.set_xticks(ticks+.5,[str(391+i) for i in ticks],rotation=45);ax.set_yticks(ticks+.5,[str(391+i) for i in ticks],rotation=45)
 ax.set(xlabel='AR residue',ylabel='AR residue',title=title)
fig,axes=plt.subplots(1,3,figsize=(21,6.5),layout='constrained')
for ax,s in zip(axes[:2],('WT','AA')):draw(ax,np.loadtxt(D/f'{s}_R2R3.csv',delimiter=','),f'{s} base pool (n=2,000)')
draw(axes[2],matrix,f'2022 apo WT MD (n={n:,} frames)')
fig.suptitle('Matched R2–R3 contact maps: AR391–446\nClosest heavy-atom distance < 1.2 nm; uniform averages',fontsize=15)
for ext in ('png','pdf'):fig.savefig(D/f'R2R3_WT_AA_MD.{ext}',dpi=160)
plt.close(fig)
fig,ax=plt.subplots(figsize=(8,6.5),layout='constrained');draw(ax,matrix,f'2022 apo WT MD — AR391–446\n{n:,} frames; closest-heavy < 1.2 nm')
for ext in ('png','pdf'):fig.savefig(D/f'MD_R2R3.{ext}',dpi=160)
plt.close(fig);print('MD map complete',n,flush=True)
