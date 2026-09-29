"""Execute the study notebook's actual plotting cell with compatibility substitutions."""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'outputs/Tau5_ASTEROIDS/Sa_Rg/repo_style';OUT.mkdir(exist_ok=True);DATA=OUT.parent
nb=json.loads((ROOT/'outputs/Tau5_joint_refinement/reference_2022/Zhu_AR_Ligands_Apo.9.27.22.ipynb').read_text())
original=''.join(nb['cells'][28]['source']);(OUT/'original_plotting_cell.py').write_text(original)
# NumPy removed normed; density=True is its replacement. Setting limits on im
# preserves the original mapping without relying on the current colorbar axes.
code=original.replace('normed=True','density=True').replace('plt.clim(vmin=0.1, vmax=3.0)','im.set_clim(vmin=0.1, vmax=3.0)')
(OUT/'compatible_plotting_cell.py').write_text(code)
metadata={}
for view,keys in [('R2R3',['WT','AA','original_apo']),('full',['WT','AA'])]:
 images=[];labels=[]
 for key in keys:
  z=np.load(DATA/('original_apo_observables.npz' if key=='original_apo' else f'{key}_observables.npz'))
  sa=z['sa'] if key=='original_apo' else z['r2r3_sa' if view=='R2R3' else 'full_sa']
  rg=z['rg'] if key=='original_apo' else z['r2r3_rg' if view=='R2R3' else 'full_rg']
  dest=OUT/f'{view}_{key}';dest.mkdir(exist_ok=True)
  plotcode=code if view=='R2R3' else code.replace('[[0.9, 2.5], [0, 25.0]]','[[1.5, 6.5], [0, 25.0]]')
  exec(plotcode,{'np':np,'plt':plt,'rg_CA':rg,'Sa_total':sa,'outdir':str(dest)+'/'})
  png=dest/'Sa_Rg.png';plt.gcf().savefig(png,dpi=160);plt.close('all');images.append(plt.imread(png))
  labels.append(('2022 apo MD' if key=='original_apo' else key+' base pool')+f' — n={len(sa):,}')
  lo,hi=(.9,2.5) if view=='R2R3' else (1.5,6.5)
  metadata[f'{view}_{key}']={'n':len(sa),'n_in_histogram_range':int(((rg>=lo)&(rg<=hi)&(sa>=0)&(sa<=25)).sum()),'Rg_range_nm':[lo,hi]}
 fig,axes=plt.subplots(1,len(images),figsize=(8*len(images),6.5))
 for ax,img,label in zip(np.atleast_1d(axes),images,labels):ax.imshow(img);ax.axis('off');ax.set_title(label,fontsize=20,pad=10)
 fig.subplots_adjust(left=0,right=1,bottom=0,top=.90,wspace=0)
 fig.savefig(OUT/f'{view}_side_by_side.png',dpi=150);fig.savefig(OUT/f'{view}_side_by_side.pdf');plt.close(fig)
metadata['method']='Exact notebook histogram/log transform: -0.001987*300*log(density+1e-6), no minimum subtraction; 30 bins; Gaussian display interpolation; jet; clim 0.1–3.0. R2R3 uses original ranges. Full chain changes only Rg histogram extent to 1.5–6.5 nm. Calculated observables retained from prior analysis.'
(OUT/'plot_metadata.json').write_text(json.dumps(metadata,indent=2));print(json.dumps(metadata,indent=2))
