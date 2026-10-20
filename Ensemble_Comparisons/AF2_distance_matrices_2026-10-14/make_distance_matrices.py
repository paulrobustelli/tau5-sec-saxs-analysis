"""Distogram distance proxies and exact short-distance probabilities from AF2 logits."""
from pathlib import Path
import json,pickle,hashlib,re
import numpy as np,pandas as pd
from scipy.special import softmax
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
ROOT=Path(__file__).resolve().parents[2];OUT=Path(__file__).resolve().parent
DATA=ROOT/'work/af2_original_dropbox';regions=[('R1',350,360),('R2 main',390,404),('R2 minor',404,410),('R3',430,440)]
results={};provenance={};labels=np.arange(330,448)
for variant in ['WT','AA']:
 files=sorted(p for p in DATA.rglob(f'*{variant}*all_rank_*.pickle') if not p.name.startswith('._'))
 assert len(files)==5
 proxies=[];contacts=[];overflow=[];rank1=None;hashes=[]
 for path in files:
  d=pickle.load(path.open('rb'))['distogram'];p=softmax(d['logits'].astype(np.float64),axis=-1);edges=d['bin_edges'].astype(float)
  centers=edges-(edges[1]-edges[0])/2
  assert p.shape==(120,120,64) and len(centers)==63
  proxy=np.einsum('ijk,k->ij',p,np.r_[centers,22.0]);over=p[:,:,-1]
  finite=np.einsum('ijk,k->ij',p[:,:,:63],centers)/(1-over)
  # Probability integral at8A includes partial boundary-bin mass, assuming
  # uniform density within that bin; exact sum of complete bins also saved.
  lower=np.r_[2.0,edges[:-1]];upper=edges
  fraction=np.clip((8.0-lower)/(upper-lower),0,1)
  contact=np.einsum('ijk,k->ij',p[:,:,:63],fraction)
  complete=np.sum(p[:,:,:63][:,:,upper<=8.0],axis=-1)
  m={'mean_distance_22A_proxy':proxy,'finite_bin_conditional_mean':finite,'probability_under_8A_interpolated':contact,'probability_complete_bins_below_8A':complete,'overflow_probability':over}
  proxies.append(proxy);contacts.append(contact);overflow.append(over)
  if 'rank_001' in path.name:rank1=m
  hashes.append({'name':path.name,'sha256':hashlib.sha256(path.read_bytes()).hexdigest()})
 assert rank1 is not None
 m=dict(rank1);m.update(five_model_mean_22A_proxy=np.mean(proxies,axis=0),five_model_probability_under_8A=np.mean(contacts,axis=0),five_model_overflow_probability=np.mean(overflow,axis=0))
 np.savez_compressed(OUT/f'{variant}_full_120res_matrices.npz',**m)
 for key,val in m.items():pd.DataFrame(val[2:,2:],index=labels,columns=labels).to_csv(OUT/f'{variant}_{key}_AR330_447.csv',index_label='AR_residue')
 results[variant]=m
 provenance[variant]={'inputs':hashes,'rank1_model':'alphafold2_ptm_model_2_seed_000','mapping':'zero-based construct index2 = AR330; GP tag excluded from plotted matrices','open_ended_bin':'last bin starts at21.6875A; assigned22A solely for mean-distance proxy, not a true full mean','finite_bin_centers':'same midpoint convention as supplied bAIes preprocessing','contact_probability':'P(d<8A); boundary bin interpolated under uniform-in-bin assumption; complete-bin sum also saved','five_model_average':'equal average of five prediction probability distributions, not physical ensemble populations'}

def decorate(ax):
 ticks=[330,350,375,400,425,447];ax.set_xticks(ticks);ax.set_yticks(ticks);ax.set_xlabel('AR residue');ax.set_ylabel('AR residue')
 for name,a,b in regions:
  ax.add_patch(Rectangle((a-.5,a-.5),b-a+1,b-a+1,fill=False,edgecolor='black',lw=1.1))
  label_x,label_y=(b+2,a-2) if name!='R2 minor' else (b+2,b+3)
  ax.text(label_x,label_y,name,fontsize=8,color='black',bbox=dict(facecolor='white',edgecolor='none',alpha=.8,pad=.6))
 ax.set_xlim(329.5,447.5);ax.set_ylim(447.5,329.5)

def draw_pair(key,title,cmap,vmin,vmax,filename,colorbar):
 fig,axs=plt.subplots(1,2,figsize=(11,5.4),layout='constrained')
 for ax,v in zip(axs,['WT','AA']):
  a=results[v][key][2:,2:].copy();np.fill_diagonal(a,np.nan)
  im=ax.imshow(a,extent=(329.5,447.5,447.5,329.5),cmap=cmap,vmin=vmin,vmax=vmax,interpolation='nearest');decorate(ax);ax.set_title(v+' · rank-1 AF2' if not key.startswith('five') else v+' · five-model average')
 cb=fig.colorbar(im,ax=axs,shrink=.8,pad=.02);cb.set_label(colorbar);fig.suptitle(title,fontsize=13)
 fig.savefig(OUT/(filename+'.png'),dpi=190);fig.savefig(OUT/(filename+'.pdf'));plt.close(fig)
draw_pair('mean_distance_22A_proxy','AF2 distogram average-distance proxy\nOpen-ended >21.69 Å bin represented by 22 Å','viridis',3,22,'WT_AA_distogram_mean_distance','Mean distance proxy (Å); lower = closer')
draw_pair('probability_under_8A_interpolated','AF2 probability of a Cβ–Cβ distance below 8 Å\nCα used for glycine; no restraint-selection filter applied','magma',0,1,'WT_AA_AF2_contact_probability','P(distance < 8 Å)')
draw_pair('overflow_probability','Probability in the unresolved >21.69 Å bin','Greys',0,1,'WT_AA_overflow_probability','Open-ended-bin probability')
draw_pair('five_model_mean_22A_proxy','AF2 distance proxy averaged over all five predictions\nOpen-ended >21.69 Å bin represented by 22 Å','viridis',3,22,'WT_AA_five_model_mean_distance','Mean distance proxy (Å)')
(OUT/'provenance.json').write_text(json.dumps(provenance,indent=2))
print('Saved matrices/figures to',OUT)
