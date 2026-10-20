from pathlib import Path
import json,numpy as np,pandas as pd,mdtraj as md
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
R=Path(__file__).resolve().parents[2];O=Path(__file__).resolve().parent
regions=list(range(371,382));pro=list(range(372,380))
def indices(top):
 def a(i,n):return next(a.index for a in top.residue(i).atoms if a.name==n)
 ix=[]
 for ar in regions:
  i=ar-328
  ix.extend([[a(i-1,'C'),a(i,'N'),a(i,'CA'),a(i,'C')],[a(i,'N'),a(i,'CA'),a(i,'C'),a(i+1,'N')],[a(i-1,'CA'),a(i-1,'C'),a(i,'N'),a(i,'CA')]])
 return np.array(ix)
def calc(t):return md.compute_dihedrals(t,indices(t.top),periodic=False).reshape(len(t),len(regions),3)*180/np.pi
datasets={};times={}
M=R/'outputs/WT_MD_300K_2026-10-15';t=md.load(str(M/'300K_skip5.xtc'),top=str(M/'strip.pdb'));datasets['MD WT']=calc(t);times['MD WT']=t.time/1000
for v in ['WT','AA']:
 p=R/'outputs/bAIes_8A_seqsep2_2us_2026-10-14'/v;datasets['bAIes '+v]=np.concatenate([calc(t) for t in md.iterload(str(p/'production.xtc'),top=str(p/'protein.pdb'),chunk=500)]);times['bAIes '+v]=np.arange(1,20001)*.1
 manifest=json.loads((R/f'outputs/IDPForge_D2D_MaxEnt_2026-10-15/{v}_manifest.json').read_text());datasets['IDPForge '+v]=np.concatenate([calc(md.load(m['file'])) for m in manifest])
 base=Path('/Users/paulrobustelli/Documents/Codex/2026-09-28/is-th/outputs/Ensemble-Selection/final_ensembles');p=base/('WT_shifts_plus_SAXS' if v=='WT' else 'AA_shifts_plus_SAXS')/'ensemble.pdb.gz';datasets['ASTEROIDS '+v]=calc(md.load(str(p)))
rows=[];summary=[]
fig,axs=plt.subplots(2,4,figsize=(15,8),layout='constrained');figo,axo=plt.subplots(2,4,figsize=(15,7),layout='constrained')
for ax,ao,(name,x) in zip(axs.flat,axo.flat,datasets.items()):
 ids=[regions.index(p) for p in pro];a=x[:,ids,:];cis=abs(a[:,:,2])<30;trans=abs(a[:,:,2])>150;pp=(a[:,:,0]>=-100)&(a[:,:,0]<=-50)&(a[:,:,1]>=100)&(a[:,:,1]<=180)&trans
 np.savez_compressed(O/(name.replace(' ','_')+'_angles.npz'),angles_degrees=x,AR_residues=regions,columns=['phi','psi','omega_incoming'],time_ns=times.get(name,np.array([])))
 for k,ar in enumerate(regions):
  z=x[:,k];rows.append({'ensemble':name,'AR_residue':ar,'N':len(x),'phi_circular_mean_deg':float(np.angle(np.mean(np.exp(1j*np.deg2rad(z[:,0]))),deg=True)),'psi_circular_mean_deg':float(np.angle(np.mean(np.exp(1j*np.deg2rad(z[:,1]))),deg=True)),'incoming_omega_cis_percent':float(np.mean(abs(z[:,2])<30)*100),'incoming_omega_trans_percent':float(np.mean(abs(z[:,2])>150)*100),'incoming_omega_other_percent':float(np.mean((abs(z[:,2])>=30)&(abs(z[:,2])<=150))*100)})
 summary.append({'ensemble':name,'frames':len(x),'P372_379_mean_cis_percent':cis.mean()*100,'P372_379_mean_trans_percent':trans.mean()*100,'P372_379_other_omega_percent':(~(cis|trans)).mean()*100,'all8_incoming_bonds_trans_percent':trans.all(1).mean()*100,'PPII_box_and_trans_percent':pp.mean()*100})
 ax.hist2d(a[:,:,0].ravel(),a[:,:,1].ravel(),bins=72,range=[[-180,180],[-180,180]],norm=LogNorm(),cmap='magma');ax.plot(-75,145,'cx',ms=8,mew=2);ax.set(title=name,xlabel='φ (°)',ylabel='ψ (°)',xlim=(-180,180),ylim=(-180,180));ax.set_xticks([-180,-90,0,90,180]);ax.set_yticks([-180,-90,0,90,180])
 # Map trans peak away from split ±180 boundary: display −90 to270.
 om=(a[:,:,2].ravel()+90)%360-90;ao.hist(om,bins=np.linspace(-90,270,91),density=True,color='#617cad');ao.axvline(0,color='red',ls='--',lw=.7);ao.axvline(180,color='black',ls='--',lw=.7);ao.set(title=f'{name}\ncis {cis.mean()*100:.2f}%, trans {trans.mean()*100:.2f}%',xlabel='Incoming peptide ω (°)',ylabel='Density',xlim=(-90,270));ao.set_xticks([-90,0,90,180,270])
axs.flat[-1].axis('off');axo.flat[-1].axis('off');axs.flat[-1].text(0,.8,'AR372–379: PPPPPPPP\nCyan ×: ideal PPII (−75°, +145°)\nLog counts within each panel\nUnweighted ensembles\nASTEROIDS: final shifts+SAXS selection',va='top');axo.flat[-1].text(0,.8,'ω belongs to the bond entering each Pro:\nCA(i−1)–C(i−1)–N(i)–CA(i)\ncis: |ω| <30°\ntrans: |ω| >150°\nOther: remaining angles\n+180° and −180° are equivalent',va='top')
fig.suptitle('Polyproline stretch Ramachandran distributions');figo.suptitle('Polyproline peptide-bond cis/trans distributions');fig.savefig(O/'Ramachandran.png',dpi=180);fig.savefig(O/'Ramachandran.pdf');figo.savefig(O/'omega.png',dpi=180);figo.savefig(O/'omega.pdf');plt.close('all')
pd.DataFrame(rows).to_csv(O/'per_residue_angles.csv',index=False);pd.DataFrame(summary).to_csv(O/'ensemble_summary.csv',index=False);print(pd.DataFrame(summary).to_string(index=False),flush=True)
# Chronological cis-state switches: observations, not equilibrium isomerization rates.
switch=[]
for name in times:
 x=datasets[name]
 for ar in pro:
  om=x[:,regions.index(ar),2];state=np.where(abs(om)<30,0,np.where(abs(om)>150,1,-1));valid=state[state>=0];switch.append({'ensemble':name,'AR_residue':ar,'observed_cis_trans_switches':int(np.sum(valid[1:]!=valid[:-1]))})
pd.DataFrame(switch).to_csv(O/'observed_switches.csv',index=False)
