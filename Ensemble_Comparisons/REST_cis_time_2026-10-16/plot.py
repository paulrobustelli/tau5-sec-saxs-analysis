from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd,mdtraj as md
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(__file__).resolve().parents[2];O=Path(__file__).resolve().parent;D=R/'outputs/WT_MD_300K_2026-10-15'
t=md.load(str(D/'300K_skip5.xtc'),top=str(D/'strip.pdb'))
ix,omega=md.compute_omega(t,periodic=False);deg=np.rad2deg(omega);cis=np.abs(deg)<30
res=[t.top.atom(int(i[2])).residue for i in ix];ar=np.array([r.index+328 for r in res]);pro=np.array([r.name=='PRO' for r in res]);poly=(ar>=372)&(ar<=379)
time=t.time/1000
assert len(t)==1866 and np.all(np.diff(time)>0) and len(ix)==120
old=np.load(R/'outputs/PPII_dihedrals_2026-10-15/MD_WT_angles.npz'); oldcis=np.abs(old['angles_degrees'][:,np.isin(old['AR_residues'],np.arange(372,380)),2])<30
assert np.array_equal(cis[:,poly],oldcis)
counts=pd.DataFrame({'time_ns':time,'all_peptide_bonds_cis':cis.sum(1),'X_Pro_bonds_cis':cis[:,pro].sum(1),'AR372_379_incoming_bonds_cis':cis[:,poly].sum(1)})
counts.to_csv(O/'REST_WT_cis_counts_vs_time.csv',index=False)
np.savez_compressed(O/'REST_WT_omega_states.npz',time_ns=time,incoming_AR_residue=ar,omega_degrees=deg,cis=cis,is_proline=pro)
fig,axs=plt.subplots(2,1,figsize=(12,6.5),sharex=True,layout='constrained')
axs[0].step(time,counts.all_peptide_bonds_cis,where='mid',color='#333333',lw=.8,label='All 120 peptide bonds')
axs[0].step(time,counts.AR372_379_incoming_bonds_cis,where='mid',color='#0072b2',lw=.65,alpha=.8,label='Eight bonds entering Pro372–379')
axs[0].set(ylabel='Number of cis bonds per frame',title='WT REST 300 K trajectory — cis peptide bonds')
axs[0].set_yticks(range(int(cis.sum(1).max())+1));axs[0].legend(loc='upper right',frameon=True);axs[0].grid(axis='y',alpha=.2)
active=np.flatnonzero(cis.any(0));labels=[]
for j,k in enumerate(active):
 axs[1].scatter(time[cis[:,k]],np.full(cis[:,k].sum(),j),s=3,marker='|',color='#0072b2',linewidths=.8)
 labels.append('Tag Pro (GP)' if res[k].index==1 else f'{res[k].name.capitalize()}{ar[k]}')
axs[1].set_yticks(range(len(active)),labels);axs[1].set(ylabel='Residue following cis bond',xlabel='Saved REST trajectory time (ns)',xlim=(time[0],time[-1]))
fig.get_layout_engine().set(rect=(0,.06,1,.94))
fig.text(.5,.015,'Cis: |ω| < 30°. Unweighted saved frames (0.4 ns spacing). REST exchanges mean frame-to-frame changes are not physical kinetics.',ha='center',fontsize=9)
for ext in ['png','pdf']:fig.savefig(O/f'REST_WT_cis_vs_time.{ext}',dpi=200)
report={'frames':len(t),'peptide_bonds':len(ix),'time_range_ns':[float(time[0]),float(time[-1])],'mean_cis_bonds':float(cis.sum(1).mean()),'max_cis_bonds':int(cis.sum(1).max()),'residues_with_cis':labels,'definition':'abs(omega)<30 degrees; bond labeled by following residue; GP tags precede AR330','polyproline_matches_previous_audit':True,'sources':{str(D/f):hashlib.sha256((D/f).read_bytes()).hexdigest() for f in ['strip.pdb','300K_skip5.xtc']}}
(O/'provenance.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
