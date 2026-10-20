from pathlib import Path
import json
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
O=Path(__file__).resolve().parent;R=O.parents[1];A=O/'comparison_analysis';OLD=R/'outputs/bAIes_8A_seqsep2_2us_2026-10-14';bases={'Jeffreys':OLD,'Harmonic':O};prefix={'Jeffreys':'Jeffreys','Harmonic':'harmonic'};colors={'Jeffreys':'#777777','Harmonic':'#0072b2'};rows=[];data={};helices={};maps={};restraints={}
for name,base in bases.items():
 for v in ['WT','AA']:
  f=base/'final_region_analysis';df=pd.read_csv(f/f'{v}_timeseries.csv').merge(pd.read_csv(f/f'{v}_Sa_timeseries.csv'),on=['frame_index','production_ns'],validate='one_to_one');assert len(df)==20000;data[name,v]=df;helices[name,v]=pd.read_csv(f/f'{v}_residue_helicity.csv');maps[name,v]=np.loadtxt(f/f'{v}_contact_map.csv',delimiter=',');rs=pd.read_csv(f/f'{v}_restraint_agreement.csv');restraints[name,v]=rs
  row={'method':name,'variant':v,'frames':len(df),'Rg_CA_A':df.Rg_CA_A.mean(),'Rg_min_CA_A':df.Rg_CA_A.min(),'fraction_Rg_lt19':(df.Rg_CA_A<19).mean(),'Sa_full':df.Sa_full.mean(),'Sa_R2R3':df.Sa_AR391_446.mean(),'restraint_occupancy_percent':rs.fraction_within_mu_sigma.mean()*100,'mean_restraint_squared_z':np.mean((rs.SD_nm**2+(rs.mean_nm-rs.mu_nm)**2)/rs.sigma_nm**2)}
  for col in ['R2_main','R2_minor','R2_historical_core','R3']:row[col+'_H_percent']=df[col+'_DSSP_H'].mean()*100
  for col in ['Rg_CA_A','R2_main_DSSP_H','R3_DSSP_H']:
   vals=df[col].to_numpy();blocks=vals.reshape(10,2000).mean(1);row[col+'_first1us']=vals[:10000].mean();row[col+'_last1us']=vals[10000:].mean();row[col+'_block_min']=blocks.min();row[col+'_block_max']=blocks.max()
  rows.append(row)
pd.DataFrame(rows).to_csv(A/'likelihood_comparison.csv',index=False)
fig,axes=plt.subplots(2,2,figsize=(13,8),layout='constrained')
for j,v in enumerate(['WT','AA']):
 for name in bases:
  h=helices[name,v];df=data[name,v];axes[0,j].plot(h.AR_residue,h.DSSP_H*100,color=colors[name],label=name,lw=1.6);axes[1,j].hist(df.Rg_CA_A,bins=np.arange(15,75,1),density=True,histtype='step',lw=1.8,color=colors[name],label=f'{name}: mean {df.Rg_CA_A.mean():.2f} Å');axes[1,j].axvline(df.Rg_CA_A.mean(),color=colors[name],ls=':')
 for lo,hi in [(390,404),(404,410),(430,440)]:axes[0,j].axvspan(lo,hi,color='grey',alpha=.12)
 axes[0,j].set(title=v,xlabel='AR residue',ylabel='DSSP H (%)',ylim=(0,100));axes[0,j].legend();axes[1,j].set(xlabel='Cα Rg (Å)',ylabel='Probability density');axes[1,j].legend()
fig.suptitle('Harmonic versus Jeffreys bAIes • full 2 μs trajectories\nSame 8 Å / seqsep2 selection; native120-residue construct; no experimental reweighting')
for ext in ['png','pdf']:fig.savefig(A/f'helicity_Rg_comparison.{ext}',dpi=180)
plt.close(fig)
fig,axes=plt.subplots(3,2,figsize=(12,9),layout='constrained')
for j,v in enumerate(['WT','AA']):
 for i,(col,scale,label) in enumerate([('Rg_CA_A',1,'Cα Rg (Å)'),('R2_main_DSSP_H',100,'R2 main H (%)'),('Sa_full',1,'Sα full')]):
  ax=axes[i,j]
  for name in bases:
   df=data[name,v];x=df[col].to_numpy()*scale;b=x.reshape(10,2000).mean(1);ax.plot(df.production_ns,pd.Series(x).rolling(200,min_periods=1).mean(),lw=.7,alpha=.7,color=colors[name]);ax.plot(np.arange(100,2000,200),b,'o-',color=colors[name],label=name,ms=4)
  ax.set(title=v,xlabel='Production time (ns)',ylabel=label);ax.legend(fontsize=8)
fig.suptitle('Within-trajectory block comparison • dots are200-ns means; thin lines20-ns rolling means\nOne trajectory per condition; these do not establish independent-replicate convergence')
for ext in ['png','pdf']:fig.savefig(A/f'time_block_comparison.{ext}',dpi=170)
plt.close(fig)
# Compare each same restraint pair directly.
fig,axes=plt.subplots(2,2,figsize=(14,9),layout='constrained');pairrows=[]
for j,v in enumerate(['WT','AA']):
 old=restraints['Jeffreys',v];new=restraints['Harmonic',v];assert np.array_equal(old[['AR_i','AR_j','mu_nm','sigma_nm']].to_numpy(),new[['AR_i','AR_j','mu_nm','sigma_nm']].to_numpy());xx=np.arange(len(old));labels=[f'{i}–{k}' for i,k in zip(old.AR_i,old.AR_j)]
 for name,rs in [('Jeffreys',old),('Harmonic',new)]:
  axes[0,j].plot(xx,rs.fraction_within_mu_sigma*100,'o-',ms=3,color=colors[name],label=name);axes[1,j].plot(xx,(rs.mean_nm-rs.mu_nm)/rs.sigma_nm,'o-',ms=3,color=colors[name],label=name)
 for i in [0,1]:axes[i,j].set_xticks(xx,labels,rotation=90,fontsize=6);axes[i,j].legend();axes[i,j].set_title(v);axes[i,j].grid(axis='y',alpha=.15)
 axes[0,j].set(ylabel='Frames within μ ± σ (%)',ylim=(0,100));axes[1,j].set(ylabel='Mean distance deviation / σ',xlabel='Restrained AR pair');axes[1,j].axhline(0,color='k',lw=.7)
 x=new.copy();x['variant']=v;x['Jeffreys_occupancy']=old.fraction_within_mu_sigma;x['Harmonic_minus_Jeffreys_occupancy']=new.fraction_within_mu_sigma-old.fraction_within_mu_sigma;pairrows.append(x)
pd.concat(pairrows).to_csv(A/'paired_restraint_comparison.csv',index=False)
fig.suptitle('Same45 WT /37 AA restraints • harmonic versus Jeffreys\nμ ± σ occupancy is descriptive agreement with soft restraints, not a pass/fail test')
for ext in ['png','pdf']:fig.savefig(A/f'restraint_comparison.{ext}',dpi=170)
plt.close(fig)
# Full and R2R3 maps, paper linear style plus log weak-contact view.
for scale in ['linear','log']:
 fig,axes=plt.subplots(2,4,figsize=(17,9),layout='constrained')
 for row,v in enumerate(['WT','AA']):
  for col,(name,crop) in enumerate([('Jeffreys',False),('Harmonic',False),('Jeffreys',True),('Harmonic',True)]):
   m=maps[name,v];m=m[63:119,63:119] if crop else m;ax=axes[row,col];kwargs={'norm':LogNorm(5e-5,1),'cmap':'viridis'} if scale=='log' else {'vmin':0,'vmax':1,'cmap':'jet'};im=ax.imshow(m,origin='upper',**kwargs)
   if crop:ax.set_xticks([0,19,39,55],[391,410,430,446]);ax.set_yticks([0,19,39,55],[391,410,430,446])
   else:ax.set_xticks([2,22,42,62,82,102,119],[330,350,370,390,410,430,447]);ax.set_yticks([2,22,42,62,82,102,119],[330,350,370,390,410,430,447])
   ax.set(title=f'{v} {name} • '+('R2R3' if crop else 'full + GP'),xlabel='AR residue',ylabel='AR residue')
 fig.colorbar(im,ax=axes.ravel().tolist(),shrink=.7,label='Closest-heavy contact probability (<12 Å)');fig.suptitle('Full-trajectory contact maps • '+scale+' scale • adjacent-residue band retained')
 for ext in ['png','pdf']:fig.savefig(A/f'contact_comparison_{scale}.{ext}',dpi=170)
 plt.close(fig)
# Reproduce prior paper-style Sα/Rg transformation with common bins.
coverage=[]
for crop in [False,True]:
 fig,axes=plt.subplots(2,2,figsize=(11,9),layout='constrained')
 for row,v in enumerate(['WT','AA']):
  for col,name in enumerate(bases):
   df=data[name,v];x=df['Rg_R2R3_CA_A' if crop else 'Rg_CA_A']/10;y=df['Sa_AR391_446' if crop else 'Sa_full'];xr=(.9,2.5) if crop else (1.5,6.5);hist,xe,ye=np.histogram2d(x,y,bins=30,range=[xr,(0,25)],density=True);z=-.001987*300*np.log(hist.T+1e-6);ax=axes[row,col];im=ax.imshow(np.flipud(z),extent=[*xr,0,25],aspect='auto',interpolation='gaussian',cmap='jet',vmin=.1,vmax=3)
   frac=float(((x>=xr[0])&(x<=xr[1])&(y>=0)&(y<=25)).mean());coverage.append({'variant':v,'method':name,'scope':'R2R3' if crop else 'full','retained_frame_fraction':frac});ax.set(title=f'{v} {name} • displayed {frac:.1%}',xlabel='Cα Rg (nm)',ylabel='Sα')
 fig.colorbar(im,ax=axes.ravel().tolist(),label='−RT ln(density + 10⁻⁶), display units',shrink=.8);fig.suptitle(('AR391–446' if crop else 'Full120-residue')+' Sα–Rg • original paper-style bins/color scale\nHistogram display transform, not an equilibrium free-energy estimate')
 for ext in ['png','pdf']:fig.savefig(A/f'Sa_Rg_{"R2R3" if crop else "full"}.{ext}',dpi=180)
 plt.close(fig)
pd.DataFrame(coverage).to_csv(A/'Sa_Rg_plot_coverage.csv',index=False)
# Packing distributions from evenly sampled native120-residue frames.
pr=[];fig,axes=plt.subplots(2,3,figsize=(14,8),layout='constrained')
for j,v in enumerate(['WT','AA']):
 for name in bases:
  d=pd.read_csv(A/f'{prefix[name]}_{v}_packing.csv');pr.append({'method':name,'variant':v,'frames':len(d),**d.select_dtypes('number').drop(columns=['frame_index','production_ns']).mean().to_dict()})
  for col,(f,label) in enumerate([('buried_HSA_nonlocal_A2','Nonlocal buried hydrophobic area (Å²)'),('hydrophobic_pairs_sep6_50','Hydrophobic contact pairs'),('aromatic_hydrophobic_pairs_sep6_50','Aromatic–hydrophobic contact pairs')]):
   bins=np.arange(0,510,10) if col==0 else np.arange(-.5,16.5,1);axes[j,col].hist(d[f],bins=bins,weights=np.ones(len(d))/len(d),histtype='step',lw=1.6,color=colors[name],label=f'{name}: mean {d[f].mean():.2f}');axes[j,col].set(title=v,xlabel=label,ylabel='Fraction of sampled frames');axes[j,col].legend(fontsize=8)
fig.suptitle('Matched native120-residue packing comparison •1,000 evenly spaced frames per trajectory\nContacts<5 Å, |i−j|>6; buried area relative to unchanged local ±6-residue fragments')
for ext in ['png','pdf']:fig.savefig(A/f'packing_comparison.{ext}',dpi=180)
plt.close(fig);pd.DataFrame(pr).to_csv(A/'packing_comparison.csv',index=False)
# Cis populations and actual cis/trans well changes (exclude intermediate angles).
cs=pd.read_csv(A/'cis_residue_summary.csv');nonzero=cs[cs.cis_percent>0];arlist=sorted(nonzero.following_residue.astype(int).unique());fig,axes=plt.subplots(2,2,figsize=(13,8),layout='constrained');well=[]
for row,v in enumerate(['WT','AA']):
 for j,name in enumerate(bases):
  z=np.load(A/f'{prefix[name]}_{v}_cis_states.npz');states=z['states'];ars=z['AR_following'];pp=(ars>=372)&(ars<=379)
  vals=[]
  for ar in arlist:
   k=np.where(ars==ar)[0][0];s=states[:,k];valid=s[s>=0];n=int(np.sum(valid[1:]!=valid[:-1]));vals.append((s==1).mean()*100);well.append({'method':name,'variant':v,'AR_residue':int(ar),'cis_percent':vals[-1],'cis_trans_well_changes':n,'note':'Intermediate angles excluded; descriptive observed changes, not rate estimates'})
  axes[row,0].bar(np.arange(len(arlist))+(j-.5)*.36,vals,width=.36,color=colors[name],label=name)
  y=(states[:,pp]==1).mean(1)*100;axes[row,1].plot(np.arange(100,2000,200),y.reshape(10,2000).mean(1),'o-',color=colors[name],label=name)
 axes[row,0].set(xticks=range(len(arlist)),xticklabels=arlist,ylabel='Incoming peptide bond cis (%)',xlabel='AR proline residue',title=v,ylim=(0,105));axes[row,0].legend();axes[row,1].set(title=v,xlabel='Production time (ns)',ylabel='Mean cis fraction across Pro372–379 (%)');axes[row,1].legend()
fig.suptitle('Persistent cis peptide bonds • cis |ω|<30°; trans |ω|>150°\nTime-block means do not establish equilibrium cis populations')
for ext in ['png','pdf']:fig.savefig(A/f'cis_comparison.{ext}',dpi=180)
plt.close(fig);pd.DataFrame(well).to_csv(A/'cis_well_change_summary.csv',index=False)
print(pd.DataFrame(rows).to_string(index=False))
