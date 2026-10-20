from analyze import *
import hashlib,html
D=O/'fit_experiments';reports=[];checks=[]
rest,rw,rr,q,y,sy=load('REST_WT')
for p in sorted(D.glob('*/complete.json')):
 a=json.loads(p.read_text());name=a['base_pool'];d,w,raw,_,_,_=load(name);met=d.copy()
 if 'source_indices' in a:raw=raw[:,a['source_indices']];met=d.iloc[a['source_indices']].copy()
 elif 'added_REST_indices' in a:raw=np.c_[raw,rr[:,a['added_REST_indices']]];met=pd.concat([d,rest.iloc[a['added_REST_indices']]],ignore_index=True)
 ww=pd.read_csv(p.parent/'weights.csv').weight.to_numpy();stored=pd.read_csv(p.parent/'curve.csv');pred=a['scale']*np.einsum('ij,j->i',raw,ww,optimize=False)
 assert np.isfinite(pred).all() and np.allclose(pred,stored.prediction,atol=1e-12,rtol=1e-10)
 neff=1/np.einsum('i,i->',ww,ww,optimize=False);assert np.isclose(neff,a['Neff'])
 aa={k:v for k,v in a.items() if k not in ['source_indices','added_REST_indices','gradient_checks']}
 aa['Rg_CA_A']=float(np.einsum('i,i->',ww,met.Rg_CA_A));aa['weight_Rg_lt19']=float(ww[met.Rg_CA_A<19].sum())
 for f in FIELDS:aa[f]=float(np.einsum('i,i->',ww,met[f]))
 reports.append(aa);checks.append({'experiment':a['label'],'max_curve_absolute_difference':float(abs(pred-stored.prediction).max()),'Neff_independent_einsum':float(neff),'finite':True})
result=pd.DataFrame(reports);result.to_csv(O/'diagnostic_fit_results.csv',index=False);(O/'independent_fit_validation.json').write_text(json.dumps(checks,indent=2))
# Clearly separate published combined fits from the new SAXS-only sensitivity fits.
fig,axes=plt.subplots(2,3,figsize=(15,8),layout='constrained')
for col,name in enumerate(['IDPForge_WT','bAIes_WT','ASTEROIDS_WT']):
 for suffix,label,color in [('_baseline','Original pool','#777777'),('_add_compact_REST','+160 REST Rg<19 Å','#0072b2'),('_add_expanded_REST','+160 REST Rg≥25 Å','#d55e00')]:
  a=result[result.label==name+suffix].iloc[0];curve=pd.read_csv(D/(name+suffix)/'curve.csv')
  axes[0,col].plot(curve.q,curve.prediction/curve.I_exp,color=color,label=f'{label}\nmid χ²={a.chi2_mid:.2f}')
 axes[0,col].axhline(1,color='k',lw=.7);axes[0,col].axvspan(.08,.25,color='grey',alpha=.1);axes[0,col].set(title=name.replace('_',' '),xlabel='q (Å⁻¹)',ylabel='Prediction / experiment',ylim=(.3,1.5));axes[0,col].legend(fontsize=8)
 for j,suffix in enumerate(['_baseline','_add_compact_REST','_add_expanded_REST']):
  a=result[result.label==name+suffix].iloc[0];axes[1,col].bar(j,a.chi2_mid,color=['#777777','#0072b2','#d55e00'][j]);axes[1,col].text(j,a.chi2_mid*1.13,f'{a.chi2_mid:.2f}',ha='center')
 axes[1,col].set(yscale='log',ylim=(.4,350),xticks=[0,1,2],xticklabels=['Original','+ compact','+ expanded'],ylabel='Intermediate-q χ² per point');axes[1,col].axhline(1,color='k',ls=':',lw=.7)
fig.suptitle('WT pool-addition experiments • SAXS-only Kish MaxEnt\nSame 160 added structures per comparison; one global intensity scale; original absolute Kish support retained',fontsize=12)
for ext in ['png','pdf']:fig.savefig(O/f'WT_pool_addition_experiments.{ext}',dpi=180)
plt.close(fig)
fig,ax=plt.subplots(1,2,figsize=(12,4.5),layout='constrained')
for label,col in [('REST_WT_baseline','#333333'),('REST_WT_remove_Rg_lt19','#0072b2'),('REST_WT_remove_Rg_lt22','#009e73'),('REST_WT_remove_Rg_lt24','#d55e00')]:
 a=result[result.label==label].iloc[0];curve=pd.read_csv(D/label/'curve.csv');name='All REST' if label.endswith('baseline') else 'Remove Rg < '+label.split('lt')[-1]+' Å';ax[0].plot(curve.q,curve.prediction/curve.I_exp,color=col,label=f'{name}: mid χ²={a.chi2_mid:.2f}')
ax[0].axhline(1,color='k',lw=.7);ax[0].axvspan(.08,.25,color='grey',alpha=.1);ax[0].legend(fontsize=8);ax[0].set(xlabel='q (Å⁻¹)',ylabel='Prediction / experiment',ylim=(.6,1.4),title='Removing compact REST frames')
labels=['baseline','low_burial_Rgmatched','high_burial_Rgmatched','low_aromatic_Rgmatched','high_aromatic_Rgmatched'];values=[result[result.label=='REST_WT_'+l].iloc[0].chi2_mid for l in labels];ax[1].bar(range(5),values,color=['#333333','#56b4e9','#0072b2','#e69f00','#d55e00']);ax[1].set(xticks=range(5),xticklabels=['All','Low burial','High burial','Low aromatic','High aromatic'],ylabel='Intermediate-q χ² per point',title='REST subsets matched in 2-Å Rg bins');ax[1].tick_params(axis='x',rotation=25);ax[1].axhline(1,color='k',ls=':')
fig.suptitle('SAXS-only removal tests • absolute Kish Neff ≥186.6 in every fit')
for ext in ['png','pdf']:fig.savefig(O/f'REST_removal_experiments.{ext}',dpi=180)
plt.close(fig)
# Coordination vs buried area, colored by size: one row per variant.
fig,ax=plt.subplots(2,3,figsize=(15,8),layout='constrained')
colors={'REST':'#cc79a7','IDPForge':'#0072b2','bAIes':'#d55e00','ASTEROIDS':'#009e73'};summary=[]
for name in NAMES:
 d,w,raw,q,y,sy=load(name);v=name.split('_')[-1];method=name.rsplit('_',1)[0];row=0 if v=='WT' else 1
 for col,f in enumerate(['hydrophobic_mean_coordination','aromatic_hydrophobic_pairs_sep6_50','all_tertiary_pairs_sep6_50']):
  if col==2:
   tc=pd.read_csv(O/name/'tertiary_contacts.csv');d[f]=tc[f]
  # Binned means to keep the all-method comparison legible; raw per-frame data retained.
  bins=np.arange(0,2601,100);xx=[];yy=[]
  for lo,hi in zip(bins[:-1],bins[1:]):
   m=(d.buried_HSA_nonlocal_A2>=lo)&(d.buried_HSA_nonlocal_A2<hi)
   if m.sum()>=5:xx.append(d.loc[m,'buried_HSA_nonlocal_A2'].mean());yy.append(d.loc[m,f].mean())
  ax[row,col].plot(xx,yy,'o-',ms=3,color=colors[method],label=method);ax[row,col].set(xlabel='Nonlocal buried hydrophobic area (Å²)',ylabel=['Mean hydrophobic coordination','Aromatic–hydrophobic contact pairs','All tertiary contact pairs'][col],title=v);ax[row,col].legend(fontsize=8)
 for stage in ['uniform','combined']:
  summary.append({'ensemble':name,'stage':stage,**{c:float(np.einsum('i,i->',w[stage],tc[c])) for c in tc.columns}})
pd.DataFrame(summary).to_csv(O/'tertiary_contacts_before_after.csv',index=False)
for ext in ['png','pdf']:fig.savefig(O/f'coordination_vs_buried_HSA.{ext}',dpi=180)
plt.close(fig)
# General tertiary-contact conditional curves (5 Å closest-heavy, |i-j|>6).
fig,axes=plt.subplots(2,4,figsize=(16,8),layout='constrained');rows=[]
for ax,name in zip(axes.flat,NAMES):
 d,w,raw,q,y,sy=load(name);tc=pd.read_csv(O/name/'tertiary_contacts.csv');f=tc.all_tertiary_pairs_sep6_50;cuts=np.unique(np.quantile(f,[1/3,2/3]));g=np.digitize(f,cuts)
 for k in sorted(set(g)):
  m=g==k;scale,met,pred=score(raw[:,m].mean(1),q,y,sy);ax.plot(q,pred/y,label=f'Group {k+1}, n={sum(m)}, mid χ²={met["chi2_mid"]:.1f}');rows.append({'ensemble':name,'group':k,'N':int(sum(m)),'mean_tertiary_pairs':float(f[m].mean()),**met})
  pd.DataFrame({'q':q,'raw_intensity':raw[:,m].mean(1),'scaled_intensity':pred}).to_csv(O/name/f'tertiary_group{k}_SAXS.csv',index=False)
 ax.axhline(1,color='k',lw=.7);ax.axvspan(.08,.25,color='grey',alpha=.1);ax.set(title=name.replace('_',' '),xlabel='q (Å⁻¹)',ylabel='Prediction / experiment',ylim=(.2,1.8));ax.legend(fontsize=7)
axes.flat[-1].axis('off');fig.suptitle('SAXS split by general tertiary-contact count\nClosest heavy atoms <5 Å, sequence separation >6; within-pool tertiles; one full-q ensemble scale per group')
for ext in ['png','pdf']:fig.savefig(O/f'tertiary_contact_SAXS.{ext}',dpi=160)
pd.DataFrame(rows).to_csv(O/'tertiary_group_scores.csv',index=False)
print(result[['label','chi2_all','chi2_mid','Rg_CA_A','buried_HSA_nonlocal_A2','weight_Rg_lt19']].to_string(index=False))
