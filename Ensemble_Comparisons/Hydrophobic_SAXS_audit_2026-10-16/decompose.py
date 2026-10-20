from analyze import *
D=O/'fit_experiments';d,w,raw,q,y,sy=load('REST_WT');tc=pd.read_csv(O/'REST_WT/tertiary_contacts.csv');rows=[]
fig,axes=plt.subplots(2,2,figsize=(13,9),layout='constrained');a=json.loads((D/'REST_WT_baseline/complete.json').read_text());ww=pd.read_csv(D/'REST_WT_baseline/weights.csv').weight.to_numpy();scale=a['scale']
for ax,f in zip(axes.flat,['Rg_CA_A','buried_HSA_nonlocal_A2','aromatic_hydrophobic_pairs_sep6_50','all_tertiary_pairs_sep6_50']):
 x=d[f] if f in d else tc[f];cuts=[19,25,35] if f=='Rg_CA_A' else np.unique(np.quantile(x,[1/3,2/3]));g=np.digitize(x,cuts)
 pieces=[]
 for k in sorted(set(g)):
  m=g==k;piece=scale*np.einsum('ij,j->i',raw[:,m],ww[m]);pieces.append(piece)
  lo='−∞' if k==0 else f'{cuts[k-1]:.0f}';hi='∞' if k==len(cuts) else f'{cuts[k]:.0f}'
  ax.plot(q,piece,label=f'{lo} to {hi}: {ww[m].sum()*100:.1f}% weight')
  rows.extend({'split':f,'group':int(k),'q':float(qi),'weighted_intensity_contribution':float(pi)} for qi,pi in zip(q,piece))
 total=np.sum(pieces,0);assert np.allclose(total,scale*np.einsum('ij,j->i',raw,ww));ax.plot(q,total,'k-',lw=2,label='Total fitted curve');ax.errorbar(q,y,yerr=sy,fmt='.',color='grey',ms=2,alpha=.5,label='Experiment');ax.axvspan(.08,.25,color='grey',alpha=.08);ax.set(yscale='log',xlabel='q (Å⁻¹)',ylabel='Intensity contribution',title=f.replace('_',' '));ax.legend(fontsize=7)
fig.suptitle('REST SAXS-only fit decomposed by conformation type\nOne shared fitted scale: colored contributions sum exactly to black total; no group rescaling')
for ext in ['png','pdf']:fig.savefig(O/f'REST_fitted_SAXS_components.{ext}',dpi=180)
pd.DataFrame(rows).to_csv(O/'REST_fitted_SAXS_components.csv',index=False)
# Replicate conditional group curves on the SAME parent-prior scale as an alternative view.
c=pd.read_csv(O/'conditional_SAXS_curves.csv');fig,axes=plt.subplots(7,4,figsize=(18,23),layout='constrained')
for row,name in enumerate(NAMES):
 dd,ww,rr,q,y,sy=load(name);scale,_,_=score(rr.mean(1),q,y,sy)
 for col,field in enumerate(['Rg_CA_A']+FIELDS):
  ax=axes[row,col]
  for group,g in c[(c.ensemble==name)&(c.split==field)].groupby('group'):
   ax.plot(q,g.raw_intensity.to_numpy()*scale/y,label=f'Group {group+1}')
  ax.axhline(1,color='k',lw=.7);ax.axvspan(.08,.25,color='grey',alpha=.1);ax.set(ylim=(0,2),title=name+' • '+field.replace('_',' '),xlabel='q (Å⁻¹)',ylabel='Prediction / experiment');ax.legend(fontsize=7)
fig.suptitle('Conditional uniform curves at one FIXED parent-prior intensity scale per pool\nGroups are not rescaled separately. Group definitions match conditional_group_scores.csv.',fontsize=12)
for ext in ['png','pdf']:fig.savefig(O/f'conditional_SAXS_fixed_parent_scale.{ext}',dpi=160)
