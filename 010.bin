from audit_component_recovery import *

def make_plots():
 plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none'})
 fig,ax=plt.subplots(3,2,figsize=(12,12));diag,dx=plt.subplots(3,2,figsize=(12,11));numbers=[]
 for col,(name,spec) in enumerate(SPECS.items()):
  d=np.load(OUT/f'{name}_audit.npz');r=json.loads((OUT/f'{name}_audit.json').read_text());q=d['q'];mask=(q>=.012)&(q<=.25);color=['#254e9e','#c65e12'][col];label=['WT','AA'][col]
  avg=d['average'];candidate=d['candidate'];fa=r['average_guinier'];fc=r['candidate_guinier'];na=normalized(q,avg);nc=normalized(q,candidate)
  bs=np.array([normalized(q,v) for v in d['bootstrap']]);ba=np.array([normalized(q,v) for v in d['bootstrap_average']]);ratios=bs/ba
  a=ax[0,col];a.loglog(q[mask],na[mask],color='0.5',lw=2,label='Frame average, same baseline');a.loglog(q[mask],nc[mask],color=color,label='EFA main-peak candidate',lw=1.5)
  from matplotlib.ticker import FixedLocator, FuncFormatter, NullFormatter
  a.xaxis.set_major_locator(FixedLocator([.02,.05,.1,.2]));a.xaxis.set_major_formatter(FuncFormatter(lambda x,pos:f'{x:g}'));a.xaxis.set_minor_formatter(NullFormatter());a.yaxis.set_minor_formatter(NullFormatter())
  a.set(title=label+' — actual component I(q)',xlabel='q (Å⁻¹)',ylabel='I / mean I(0.06–0.10 Å⁻¹)');a.legend(fontsize=8)
  a=ax[1,col];sel=mask&(q<=.12);a.axhline(1,color='0.5',ls='--');a.plot(q[sel],nc[sel]/na[sel],color=color)
  lo,hi=np.quantile(ratios,[.025,.975],axis=0);a.fill_between(q[sel],lo[sel],hi[sel],color=color,alpha=.2)
  a.set(xlabel='q (Å⁻¹)',ylabel='Normalized candidate / average',title='Paired bootstrap band; fixed model/supports')
  a=ax[2,col]
  for tag,c,l in [('average','0.5','Average'),('candidate',color,'EFA candidate')]:
   p=np.load(OUT/f'{name}_{tag}_bift.npz');j=json.loads((OUT/f'{name}_{tag}_bift.json').read_text());norm=4*np.pi*np.trapezoid(p['p'],p['r']);a.plot(p['r'],p['p']/norm,c=c,label=f"{l}: BIFT Rg {j['rg']:.1f} Å")
  p=np.load(OUT/f'{name}_average_matched_candidate_errors_bift.npz');norm=4*np.pi*np.trapezoid(p['p'],p['r']);a.plot(p['r'],p['p']/norm,c='0.3',ls='--',label='Average, candidate error weights')
  a.set(xlabel='Pair distance r (Å)',ylabel='P(r) / I(0) (Å⁻¹)',title='Conditional P(r); not proof of purification');a.legend(fontsize=8)
  a=dx[0,col]
  for k,lab in [(0,'Other fitted component'),(1,'Main-peak candidate')]:a.plot(d['frames'],d['C'][:,k],label=lab,color='0.5' if k==0 else color)
  a.set(title=label+' — fitted elution profiles',xlabel='Frame',ylabel='Coefficient (sum = 1)');a.legend(fontsize=8)
  a=dx[1,col];a.plot(q[mask],d['other'][mask],c='0.5',label='Other fitted curve');a.fill_between(q[mask],d['other'][mask]-d['other_error'][mask],d['other'][mask]+d['other_error'][mask],color='0.5',alpha=.2);a.axhline(0,c='black',lw=.7)
  a.set(xlabel='q (Å⁻¹)',ylabel='I in average-window contribution scale',title='Inspect removed signal, including negatives')
  a=dx[2,col];im=a.imshow(d['residual'],origin='lower',aspect='auto',extent=[d['frames'][0],d['frames'][-1],0,mask.sum()],vmin=-3,vmax=3,cmap='RdBu_r');a.set(xlabel='Frame',ylabel='Increasing q-bin index',title=f"Residual / σ; mean square {r['residual_ms']:.2f}");diag.colorbar(im,ax=a)
  numbers.append(dict(sample=name,average_Rg=fa['Rg'],candidate_Rg=fc['Rg'],bootstrap_Rg_95=r['bootstrap_Rg_95'],bootstrap_delta_95=r['bootstrap_delta_Rg_vs_average_95'],s2=r['singular_values'][1],s2_null99=r['null_s2_99'],candidate_bift_Rg=json.loads((OUT/f'{name}_candidate_bift.json').read_text())['rg']))
 fig.suptitle('EFA component curves compared with the measured frame averages\nConditional candidates — monomer identity / oligomer removal not validated',fontsize=14);fig.tight_layout(rect=[0,0,1,.94])
 diag.suptitle('What EFA is separating: elution, other component and reconstruction residuals',fontsize=13);diag.tight_layout(rect=[0,0,1,.96])
 for f,stem in [(fig,'component_Iq_comparison'),(diag,'decomposition_diagnostics')]:
  f.savefig(OUT/(stem+'.png'),dpi=200);f.savefig(OUT/(stem+'.svg'));plt.close(f)
 # Common-range Guinier comparison: no change in q range between curves/samples.
 fig,ax=plt.subplots(2,2,figsize=(11,8),gridspec_kw={'height_ratios':[3,1]})
 for col,name in enumerate(SPECS):
  d=np.load(OUT/f'{name}_audit.npz');r=json.loads((OUT/f'{name}_audit.json').read_text());q=d['q'];m=(q>=.012-1e-10)&(q<=.030+1e-10)
  for key,color,lab in [('average','0.5','Average'),('candidate',['#254e9e','#c65e12'][col],'EFA candidate')]:
   f=r[key+'_guinier'];I=d[key];e=d['average_error'] if key=='average' else d['candidate_bootstrap_error'];y=np.log(I[m]/f['I0']);pred=-f['Rg']**2*q[m]**2/3
   ax[0,col].errorbar(q[m]**2,y,yerr=e[m]/I[m],fmt='.',color=color,label=f"{lab}: {f['Rg']:.2f} Å")
   ax[0,col].plot(q[m]**2,pred,color=color);ax[1,col].plot(q[m]**2,(y-pred)/(e[m]/I[m]),'.-',color=color)
  ax[0,col].set(title=name,xlabel='q² (Å⁻²)',ylabel='ln[I / fitted I(0)]');ax[0,col].legend(fontsize=8);ax[1,col].axhline(0,c='black',lw=.7);ax[1,col].set(xlabel='q² (Å⁻²)',ylabel='Residual / σ')
 fig.suptitle('Identical Guinier range: q = 0.012–0.030 Å⁻¹\nCandidate bars use bootstrap σ; full bootstrap intervals reported separately');fig.tight_layout(rect=[0,0,1,.94]);fig.savefig(OUT/'component_Guinier.png',dpi=200);plt.close(fig)
 # Sensitivity, grouped by baseline. Filter only numerical validity and elution identity.
 rows=json.loads((OUT/'rotation_grid.json').read_text());fig,ax=plt.subplots(1,2,figsize=(11,4));sensitivity=[]
 for col,name in enumerate(SPECS):
  peak=167 if col==0 else 176
  for j,b in enumerate(['pre','linear','post']):
   rr=[r for r in rows if r['sample']==name and r['baseline']==b and r.get('main_screen') and abs(r['peaks'][1]-peak)<=5 and (col==1 or r['peaks'][0]<155)]
   ys=[r['main_guinier']['Rg'] for r in rr];ax[col].scatter(j+np.linspace(-.14,.14,len(ys)),ys,alpha=.7,s=22)
   sensitivity.append(dict(sample=name,baseline=b,n=len(ys),Rg_range=[min(ys),max(ys)] if ys else None))
  ax[col].set(xticks=[0,1,2],xticklabels=['Pre only','Interpolated','Post only'],ylabel='Candidate Guinier Rg (Å)',title=name)
 fig.suptitle('Component-support and q-range sensitivity within each baseline\nSpread is model ambiguity, not a statistical confidence interval');fig.tight_layout(rect=[0,0,1,.88]);fig.savefig(OUT/'component_model_sensitivity.png',dpi=200);plt.close(fig)
 (OUT/'summary.json').write_text(json.dumps(dict(primary=numbers,sensitivity=sensitivity),indent=2))
 # Leading/central/trailing AA check with same baseline and identical fixed q interval.
 rows=[]
 for name,windows in [('WT_092025',[(155,162),(163,170),(171,178),(179,186)]),('AA_052026',[(163,171),(172,180),(181,189),(190,198)])]:
  q,Y,E,frames,A,I,err=transform(name)
  for a,b in windows:
   ix=(frames>=a)&(frames<=b);w=A[:,ix].mean(1);f=fixed_guinier(q,Y@w,np.sqrt((E**2)@(w*w)));rows.append(dict(sample=name,frames=[a,b],**f))
 (OUT/'elution_window_check.json').write_text(json.dumps(rows,indent=2))
if __name__=='__main__':make_plots()
