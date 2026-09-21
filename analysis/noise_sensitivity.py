from run_efa import *
cal=json.loads((OUT/'rg_bootstrap.json').read_text());rows=[]
for name,a,b in [('WT_092025',155,180),('AA_052026',163,204)]:
 result,data=analyze(name,a,b);q,y,e,m,D,u,s,vt=data;q,Y,E=load(name);w=e[m,a:b+1].mean(1);fit1=np.outer(u[:,0]*s[0],vt[0]);scale=next(x['background_difference_noise_scale'] for x in cal if x['sample']==name)
 rng=np.random.default_rng(20260922);vals=[]
 for j in range(300):
  n=rng.normal(size=Y.shape)*E*scale;n,_=correct(n,E,name);vals.append(np.linalg.svd(fit1+n[m,a:b+1]/w[:,None],compute_uv=False)[1])
 r=dict(sample=name,assumed_noise_multiplier=scale,observed_s2=float(s[1]),null_s2_99=float(np.quantile(vals,.99)),p_rank1=float((1+sum(x>=s[1] for x in vals))/301),interpretation='Exploratory: noise calibrated from pre-peak frames may not transfer to protein peak or account for correlations.')
 rows.append(r);print(r)
(OUT/'noise_scale_sensitivity.json').write_text(json.dumps(rows,indent=2))
